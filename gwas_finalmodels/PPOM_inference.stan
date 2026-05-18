// Partial proportional-odds (PPOM) variant of POM_inference.stan.

// Differs from POM_inference.stan in exactly one place: beta_variant_std is a
// V x (K-1) matrix with an independent horseshoe scale per (variant, cutpoint),
// so proportional odds is relaxed for variant effects only. Lineage and
// sublineage effects remain proportional across cutpoints.

// All other structure is inherited from POM_inference.stan: data-informed
// alpha_prior_mean used to shift the cutpoint anchor, treatment-contrast
// sublineage encoding with within-lineage sum-to-zero, etc.

data {
  int<lower=1> N;                                  // samples
  int<lower=1> V;                                  // variants
  int<lower=1> L;                                  // lineage clusters
  int<lower=1> S;                                  // lineage subclusters
  int<lower=1> K;                                  // ordered categories

  array[N] int<lower=1, upper=K> phenotype;
  vector<lower=1e-12>[K-1] mic_breakpoints;

  matrix[N, V] variant_matrix;
  matrix[N, S] sublineage_matrix;
  array[S] int<lower=1, upper=L> parent_lineage;

  int<lower=0, upper=N> N_ppc;                     // PPC subset size (20% of N, capped at 500)
  array[N_ppc] int<lower=1, upper=N> ppc_idx;      // deterministic indices into 1:N for PPC
}

transformed data {
  // Treatment-contrast encoding: drop the reference sublineage (column 1)
  matrix[N, S-1] X_sublineage = block(sublineage_matrix, 1, 2, N, S-1);

  array[S-1] int parent_lineage_sub;
  for (k in 1:(S-1))
    parent_lineage_sub[k] = parent_lineage[k+1];

  // Center and scale the genotype matrix (sd=0 columns pinned to 0)
  matrix[N, V] X_std;
  vector[V] mean_variant;
  vector[V] sd_variant;
  for (v in 1:V) {
    mean_variant[v] = mean(variant_matrix[, v]);
    vector[N] centred = variant_matrix[, v] - mean_variant[v];
    sd_variant[v] = sd(centred);
    for (n in 1:N)
      X_std[n, v] = (sd_variant[v] > 0) ? centred[n] / sd_variant[v] : 0;
  }

  // MIC breakpoints on the log2 scale (used to build the cutpoint prior)
  vector[K-1] mic_cutpoints = log2(mic_breakpoints);

  // Empirical baseline: Laplace-smoothed fraction of reference-sublineage
  // samples in MIC category 1. You should pick the sublineage with the lowest
  // average phenotype (e.g. most susceptible) as the reference during data preprocessing.
  // The [0.5, 0.995] clamp just guards against logit(1) = inf
  int n_ref = 0;
  int n_ref_cat1 = 0;
  for (n in 1:N) {
    if (sublineage_matrix[n, 1] > 0.5) {
      n_ref += 1;
      if (phenotype[n] == 1) n_ref_cat1 += 1;
    }
  }
  real p_baseline_emp = (n_ref_cat1 + 0.5) / (n_ref + 1.0);
  p_baseline_emp = fmin(fmax(p_baseline_emp, 0.5), 0.995);

  real alpha_prior_mean = mic_cutpoints[1] - logit(p_baseline_emp);

  // Shift the MIC anchors down by alpha_prior_mean so cutpoint 1 sits
  // near logit(p_baseline_emp)
  vector[K-1] cutpoint_prior_mean;
  for (k in 1:(K-1))
    cutpoint_prior_mean[k] = mic_cutpoints[k] - alpha_prior_mean;

  // The wide interval prevents misspecification from a non-uniform latent-scale
  real cutpoint_prior_sd = 1.5;
}

parameters {
  // Estimated cutpoints
  // ordered[K-1] type enforces strict ordering and cutpoint_prior_mean /
  // cutpoint_prior_sd move cutpoint location within that constraint
  ordered[K-1] cutpoints;

  vector[L] beta_lineage;
  real<lower=0> sigma_sublineage;
  vector[S-1] z_sub;

  // Regularized horseshoe, independent per (variant, cutpoint)
  real<lower=0> tau;
  real<lower=0> c2;
  matrix[V, K-1] z_variant;
  matrix<lower=0>[V, K-1] lambda_variant;
}

transformed parameters {
  vector[S-1] beta_sublineage_raw =
    beta_lineage[parent_lineage_sub] + sigma_sublineage * z_sub;

  // Within-lineage sum-to-zero centering (global centering conflicts with
  // treatment-contrast encoding and makes the cutpoint anchor uninterpretable)
  vector[S-1] beta_sublineage = beta_sublineage_raw;
  for (l in 1:L) {
    real m = 0;
    int cnt = 0;
    for (s in 1:(S-1)) if (parent_lineage_sub[s] == l) { m += beta_sublineage_raw[s]; cnt += 1; }
    if (cnt > 0) {
      m /= cnt;
      for (s in 1:(S-1)) if (parent_lineage_sub[s] == l)
        beta_sublineage[s] = beta_sublineage_raw[s] - m;
    }
  }

  // Horseshoe hyperparameters (Piironen & Vehtari 2017, eq. 3.12).
  real slab_scale = 5;
  real nu         = 4;
  real tau_0      = 1;

  matrix[V, K-1] beta_variant_std;
  matrix<lower=0>[V, K-1] lambda_tilde_variant;
  for (k in 1:(K-1)) {
    for (v in 1:V) {
      lambda_tilde_variant[v, k] =
        sqrt( (c2 * square(lambda_variant[v, k])) /
              (c2 + square(tau) * square(lambda_variant[v, k])) );
      beta_variant_std[v, k] = z_variant[v, k] * tau * lambda_tilde_variant[v, k];
    }
  }
}

model {
  c2 ~ inv_gamma(0.5 * nu, 0.5 * nu * square(slab_scale));
  tau ~ cauchy(0, tau_0);
  to_vector(z_variant) ~ normal(0, 1);
  to_vector(lambda_variant) ~ cauchy(0, 2);

  // Wide normal prior on each cutpoint around its shifted MIC breakpoint anchor
  cutpoints ~ normal(cutpoint_prior_mean, cutpoint_prior_sd);

  beta_lineage ~ normal(0, 0.1);
  sigma_sublineage ~ normal(0, 0.1);
  z_sub ~ normal(0, 1);

  // PPOM likelihood -- one bernoulli per (sample, cutpoint)
  {
    vector[N] sub_eta = X_sublineage * beta_sublineage;
    for (n in 1:N) {
      for (k in 1:(K-1)) {
        real mu_nk = dot_product(X_std[n], beta_variant_std[, k]) + sub_eta[n];
        if (phenotype[n] <= k) {
          target += bernoulli_logit_lpmf(1 | cutpoints[k] - mu_nk);
        } else {
          target += bernoulli_logit_lpmf(0 | cutpoints[k] - mu_nk);
        }
      }
    }
  }
}

generated quantities {
  // Unstandardized variant effects (per 0->1 allele) and cumulative odds ratios,
  // indexed by cutpoint k
  matrix[V, K-1] beta_variant;
  matrix[V, K-1] OR_variant_allele;
  for (k in 1:(K-1)) {
    for (v in 1:V) {
      if (sd_variant[v] > 0) {
        beta_variant[v, k] = beta_variant_std[v, k] / sd_variant[v];
        OR_variant_allele[v, k] = exp(beta_variant[v, k]);
      } else {
        beta_variant[v, k] = 0;
        OR_variant_allele[v, k] = 1;
      }
    }
  }

  // Per-cutpoint drift vs the shifted MIC-grid anchor. Interpretable as how
  // much each cutpoint has moved from its "known" location after accounting
  // for the level shift in cutpoint_prior_mean; large drifts at specific k
  // pinpoint categories where the MIC anchor may be uninformative/misspecified
  vector[K-1] cutpoint_drift;
  for (k in 1:(K-1))
    cutpoint_drift[k] = cutpoints[k] - cutpoint_prior_mean[k];

  // Posterior predictive check on the specified subset of N
  // y_true_ppc just included to make downstream analysis less annoying
  array[N_ppc] int y_rep_ppc;
  array[N_ppc] int y_true_ppc;
  {
    vector[N] sub_eta_gen = X_sublineage * beta_sublineage;
    for (i in 1:N_ppc) {
      int n = ppc_idx[i];
      vector[K-1] cdf_n;
      for (k in 1:(K-1)) {
        real mu_nk = dot_product(X_std[n], beta_variant_std[, k]) + sub_eta_gen[n];
        cdf_n[k] = inv_logit(cutpoints[k] - mu_nk);
      }
      vector[K] probs;
      probs[1] = cdf_n[1];
      for (k in 2:(K-1)) probs[k] = cdf_n[k] - cdf_n[k-1];
      probs[K] = 1 - cdf_n[K-1];
      for (k in 1:K) if (probs[k] < 0) probs[k] = 0;
      probs /= sum(probs);
      y_rep_ppc[i]  = categorical_rng(probs);
      y_true_ppc[i] = phenotype[n];
    }
  }

  // Expose priors for post-hoc verification
  real alpha_prior_mean_out = alpha_prior_mean;
  real p_baseline_emp_out   = p_baseline_emp;

  matrix[V, K-1] beta_variant_std_prior;
  for (k in 1:(K-1))
    for (v in 1:V)
      beta_variant_std_prior[v, k] = z_variant[v, k] * tau * lambda_tilde_variant[v, k];

  // Heritability, on the liability (logistic-latent) scale, per cutpoint k
  
  // Some tips on interpreting these:
  
    // The reason we get per-threshold (K-1) additive genetic score variances/heritabilities
    // is because beta_variant_std is V x (K-1) under partial proportional odds
  
    // h2_narrow_k[k] / h2_broad_k[k] are per-threshold heritability
    // h2_*_median_k is a scalar summary for tables etc.
  
    // h2 differing across K can help you understand which resistance breakpoints are better
    // predicted by the model (and therefore possibly more genetically determined/meaningful)

    // Horseshoe shrinkage biases V_A_k downward, so h2_narrow_k is a lower bound
  
  vector<lower=0>[K-1] V_A_k;
  real<lower=0> V_pop;
  real<lower=0> V_E = pi()^2 / 3;
  vector<lower=0, upper=1>[K-1] h2_narrow_k;
  vector<lower=0, upper=1>[K-1] h2_broad_k;
  real<lower=0, upper=1> h2_narrow_median_k;
  real<lower=0, upper=1> h2_broad_median_k;
  {
    vector[N] g_pop = X_sublineage * beta_sublineage;
    V_pop = variance(g_pop);
    for (k in 1:(K-1)) {
      vector[N] g_variant_k = X_std * beta_variant_std[, k];
      V_A_k[k] = variance(g_variant_k);
      real V_tot_k = V_A_k[k] + V_pop + V_E;
      h2_narrow_k[k] = V_A_k[k] / V_tot_k;
      h2_broad_k[k]  = (V_A_k[k] + V_pop) / V_tot_k;
    }
    h2_narrow_median_k = quantile(h2_narrow_k, 0.5);
    h2_broad_median_k  = quantile(h2_broad_k,  0.5);
  }
}
