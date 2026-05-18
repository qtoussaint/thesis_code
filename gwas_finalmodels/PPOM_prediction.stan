// Partial proportional-odds model (PPOM) for phenotype prediction.

// Fits on N_train samples with the same likelihood, priors, free cutpoints, and
// generated summaries as the inference model, then generates predicted_phenotype on
// N_test held-out samples as an N_test x K probability matrix

// Needed to use the per-category probability draws for the predicted phenotype
// as there's no categorical RNG that supports PPOMs

// Differs from POM_prediction.stan in exactly one place: beta_variant_std is a
// V x (K-1) matrix with an independent horseshoe scale per (variant, cutpoint),
// so proportional odds is relaxed for variant effects only. Lineage and
// sublineage effects remain proportional across cutpoints.

data {
  int<lower=1> N_train;                                          // training samples
  int<lower=1> N_test;                                           // test samples
  int<lower=1> V;                                                // variants
  int<lower=1> L;                                                // lineage clusters
  int<lower=1> S;                                                // lineage subclusters
  int<lower=1> K;                                                // ordered categories

  array[N_train] int<lower=1, upper=K> training_phenotype;
  vector<lower=1e-12>[K-1] mic_breakpoints;

  matrix[N_train, V] training_variants;
  matrix[N_test,  V] test_variants;
  matrix[N_train, S] training_sublineages;
  matrix[N_test,  S] test_sublineages;
  array[S] int<lower=1, upper=L> parent_lineage;
}

transformed data {
  // Treatment-contrast encoding: drop the reference sublineage (column 1)
  matrix[N_train, S-1] X_sublineage_train =
    block(training_sublineages, 1, 2, N_train, S-1);
  matrix[N_test,  S-1] X_sublineage_test  =
    block(test_sublineages,     1, 2, N_test,  S-1);

  array[S-1] int parent_lineage_sub;
  for (k in 1:(S-1))
    parent_lineage_sub[k] = parent_lineage[k+1];

  // Standardize genotype using TRAINING column mean/SD; apply the same
  // transform to test. (sd=0 columns pinned to 0)
  matrix[N_train, V] X_std_train;
  matrix[N_test,  V] X_std_test;
  vector[V] mean_variant;
  vector[V] sd_variant;
  for (v in 1:V) {
    mean_variant[v] = mean(training_variants[, v]);
    vector[N_train] centred_train = training_variants[, v] - mean_variant[v];
    sd_variant[v] = sd(centred_train);
    for (n in 1:N_train)
      X_std_train[n, v] = (sd_variant[v] > 0) ? centred_train[n] / sd_variant[v] : 0;
    for (n in 1:N_test)
      X_std_test[n, v] = (sd_variant[v] > 0)
        ? (test_variants[n, v] - mean_variant[v]) / sd_variant[v] : 0;
  }

  // MIC breakpoints on the log2 scale (used to build the cutpoint prior)
  vector[K-1] mic_cutpoints = log2(mic_breakpoints);

  // Empirical baseline: Laplace-smoothed fraction of reference-sublineage
  // samples in MIC category 1. You should pick the sublineage with the lowest
  // average phenotype (e.g. most susceptible) as the reference during data preprocessing.
  // The [0.5, 0.995] clamp just guards against logit(1) = inf
  int n_ref = 0;
  int n_ref_cat1 = 0;
  for (n in 1:N_train) {
    if (training_sublineages[n, 1] > 0.5) {
      n_ref += 1;
      if (training_phenotype[n] == 1) n_ref_cat1 += 1;
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

  // PPOM likelihood -- one bernoulli per (sample, cutpoint), training only
  {
    vector[N_train] sub_eta = X_sublineage_train * beta_sublineage;
    for (n in 1:N_train) {
      for (k in 1:(K-1)) {
        real mu_nk = dot_product(X_std_train[n], beta_variant_std[, k]) + sub_eta[n];
        if (training_phenotype[n] <= k) {
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

  // Held-out per-category probabilities (since no categorical RNG exists for PPOMs in Stan)
  matrix[N_test, K] predicted_phenotype;
  {
    vector[N_test] sub_eta_test = X_sublineage_test * beta_sublineage;
    for (n in 1:N_test) {
      vector[K-1] cdf_n;
      for (k in 1:(K-1)) {
        real mu_nk = dot_product(X_std_test[n], beta_variant_std[, k]) + sub_eta_test[n];
        cdf_n[k] = inv_logit(cutpoints[k] - mu_nk);
      }
      vector[K] probs;
      probs[1] = cdf_n[1];
      for (k in 2:(K-1)) probs[k] = cdf_n[k] - cdf_n[k-1];
      probs[K] = 1 - cdf_n[K-1];
      // Prevent tiny negative probabilities due to rounding
      for (k in 1:K) if (probs[k] < 1e-12) probs[k] = 1e-12;
      probs /= sum(probs);
      for (k in 1:K) predicted_phenotype[n, k] = probs[k];
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
    vector[N_train] g_pop = X_sublineage_train * beta_sublineage;
    V_pop = variance(g_pop);
    for (k in 1:(K-1)) {
      vector[N_train] g_variant_k = X_std_train * beta_variant_std[, k];
      V_A_k[k] = variance(g_variant_k);
      real V_tot_k = V_A_k[k] + V_pop + V_E;
      h2_narrow_k[k] = V_A_k[k] / V_tot_k;
      h2_broad_k[k]  = (V_A_k[k] + V_pop) / V_tot_k;
    }
    h2_narrow_median_k = quantile(h2_narrow_k, 0.5);
    h2_broad_median_k  = quantile(h2_broad_k,  0.5);
  }
}
