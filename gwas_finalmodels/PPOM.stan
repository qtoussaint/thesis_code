// Partial proportional-odds (PPOM) variant of POM.stan.
//
// Differs from POM.stan in exactly one place: beta_variant_std is a V x (K-1)
// matrix with an independent horseshoe scale per (variant, cutpoint), so
// proportional odds is relaxed for variant effects only. Alpha, lineage, and
// sublineage effects remain proportional across cutpoints.
//
// All other structure is inherited from POM.stan: fixed cutpoints at
// log2(mic_breakpoints), data-informed tight alpha prior anchored on the
// reference sublineage, treatment-contrast sublineage encoding with
// within-lineage sum-to-zero, regularized horseshoe on variant effects.

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
}

transformed data {
  // Treatment-contrast encoding: drop the reference sublineage (column 1)
  matrix[N, S-1] X_sublineage = block(sublineage_matrix, 1, 2, N, S-1);

  array[S-1] int parent_lineage_sub;
  for (k in 1:(S-1))
    parent_lineage_sub[k] = parent_lineage[k+1];

  // Centre and scale the genotype matrix (sd=0 columns pinned to 0)
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

  // Fixed cutpoints on the log2-MIC scale
  ordered[K-1] cutpoints = log2(mic_breakpoints);

  // Empirical baseline for the alpha prior: fraction of reference-sublineage
  // samples in MIC category 1. Laplace smoothing handles sparse reference
  // clusters & a 0.9 fallback applies when the reference cluster has <5 samples
  int n_ref = 0;
  int n_ref_cat1 = 0;
  for (n in 1:N) {
    if (sublineage_matrix[n, 1] > 0.5) {
      n_ref += 1;
      if (phenotype[n] == 1) n_ref_cat1 += 1;
    }
  }
  real p_baseline_emp;
  if (n_ref >= 5)
    p_baseline_emp = (n_ref_cat1 + 0.5) / (n_ref + 1.0);
  else
    p_baseline_emp = 0.9;
  // clamp to avoid a degenerate logit at 0 or 1
  p_baseline_emp = fmin(fmax(p_baseline_emp, 0.5), 0.995);

  real alpha_prior_mean = cutpoints[1] - logit(p_baseline_emp);
  real alpha_prior_sd   = 0.5;  // roughly one log2-dilution at 95%, matches MIC rounding
}

parameters {
  real alpha;

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

  // Within-lineage sum-to-zero centering (global centering conflicts
  // with treatment-contrast encoding and make alpha uninterpretable)
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

  // Horseshoe hyperparameters (Piironen & Vehtari 2017, eq. 3.12)
  real slab_scale = 200;
  real nu         = 4;
  real tau_0      = 0.1;

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

  // Tight and data-informed intercept prior to prevent alpha from absorbing
  // variant-level signal while still accounting for reference-cluster noise
  alpha ~ normal(alpha_prior_mean, alpha_prior_sd);

  beta_lineage ~ normal(0, 0.1);
  sigma_sublineage ~ normal(0, 0.1);
  z_sub ~ normal(0, 1);

  // PPOM likelihood -- one bernoulli per (sample, cutpoint)
  {
    vector[N] sub_eta = X_sublineage * beta_sublineage;
    for (n in 1:N) {
      for (k in 1:(K-1)) {
        real mu_nk = alpha + dot_product(X_std[n], beta_variant_std[, k]) + sub_eta[n];
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
  // indexed by cutpoint k.
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

  // Posterior predictive MIC histogram 
  array[N] int y_rep;
  vector[K] cat_freq_rep = rep_vector(0, K);
  {
    vector[N] sub_eta_gen = X_sublineage * beta_sublineage;
    for (n in 1:N) {
      vector[K-1] cdf_n;
      for (k in 1:(K-1)) {
        real mu_nk = alpha + dot_product(X_std[n], beta_variant_std[, k]) + sub_eta_gen[n];
        cdf_n[k] = inv_logit(cutpoints[k] - mu_nk);
      }
      vector[K] probs;
      probs[1] = cdf_n[1];
      for (k in 2:(K-1)) probs[k] = cdf_n[k] - cdf_n[k-1];
      probs[K] = 1 - cdf_n[K-1];
      for (k in 1:K) if (probs[k] < 0) probs[k] = 0;
      probs /= sum(probs);
      y_rep[n] = categorical_rng(probs);
    }
    for (n in 1:N) cat_freq_rep[y_rep[n]] += 1;
    cat_freq_rep /= N;
  }

  // Expose priors for post-hoc verification
  real alpha_prior_mean_out = alpha_prior_mean;
  real p_baseline_emp_out   = p_baseline_emp;

  matrix[V, K-1] beta_variant_std_prior;
  for (k in 1:(K-1))
    for (v in 1:V)
      beta_variant_std_prior[v, k] = z_variant[v, k] * tau * lambda_tilde_variant[v, k];
}
