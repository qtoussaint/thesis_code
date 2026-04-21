// Partial proportional odds (PPOM) version of final_ordered_categorical_POM_tight_alpha.stan.
// Relaxes proportional odds ONLY for variant effects (beta_variant_std is now a
// V x (K-1) matrix with independent horseshoe local scales per (v,k)).
// Sublineage effects and the free alpha intercept remain proportional across cutpoints.
// Retained from the tight-alpha POM base:
//   - Fixed cutpoints at log2(mic_breakpoints)
//   - Data-informed, tight alpha prior anchored on empirical reference-cluster baseline
//   - Within-lineage centering only (no global sublineage centering)
//   - Regularized horseshoe on variant coefficients
//   - Diagnostic outputs (alpha_prior_mean_out, p_baseline_emp_out)

data {
  int<lower=1> N; // number of samples
  int<lower=1> V; // number of variants
  int<lower=1> L; // number of lineage clusters
  int<lower=1> S; // number of lineage subclusters
  int<lower=1> K; // number of ordered categories

  array[N] int<lower = 1, upper = K> phenotype;
  vector<lower=1e-12>[K-1] mic_breakpoints;

  matrix[N,V] variant_matrix;

  matrix[N,S] sublineage_matrix;
  array[S] int<lower = 1, upper = L> parent_lineage;

}

transformed data {

  matrix[N, S-1] sublineages_treatmentcontrast = block(sublineage_matrix, 1, 2, N, S-1);

  array[S-1] int parent_lineage_treatmentcontrast;
  for (k in 1:(S-1))
    parent_lineage_treatmentcontrast[k] = parent_lineage[k+1];

  matrix[N,V] X_ctr;
  matrix[N,V] X_std;
  vector[V] p;
  vector[V] sss;

  for (v in 1:V) {
    p[v] = mean(variant_matrix[, v]);
    for (n in 1:N)
      X_ctr[n, v] = variant_matrix[n, v] - p[v];

    sss[v] = sd(col(X_ctr, v));
    for (n in 1:N)
      X_std[n, v] = (sss[v] > 0) ? X_ctr[n, v] / sss[v] : 0;
  }

  ordered[K-1] cutpoints;
  cutpoints = log2(mic_breakpoints);

  // ---- Empirical baseline for alpha prior ---------------------------
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
  p_baseline_emp = fmin(fmax(p_baseline_emp, 0.5), 0.995);

  real alpha_prior_mean = cutpoints[1] - logit(p_baseline_emp);
  real alpha_prior_sd   = 0.5;
}

parameters {

  real alpha;

  vector[L] beta_lineage;
  real<lower=0> sigma_sublineage;
  vector[S-1] z_sub;

  real<lower=0> tau;
  real<lower=0> c2;
  matrix[V, K-1] z;
  matrix<lower=0>[V, K-1] lambda;

}

transformed parameters {

  vector[S-1] beta_sublineage_raw = beta_lineage[parent_lineage_treatmentcontrast] + sigma_sublineage * z_sub;

  vector[S-1] beta_sublineage = beta_sublineage_raw;

  for (l in 1:L) {
    real m = 0;
    int cnt = 0;
    for (s in 1:(S-1)) if (parent_lineage_treatmentcontrast[s] == l) { m += beta_sublineage_raw[s]; cnt += 1; }
    if (cnt > 0) {
      m /= cnt;
      for (s in 1:(S-1)) if (parent_lineage_treatmentcontrast[s] == l)
        beta_sublineage[s] = beta_sublineage_raw[s] - m;
    }
  }

  real sigma_latent = 5;
  real slab_scale   = 200;
  real nu           = 4;

  real m_0 = 3000;
  real tau_scale = 50;
  real tau_0;
  tau_0 = 1;

  matrix[V, K-1] beta_variant_std;
  matrix<lower=0>[V, K-1] lambda_tilde;

  for (k in 1:(K-1)) {
    for (v in 1:V) {
      lambda_tilde[v, k] =
        sqrt( (c2 * square(lambda[v, k])) /
              (c2 + square(tau) * square(lambda[v, k])) );

      beta_variant_std[v, k] = z[v, k] * tau * lambda_tilde[v, k];
    }
  }

}

model {

  c2 ~ inv_gamma(0.5 * nu, 0.5 * nu * square(slab_scale));
  tau ~ cauchy(0, tau_0);
  to_vector(z) ~ normal(0, 1);
  to_vector(lambda) ~ cauchy(0, 2);

  alpha ~ normal(alpha_prior_mean, alpha_prior_sd);

  beta_lineage ~ normal(0, 0.1);
  sigma_sublineage ~ normal(0, 0.1);
  z_sub ~ normal(0, 1);

  {
    vector[N] sub_eta = sublineages_treatmentcontrast * beta_sublineage;
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
  matrix[V, K-1] beta_variant;
  matrix[V, K-1] OR_variant_allele;

  for (k in 1:(K-1)) {
    for (v in 1:V) {
      if (sss[v] > 0) {
        beta_variant[v, k] = beta_variant_std[v, k] / sss[v];
        OR_variant_allele[v, k] = exp(beta_variant[v, k]);
      } else {
        beta_variant[v, k] = 0;
        OR_variant_allele[v, k] = 1;
      }
    }
  }

  // Posterior predictive MIC category frequencies
  array[N] int y_rep;
  vector[K] cat_freq_rep = rep_vector(0, K);
  {
    vector[N] sub_eta_gen = sublineages_treatmentcontrast * beta_sublineage;
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

  real alpha_prior_mean_out = alpha_prior_mean;
  real p_baseline_emp_out   = p_baseline_emp;

  matrix[V, K-1] beta_variant_std_prior;
  for (k in 1:(K-1))
    for (v in 1:V)
      beta_variant_std_prior[v, k] = z[v, k] * tau * lambda_tilde[v, k];
}
