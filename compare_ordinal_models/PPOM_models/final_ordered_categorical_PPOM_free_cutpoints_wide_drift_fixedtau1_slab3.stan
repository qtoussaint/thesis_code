// DIAGNOSTIC PPOM variant: free_cutpoints + wide cutpoint drift + FIXED-tau horseshoe.
//
// Identical to final_ordered_categorical_PPOM_free_cutpoints_wide_drift.stan EXCEPT
// the global horseshoe scale tau is FIXED at 0.05 instead of being sampled,
// slab_scale is reduced from 5 -> 3, and lambda's Cauchy prior is tightened
// from scale 2 -> 1.
//
// Motivation: in the wide_drift model (and tau5/tau5_slab50 cousins), tau ~
// cauchy(0, tau_0) collapses to ~5e-4 because data prefer near-universal
// shrinkage. With tau that small, the slab cap never engages
// (tau * lambda << c) so beta = z * tau * lambda for everyone -- the slab is
// dead and every effect is in the spike. Loosening the half-Cauchy prior on
// tau cannot rescue this: heavy tails in BOTH directions mean the data still
// pulls tau down. Anchoring tau at a target spike width fixes the spike at
// 2 sigma ~= [-0.1, 0.1] for "off" variants while letting Cauchy(0,1) lambda
// outliers drive their corresponding beta toward the slab cap c. slab_scale=3
// caps max|beta| near ~3, matching the no-HS comparison model's empirical max.

data {
  int<lower=1> N;
  int<lower=1> V;
  int<lower=1> L;
  int<lower=1> S;
  int<lower=1> K;

  array[N] int<lower=1, upper=K> phenotype;
  vector<lower=1e-12>[K-1] mic_breakpoints;

  matrix[N, V] variant_matrix;
  matrix[N, S] sublineage_matrix;
  array[S] int<lower=1, upper=L> parent_lineage;

  int<lower=0, upper=N> N_ppc;
  array[N_ppc] int<lower=1, upper=N> ppc_idx;
}

transformed data {
  matrix[N, S-1] X_sublineage = block(sublineage_matrix, 1, 2, N, S-1);

  array[S-1] int parent_lineage_sub;
  for (k in 1:(S-1))
    parent_lineage_sub[k] = parent_lineage[k+1];

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

  // MIC-grid cutpoint anchors (no longer the cutpoints themselves; used for prior)
  vector[K-1] mic_cutpoints = log2(mic_breakpoints);

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

  // Shift MIC-grid anchors down by alpha_prior_mean so cutpoints[1] is centred
  // near logit(p_baseline_emp). Preserves MIC-ladder spacing.
  vector[K-1] cutpoint_prior_mean;
  for (k in 1:(K-1))
    cutpoint_prior_mean[k] = mic_cutpoints[k] - alpha_prior_mean;

  // Inherits wide_drift's 1.5 cutpoint prior SD (admits implied non-uniform misfit)
  real cutpoint_prior_sd = 1.5;
}

parameters {
  // Estimated cutpoints (replacing fixed cutpoints + free alpha)
  ordered[K-1] cutpoints;

  vector[L] beta_lineage;
  real<lower=0> sigma_sublineage;
  vector[S-1] z_sub;

  // Horseshoe with FIXED tau: tau is no longer a parameter.
  real<lower=0> c2;
  matrix[V, K-1] z_variant;
  matrix<lower=0>[V, K-1] lambda_variant;
}

transformed parameters {
  vector[S-1] beta_sublineage_raw =
    beta_lineage[parent_lineage_sub] + sigma_sublineage * z_sub;

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

  real tau        = 1;  // FIXED -- sets spike width [-0.1, 0.1] at 2 sigma
  real slab_scale = 3;     // caps max|beta| near ~3 (matches no-HS empirical max)
  real nu         = 4;

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
  // no prior on tau -- fixed at 0.05 in transformed parameters
  to_vector(z_variant) ~ normal(0, 1);
  to_vector(lambda_variant) ~ cauchy(0, 1);

  // Wider prior on each cutpoint around its MIC-grid anchor. The ordered[K-1]
  // type still enforces ordering; this prior shapes each cutpoint's location
  // within that constraint, now with enough room to absorb the implied
  // non-uniform latent-scale misfit (~log2(5) ~ 2.3 units).
  cutpoints ~ normal(cutpoint_prior_mean, cutpoint_prior_sd);

  beta_lineage ~ normal(0, 0.1);
  sigma_sublineage ~ normal(0, 0.1);
  z_sub ~ normal(0, 1);

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

  // Per-cutpoint drift vs MIC-grid anchor (interpretable: how much has each
  // cutpoint moved from its "known" location, after accounting for the
  // alpha_prior_mean level shift). With the wider prior, large drifts at
  // specific k pinpoint exactly which categories the MIC anchor is wrong about.
  vector[K-1] cutpoint_drift;
  for (k in 1:(K-1))
    cutpoint_drift[k] = cutpoints[k] - cutpoint_prior_mean[k];

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

  real alpha_prior_mean_out = alpha_prior_mean;
  real p_baseline_emp_out   = p_baseline_emp;

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
