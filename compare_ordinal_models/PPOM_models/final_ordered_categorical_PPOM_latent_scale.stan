// DIAGNOSTIC PPOM variant: free latent-scale parameter.
//
// Baseline tight_alpha_tau1 fixes the logistic-latent scale at 1 and pins
// cutpoints at log2(mic_breakpoints). If the empirical MIC dispersion in a
// given dataset does not match a standard logistic (scale 1) at those fixed
// cutpoints, nothing in the baseline model can absorb the spacing mismatch
// except the per-(v, k) variant effects -- which is the "late-k inflation"
// failure mode described in mechanism (3).
//
// This variant adds a single free parameter sigma_latent that rescales the
// latent (logistic) scale relative to the log2-MIC cutpoint spacing:
//
//   P(y <= k | n) = logit^-1( (cutpoints[k] - mu_nk) / sigma_latent )
//
// sigma_latent = 1 recovers the baseline. sigma_latent > 1 means the cutpoints
// on the logit scale are effectively compressed (wide MIC distribution); < 1
// means they are stretched. Variant effects are interpreted as shifts on the
// same rescaled latent scale.
//
// Prior: sigma_latent ~ normal+(1, 0.3). Tight at 1 to prevent identifiability
// drift against alpha's tight prior, but wide enough to absorb the factor-of-2
// mismatches typical of MIC datasets.
//
// If switching to this variant collapses late-k beta magnitudes and posterior
// sigma_latent clearly exceeds 1, mechanism (3) is present. If sigma_latent
// sits near 1 and effects still inflate, mechanism (3) can be ruled out and
// attention goes to (1) and (2).

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

  ordered[K-1] cutpoints = log2(mic_breakpoints);

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
  matrix[V, K-1] z_variant;
  matrix<lower=0>[V, K-1] lambda_variant;

  // Free latent-scale parameter (the diagnostic addition)
  real<lower=0.1> sigma_latent;
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

  sigma_latent ~ normal(1, 0.3);  // tight-near-1, prevents drift vs alpha

  alpha ~ normal(alpha_prior_mean, alpha_prior_sd);

  beta_lineage ~ normal(0, 0.1);
  sigma_sublineage ~ normal(0, 0.1);
  z_sub ~ normal(0, 1);

  {
    vector[N] sub_eta = X_sublineage * beta_sublineage;
    for (n in 1:N) {
      for (k in 1:(K-1)) {
        real mu_nk = alpha + dot_product(X_std[n], beta_variant_std[, k]) + sub_eta[n];
        real logit_arg = (cutpoints[k] - mu_nk) / sigma_latent;
        if (phenotype[n] <= k) {
          target += bernoulli_logit_lpmf(1 | logit_arg);
        } else {
          target += bernoulli_logit_lpmf(0 | logit_arg);
        }
      }
    }
  }
}

generated quantities {
  matrix[V, K-1] beta_variant;
  matrix[V, K-1] OR_variant_allele;
  // Latent-scale-adjusted effects: beta on the 'standard logistic' scale, directly
  // comparable to baseline tight_alpha_tau1's beta_variant.
  matrix[V, K-1] beta_variant_on_unit_latent;
  for (k in 1:(K-1)) {
    for (v in 1:V) {
      if (sd_variant[v] > 0) {
        beta_variant[v, k] = beta_variant_std[v, k] / sd_variant[v];
        OR_variant_allele[v, k] = exp(beta_variant[v, k]);
        beta_variant_on_unit_latent[v, k] = beta_variant[v, k] / sigma_latent;
      } else {
        beta_variant[v, k] = 0;
        OR_variant_allele[v, k] = 1;
        beta_variant_on_unit_latent[v, k] = 0;
      }
    }
  }

  array[N_ppc] int y_rep_ppc;
  array[N_ppc] int y_true_ppc;
  {
    vector[N] sub_eta_gen = X_sublineage * beta_sublineage;
    for (i in 1:N_ppc) {
      int n = ppc_idx[i];
      vector[K-1] cdf_n;
      for (k in 1:(K-1)) {
        real mu_nk = alpha + dot_product(X_std[n], beta_variant_std[, k]) + sub_eta_gen[n];
        cdf_n[k] = inv_logit((cutpoints[k] - mu_nk) / sigma_latent);
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
  real sigma_latent_out     = sigma_latent;

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
