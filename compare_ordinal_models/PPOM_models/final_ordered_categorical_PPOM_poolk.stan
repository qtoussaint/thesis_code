// DIAGNOSTIC PPOM variant: partial pooling of variant effects across cutpoints.
//
// Baseline tight_alpha_tau1 treats z_variant[v, k] as independent across k so
// per-variant threshold-specific effects inherit no regularisation from each
// other. This variant decomposes
//
//   z_variant[v, k] = z_base[v] + sigma_dev * z_dev[v, k]
//
// with sigma_dev ~ normal+(0, 0.3). sigma_dev -> 0 reduces to proportional
// odds (POM); sigma_dev unbounded reduces to the independent tight_alpha_tau1.
// A half-normal(0, 0.3) prior corresponds to "threshold deviations of at most
// a few standardised units around the variant-level anchor", which is a
// reasonable a-priori cap given MIC data rarely supports genuine dose-response
// gradients that large.
//
// Horseshoe is kept on the per-variant anchor z_base[v] (not on the deviations
// z_dev[v, k]) so sparsity is enforced at the variant level -- a variant that
// is zero at one cutpoint but huge at the next is unlikely biology and should
// pay a prior cost.
//
// If switching to this variant collapses the late-k inflation and the POM and
// PPOM beta estimates line up closely, mechanism (2) (no borrowing across k)
// is dominant. If effects still inflate at late k despite the pooling, the
// problem is likely mechanism (3) (cutpoint scale misfit) instead.

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

  // Horseshoe on the per-variant anchor (shared across k)
  real<lower=0> tau;
  real<lower=0> c2;
  vector[V] z_base;
  vector<lower=0>[V] lambda_variant;

  // Per-(v, k) deviations from the anchor, partial-pooled
  matrix[V, K-1] z_dev;
  real<lower=0> sigma_dev;
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

  real slab_scale = 5;   // tight slab -- redundant with sigma_dev cap but keeps the anchor sane
  real nu         = 4;
  real tau_0      = 1;

  vector<lower=0>[V] lambda_tilde_variant;
  vector[V] beta_anchor_std;   // per-variant (anchor) effect, shared across k
  for (v in 1:V) {
    lambda_tilde_variant[v] =
      sqrt( (c2 * square(lambda_variant[v])) /
            (c2 + square(tau) * square(lambda_variant[v])) );
    beta_anchor_std[v] = z_base[v] * tau * lambda_tilde_variant[v];
  }

  matrix[V, K-1] beta_variant_std;
  for (k in 1:(K-1)) {
    for (v in 1:V) {
      beta_variant_std[v, k] = beta_anchor_std[v] + sigma_dev * z_dev[v, k];
    }
  }
}

model {
  c2 ~ inv_gamma(0.5 * nu, 0.5 * nu * square(slab_scale));
  tau ~ cauchy(0, tau_0);
  z_base ~ normal(0, 1);
  lambda_variant ~ cauchy(0, 2);

  sigma_dev ~ normal(0, 0.3);
  to_vector(z_dev) ~ normal(0, 1);

  alpha ~ normal(alpha_prior_mean, alpha_prior_sd);

  beta_lineage ~ normal(0, 0.1);
  sigma_sublineage ~ normal(0, 0.1);
  z_sub ~ normal(0, 1);

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
  matrix[V, K-1] beta_variant;
  matrix[V, K-1] OR_variant_allele;
  vector[V]       beta_variant_anchor;   // unstandardised per-variant anchor
  for (v in 1:V) {
    if (sd_variant[v] > 0) {
      beta_variant_anchor[v] = beta_anchor_std[v] / sd_variant[v];
    } else {
      beta_variant_anchor[v] = 0;
    }
  }
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

  array[N_ppc] int y_rep_ppc;
  array[N_ppc] int y_true_ppc;
  {
    vector[N] sub_eta_gen = X_sublineage * beta_sublineage;
    for (i in 1:N_ppc) {
      int n = ppc_idx[i];
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
      y_rep_ppc[i]  = categorical_rng(probs);
      y_true_ppc[i] = phenotype[n];
    }
  }

  real alpha_prior_mean_out = alpha_prior_mean;
  real p_baseline_emp_out   = p_baseline_emp;
  real sigma_dev_out        = sigma_dev;

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
