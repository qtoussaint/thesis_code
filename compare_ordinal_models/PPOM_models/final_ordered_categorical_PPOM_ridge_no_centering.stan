// Ridge (Gaussian, sigma=5) variant of final_ordered_categorical_PPOM, using
// raw genotype X (no centering/scaling). Replaces the regularized horseshoe
// on beta_variant with a normal(0, 5) prior. Architecture: fixed cutpoints
// from log2(mic_breakpoints), fixed alpha_mean anchored intercept, raw
// variant_matrix in linear predictor, deterministic-subset PPC, per-cutpoint
// heritability metrics.

data {
  int<lower=1> N;
  int<lower=1> V;
  int<lower=1> L;
  int<lower=1> S;
  int<lower=1> K;

  array[N] int<lower = 1, upper = K> phenotype;
  vector<lower=1e-12>[K-1] mic_breakpoints;

  matrix[N,V] variant_matrix;

  matrix[N,S] sublineage_matrix;
  array[S] int<lower = 1, upper = L> parent_lineage;

  int<lower=0, upper=N> N_ppc;
  array[N_ppc] int<lower=1, upper=N> ppc_idx;
}

transformed data {

  matrix[N, S-1] sublineages_treatmentcontrast = block(sublineage_matrix, 1, 2, N, S-1);

  array[S-1] int parent_lineage_treatmentcontrast;
  for (k in 1:(S-1))
    parent_lineage_treatmentcontrast[k] = parent_lineage[k+1];

  ordered[K-1] cutpoints;
  cutpoints = log2(mic_breakpoints);

  real p_baseline = 0.99;
  real alpha_mean = cutpoints[1] - logit(p_baseline);

}

parameters {

  vector[L] beta_lineage;
  real<lower=0> sigma_sublineage;
  vector[S-1] z_sub;

  matrix[V, K-1] beta_variant_raw;

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

  beta_sublineage -= mean(beta_sublineage);

}

model {

  to_vector(beta_variant_raw) ~ normal(0, 5);

  beta_lineage ~ normal(0, 0.1);
  sigma_sublineage ~ normal(0, 0.1);
  z_sub ~ normal(0, 1);

  {
    vector[N] sub_eta = sublineages_treatmentcontrast * beta_sublineage;
    for (n in 1:N) {
      for (k in 1:(K-1)) {
        real mu_nk = dot_product(variant_matrix[n], beta_variant_raw[, k]) + sub_eta[n] + alpha_mean;
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
      beta_variant[v, k] = beta_variant_raw[v, k];
      OR_variant_allele[v, k] = exp(beta_variant[v, k]);
    }
  }

  // Deterministic-subset PPC
  array[N_ppc] int y_rep_ppc;
  array[N_ppc] int y_true_ppc;
  {
    vector[N] sub_eta_gen = sublineages_treatmentcontrast * beta_sublineage;
    for (i in 1:N_ppc) {
      int n = ppc_idx[i];
      vector[K-1] cdf_n;
      for (k in 1:(K-1)) {
        real mu_nk = dot_product(variant_matrix[n], beta_variant_raw[, k]) + sub_eta_gen[n] + alpha_mean;
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

  // Per-cutpoint heritability on the liability (logistic-latent) scale
  vector<lower=0>[K-1] V_A_k;
  real<lower=0> V_pop;
  real<lower=0> V_E = pi()^2 / 3;
  vector<lower=0, upper=1>[K-1] h2_narrow_k;
  vector<lower=0, upper=1>[K-1] h2_broad_k;
  real<lower=0, upper=1> h2_narrow_median_k;
  real<lower=0, upper=1> h2_broad_median_k;
  {
    vector[N] g_pop = sublineages_treatmentcontrast * beta_sublineage;
    V_pop = variance(g_pop);
    for (k in 1:(K-1)) {
      vector[N] g_variant_k = variant_matrix * beta_variant_raw[, k];
      V_A_k[k] = variance(g_variant_k);
      real V_tot_k = V_A_k[k] + V_pop + V_E;
      h2_narrow_k[k] = V_A_k[k] / V_tot_k;
      h2_broad_k[k]  = (V_A_k[k] + V_pop) / V_tot_k;
    }
    h2_narrow_median_k = quantile(h2_narrow_k, 0.5);
    h2_broad_median_k  = quantile(h2_broad_k,  0.5);
  }

}
