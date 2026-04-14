// Partial proportional odds (PPOM) version of final_no_horseshoe_with_centering.stan.
// Relaxes proportional odds ONLY for variant effects (beta_variant_std is now a
// V x (K-1) matrix). Sublineage effects and the fixed intercept remain proportional.
// PPOM coding follows PPOM_models/paper_model_PPOM.stan (per-threshold bernoulli_logit).

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

  real p_baseline = 0.99;
  real alpha_mean = cutpoints[1] - logit(p_baseline);

}

parameters {

  vector[L] beta_lineage;
  real<lower=0> sigma_sublineage;
  vector[S-1] z_sub;

  // Per-cutpoint standardized variant effects (PPOM).
  matrix[V, K-1] beta_variant_std;

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

  to_vector(beta_variant_std) ~ normal(0, 1);

  beta_lineage ~ normal(0, 0.1);
  sigma_sublineage ~ normal(0, 0.1);
  z_sub ~ normal(0, 1);

  {
    vector[N] sub_eta = sublineages_treatmentcontrast * beta_sublineage;
    for (n in 1:N) {
      for (k in 1:(K-1)) {
        real mu_nk = dot_product(X_std[n], beta_variant_std[, k]) + sub_eta[n] + alpha_mean;
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

}
