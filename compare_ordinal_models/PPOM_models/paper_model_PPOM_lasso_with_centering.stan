// Lasso (Laplace) variant of paper_model_PPOM, using centered+scaled X_std.
// Replaces the normal(0, 5) prior on beta_variant with double_exponential(0,
// 1) AND swaps raw variant_matrix for centered+scaled X_std in the linear
// predictor. Architecture otherwise: cutpoints estimated as parameters with
// N(log2(MIC), 0.5) prior, scalar lineage-level heritability, no alpha_mean.
// Adds N_ppc/ppc_idx data + deterministic-subset PPC + per-cutpoint
// liability-scale heritability in generated quantities. beta_variant is
// reported on both standardized and natural (allele-scale) scales.

data {
  int<lower=1> K; // K is number of ordinal categories
  int<lower=1> N;
  int<lower=1> V;
  int<lower=1> L; // number of lineage clusters
  int<lower=1> S; // number of sublineage clusters
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
}

parameters {
  ordered[K-1] cutpoints;
  real<lower=0, upper=1> heritability;
  real<lower=0> h_env_var;
  matrix[V, K-1] beta_variant_std;
  vector[L] beta_lineage;
  vector[S-1] beta_sublineage;
  real<lower=0> sigma_sublineage;
}

transformed parameters {
  real<lower=0> h_gene_var = (h_env_var * heritability) / (1 - heritability);
}

model {

  to_vector(beta_variant_std) ~ double_exponential(0, 1);
  beta_lineage ~ normal(0, h_gene_var);
  sigma_sublineage ~ normal(0,1);
  heritability ~ beta(1,1);
  cutpoints ~ normal(log2(mic_breakpoints), 0.5);

  for (k in 1:(S-1)) {
    beta_sublineage[k] ~ normal(beta_lineage[parent_lineage_treatmentcontrast[k]], sigma_sublineage);
  }

  for (n in 1:N) {
    for (k in 1:(K-1)) {
      real mu_k = dot_product(X_std[n], beta_variant_std[,k]) +
                  dot_product(sublineages_treatmentcontrast[n], beta_sublineage);

      if (phenotype[n] <= k) {
        target += bernoulli_logit_lpmf(1 | cutpoints[k] - mu_k);
      } else {
        target += bernoulli_logit_lpmf(0 | cutpoints[k] - mu_k);
      }
    }
  }

}

generated quantities {
  // Beta on natural (allele) scale, plus odds ratio
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

  // Deterministic-subset PPC
  array[N_ppc] int y_rep_ppc;
  array[N_ppc] int y_true_ppc;
  {
    for (i in 1:N_ppc) {
      int n = ppc_idx[i];
      vector[K-1] cdf_n;
      for (k in 1:(K-1)) {
        real mu_k = dot_product(X_std[n], beta_variant_std[, k]) +
                    dot_product(sublineages_treatmentcontrast[n], beta_sublineage);
        cdf_n[k] = inv_logit(cutpoints[k] - mu_k);
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

  // Per-cutpoint variant-level heritability on the liability (logistic-latent)
  // scale. Distinct from the scalar `heritability` parameter (lineage-level).
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
