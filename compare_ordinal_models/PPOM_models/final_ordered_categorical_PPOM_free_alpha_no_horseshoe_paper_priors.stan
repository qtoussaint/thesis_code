// Partial proportional odds (PPOM) version of
// final_ordered_categorical_POM_free_alpha_no_horseshoe_paper_priors.stan.
// Relaxes proportional odds ONLY for variant effects (beta_variant is now a
// V x (K-1) matrix). Sublineage effects and the free alpha intercept remain proportional.
// Retained from the POM base:
//   - Fixed cutpoints at log2(mic_breakpoints)
//   - Free alpha with informative prior
//   - Raw variant_matrix (no centering/standardization), N(0,1) prior on coefficients
//   - Paper-style heritability parameterization of beta_lineage prior scale
//   - sigma_sublineage ~ normal(0, 1)  (matches paper model)
//   - Within-lineage sublineage centering, no global centering
//   - No horseshoe shrinkage

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

  ordered[K-1] cutpoints;
  cutpoints = log2(mic_breakpoints);

}

parameters {

  real alpha;

  // Heritability parameterization of beta_lineage prior scale (from paper_model_POM)
  real<lower=0, upper=1> heritability;
  real<lower=0> h_env_var;

  vector[L] beta_lineage;
  real<lower=0> sigma_sublineage;
  vector[S-1] z_sub;

  // Per-cutpoint variant effects (PPOM), simple N(0,1) prior, raw genotype.
  matrix[V, K-1] beta_variant;

}

transformed parameters {

  real<lower=0> h_gene_var = (h_env_var * heritability) / (1 - heritability);

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

}

model {

  to_vector(beta_variant) ~ normal(0, 1);

  alpha ~ normal(cutpoints[1] - 2, 1.5);

  heritability ~ beta(1, 1);
  beta_lineage ~ normal(0, h_gene_var);
  sigma_sublineage ~ normal(0, 1);
  z_sub ~ normal(0, 1);

  {
    vector[N] sub_eta = sublineages_treatmentcontrast * beta_sublineage;
    for (n in 1:N) {
      for (k in 1:(K-1)) {
        real mu_nk = alpha + dot_product(variant_matrix[n], beta_variant[, k]) + sub_eta[n];
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
  matrix[V, K-1] OR_variant_allele;
  for (k in 1:(K-1))
    for (v in 1:V)
      OR_variant_allele[v, k] = exp(beta_variant[v, k]);
}
