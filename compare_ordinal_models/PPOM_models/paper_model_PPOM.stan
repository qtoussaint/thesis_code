// Ordered association model with lineage subclusters

data {
  int<lower=1> K; // K is number of ordinal categories
  int<lower=1> N;
  int<lower=1> V;
  int<lower=1> L; // number of lineage clusters
  int<lower=1> S; // number of sublineage clusters
  array[N] int<lower = 1, upper = K> phenotype;
  vector<lower=1e-12>[K-1] mic_breakpoints;  // concentration breakpoints used to define K (must be >0)
  matrix[N,V] variant_matrix;
  matrix[N,S] sublineage_matrix; // lineage subclusters (full one-hot, reference in column 1)
  array[S] int<lower = 1, upper = L> parent_lineage; // vector in same order as subcluster matrix, containing parent lineage of that cluster
}

transformed data {
  // Drop reference sublineage (column 1) — treatment-contrast encoding
  matrix[N, S-1] sublineages_treatmentcontrast = block(sublineage_matrix, 1, 2, N, S-1);

  // Drop reference sublineage from parent_lineage mapping
  array[S-1] int parent_lineage_treatmentcontrast;
  for (k in 1:(S-1))
    parent_lineage_treatmentcontrast[k] = parent_lineage[k+1];
}

parameters {
  ordered[K-1] cutpoints; // for ordinal regression
  real<lower=0, upper=1> heritability;
  real<lower=0> h_env_var;
  matrix[V, K-1] beta_variant; // made into matrix of variants x cutpoints
  vector[L] beta_lineage; // effect sizes of lineage clusters
  vector[S-1] beta_sublineage; // effect sizes of lineage subclusters
  real<lower=0> sigma_sublineage; // variance in sublineage effect sizes
}

transformed parameters {
  real<lower=0> h_gene_var = (h_env_var * heritability) / (1 - heritability);
}

model {

  // define priors
  to_vector(beta_variant) ~ normal(0, 5);  // flattened to vector
  beta_lineage ~ normal(0, h_gene_var);
  sigma_sublineage ~ normal(0,1); // variance from parent lineage mean
  heritability ~ beta(1,1);
  cutpoints ~ normal(log2(mic_breakpoints), 0.5); // prior centered on log2(MIC breakpoints), reduces VI instability

  // each sublineage EF is centered around its parent lineage dist with variance of sigma_sublineage
  for (k in 1:(S-1)) {
    beta_sublineage[k] ~ normal(beta_lineage[parent_lineage_treatmentcontrast[k]], sigma_sublineage);
  }

  // (nonproportional odds for variant effects, proportional effects for sublins)
  for (n in 1:N) {
    for (k in 1:(K-1)) {
      // for each cutpoint, compute linear predictor with that cutpoint's variant effects
      real mu_k = dot_product(variant_matrix[n], beta_variant[,k]) +
                  dot_product(sublineages_treatmentcontrast[n], beta_sublineage);

      // cumulative logit for this cutpoint
      // the ordered_logistic_lpmf expects a vector of cutpoints,
      // so instead we add up binary logit likelihoods for each threshold
      if (phenotype[n] <= k) {
        target += bernoulli_logit_lpmf(1 | cutpoints[k] - mu_k);
      } else {
        target += bernoulli_logit_lpmf(0 | cutpoints[k] - mu_k);
      }
    }
  }

}


