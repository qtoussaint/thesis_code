// Ordered association model with lineage subclusters

data {
  int<lower=1> K; // K is number of ordinal categories
  int<lower=1> N;
  int<lower=1> V;
  int<lower=1> L; // number of lineage clusters
  int<lower=1> S; // number of sublineage clusters
  array[N] int<lower = 1, upper = K> phenotype;
  matrix[N,V] variant_matrix;
  matrix[N,L-1] lineage_matrix; // lineage clusters
  matrix[N,S-1] sublineage_matrix; // lineage subclusters
  array[S-1] int<lower = 1, upper = L-1> parent_lineage; // vector in same order as subcluster matrix, containing parent lineage of that cluster
}

parameters {
  ordered[K-1] cutpoints; // for ordinal regression
  real<lower=0, upper=1> heritability;
  real<lower=0> h_env_var;
  vector[V] beta_variant;
  vector[L-1] beta_lineage; // effect sizes of lineage clusters
  vector[S-1] beta_sublineage; // effect sizes of lineage subclusters
  real<lower=0> sigma_sublineage; // variance in sublineage effect sizes
}

transformed parameters {
  real<lower=0> h_gene_var = (h_env_var * heritability) / (1 - heritability);
}

model {
  //create the linear predictor
  vector[N] mu =  variant_matrix * beta_variant + sublineage_matrix * beta_sublineage;

  //define priors
  beta_variant ~ normal(0, 5);
  beta_lineage ~ normal(0, h_gene_var);
  sigma_sublineage ~ normal(0,1); // variance from parent lineage mean
  heritability ~ beta(1,1);
  cutpoints ~ normal([-4.97, -3.94, -2.32, 1], 0.5); // reduces VI instability (cutpoints are not on the same scale as ordcats -- this is a good prior given my integer categories)
  // cutpoint priors are log2(MIC breakpoints)

  // each sublineage EF is centered around its parent lineage dist with variance of sigma_sublineage
  for (k in 1:(S-1)) {
    beta_sublineage[k] ~ normal(beta_lineage[parent_lineage[k]], sigma_sublineage);
  }

  //write the likelihood
  target += ordered_logistic_lpmf(phenotype | mu, cutpoints);
}


