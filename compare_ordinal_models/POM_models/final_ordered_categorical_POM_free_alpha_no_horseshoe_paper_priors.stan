// Ordered association model with lineage subclusters
//
// Changes from final_ordered_categorical_POM_free_alpha_no_horseshoe.stan:
//
//   Lineage and sublineage prior *distributions* now match paper_model_POM.stan:
//     - beta_lineage ~ normal(0, h_gene_var), where
//         h_gene_var = (h_env_var * heritability) / (1 - heritability)
//       with free parameters heritability ~ beta(1,1) and h_env_var > 0.
//     - sigma_sublineage ~ normal(0, 1)   (half-normal with sd 1, was sd 0.1)
//
//   The structural choices from the free_alpha_no_horseshoe model are retained:
//     - Non-centered parameterization for sublineage effects (z_sub).
//     - Within-lineage sum-to-zero centering of beta_sublineage.
//     - Treatment-contrast sublineage encoding (reference column dropped).
//     - Fixed cutpoints = log2(mic_breakpoints).
//     - Free alpha intercept.
//     - N(0,1) prior on beta_variant; variant_matrix used directly (no centering/standardization).

data {
  int<lower=1> N; // number of samples
  int<lower=1> V; // number of variants
  int<lower=1> L; // number of lineage clusters
  int<lower=1> S; // number of lineage subclusters
  int<lower=1> K; // number of ordered categories

  array[N] int<lower = 1, upper = K> phenotype; // MICs expressed as ordered categories
  vector<lower=1e-12>[K-1] mic_breakpoints;  // concentration breakpoints used to define K (must be >0)

  matrix[N,V] variant_matrix; // genotype

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

  // Fixed cutpoints — experimental design constants
  ordered[K-1] cutpoints;
  cutpoints = log2(mic_breakpoints);

}

parameters {

  // Intercept: latent position of the reference sublineage at zero variant effects.
  real alpha;

  // Heritability parameterization of beta_lineage prior scale (from paper_model_POM)
  real<lower=0, upper=1> heritability;
  real<lower=0> h_env_var;

  // basic parameters
  vector[L] beta_lineage; // effect sizes of lineage clusters
  real<lower=0> sigma_sublineage; // variance in sublineage effect sizes
  vector[S-1] z_sub;              // std normals, noncentered

  // Variant effects — simple N(0,1) prior, no horseshoe, no standardization
  vector[V] beta_variant;

}

transformed parameters {

  real<lower=0> h_gene_var = (h_env_var * heritability) / (1 - heritability);

  vector[S-1] beta_sublineage_raw = beta_lineage[parent_lineage_treatmentcontrast] + sigma_sublineage * z_sub;

  // Sum-to-zero within each parent lineage (partial pooling / identifiability within lineage groups).
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

  vector[N] mu;

  // VARIANT EFFECTS — simple N(0,1) prior
  beta_variant ~ normal(0, 1);

  // INTERCEPT
  alpha ~ normal(cutpoints[1] - 2, 1.5);

  // POPULATION STRUCTURE CORRECTION — prior distributions from paper_model_POM.stan
  heritability ~ beta(1, 1);
  beta_lineage ~ normal(0, h_gene_var);
  sigma_sublineage ~ normal(0, 1);  // half-normal with sd 1
  z_sub ~ normal(0, 1);

  // LINEAR PREDICTOR
  mu = (variant_matrix * beta_variant) + (sublineages_treatmentcontrast * beta_sublineage) + alpha;

  // fit with ordered logistic
  phenotype ~ ordered_logistic(mu, cutpoints);

}

generated quantities {
  // beta_variant is already on the per 0->1 allele scale since variant_matrix
  // is used directly (no centering/standardization).
  vector[V] OR_variant_allele = exp(beta_variant);
}
