// Ordered association model with lineage subclusters
// Variant effects use a simple weakly informative normal prior (no horseshoe)

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

  // log2-scale the MIC breakpoints
  // use the transformed breakpoints as fixed cutpoints
  ordered[K-1] cutpoints;
  cutpoints = log2(mic_breakpoints);

  real p_baseline = 0.99;
  real alpha_mean = cutpoints[1] - logit(p_baseline);

}

parameters {

  // basic parameters
  vector[L] beta_lineage; // effect sizes of lineage clusters
  real<lower=0> sigma_sublineage; // variance in sublineage effect sizes
  vector[S-1] z_sub;              // std normals, noncentered

  // variant effects with weakly informative normal prior
  vector[V] beta_variant_raw;

}

transformed parameters {

  vector[S-1] beta_sublineage_raw = beta_lineage[parent_lineage_treatmentcontrast] + sigma_sublineage * z_sub;

  // Optional but recommended: sum-to-zero within each parent lineage
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

  vector[N] mu; // likelihood

  // weakly informative normal prior on variant effects
  beta_variant_raw ~ normal(0, 1);

  // POPULATION STRUCTURE CORRECTION

  // lineage prior
  beta_lineage ~ normal(0, 0.1);

  // variance of lineage subclusters from their parent lineages' mean
  sigma_sublineage ~ normal(0, 0.1);  // half-normal with sd 0.1
  z_sub ~ normal(0, 1);

  // LINEAR PREDICTOR

  // likelihood statement
  mu = (variant_matrix * beta_variant_raw) + (sublineages_treatmentcontrast * beta_sublineage) + alpha_mean;

  // fit with ordered logistic
  phenotype ~ ordered_logistic(mu, cutpoints);

}

generated quantities {
  vector[V] beta_variant;
  vector[V] OR_variant_allele;

  for (v in 1:V) {
    beta_variant[v] = beta_variant_raw[v];          // already on 0->1 allele scale
    OR_variant_allele[v]   = exp(beta_variant[v]); // cumulative OR
  }

}
