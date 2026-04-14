// Ordered association model with lineage subclusters
//
// Changes from final_ordered_categorical_POM_free_alpha.stan:
//
//   1. HORSESHOE REMOVED: The regularized horseshoe prior on variant effects has been
//      replaced with a simple N(0,1) prior on beta_variant.
//
//   2. VARIANT CENTERING/STANDARDIZATION REMOVED: variant_matrix is used directly in
//      the linear predictor (no X_ctr, no X_std). beta_variant is therefore already
//      on the per-allele (0->1) scale on the raw genotype, so OR_variant_allele is
//      just exp(beta_variant) with no rescaling.
//
// Unchanged aspects:
//
//   1. CUTPOINTS: Fixed at log2(mic_breakpoints). These are experimental design
//      constants (dilution series thresholds), not biological parameters.
//
//   2. INTERCEPT: Free parameter `alpha` with an informative prior. Identification
//      relies on fixed cutpoints + free intercept (canonical parameterization).
//
//   3. SUBLINEAGE ENCODING: treatment-contrast (drop reference column). Within-lineage
//      centering is retained; global centering is not applied (it would conflict with
//      treatment-contrast encoding and make alpha's prior uninterpretable).

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
  // These are experimental design constants — do not estimate them
  ordered[K-1] cutpoints;
  cutpoints = log2(mic_breakpoints);

}

parameters {

  // Intercept: latent position of the reference sublineage at zero variant effects.
  real alpha;

  // basic parameters
  vector[L] beta_lineage; // effect sizes of lineage clusters
  real<lower=0> sigma_sublineage; // variance in sublineage effect sizes
  vector[S-1] z_sub;              // std normals, noncentered

  // Variant effects — simple N(0,1) prior, no horseshoe, no standardization
  vector[V] beta_variant;

}

transformed parameters {

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

  // beta_sublineage[s] now represents the deviation of sublineage s from its
  // parent lineage's mean. The reference sublineage's effect is captured by alpha.

}

model {

  vector[N] mu; // likelihood

  // VARIANT EFFECTS — simple N(0,1) prior (no horseshoe shrinkage)
  beta_variant ~ normal(0, 1);

  // INTERCEPT
  alpha ~ normal(cutpoints[1] - 2, 1.5);

  // POPULATION STRUCTURE CORRECTION

  // lineage prior
  beta_lineage ~ normal(0, 0.1);

  // variance of lineage subclusters from their parent lineages' mean
  sigma_sublineage ~ normal(0, 0.1);  // half-normal with sd 0.1
  z_sub ~ normal(0, 1);


  // LINEAR PREDICTOR
  // mu = variant effects (raw genotype) + sublineage population correction + intercept
  mu = (variant_matrix * beta_variant) + (sublineages_treatmentcontrast * beta_sublineage) + alpha;

  // fit with ordered logistic
  phenotype ~ ordered_logistic(mu, cutpoints);

}

generated quantities {
  // beta_variant is already on the per 0->1 allele scale since variant_matrix
  // is used directly (no centering/standardization).
  vector[V] OR_variant_allele = exp(beta_variant);
}
