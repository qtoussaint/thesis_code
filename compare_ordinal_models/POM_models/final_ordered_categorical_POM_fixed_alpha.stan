// Ordered association model with lineage subclusters
//
// Variant of final_ordered_categorical_POM_free_alpha.stan in which the
// intercept is NOT estimated. Instead, alpha is fixed at
//     alpha_mean = cutpoints[1] - logit(0.99)
// (the heuristic from the original final_ordered_categorical_POM.stan).
// All other improvements from _free_alpha are retained:
//   - cutpoints fixed at log2(mic_breakpoints)
//   - global sublineage mean-centering REMOVED (only within-lineage centering kept)
//   - treatment-contrast encoding of sublineage matrix
// Use this variant when you want to compare a fixed-alpha baseline against
// the free-alpha model on otherwise-identical machinery.

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

  // log2-scale the MIC breakpoints
  // use the transformed breakpoints as fixed cutpoints
  // These are experimental design constants — do not estimate them
  ordered[K-1] cutpoints;
  cutpoints = log2(mic_breakpoints);

  // Fixed intercept: places baseline at ~0.99 probability below the lowest cutpoint
  real p_baseline = 0.99;
  real alpha_mean = cutpoints[1] - logit(p_baseline);

}

parameters {

  // basic parameters
  vector[L] beta_lineage; // effect sizes of lineage clusters
  real<lower=0> sigma_sublineage; // variance in sublineage effect sizes
  vector[S-1] z_sub;              // std normals, noncentered

  // REGULARIZED HORSESHOE FOR BETA_VARIANTS

  // the following is a widely used calibration (from Piironen & Vehtari's horseshoe prior work)

  // following their recommendation, because we only known the best distribution for a GLM not other models like an ordcat,
  // start with half-cauchy on (0, tau_0) and then manually change values to optimize prior predictive check
  // see eq. 3.12

  // number of predictors (V)
  // number of samples (N)


  // global shrinkage (half-Cauchy)
  real<lower=0> tau;

  real<lower=0> c2; // slab variance

  vector[V] z;
  vector<lower=0>[V] lambda;


}

transformed parameters {

  vector[S-1] beta_sublineage_raw = beta_lineage[parent_lineage_treatmentcontrast] + sigma_sublineage * z_sub;

  // Sum-to-zero within each parent lineage (partial pooling / identifiability within lineage groups).
  // NOTE: global mean-centering (beta_sublineage -= mean(beta_sublineage)) has been removed.
  // That step conflicted with treatment-contrast encoding: after global centering, the reference
  // sublineage (all-zeros rows) acquired an implicit effect of -mean(beta_sublineage), making
  // the intercept uninterpretable. Within-lineage centering is retained because it is
  // hierarchically meaningful and does not affect the reference sublineage.
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
  // parent lineage's mean. The reference sublineage's effect is absorbed into alpha_mean.


    // ---- calibration constants (as in your original transformed parameters) ----
  real sigma_latent = 5;     // recommended latent scale for logistic-type models
  real slab_scale   = 200;
  real nu           = 4;

  // prior guess of relevant predictors
  real m_0 = 3000;
  real tau_scale=50;
  real tau_0;
  tau_0 = 0.1;



  vector[V] beta_variant_std; // variant effects

  // local shrinkage parameters

  vector<lower=0>[V] lambda_tilde;

  for (v in 1:V) {

    // regularized local scale (Piironen & Vehtari)
    lambda_tilde[v] =
      sqrt( (c2 * square(lambda[v])) /
            (c2 + square(tau) * square(lambda[v])) );

    // regularized horseshoe prior draw for coefficient
    beta_variant_std[v] = z[v] * tau * lambda_tilde[v];
  }

}

model {

  vector[N] mu; // likelihood

  c2 ~ inv_gamma(0.5 * nu, 0.5 * nu * square(slab_scale));

  tau ~ cauchy(0, tau_0);

  to_vector(z) ~ normal(0,1);
  to_vector(lambda) ~ cauchy(0, 2); // half-Cauchy(0,2)


  // POPULATION STRUCTURE CORRECTION

  // lineage prior
  beta_lineage ~ normal(0, 0.1);

  // variance of lineage subclusters from their parent lineages' mean
  sigma_sublineage ~ normal(0, 0.1);  // half-normal with sd 0.1
  z_sub ~ normal(0, 1);


  // LINEAR PREDICTOR
  // mu = variant effects + sublineage population correction + fixed intercept
  // alpha_mean: cutpoints[1] - logit(0.99), anchors baseline at ~99% susceptibility
  // beta_sublineage: within-lineage-centered deviations; reference sublineage effect = 0 by encoding
  mu = (X_std * beta_variant_std) + (sublineages_treatmentcontrast * beta_sublineage) + alpha_mean;

  // fit with ordered logistic
  phenotype ~ ordered_logistic(mu, cutpoints);

}

generated quantities {
  vector[V] beta_variant;
  vector[V] OR_variant_allele;

  for (v in 1:V) {
    if (sss[v] > 0) {
      beta_variant[v] = beta_variant_std[v] / sss[v];   // per 0->1 allele
      OR_variant_allele[v]   = exp(beta_variant[v]); // cumulative OR
    } else {
      beta_variant[v] = 0;
      OR_variant_allele[v]   = 1;
    }
  }

  // HERITABILITY

  // latent residual variance for logistic link is fixed at pi^2 / 3
  // genetic variance Vg computed from the variant component of the linear predictor
  // so: h2 = Vg / (Vg + Ve)

  //real Ve = pi()^2 / 3;
  //vector[N] eta_var;
  //real Vg;
  //real h2_latent;

  // for (n in 1:N) {
    //   eta_var[n] = (variant_matrix_ctr[n] * beta_variant_std);
    //}
  // Vg = variance(eta_var);

  //h2_latent = Vg / (Vg + Ve);


  // POSTERIOR PREDICTIVE CHECKS

  // for WAIC/LOO
  //array[N_train] real log_lik;

  // simulate phenotypes from likelihood for prior/posterior predictive checks
  //array[N] real y_rep;

  //for (n in 1:N) {
    //log_lik[n] = ordered_logistic_lpmf(training_phenotypes[n] | mu[n], cutpoints);
    //y_rep[n] = ordered_logistic_rng(mu[n], cutpoints);
    //}

  // -- ADDITIONAL DIAGNOSTICS FOR PRIOR PREDICTION --


  vector[V] beta_variant_std_prior;

  for (v in 1:V) {
      beta_variant_std_prior[v] = z[v] * tau * lambda_tilde[v];
  }

  //vector[N] eta_variant;
  //vector[N] eta_sub;

  //for (n in 1:N) {
    // eta_variant[n] = (variant_matrix_std[n] * beta_variant_std);
    // eta_sub[n]     = sublineage_matrix[n] * beta_sublineage;
    //}

  //real sd_eta_variant = sd(eta_variant);
  //real sd_eta_sub     = sd(eta_sub);
  //real sd_mu          = sd(to_vector(mu));

}
