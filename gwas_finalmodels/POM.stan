// ORIGINALLY: final_ordered_categorical_POM_tight_alpha.stan

// Ordered association model with lineage subclusters
//
// Changes from final_ordered_categorical_POM_free_alpha.stan:
//
//   1. DATA-INFORMED ALPHA PRIOR: The previous `alpha ~ normal(cutpoints[1]-2, 1.5)`
//      was too loose (±3 log2-dilutions of flex), letting alpha absorb variant-level
//      signal. Here alpha is anchored on the *empirical* fraction of reference-sublineage
//      samples that sit in MIC category 1:
//          p_baseline_emp    = #(reference rows in cat 1) / #(reference rows)
//          alpha_prior_mean  = cutpoints[1] - logit(p_baseline_emp)
//          alpha_prior_sd    = 0.5   (~1 log2-dilution at 95%, matches MIC rounding)
//      This works because `gwas_datasets/utils.R:select_reference_sublineage()` picks
//      the minimum-mean-phenotype sublineage upstream, so column 1 of sublineage_matrix
//      is guaranteed to be a deliberately-susceptible reference. Laplace smoothing
//      (+0.5/+1) handles sparse reference clusters; a small fallback applies if the
//      reference cluster has <5 samples.
//
//   2. POSTERIOR PREDICTIVE CATEGORY FREQUENCIES: generated quantities now produces
//      `y_rep` and `cat_freq_rep` so the posterior predictive MIC histogram can be
//      compared directly to the observed histogram (the direct check on whether
//      alpha is in the right place).
//
//   3. DIAGNOSTIC OUTPUTS: `alpha_prior_mean_out` and `p_baseline_emp_out` are exposed
//      in generated quantities for post-hoc verification against an R-side hand
//      computation over the reference-sublineage subset.
//
// Unchanged from free_alpha: fixed cutpoints at log2(mic_breakpoints), within-lineage
// sum-to-zero sublineage centering (no global centering), regularized horseshoe on
// beta_variant with tau_0 = 0.1.

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

  // ---- Empirical baseline for alpha prior ---------------------------
  // The reference sublineage (column 1 of sublineage_matrix) is chosen upstream
  // to be the minimum-mean-phenotype sublineage (see gwas_datasets/utils.R
  // select_reference_sublineage), so it is deliberately susceptible. We compute
  // the empirical fraction of reference-sublineage samples in MIC category 1 and
  // use it to anchor alpha.
  int n_ref = 0;
  int n_ref_cat1 = 0;
  for (n in 1:N) {
    if (sublineage_matrix[n, 1] > 0.5) {
      n_ref += 1;
      if (phenotype[n] == 1) n_ref_cat1 += 1;
    }
  }
  real p_baseline_emp;
  if (n_ref >= 5)
    p_baseline_emp = (n_ref_cat1 + 0.5) / (n_ref + 1.0);  // Laplace smoothing
  else
    p_baseline_emp = 0.9;                                 // fallback for tiny ref clusters
  // Clamp to avoid degenerate logit at 0 or 1
  p_baseline_emp = fmin(fmax(p_baseline_emp, 0.5), 0.995);

  real alpha_prior_mean = cutpoints[1] - logit(p_baseline_emp);
  real alpha_prior_sd   = 0.5;  // ≈ 1 log2-dilution at 95% — matches MIC rounding
}

parameters {

  // Intercept: latent position of the reference sublineage at zero variant effects.
  // Prior is tight and anchored on the empirical reference-cluster MIC distribution
  // (computed in transformed data). Tight enough to block alpha from absorbing
  // variant-level signal; loose enough to accommodate sampling noise in the
  // reference cluster.
  real alpha;

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
  // NOTE: global mean-centering is intentionally absent (it would conflict with the treatment-contrast
  // encoding and make alpha's prior uninterpretable).
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


  // INTERCEPT — tight, data-informed prior anchored on the empirical reference-cluster baseline
  alpha ~ normal(alpha_prior_mean, alpha_prior_sd);

  // POPULATION STRUCTURE CORRECTION

  // lineage prior
  beta_lineage ~ normal(0, 0.1);

  // variance of lineage subclusters from their parent lineages' mean
  sigma_sublineage ~ normal(0, 0.1);  // half-normal with sd 0.1
  z_sub ~ normal(0, 1);


  // LINEAR PREDICTOR
  mu = (X_std * beta_variant_std) + (sublineages_treatmentcontrast * beta_sublineage) + alpha;

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

  // Posterior predictive MIC category frequencies — direct check on alpha placement.
  array[N] int y_rep;
  vector[K] cat_freq_rep = rep_vector(0, K);
  {
    vector[N] mu_gen = (X_std * beta_variant_std)
                       + (sublineages_treatmentcontrast * beta_sublineage)
                       + alpha;
    for (n in 1:N) y_rep[n] = ordered_logistic_rng(mu_gen[n], cutpoints);
    for (n in 1:N) cat_freq_rep[y_rep[n]] += 1;
    cat_freq_rep /= N;
  }

  // Expose the data-informed prior inputs for post-hoc verification in R.
  real alpha_prior_mean_out = alpha_prior_mean;
  real p_baseline_emp_out   = p_baseline_emp;

  vector[V] beta_variant_std_prior;
  for (v in 1:V) {
      beta_variant_std_prior[v] = z[v] * tau * lambda_tilde[v];
  }

}
