// Proportional-odds ordinal logistic model (POM) for phenotype prediction 
//
// Fits on N_train samples with the same likelihood, priors, free cutpoints, and
// generated summaries as the inference model, then generates predicted_phenotype on
// N_test held-out samples as integer ordinal draws.

data {
  int<lower=1> N_train;                                          // training samples
  int<lower=1> N_test;                                           // test samples
  int<lower=1> V;                                                // variants
  int<lower=1> L;                                                // lineage clusters
  int<lower=1> S;                                                // lineage subclusters
  int<lower=1> K;                                                // ordered categories

  array[N_train] int<lower=1, upper=K> training_phenotype;       // MIC category (train)
  vector<lower=1e-12>[K-1] mic_breakpoints;                      // concentration breakpoints

  matrix[N_train, V] training_variants;
  matrix[N_test,  V] test_variants;
  matrix[N_train, S] training_sublineages;
  matrix[N_test,  S] test_sublineages;
  array[S] int<lower=1, upper=L> parent_lineage;
}

transformed data {
  // Treatment-contrast encoding: drop the reference sublineage (column 1)
  matrix[N_train, S-1] X_sublineage_train =
    block(training_sublineages, 1, 2, N_train, S-1);
  matrix[N_test,  S-1] X_sublineage_test  =
    block(test_sublineages,     1, 2, N_test,  S-1);

  array[S-1] int parent_lineage_sub;
  for (k in 1:(S-1))
    parent_lineage_sub[k] = parent_lineage[k+1];

  // Standardize genotype using TRAINING column mean/SD; apply the same
  // transform to test. (sd=0 columns pinned to 0)
  matrix[N_train, V] X_std_train;
  matrix[N_test,  V] X_std_test;
  vector[V] mean_variant;
  vector[V] sd_variant;
  for (v in 1:V) {
    mean_variant[v] = mean(training_variants[, v]);
    vector[N_train] centred_train = training_variants[, v] - mean_variant[v];
    sd_variant[v] = sd(centred_train);
    for (n in 1:N_train)
      X_std_train[n, v] = (sd_variant[v] > 0) ? centred_train[n] / sd_variant[v] : 0;
    for (n in 1:N_test)
      X_std_test[n, v] = (sd_variant[v] > 0)
        ? (test_variants[n, v] - mean_variant[v]) / sd_variant[v] : 0;
  }

  // MIC breakpoints on the log2 scale (used to build the cutpoint prior)
  vector[K-1] mic_cutpoints = log2(mic_breakpoints);

  // Empirical baseline: Laplace-smoothed fraction of reference-sublineage
  // samples in MIC category 1. You should pick the sublineage with the lowest
  // average phenotype (e.g. most susceptible) as the reference during data preprocessing.
  // The [0.5, 0.995] clamp just guards against logit(1) = inf
  int n_ref = 0;
  int n_ref_cat1 = 0;
  for (n in 1:N_train) {
    if (training_sublineages[n, 1] > 0.5) {
      n_ref += 1;
      if (training_phenotype[n] == 1) n_ref_cat1 += 1;
    }
  }
  real p_baseline_emp = (n_ref_cat1 + 0.5) / (n_ref + 1.0);
  p_baseline_emp = fmin(fmax(p_baseline_emp, 0.5), 0.995);

  real alpha_prior_mean = mic_cutpoints[1] - logit(p_baseline_emp);

  // Shift the MIC anchors down by alpha_prior_mean so cutpoint 1 sits
  // near logit(p_baseline_emp)
  vector[K-1] cutpoint_prior_mean;
  for (k in 1:(K-1))
    cutpoint_prior_mean[k] = mic_cutpoints[k] - alpha_prior_mean;

  // The wide interval prevents misspecification from a non-uniform latent-scale
  real cutpoint_prior_sd = 1.5;
}

parameters {
  // Estimated cutpoints
  // ordered[K-1] type enforces strict ordering and cutpoint_prior_mean /
  // cutpoint_prior_sd move cutpoint location within that constraint
  ordered[K-1] cutpoints;

  vector[L] beta_lineage;
  real<lower=0> sigma_sublineage;
  vector[S-1] z_sub;

  // Regularized horseshoe on variant effects (per-variant penalty is pooled across cutpoints)
  real<lower=0> tau;
  real<lower=0> c2;
  vector[V] z_variant;
  vector<lower=0>[V] lambda_variant;
}

transformed parameters {
  vector[S-1] beta_sublineage_raw =
    beta_lineage[parent_lineage_sub] + sigma_sublineage * z_sub;

  // Within-lineage sum-to-zero centering (global centering conflicts with
  // treatment-contrast encoding and makes the cutpoint anchor uninterpretable)
  vector[S-1] beta_sublineage = beta_sublineage_raw;
  for (l in 1:L) {
    real m = 0;
    int cnt = 0;
    for (s in 1:(S-1)) if (parent_lineage_sub[s] == l) { m += beta_sublineage_raw[s]; cnt += 1; }
    if (cnt > 0) {
      m /= cnt;
      for (s in 1:(S-1)) if (parent_lineage_sub[s] == l)
        beta_sublineage[s] = beta_sublineage_raw[s] - m;
    }
  }

  // Horseshoe hyperparameters (Piironen & Vehtari 2017, eq. 3.12).
  real slab_scale = 5;
  real nu         = 4;
  real tau_0      = 1;

  vector[V] beta_variant_std;
  vector<lower=0>[V] lambda_tilde_variant;
  for (v in 1:V) {
    lambda_tilde_variant[v] =
      sqrt( (c2 * square(lambda_variant[v])) /
            (c2 + square(tau) * square(lambda_variant[v])) );
    beta_variant_std[v] = z_variant[v] * tau * lambda_tilde_variant[v];
  }
}

model {
  c2 ~ inv_gamma(0.5 * nu, 0.5 * nu * square(slab_scale));
  tau ~ cauchy(0, tau_0);
  z_variant ~ normal(0, 1);
  lambda_variant ~ cauchy(0, 2);

  // Wide normal prior on each cutpoint around its shifted MIC breakpoint anchor
  cutpoints ~ normal(cutpoint_prior_mean, cutpoint_prior_sd);

  beta_lineage ~ normal(0, 0.1);
  sigma_sublineage ~ normal(0, 0.1);
  z_sub ~ normal(0, 1);

  // POM likelihood -- single effect per variant shared across cutpoints
  vector[N_train] mu_train = X_std_train * beta_variant_std
                           + X_sublineage_train * beta_sublineage;
  training_phenotype ~ ordered_logistic(mu_train, cutpoints);
}

generated quantities {
  // Unstandardized variant effects (per 0->1 allele) and cumulative odds ratios
  vector[V] beta_variant;
  vector[V] OR_variant_allele;
  for (v in 1:V) {
    if (sd_variant[v] > 0) {
      beta_variant[v] = beta_variant_std[v] / sd_variant[v];
      OR_variant_allele[v] = exp(beta_variant[v]);
    } else {
      beta_variant[v] = 0;
      OR_variant_allele[v] = 1;
    }
  }

  // Predicted phenotypes (draws)
  vector[N_test] mu_test = X_std_test * beta_variant_std
                         + X_sublineage_test * beta_sublineage;
  array[N_test] int predicted_phenotype;
  for (n in 1:N_test)
    predicted_phenotype[n] = ordered_logistic_rng(mu_test[n], cutpoints);

  // Expose priors for post-hoc verification
  real alpha_prior_mean_out = alpha_prior_mean;
  real p_baseline_emp_out   = p_baseline_emp;

  vector[V] beta_variant_std_prior;
  for (v in 1:V)
    beta_variant_std_prior[v] = z_variant[v] * tau * lambda_tilde_variant[v];


  // Heritability, on the liability (logistic-latent) scale

  // Some tips on interpreting these:

    // h2_narrow counts only measured variant effects
    // h2_broad counts variant effects + lineage/sublineage effects

    // Horseshoe shrinkage biases V_A downward, so h2_narrow is a lower bound

    // Residual latent variance is calculated as pi^2/3 as in bernoulli_logit

    // For info on the per-k PPO model heritability scores, see notes on that model

  real<lower=0> V_A;
  real<lower=0> V_pop;
  real<lower=0> V_E = pi()^2 / 3;
  real<lower=0, upper=1> h2_narrow;
  real<lower=0, upper=1> h2_broad;
  {
    vector[N_train] g_variant = X_std_train * beta_variant_std;
    vector[N_train] g_pop     = X_sublineage_train * beta_sublineage;
    V_A   = variance(g_variant);
    V_pop = variance(g_pop);
    real V_tot = V_A + V_pop + V_E;
    h2_narrow = V_A / V_tot;
    h2_broad  = (V_A + V_pop) / V_tot;
  }
}
