// Continuous (log2-MIC) prediction model -- training/test split counterpart of
// continuous_inference.stan.
//
// Fits on N_train samples with the same likelihood, priors, and generated
// summaries as the inference model, and emits predicted_phenotype on N_test
// held-out samples for gwas_workflow prediction-accuracy metrics. No PPC
// outputs: y_rep / y_true are absent by design (downstream consumes only
// predicted_phenotype).

data {
  int<lower=1> N_train;                                          // training samples
  int<lower=1> N_test;                                           // held-out samples
  int<lower=1> V;                                                // variants
  int<lower=1> L;                                                // lineage clusters
  int<lower=1> S;                                                // lineage subclusters

  vector[N_train] training_phenotype;                            // log2 MIC (training only)

  matrix[N_train, V] training_variants;                          // genotype (train)
  matrix[N_test,  V] test_variants;                              // genotype (test)
  matrix[N_train, S] training_sublineages;                       // full one-hot, ref in col 1
  matrix[N_test,  S] test_sublineages;                           // full one-hot, ref in col 1
  array[S] int<lower=1, upper=L> parent_lineage;                 // parent lineage of each subcluster
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
  // transform to test so coefficients on X_std_train transfer to X_std_test.
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

  // Empirical baseline from training only: mean log2 MIC of the reference
  // sublineage. Matches continuous_inference.stan exactly.
  int n_ref = 0;
  real ref_sum = 0;
  for (n in 1:N_train) {
    if (training_sublineages[n, 1] > 0.5) {
      n_ref += 1;
      ref_sum += training_phenotype[n];
    }
  }
  real alpha_prior_mean = ref_sum / n_ref;
  real alpha_prior_sd = 0.5;
}

parameters {
  real alpha;
  real<lower=0> sigma;  // residual std dev on log2-MIC scale

  vector[L] beta_lineage;
  real<lower=0> sigma_sublineage;
  vector[S-1] z_sub;

  // Regularized horseshoe on variant effects
  real<lower=0> tau;
  real<lower=0> c2;
  vector[V] z_variant;
  vector<lower=0>[V] lambda_variant;
}

transformed parameters {
  vector[S-1] beta_sublineage_raw =
    beta_lineage[parent_lineage_sub] + sigma_sublineage * z_sub;

  // Within-lineage sum-to-zero centering (global centering conflicts
  // with treatment-contrast encoding and makes alpha uninterpretable)
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

  // Horseshoe hyperparameters (Piironen & Vehtari 2017, eq. 3.12)
  real slab_scale = 200;
  real nu         = 4;
  real tau_0      = 0.1;

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

  alpha ~ normal(alpha_prior_mean, alpha_prior_sd);
  sigma ~ normal(0, 1);  // half-normal, one log2-dilution is the MIC-rounding scale

  beta_lineage ~ normal(0, 0.1);
  sigma_sublineage ~ normal(0, 0.1);
  z_sub ~ normal(0, 1);

  vector[N_train] mu_train = X_std_train * beta_variant_std
                           + X_sublineage_train * beta_sublineage
                           + alpha;
  training_phenotype ~ normal(mu_train, sigma);
}

generated quantities {
  // Unstandardized variant effects (log2-MIC change per 0->1 allele)
  vector[V] beta_variant;
  for (v in 1:V) {
    if (sd_variant[v] > 0)
      beta_variant[v] = beta_variant_std[v] / sd_variant[v];
    else
      beta_variant[v] = 0;
  }

  // Held-out posterior predictive draws consumed by
  // gwas_workflow/R/prediction.R -> compute_continuous_accuracy().
  vector[N_test] mu_test = X_std_test * beta_variant_std
                         + X_sublineage_test * beta_sublineage
                         + alpha;
  array[N_test] real predicted_phenotype;
  for (n in 1:N_test)
    predicted_phenotype[n] = normal_rng(mu_test[n], sigma);

  // Expose priors for post-hoc verification
  real alpha_prior_mean_out = alpha_prior_mean;

  vector[V] beta_variant_std_prior;
  for (v in 1:V)
    beta_variant_std_prior[v] = z_variant[v] * tau * lambda_tilde_variant[v];

  // Heritability (observed-scale, Gaussian residual), computed on the training
  // partition to match the inference model's V_A / V_pop definitions.
  real<lower=0> V_A;
  real<lower=0> V_pop;
  real<lower=0> V_E = square(sigma);
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
