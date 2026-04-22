// Binary logistic prediction model -- training/test split counterpart of
// logistic_inference.stan.
//
// Fits on N_train samples with the same likelihood, priors, and generated
// summaries as the inference model, and emits predicted_phenotype on N_test
// held-out samples (Bernoulli draws) for gwas_workflow prediction-accuracy
// metrics. No PPC outputs.

data {
  int<lower=1> N_train;                                          // training samples
  int<lower=1> N_test;                                           // held-out samples
  int<lower=1> V;                                                // variants
  int<lower=1> L;                                                // lineage clusters
  int<lower=1> S;                                                // lineage subclusters

  array[N_train] int<lower=0, upper=1> training_phenotype;       // binary outcome (train)

  matrix[N_train, V] training_variants;
  matrix[N_test,  V] test_variants;
  matrix[N_train, S] training_sublineages;
  matrix[N_test,  S] test_sublineages;
  array[S] int<lower=1, upper=L> parent_lineage;
}

transformed data {
  matrix[N_train, S-1] X_sublineage_train =
    block(training_sublineages, 1, 2, N_train, S-1);
  matrix[N_test,  S-1] X_sublineage_test  =
    block(test_sublineages,     1, 2, N_test,  S-1);

  array[S-1] int parent_lineage_sub;
  for (k in 1:(S-1))
    parent_lineage_sub[k] = parent_lineage[k+1];

  // Standardize genotype using TRAINING column mean/SD; apply the same
  // transform to test.
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

  // Empirical baseline from training only: Laplace-smoothed susceptible
  // fraction of the reference sublineage. Matches logistic_inference.stan.
  int n_ref = 0;
  int n_ref_res = 0;
  for (n in 1:N_train) {
    if (training_sublineages[n, 1] > 0.5) {
      n_ref += 1;
      if (training_phenotype[n] == 1) n_ref_res += 1;
    }
  }
  real p_baseline_emp = (n_ref - n_ref_res + 0.5) / (n_ref + 1.0);
  p_baseline_emp = fmin(fmax(p_baseline_emp, 0.5), 0.995);

  real alpha_prior_mean = -logit(p_baseline_emp);
  real alpha_prior_sd   = 0.5;
}

parameters {
  real alpha;

  vector[L] beta_lineage;
  real<lower=0> sigma_sublineage;
  vector[S-1] z_sub;

  real<lower=0> tau;
  real<lower=0> c2;
  vector[V] z_variant;
  vector<lower=0>[V] lambda_variant;
}

transformed parameters {
  vector[S-1] beta_sublineage_raw =
    beta_lineage[parent_lineage_sub] + sigma_sublineage * z_sub;

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

  beta_lineage ~ normal(0, 0.1);
  sigma_sublineage ~ normal(0, 0.1);
  z_sub ~ normal(0, 1);

  vector[N_train] mu_train = X_std_train * beta_variant_std
                           + X_sublineage_train * beta_sublineage
                           + alpha;
  training_phenotype ~ bernoulli_logit(mu_train);
}

generated quantities {
  // Unstandardized variant effects (per 0->1 allele) and odds ratios
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

  // Held-out posterior predictive draws consumed by
  // gwas_workflow/R/prediction.R -> compute_binary_accuracy().
  vector[N_test] mu_test = X_std_test * beta_variant_std
                         + X_sublineage_test * beta_sublineage
                         + alpha;
  array[N_test] int predicted_phenotype;
  for (n in 1:N_test)
    predicted_phenotype[n] = bernoulli_logit_rng(mu_test[n]);

  // Expose priors for post-hoc verification
  real alpha_prior_mean_out = alpha_prior_mean;
  real p_baseline_emp_out   = p_baseline_emp;

  vector[V] beta_variant_std_prior;
  for (v in 1:V)
    beta_variant_std_prior[v] = z_variant[v] * tau * lambda_tilde_variant[v];

  // Liability-scale heritability, computed on the training partition.
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
