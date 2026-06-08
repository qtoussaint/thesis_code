// AUTO-GENERATED prediction (training/test split) twin of
// final_ordered_categorical_PPOM_free_cutpoints_wide_drift_ridge_estscale_ncp.stan
// Produced by parameter_sweeps_prediction/generate_prediction_stan.R; do not
// edit by hand. Fits on N_train and emits an N_test x K predicted_phenotype
// matrix (Rao-Blackwellized) for out-of-sample prediction scoring.

data {
  int<lower=1> N_train;
  int<lower=1> N_test;
  int<lower=1> V;
  int<lower=1> L;
  int<lower=1> S;
  int<lower=1> K;

  array[N_train] int<lower=1, upper=K> training_phenotype;
  vector<lower=1e-12>[K-1] mic_breakpoints;

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

  vector[K-1] mic_cutpoints = log2(mic_breakpoints);

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

  vector[K-1] cutpoint_prior_mean;
  for (k in 1:(K-1))
    cutpoint_prior_mean[k] = mic_cutpoints[k] - alpha_prior_mean;

  real cutpoint_prior_sd = 1.5;
}

parameters {
  ordered[K-1] cutpoints;

  vector[L] beta_lineage;
  real<lower=0> sigma_sublineage;
  vector[S-1] z_sub;

  matrix[V, K-1] z_variant;
  real<lower=0> sigma_beta_variant;
}

transformed parameters {
  matrix[V, K-1] beta_variant_std = z_variant * sigma_beta_variant;

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
}

model {
  to_vector(z_variant) ~ std_normal();
  sigma_beta_variant ~ student_t(3, 0, 2);

  cutpoints ~ normal(cutpoint_prior_mean, cutpoint_prior_sd);

  beta_lineage ~ normal(0, 0.1);
  sigma_sublineage ~ normal(0, 0.1);
  z_sub ~ normal(0, 1);

  {
    vector[N_train] sub_eta = X_sublineage_train * beta_sublineage;
    for (n in 1:N_train) {
      for (k in 1:(K-1)) {
        real mu_nk = dot_product(X_std_train[n], beta_variant_std[, k]) + sub_eta[n];
        if (training_phenotype[n] <= k) {
          target += bernoulli_logit_lpmf(1 | cutpoints[k] - mu_nk);
        } else {
          target += bernoulli_logit_lpmf(0 | cutpoints[k] - mu_nk);
        }
      }
    }
  }
}

generated quantities {
  matrix[V, K-1] beta_variant;
  matrix[V, K-1] OR_variant_allele;
  for (k in 1:(K-1)) {
    for (v in 1:V) {
      if (sd_variant[v] > 0) {
        beta_variant[v, k] = beta_variant_std[v, k] / sd_variant[v];
        OR_variant_allele[v, k] = exp(beta_variant[v, k]);
      } else {
        beta_variant[v, k] = 0;
        OR_variant_allele[v, k] = 1;
      }
    }
  }

  vector[K-1] cutpoint_drift;
  for (k in 1:(K-1))
    cutpoint_drift[k] = cutpoints[k] - cutpoint_prior_mean[k];

  matrix[N_test, K] predicted_phenotype;
  {
    vector[N_test] sub_eta_test = X_sublineage_test * beta_sublineage;
    for (n in 1:N_test) {
      vector[K-1] cdf_n;
      for (k in 1:(K-1)) {
        real mu_nk = dot_product(X_std_test[n], beta_variant_std[, k]) + sub_eta_test[n];
        cdf_n[k] = inv_logit(cutpoints[k] - mu_nk);
      }
      vector[K] probs;
      probs[1] = cdf_n[1];
      for (k in 2:(K-1)) probs[k] = cdf_n[k] - cdf_n[k-1];
      probs[K] = 1 - cdf_n[K-1];
      for (k in 1:K) if (probs[k] < 1e-12) probs[k] = 1e-12;
      probs /= sum(probs);
      for (k in 1:K) predicted_phenotype[n, k] = probs[k];
    }
  }

  real alpha_prior_mean_out = alpha_prior_mean;
  real p_baseline_emp_out   = p_baseline_emp;

  vector<lower=0>[K-1] V_A_k;
  real<lower=0> V_pop;
  real<lower=0> V_E = pi()^2 / 3;
  vector<lower=0, upper=1>[K-1] h2_narrow_k;
  vector<lower=0, upper=1>[K-1] h2_broad_k;
  real<lower=0, upper=1> h2_narrow_median_k;
  real<lower=0, upper=1> h2_broad_median_k;
  {
    vector[N_train] g_pop = X_sublineage_train * beta_sublineage;
    V_pop = variance(g_pop);
    for (k in 1:(K-1)) {
      vector[N_train] g_variant_k = X_std_train * beta_variant_std[, k];
      V_A_k[k] = variance(g_variant_k);
      real V_tot_k = V_A_k[k] + V_pop + V_E;
      h2_narrow_k[k] = V_A_k[k] / V_tot_k;
      h2_broad_k[k]  = (V_A_k[k] + V_pop) / V_tot_k;
    }
    h2_narrow_median_k = quantile(h2_narrow_k, 0.5);
    h2_broad_median_k  = quantile(h2_broad_k,  0.5);
  }
}
