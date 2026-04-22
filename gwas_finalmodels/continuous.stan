// Continuous (log2-MIC) association model

// Standardized genotype with regularized horseshoe on variant effects
// (Piironen & Vehtari 2017, eq. 3.12), lineage/sublineage hierarchy with
// treatment-contrast encoding (reference sublineage = column 1) and
// within-lineage sum-to-zero centering, and a tight, data-informed alpha prior
// anchored on the reference sublineage.
//
// Alpha prior is anchored on the empirical mean phenotype of the reference
// sublineage, which gwas_datasets/utils.R:select_reference_sublineage picks
// upstream as the minimum-mean-phenotype sublineage (lowest log2 MIC). This
// stops alpha from absorbing variant-level signal.

data {
  int<lower=1> N;                                  // samples
  int<lower=1> V;                                  // variants
  int<lower=1> L;                                  // lineage clusters
  int<lower=1> S;                                  // lineage subclusters

  vector[N] phenotype;                             // log2 MIC

  matrix[N, V] variant_matrix;                     // genotype
  matrix[N, S] sublineage_matrix;                  // full one-hot, reference in column 1
  array[S] int<lower=1, upper=L> parent_lineage;   // parent lineage of each subcluster

  int<lower=0, upper=N> N_ppc;                     // PPC subset size (20% of N, capped at 500)
  array[N_ppc] int<lower=1, upper=N> ppc_idx;      // deterministic indices into 1:N for PPC
}

transformed data {
  // Treatment-contrast encoding: drop the reference sublineage (column 1)
  matrix[N, S-1] X_sublineage = block(sublineage_matrix, 1, 2, N, S-1);

  array[S-1] int parent_lineage_sub;
  for (k in 1:(S-1))
    parent_lineage_sub[k] = parent_lineage[k+1];

  // Centre and scale the genotype matrix (sd=0 columns pinned to 0)
  matrix[N, V] X_std;
  vector[V] mean_variant;
  vector[V] sd_variant;
  for (v in 1:V) {
    mean_variant[v] = mean(variant_matrix[, v]);
    vector[N] centred = variant_matrix[, v] - mean_variant[v];
    sd_variant[v] = sd(centred);
    for (n in 1:N)
      X_std[n, v] = (sd_variant[v] > 0) ? centred[n] / sd_variant[v] : 0;
  }

  // Empirical baseline: mean log2 MIC of the reference sublineage. Upstream
  // (utils.R:select_reference_sublineage) guarantees n_ref >= 1, so no
  // small-n fallback is needed. alpha_prior_sd = 0.5 ≈ one log2-dilution at
  // 95%, matches MIC rounding.
  int n_ref = 0;
  real ref_sum = 0;
  for (n in 1:N) {
    if (sublineage_matrix[n, 1] > 0.5) {
      n_ref += 1;
      ref_sum += phenotype[n];
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

  // Tight and data-informed intercept prior to prevent alpha from absorbing
  // variant-level signal while still accounting for reference-cluster noise
  alpha ~ normal(alpha_prior_mean, alpha_prior_sd);
  sigma ~ normal(0, 1);  // half-normal, one log2-dilution is the MIC-rounding scale

  beta_lineage ~ normal(0, 0.1);
  sigma_sublineage ~ normal(0, 0.1);
  z_sub ~ normal(0, 1);

  vector[N] mu = X_std * beta_variant_std
               + X_sublineage * beta_sublineage
               + alpha;
  phenotype ~ normal(mu, sigma);
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

  // Posterior predictive check on a deterministic 20%-of-N subset (capped at 500,
  // chosen upstream in build_stan_inference). y_true_ppc travels alongside so
  // downstream metrics pair predicted and true values without reloading phenotype.
  vector[N_ppc] y_rep_ppc;
  vector[N_ppc] y_true_ppc;
  {
    for (i in 1:N_ppc) {
      int n = ppc_idx[i];
      real mu_n = X_std[n] * beta_variant_std
                + X_sublineage[n] * beta_sublineage
                + alpha;
      y_rep_ppc[i]  = normal_rng(mu_n, sigma);
      y_true_ppc[i] = phenotype[n];
    }
  }

  // Expose priors for post-hoc verification
  real alpha_prior_mean_out = alpha_prior_mean;

  vector[V] beta_variant_std_prior;
  for (v in 1:V)
    beta_variant_std_prior[v] = z_variant[v] * tau * lambda_tilde_variant[v];
}
