#!/usr/bin/env Rscript
# Posterior predictive check for a single (dataset, model) POM run.
#
# Usage (activate gwas_pipeline conda env first):
#   Rscript ppc_compare_alpha_models.R --dataset 02_spn_penicillin_MIC --model tight_alpha
#
# All PPC quantities (y_rep, cat_freq_rep) are read directly from the Stan fit's
# generated quantities block — no R-side reconstruction. The three POM Stan files
# all emit y_rep and cat_freq_rep. If your fit does not contain them, re-run the
# model after pulling the updated Stan file.
#
# Known models (extend MODEL_REGISTRY below to add more):
#   fixed_alpha  - final_ordered_categorical_POM             (alpha hard-coded)
#   free_alpha   - final_ordered_categorical_POM_free_alpha  (alpha ~ N(c1-2, 1.5))
#   tight_alpha  - final_ordered_categorical_POM_tight_alpha (alpha ~ N(c1-logit(p_emp), 0.5))
#
# Outputs (written to <results_root>/<model_subdir>/<dataset>/ppc/):
#   ppc_category_frequencies.png   observed MIC histogram vs cat_freq_rep posterior
#   alpha_posterior.png            density of alpha (or a note if alpha is fixed)
#   effect_size_summary.png        |beta_variant| histogram + top-20 ranked variants
#   metrics_summary.csv            one-row metrics for this (dataset, model) run
#   cat_freq_posterior.csv         per-draw cat_freq_rep (for downstream aggregation)

suppressPackageStartupMessages({
  library(optparse)
  library(posterior)
  library(jsonlite)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

# ---- MODEL REGISTRY --------------------------------------------------------
MODEL_REGISTRY <- list(
  fixed_alpha = list(
    label      = "fixed alpha (logit 0.99)",
    subdir     = "final_ordered_categorical_POM",
    alpha_type = "fixed"
  ),
  free_alpha = list(
    label      = "free alpha (loose)",
    subdir     = "final_ordered_categorical_POM_free_alpha",
    alpha_type = "free"
  ),
  tight_alpha = list(
    label      = "tight data-informed alpha",
    subdir     = "final_ordered_categorical_POM_tight_alpha",
    alpha_type = "free"
  )
)

# ---- CLI ARGS --------------------------------------------------------------
option_list <- list(
  make_option("--dataset", type = "character", default = NULL,
              help = "Dataset nickname, e.g. 02_spn_penicillin_MIC"),
  make_option("--model", type = "character", default = NULL,
              help = paste0("Model key: ", paste(names(MODEL_REGISTRY), collapse = ", "))),
  make_option("--results_root", type = "character",
              default = "/nfs/research/jlees/jacqueline/thesis_results/compare_ordinal_models",
              help = "Parent dir holding per-model result subdirs"),
  make_option("--data_root", type = "character",
              default = "/nfs/research/jlees/jacqueline/thesis_results/gwas_datasets/inference",
              help = "Parent dir holding per-dataset JSON input (for observed phenotype)"),
  make_option("--n_draws", type = "integer", default = 1000L,
              help = "Max posterior draws to load for plots/metrics [1000]"),
  make_option("--seed", type = "integer", default = 42L,
              help = "RNG seed for draw subsampling [42]")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$dataset) || is.null(opt$model))
  stop("Both --dataset and --model are required")
if (!(opt$model %in% names(MODEL_REGISTRY)))
  stop("Unknown --model '", opt$model, "'; choose from: ",
       paste(names(MODEL_REGISTRY), collapse = ", "))

MODEL  <- MODEL_REGISTRY[[opt$model]]
DSET   <- opt$dataset
RUNDIR <- file.path(opt$results_root, MODEL$subdir, DSET)
JSON   <- file.path(opt$data_root, DSET, paste0(DSET, ".json"))
OUTDIR <- file.path(RUNDIR, "ppc")

if (!dir.exists(RUNDIR)) stop("Run directory not found: ", RUNDIR)
if (!file.exists(JSON))  stop("Data JSON not found: ", JSON)
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
set.seed(opt$seed)

message("dataset : ", DSET)
message("model   : ", opt$model, "  (", MODEL$label, ")")
message("run dir : ", RUNDIR)
message("out dir : ", OUTDIR)

# ---- OBSERVED MIC HISTOGRAM (from input JSON) ------------------------------
dat <- fromJSON(JSON, simplifyVector = TRUE)
K <- as.integer(dat$K)
N <- as.integer(dat$N)
phenotype <- as.integer(dat$phenotype)
cutpoints <- log2(as.numeric(dat$mic_breakpoints))
obs_freq  <- tabulate(phenotype, nbins = K) / N

# ---- LOAD FIT DRAWS --------------------------------------------------------
rds_candidates <- list.files(file.path(RUNDIR, "fitted_model"),
                             pattern = "\\.RDS$", ignore.case = TRUE, full.names = TRUE)
rds_candidates <- rds_candidates[!grepl("depruned", rds_candidates)]
if (length(rds_candidates) == 0)
  stop("No fit RDS in ", file.path(RUNDIR, "fitted_model"))
rds_path <- rds_candidates[which.max(file.mtime(rds_candidates))]
message("loading fit: ", basename(rds_path))

fit_obj  <- readRDS(rds_path)
draws_df <- as_draws_df(tryCatch(fit_obj$draws(), error = function(e) fit_obj))

D_all   <- nrow(draws_df)
if (opt$n_draws < D_all) {
  use_idx <- sort(sample.int(D_all, opt$n_draws))
  draws_df <- draws_df[use_idx, , drop = FALSE]
}
message("using ", nrow(draws_df), " of ", D_all, " posterior draws")

# ---- READ cat_freq_rep, beta_variant, alpha DIRECTLY -----------------------
subset_param_matrix <- function(df, base) {
  vars <- grep(paste0("^", base, "\\["), names(df), value = TRUE)
  if (length(vars) == 0) return(NULL)
  ord_idx <- as.integer(sub(".*\\[(\\d+)\\]$", "\\1", vars))
  as.matrix(df[, vars[order(ord_idx)]])
}

cat_freq_rep <- subset_param_matrix(draws_df, "cat_freq_rep")
if (is.null(cat_freq_rep))
  stop("cat_freq_rep not found in draws — the fit was produced by an older Stan file.\n",
       "       Re-run the model with the updated *.stan (all three emit cat_freq_rep/y_rep).")
if (ncol(cat_freq_rep) != K)
  stop("cat_freq_rep has ", ncol(cat_freq_rep), " cols; expected K = ", K)

beta_variant <- subset_param_matrix(draws_df, "beta_variant")
if (is.null(beta_variant))
  stop("beta_variant not found in draws")
V <- ncol(beta_variant)
beta_mean <- colMeans(beta_variant)

if (MODEL$alpha_type == "free") {
  if (!("alpha" %in% names(draws_df)))
    stop("`alpha` not found in draws for a free-alpha model")
  alpha_vec <- as.numeric(draws_df$alpha)
} else {
  alpha_vec <- NULL  # intercept is alpha_mean in transformed data, not sampled
}

# ---- PLOT 1: CATEGORY FREQUENCY PPC ----------------------------------------
cat_df <- data.frame(
  category = factor(seq_len(K)),
  mean = colMeans(cat_freq_rep),
  lo   = apply(cat_freq_rep, 2, quantile, 0.05),
  hi   = apply(cat_freq_rep, 2, quantile, 0.95)
)
obs_df <- data.frame(category = factor(seq_len(K)), freq = obs_freq)

p_cat <- ggplot() +
  geom_col(data = obs_df, aes(x = category, y = freq),
           fill = "grey85", colour = "grey40", width = 0.85) +
  geom_pointrange(data = cat_df,
                  aes(x = category, y = mean, ymin = lo, ymax = hi),
                  colour = "steelblue", size = 0.5) +
  labs(title = sprintf("PPC category frequencies — %s (%s)", DSET, MODEL$label),
       subtitle = "grey bars = observed; points = cat_freq_rep posterior mean, bars = 90% CI",
       x = "MIC category (1 = lowest)", y = "frequency") +
  theme_bw(base_size = 11)
ggsave(file.path(OUTDIR, "ppc_category_frequencies.png"),
       p_cat, width = 7, height = 4.5, dpi = 150)

# ---- PLOT 2: ALPHA POSTERIOR -----------------------------------------------
fixed_alpha_ref <- cutpoints[1] - qlogis(0.99)
if (MODEL$alpha_type == "free") {
  p_alpha <- ggplot(data.frame(alpha = alpha_vec), aes(x = alpha)) +
    geom_density(fill = "steelblue", alpha = 0.4) +
    geom_vline(xintercept = fixed_alpha_ref, linetype = "dashed", colour = "black") +
    annotate("text", x = fixed_alpha_ref, y = 0,
             label = " fixed-alpha reference", hjust = 0, vjust = -0.5, size = 3) +
    geom_vline(xintercept = cutpoints[1], linetype = "dotted", colour = "grey40") +
    annotate("text", x = cutpoints[1], y = 0, label = " cutpoints[1]",
             hjust = 0, vjust = -0.5, size = 3, colour = "grey40") +
    labs(title = sprintf("Posterior of alpha — %s (%s)", DSET, MODEL$label),
         x = "alpha (log2-MIC latent scale)", y = "density") +
    theme_bw(base_size = 11)
} else {
  p_alpha <- ggplot() +
    annotate("text", x = 0, y = 0,
             label = sprintf("alpha is FIXED at cutpoints[1] - logit(0.99) = %.3f",
                             fixed_alpha_ref)) +
    theme_void()
}
ggsave(file.path(OUTDIR, "alpha_posterior.png"),
       p_alpha, width = 7, height = 4.5, dpi = 150)

# ---- PLOT 3: EFFECT SIZE SUMMARY -------------------------------------------
abs_beta  <- abs(beta_mean)
order_idx <- order(abs_beta, decreasing = TRUE)
top_n     <- min(20, V)
top_df <- data.frame(
  variant_idx = order_idx[seq_len(top_n)],
  beta_mean   = beta_mean[order_idx[seq_len(top_n)]]
)

p_hist <- ggplot(data.frame(b = beta_mean), aes(x = b)) +
  geom_histogram(bins = 60, fill = "steelblue", colour = "grey30") +
  labs(title = "Posterior-mean beta_variant (allele scale)",
       x = "beta_variant", y = "count") +
  theme_bw(base_size = 10)

p_top <- ggplot(top_df, aes(x = beta_mean, y = reorder(factor(variant_idx), beta_mean))) +
  geom_col(fill = "steelblue") +
  labs(title = sprintf("Top %d variants by |beta|", top_n),
       x = "posterior mean beta_variant", y = "variant index") +
  theme_bw(base_size = 10)

if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  ggsave(file.path(OUTDIR, "effect_size_summary.png"),
         p_hist + p_top, width = 11, height = 5, dpi = 150)
} else {
  ggsave(file.path(OUTDIR, "effect_size_histogram.png"), p_hist,
         width = 6, height = 4.5, dpi = 150)
  ggsave(file.path(OUTDIR, "effect_size_top20.png"), p_top,
         width = 6, height = 5, dpi = 150)
}

# ---- METRICS ---------------------------------------------------------------
cf_mean <- colMeans(cat_freq_rep)
diff    <- cf_mean - obs_freq
chi_sq  <- sum(ifelse(obs_freq > 0, (diff^2) / obs_freq, 0)) * N

metrics <- data.frame(
  dataset              = DSET,
  model                = opt$model,
  model_label          = MODEL$label,
  n_samples            = N,
  n_variants           = V,
  rmse_cat_freq        = sqrt(mean(diff^2)),
  max_abs_dev_cat_freq = max(abs(diff)),
  chi_square_stat      = chi_sq,
  alpha_post_mean      = if (MODEL$alpha_type == "free") mean(alpha_vec) else fixed_alpha_ref,
  alpha_post_sd        = if (MODEL$alpha_type == "free") sd(alpha_vec)   else 0,
  n_nonzero_beta_mean  = sum(abs(beta_mean) > 1e-4),
  max_abs_beta         = max(abs(beta_mean)),
  stringsAsFactors     = FALSE
)
write.csv(metrics, file.path(OUTDIR, "metrics_summary.csv"), row.names = FALSE)

# Per-draw cat_freq_rep for downstream aggregation across (dataset, model) runs
cat_out <- as.data.frame(cat_freq_rep)
names(cat_out) <- paste0("cat_", seq_len(K))
cat_out$draw_idx <- seq_len(nrow(cat_out))
write.csv(cat_out, file.path(OUTDIR, "cat_freq_posterior.csv"), row.names = FALSE)

# ---- CONSOLE SUMMARY -------------------------------------------------------
message("\n=== metrics ===")
print(metrics, row.names = FALSE)
message("\nobserved cat freq : ", paste(round(obs_freq, 4), collapse = "  "))
message("expected cat freq : ", paste(round(cf_mean,  4), collapse = "  "))
message("\noutputs written to: ", OUTDIR)
