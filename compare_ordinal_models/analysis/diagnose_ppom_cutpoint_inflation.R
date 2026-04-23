#!/usr/bin/env Rscript
# Per-cutpoint diagnostic for PPOM fits.
#
# Diagnoses the "beta_variant[v, k] inflates with k" failure mode by producing:
#   1. Per-cutpoint event counts P(phen > k) -- mechanism (1) rare-event check
#   2. Posterior |beta| median + 95% width by k           -- identifiability check
#   3. Top-N variant beta[v, k] trajectories across k     -- inflation fingerprint
#   4. Reference-sublineage empirical CDF vs cutpoints    -- mechanism (3) sanity
#   5. Optional: sigma_latent / sigma_dev posterior if the fit has them
#   6. y_rep_ppc category frequencies (Stan-side PPC) vs observed frequencies
#
# Usage (activate gwas_pipeline conda env first):
#   Rscript diagnose_ppom_cutpoint_inflation.R \
#     --model_subdir final_ordered_categorical_PPOM_tight_alpha_tau1 \
#     --dataset      02_spn_penicillin_MIC
#
# Writes everything to
#   <results_root>/<model_subdir>/<dataset>/ppc/diagnose/
#
# The script is model-agnostic: it probes the draws for optional parameters
# (sigma_latent, sigma_dev, cutpoints, cutpoint_drift) and adapts accordingly.
# It is intended to be run identically against every PPOM variant
# (baseline tight_alpha_tau1, tight_slab, poolk, latent_scale, free_cutpoints).
#
# Handling of free_cutpoints: that variant promotes `cutpoints` from
# transformed data to a parameter and drops `alpha`. When the script detects
# `cutpoints[k]` in the draws it uses per-draw cutpoints (no alpha) to build
# the reference-CDF diagnostic, and adds diagnostic 5b (cutpoint_drift) to
# show how far each cutpoint has moved from its MIC-grid anchor.

suppressPackageStartupMessages({
  library(optparse)
  library(posterior)
  library(jsonlite)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

# ---- CLI ARGS --------------------------------------------------------------
option_list <- list(
  make_option("--model_subdir", type = "character", default = NULL,
              help = "Model subdir under results_root, e.g. final_ordered_categorical_PPOM_tight_alpha_tau1"),
  make_option("--dataset", type = "character", default = NULL,
              help = "Dataset nickname, e.g. 02_spn_penicillin_MIC"),
  make_option("--results_root", type = "character",
              default = "/nfs/research/jlees/jacqueline/thesis_results/compare_ordinal_models",
              help = "Parent dir holding per-model result subdirs"),
  make_option("--data_root", type = "character",
              default = "/nfs/research/jlees/jacqueline/thesis_results/gwas_datasets/inference",
              help = "Parent dir holding per-dataset JSON input"),
  make_option("--n_draws", type = "integer", default = 1000L,
              help = "Max posterior draws to load [1000]"),
  make_option("--top_n", type = "integer", default = 12L,
              help = "Number of top variants (by max |beta| across k) to plot [12]"),
  make_option("--seed", type = "integer", default = 42L,
              help = "RNG seed for draw subsampling [42]")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$model_subdir) || is.null(opt$dataset))
  stop("Both --model_subdir and --dataset are required")

RUNDIR <- file.path(opt$results_root, opt$model_subdir, opt$dataset)
JSON   <- file.path(opt$data_root,    opt$dataset,       paste0(opt$dataset, ".json"))
OUTDIR <- file.path(RUNDIR, "ppc", "diagnose")
if (!dir.exists(RUNDIR)) stop("Run directory not found: ", RUNDIR)
if (!file.exists(JSON))  stop("Data JSON not found: ", JSON)
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
set.seed(opt$seed)

message("model   : ", opt$model_subdir)
message("dataset : ", opt$dataset)
message("out dir : ", OUTDIR)

# ---- LOAD DATA JSON -------------------------------------------------------
dat <- fromJSON(JSON, simplifyVector = TRUE)
K <- as.integer(dat$K)
N <- as.integer(dat$N)
phenotype  <- as.integer(dat$phenotype)
cutpoints  <- log2(as.numeric(dat$mic_breakpoints))
sub_matrix <- as.matrix(dat$sublineage_matrix)    # N x S, reference = col 1
ref_mask   <- sub_matrix[, 1] > 0.5
n_ref      <- sum(ref_mask)

# ---- LOAD FIT DRAWS -------------------------------------------------------
rds_candidates <- list.files(file.path(RUNDIR, "fitted_model"),
                             pattern = "\\.RDS$", ignore.case = TRUE, full.names = TRUE)
rds_candidates <- rds_candidates[!grepl("depruned", rds_candidates)]
if (length(rds_candidates) == 0)
  stop("No fit RDS in ", file.path(RUNDIR, "fitted_model"))
rds_path <- rds_candidates[which.max(file.mtime(rds_candidates))]
message("loading fit: ", basename(rds_path))

fit_obj  <- readRDS(rds_path)
draws_df <- as_draws_df(tryCatch(fit_obj$draws(), error = function(e) fit_obj))

D_all <- nrow(draws_df)
if (opt$n_draws < D_all) {
  use_idx <- sort(sample.int(D_all, opt$n_draws))
  draws_df <- draws_df[use_idx, , drop = FALSE]
}
message("using ", nrow(draws_df), " of ", D_all, " posterior draws")

# ---- PARAMETER EXTRACTORS -------------------------------------------------
# beta_variant is stored with names beta_variant[v,k]. We want a
# (draws, V, K-1) array. Column order in as_draws_df is unsorted so we
# rebuild with explicit (v, k) parsing.
pull_matrix_param <- function(df, base, expected_nrow, expected_ncol) {
  nm <- grep(paste0("^", base, "\\["), names(df), value = TRUE)
  if (length(nm) == 0) return(NULL)
  vk <- do.call(rbind, strsplit(sub(".*\\[(.*)\\]$", "\\1", nm), ","))
  if (ncol(vk) < 2) return(NULL)  # 1-D array; caller will handle separately
  v_idx <- as.integer(vk[, 1]); k_idx <- as.integer(vk[, 2])
  V <- max(v_idx); K_1 <- max(k_idx)
  if (!is.null(expected_nrow) && V   != expected_nrow) stop(base, ": V mismatch")
  if (!is.null(expected_ncol) && K_1 != expected_ncol) stop(base, ": K-1 mismatch")
  arr <- array(NA_real_, dim = c(nrow(df), V, K_1))
  for (i in seq_along(nm))
    arr[, v_idx[i], k_idx[i]] <- df[[nm[i]]]
  arr
}

beta_variant_arr <- pull_matrix_param(draws_df, "beta_variant", NULL, K - 1L)
if (is.null(beta_variant_arr))
  stop("beta_variant not found in draws -- did the fit complete its generated quantities?")
V <- dim(beta_variant_arr)[2]

# Posterior summaries: median + sd + 2.5/97.5 per (v, k)
beta_med <- apply(beta_variant_arr, c(2, 3), median)
beta_sd  <- apply(beta_variant_arr, c(2, 3), sd)
beta_lo  <- apply(beta_variant_arr, c(2, 3), quantile, 0.025)
beta_hi  <- apply(beta_variant_arr, c(2, 3), quantile, 0.975)

# y_rep_ppc (optional; present for all 2026-era PPOM variants)
y_rep_arr <- pull_matrix_param(draws_df, "y_rep_ppc", NULL, NULL)
# (It is emitted as an array not a matrix, but pull_matrix_param fails gracefully
# if names don't match [v,k]; we fall back to 1-D extraction for y_rep_ppc.)
if (is.null(y_rep_arr)) {
  y_rep_cols <- grep("^y_rep_ppc\\[", names(draws_df), value = TRUE)
  if (length(y_rep_cols) > 0) {
    y_rep_mat <- as.matrix(draws_df[, y_rep_cols])
    # columns are y_rep_ppc[1], y_rep_ppc[2], ... (1-D index)
    ord <- order(as.integer(sub(".*\\[(\\d+)\\]$", "\\1", y_rep_cols)))
    y_rep_mat <- y_rep_mat[, ord]
  } else {
    y_rep_mat <- NULL
  }
} else {
  y_rep_mat <- NULL
}

# ---- DIAGNOSTIC 1: PER-K EVENT COUNTS + CUMULATIVE FRACTION ---------------
n_above_k    <- integer(K - 1L)
frac_above_k <- numeric(K - 1L)
for (k in seq_len(K - 1L)) {
  n_above_k[k]    <- sum(phenotype > k)
  frac_above_k[k] <- n_above_k[k] / N
}
n_above_k_ref <- integer(K - 1L)
frac_above_k_ref <- numeric(K - 1L)
for (k in seq_len(K - 1L)) {
  n_above_k_ref[k]   <- sum(phenotype[ref_mask] > k)
  frac_above_k_ref[k] <- n_above_k_ref[k] / max(n_ref, 1L)
}

events_df <- data.frame(
  k            = seq_len(K - 1L),
  cutpoint     = cutpoints,
  n_above_full = n_above_k,
  frac_full    = frac_above_k,
  n_above_ref  = n_above_k_ref,
  frac_ref     = frac_above_k_ref
)
write.csv(events_df, file.path(OUTDIR, "per_k_event_counts.csv"), row.names = FALSE)

p_events <- ggplot(events_df, aes(x = factor(k), y = frac_full)) +
  geom_col(fill = "steelblue") +
  geom_text(aes(label = n_above_full), vjust = -0.4, size = 3) +
  labs(title = sprintf("Event rate P(phen > k) per cutpoint -- %s / %s",
                       opt$dataset, opt$model_subdir),
       subtitle = "numbers above bars = count; tiny counts at late k => mechanism (1) plausible",
       x = "cutpoint index k", y = "fraction of samples > k") +
  theme_bw(base_size = 10)
ggsave(file.path(OUTDIR, "01_per_k_event_rate.png"),
       p_events, width = 7, height = 4, dpi = 150)

# ---- DIAGNOSTIC 2: |BETA| DISTRIBUTION BY K -------------------------------
abs_beta_med <- abs(beta_med)
beta_by_k <- do.call(rbind, lapply(seq_len(K - 1L), function(k)
  data.frame(k = k, variant_idx = seq_len(V),
             beta_median = beta_med[, k], beta_abs = abs_beta_med[, k],
             beta_sd = beta_sd[, k])
))

summary_by_k <- beta_by_k %>% group_by(k) %>%
  summarise(median_abs_beta = median(beta_abs),
            q95_abs_beta    = quantile(beta_abs, 0.95),
            max_abs_beta    = max(beta_abs),
            median_post_sd  = median(beta_sd),
            q95_post_sd     = quantile(beta_sd, 0.95),
            .groups = "drop")
write.csv(summary_by_k, file.path(OUTDIR, "summary_by_k.csv"), row.names = FALSE)

p_abs <- ggplot(summary_by_k, aes(x = factor(k))) +
  geom_col(aes(y = q95_abs_beta), fill = "steelblue", alpha = 0.4) +
  geom_col(aes(y = median_abs_beta), fill = "steelblue") +
  labs(title = sprintf("|beta_variant| shrinkage by cutpoint -- %s / %s",
                       opt$dataset, opt$model_subdir),
       subtitle = "solid = median across variants; light = 95th percentile. Monotone growth => inflation",
       x = "cutpoint index k", y = "|beta_variant| (allele scale)") +
  theme_bw(base_size = 10)
ggsave(file.path(OUTDIR, "02_abs_beta_by_k.png"),
       p_abs, width = 7, height = 4, dpi = 150)

p_sd <- ggplot(summary_by_k, aes(x = factor(k), y = median_post_sd)) +
  geom_col(fill = "firebrick") +
  geom_errorbar(aes(ymin = median_post_sd, ymax = q95_post_sd),
                width = 0.25, colour = "firebrick4") +
  labs(title = sprintf("Posterior SD of beta_variant by cutpoint -- %s / %s",
                       opt$dataset, opt$model_subdir),
       subtitle = "Rising SD with k => weakly-identified late-k effects (mechanism 1 + 2)",
       x = "cutpoint index k", y = "posterior sd(beta_variant)") +
  theme_bw(base_size = 10)
ggsave(file.path(OUTDIR, "03_post_sd_by_k.png"),
       p_sd, width = 7, height = 4, dpi = 150)

# ---- DIAGNOSTIC 3: TOP-N VARIANT TRAJECTORIES -----------------------------
max_abs_by_v <- apply(abs_beta_med, 1, max)
top_v <- order(max_abs_by_v, decreasing = TRUE)[seq_len(min(opt$top_n, V))]

traj_df <- do.call(rbind, lapply(top_v, function(v)
  data.frame(variant_idx = v,
             k           = seq_len(K - 1L),
             beta_median = beta_med[v, ],
             beta_lo     = beta_lo[v, ],
             beta_hi     = beta_hi[v, ])
))
traj_df$variant_label <- factor(traj_df$variant_idx,
  levels = top_v,
  labels = sprintf("v%d (max|beta|=%.2f)", top_v, max_abs_by_v[top_v]))
write.csv(traj_df, file.path(OUTDIR, "top_variant_trajectories.csv"), row.names = FALSE)

p_traj <- ggplot(traj_df, aes(x = k, y = beta_median)) +
  geom_hline(yintercept = 0, colour = "grey70") +
  geom_ribbon(aes(ymin = beta_lo, ymax = beta_hi), fill = "steelblue", alpha = 0.25) +
  geom_line(colour = "steelblue", linewidth = 0.8) +
  geom_point(colour = "steelblue", size = 1.8) +
  facet_wrap(~ variant_label, scales = "free_y") +
  scale_x_continuous(breaks = seq_len(K - 1L)) +
  labs(title = sprintf("Top-%d variant trajectories beta[v, k] across k -- %s / %s",
                       length(top_v), opt$dataset, opt$model_subdir),
       subtitle = "Systematic climb with k across many panels = the inflation you asked about",
       x = "cutpoint index k", y = "beta_variant (allele scale), 95% CI") +
  theme_bw(base_size = 9)
ggsave(file.path(OUTDIR, "04_top_variant_trajectories.png"),
       p_traj, width = 10, height = 7, dpi = 150)

# ---- DIAGNOSTIC 4: REFERENCE-SUBLINEAGE CDF vs LOGISTIC -------------------
# Empirical cumulative probability in the reference cluster:
#   P_emp_k = mean_{n in ref}(phen[n] <= k)   for k = 1 .. K-1
# Model-implied reference CDF (reference samples have X_sublineage = 0 by
# treatment-contrast encoding, so sub_eta_ref = 0, and variant_effect = 0):
#   baseline / tight_slab / poolk / latent_scale:  P_mod_k = E[inv_logit(c_data[k] - alpha)]
#   free_cutpoints:                                P_mod_k = E[inv_logit(cutpoints[k])]
# The free_cutpoints variant drops alpha and estimates cutpoints directly, so
# we pull per-draw cutpoints from the fit and evaluate the logistic without alpha.
pull_cutpoints_mat <- function(df, K_1) {
  nm <- sprintf("cutpoints[%d]", seq_len(K_1))
  if (!all(nm %in% names(df))) return(NULL)
  as.matrix(df[, nm])
}
cutpoints_mat <- pull_cutpoints_mat(draws_df, K - 1L)
alpha_vec     <- if ("alpha" %in% names(draws_df)) as.numeric(draws_df$alpha) else NULL

ref_cdf_df <- data.frame()
P_emp <- sapply(seq_len(K - 1L), function(k) mean(phenotype[ref_mask] <= k))

if (!is.null(cutpoints_mat)) {
  message("  cutpoints are estimated parameters -- using per-draw cutpoints (no alpha) for ref-CDF")
  lat_draws <- cutpoints_mat                  # D x (K-1)
  P_mat <- plogis(lat_draws)                  # D x (K-1)
  P_mod_draws <- colMeans(P_mat)
  P_mod_lo    <- apply(P_mat, 2, quantile, 0.025)
  P_mod_hi    <- apply(P_mat, 2, quantile, 0.975)
  cut_anchor  <- colMeans(cutpoints_mat)       # posterior-mean cutpoint locations
  cdf_source  <- "estimated cutpoints (no alpha)"
} else if (!is.null(alpha_vec) && !all(is.na(alpha_vec))) {
  lat_draws <- sapply(seq_len(K - 1L), function(k) cutpoints[k] - alpha_vec)  # D x (K-1)
  P_mat <- plogis(lat_draws)
  P_mod_draws <- colMeans(P_mat)
  P_mod_lo    <- apply(P_mat, 2, quantile, 0.025)
  P_mod_hi    <- apply(P_mat, 2, quantile, 0.975)
  cut_anchor  <- cutpoints                     # data-fixed anchor
  cdf_source  <- "fixed cutpoints - alpha"
} else {
  message("  warn: neither 'alpha' nor 'cutpoints[k]' in draws -- skipping model-implied reference CDF")
  P_mod_draws <- P_mod_lo <- P_mod_hi <- NULL
  cdf_source  <- NA
}

if (!is.null(P_mod_draws)) {
  ref_cdf_df <- data.frame(
    k         = seq_len(K - 1L),
    cutpoint  = cut_anchor,
    P_emp_ref = P_emp,
    P_mod_ref = P_mod_draws,
    P_mod_lo  = P_mod_lo,
    P_mod_hi  = P_mod_hi,
    n_ref     = n_ref,
    source    = cdf_source
  )
  write.csv(ref_cdf_df, file.path(OUTDIR, "ref_sublineage_cdf.csv"), row.names = FALSE)

  p_ref <- ggplot(ref_cdf_df, aes(x = k)) +
    geom_ribbon(aes(ymin = P_mod_lo, ymax = P_mod_hi),
                fill = "steelblue", alpha = 0.25) +
    geom_line(aes(y = P_mod_ref, colour = "model-implied"), linewidth = 0.9) +
    geom_point(aes(y = P_mod_ref, colour = "model-implied"), size = 2) +
    geom_line(aes(y = P_emp_ref, colour = "empirical"), linewidth = 0.9) +
    geom_point(aes(y = P_emp_ref, colour = "empirical"), size = 2) +
    scale_colour_manual(values = c("model-implied" = "steelblue",
                                   "empirical"     = "firebrick")) +
    scale_x_continuous(breaks = seq_len(K - 1L)) +
    labs(title = sprintf("Reference-sublineage cumulative CDF -- %s / %s",
                         opt$dataset, opt$model_subdir),
         subtitle = sprintf("P(phen <= k | reference), n_ref = %d, source: %s",
                            n_ref, cdf_source),
         x = "cutpoint index k", y = "P(phen <= k)", colour = "") +
    theme_bw(base_size = 10) + theme(legend.position = "bottom")
  ggsave(file.path(OUTDIR, "05_ref_cdf.png"),
         p_ref, width = 7, height = 4, dpi = 150)
}

# ---- DIAGNOSTIC 4b: CUTPOINT POSTERIOR vs MIC-GRID ANCHOR (free_cutpoints) -
# Only meaningful for free_cutpoints, which estimates `cutpoints` with a tight
# Gaussian prior around the MIC-grid anchor. cutpoint_drift[k] = cutpoints[k] -
# cutpoint_prior_mean[k]. Non-zero drift at late k co-occurring with variant
# collapse would evidence non-uniform latent-scale misfit (mechanism 3) being
# absorbed by cutpoints rather than variants.
if (!is.null(cutpoints_mat)) {
  drift_cols <- sprintf("cutpoint_drift[%d]", seq_len(K - 1L))
  if (all(drift_cols %in% names(draws_df))) {
    drift_mat <- as.matrix(draws_df[, drift_cols])
  } else {
    # reconstruct drift from cutpoints and the anchor stored as alpha_prior_mean_out
    apm <- if ("alpha_prior_mean_out" %in% names(draws_df))
      median(as.numeric(draws_df$alpha_prior_mean_out)) else 0
    drift_mat <- sweep(cutpoints_mat, 2, cutpoints - apm, `-`)
  }
  cutpoint_df <- data.frame(
    k            = seq_len(K - 1L),
    mic_anchor   = cutpoints,
    cut_median   = apply(cutpoints_mat, 2, median),
    cut_lo       = apply(cutpoints_mat, 2, quantile, 0.025),
    cut_hi       = apply(cutpoints_mat, 2, quantile, 0.975),
    drift_median = apply(drift_mat, 2, median),
    drift_lo     = apply(drift_mat, 2, quantile, 0.025),
    drift_hi     = apply(drift_mat, 2, quantile, 0.975)
  )
  write.csv(cutpoint_df, file.path(OUTDIR, "cutpoint_posterior.csv"), row.names = FALSE)

  p_drift <- ggplot(cutpoint_df, aes(x = k, y = drift_median)) +
    geom_hline(yintercept = 0, colour = "grey60") +
    geom_pointrange(aes(ymin = drift_lo, ymax = drift_hi),
                    colour = "darkorange3", size = 0.5) +
    scale_x_continuous(breaks = seq_len(K - 1L)) +
    labs(title = sprintf("Cutpoint drift from MIC-grid anchor -- %s / %s",
                         opt$dataset, opt$model_subdir),
         subtitle = "drift[k] = cutpoints[k] - (log2(mic_bp[k]) - alpha_prior_mean). Non-uniform drift => mechanism (3) absorbed by cutpoints",
         x = "cutpoint index k", y = "cutpoint_drift (logit units), 95% CI") +
    theme_bw(base_size = 10)
  ggsave(file.path(OUTDIR, "05b_cutpoint_drift.png"),
         p_drift, width = 7, height = 4, dpi = 150)
}

# ---- DIAGNOSTIC 5: OPTIONAL SCALAR HYPERPARAMETERS ------------------------
`%||%` <- function(a, b) if (is.null(a)) b else a
scalar_hyp <- list(
  sigma_latent     = draws_df$sigma_latent_out     %||% draws_df$sigma_latent,
  sigma_dev        = draws_df$sigma_dev_out        %||% draws_df$sigma_dev,
  tau              = draws_df$tau,
  h2_narrow_median = draws_df$h2_narrow_median_k,
  h2_broad_median  = draws_df$h2_broad_median_k
)
hyp_rows <- list()
for (nm in names(scalar_hyp)) {
  v <- scalar_hyp[[nm]]
  if (is.null(v) || all(is.na(v))) next
  hyp_rows[[nm]] <- data.frame(
    parameter = nm,
    n_draws   = length(v),
    mean      = mean(v),
    sd        = sd(v),
    q025      = quantile(v, 0.025),
    q500      = median(v),
    q975      = quantile(v, 0.975),
    row.names = NULL
  )
}
if (length(hyp_rows) > 0) {
  hyp_df <- do.call(rbind, hyp_rows)
  write.csv(hyp_df, file.path(OUTDIR, "scalar_hyperparameters.csv"), row.names = FALSE)
  message("\n=== scalar hyperparameters present in fit ===")
  print(hyp_df, row.names = FALSE)
}

# ---- DIAGNOSTIC 6: y_rep_ppc vs OBSERVED CATEGORY FREQUENCIES ------------
if (!is.null(y_rep_mat)) {
  # y_rep_mat is D x N_ppc. Compute per-draw category freqs then summarise.
  cat_freq_rep <- t(apply(y_rep_mat, 1, function(row)
    tabulate(row, nbins = K) / length(row)))
  obs_freq <- tabulate(phenotype, nbins = K) / N
  cat_df <- data.frame(
    category = factor(seq_len(K)),
    observed = obs_freq,
    rep_mean = colMeans(cat_freq_rep),
    rep_lo   = apply(cat_freq_rep, 2, quantile, 0.05),
    rep_hi   = apply(cat_freq_rep, 2, quantile, 0.95)
  )
  write.csv(cat_df, file.path(OUTDIR, "ppc_category_frequencies.csv"), row.names = FALSE)

  p_cat <- ggplot(cat_df, aes(x = category)) +
    geom_col(aes(y = observed), fill = "grey85", colour = "grey40", width = 0.85) +
    geom_pointrange(aes(y = rep_mean, ymin = rep_lo, ymax = rep_hi),
                    colour = "steelblue", size = 0.5) +
    labs(title = sprintf("PPC category frequencies -- %s / %s",
                         opt$dataset, opt$model_subdir),
         subtitle = "grey = observed; points = y_rep_ppc posterior mean, bars = 90% CI",
         x = "MIC category (1 = lowest)", y = "frequency") +
    theme_bw(base_size = 10)
  ggsave(file.path(OUTDIR, "06_ppc_category_frequencies.png"),
         p_cat, width = 7, height = 4, dpi = 150)
}

# ---- ONE-ROW SUMMARY METRICS ---------------------------------------------
metrics <- data.frame(
  model_subdir                  = opt$model_subdir,
  dataset                       = opt$dataset,
  N                             = N,
  V                             = V,
  K                             = K,
  n_ref_samples                 = n_ref,
  median_abs_beta_k1            = summary_by_k$median_abs_beta[1],
  median_abs_beta_kKminus1      = summary_by_k$median_abs_beta[K - 1L],
  max_abs_beta_k1               = summary_by_k$max_abs_beta[1],
  max_abs_beta_kKminus1         = summary_by_k$max_abs_beta[K - 1L],
  inflation_ratio_median_abs    = summary_by_k$median_abs_beta[K - 1L] /
                                  max(summary_by_k$median_abs_beta[1], 1e-9),
  median_post_sd_k1             = summary_by_k$median_post_sd[1],
  median_post_sd_kKminus1       = summary_by_k$median_post_sd[K - 1L],
  frac_above_k1                 = frac_above_k[1],
  frac_above_kKminus1           = frac_above_k[K - 1L],
  stringsAsFactors              = FALSE
)
write.csv(metrics, file.path(OUTDIR, "metrics_summary.csv"), row.names = FALSE)

message("\n=== per-k summary ===")
print(summary_by_k, row.names = FALSE, digits = 3)
message("\n=== inflation metric ===")
message(sprintf("median |beta| at k=1:     %.4f", summary_by_k$median_abs_beta[1]))
message(sprintf("median |beta| at k=%d:    %.4f",
                K - 1L, summary_by_k$median_abs_beta[K - 1L]))
message(sprintf("ratio (late / early):     %.2fx",
                metrics$inflation_ratio_median_abs))
message("\noutputs written to: ", OUTDIR)
