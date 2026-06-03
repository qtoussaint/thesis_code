#!/usr/bin/env Rscript
# Intercept-prior-SD sweep summary figures for the logistic and continuous models.
#
# The production logistic_inference.stan and continuous_inference.stan both pin
# the intercept prior at alpha_prior_sd = 0.5. These sweeps widen it to 1.5 / 3 / 5
# (one copy of the model each) to check how sensitive the fit is to that width.
#
# Each figure mirrors tau_sweep_summary.R:
#   Panel A: out-of-sample prediction accuracy vs alpha_prior_sd (one row per metric)
#   Panel B: variant sparsity (# significant variants) vs alpha_prior_sd
# The production default (SD = 0.5) is drawn as a dashed reference line in Panel A.
# Its variant count is unavailable (no depruned CSV; RATE-based, not comparable to
# these --norate runs), so it is omitted from Panel B.
#
# Reads the same per-run files as the prediction pipeline:
#   inference_ppc/prediction_accuracy_metrics.csv   (accuracy metrics)
#   fitted_model/depruned_variant_effects.csv       (per-variant signif flags)
#
# Run with:
#   mamba activate gwas_pipeline
#   Rscript paper_figures/alpha_sweep_summary.R
#
# Output: <output_dir>/{logistic,continuous}_alpha_sweep_summary.{png,csv}

suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
  library(viridis)
})

RESULTS_ROOT <- "/nfs/research/jlees/jacqueline/thesis_results"
SWEEP_ROOT   <- file.path(RESULTS_ROOT, "compare_ordinal_models")
OUTPUT_DIR   <- file.path(RESULTS_ROOT, "paper_figures")

BASE_SD  <- 0.5            # production default, drawn as a reference line

SWEEP_LABEL <- "swept α-prior SD"
REF_LABEL   <- "production default (SD = 0.5)"

# 1.5 -> "1p5", 3 -> "3", 5 -> "5" (matches the *_alphasd<token> dir names)
sd_token <- function(sd) gsub("\\.", "p", as.character(sd))

# -----------------------------------------------------------------------------
# Per-family config: the two model families plotted, with their accuracy metrics
# (named vector: csv_column -> facet label) and the SD = 0.5 baseline run.
# -----------------------------------------------------------------------------

MODELS <- list(
  list(
    key      = "logistic",
    prefix   = "logistic_inference_alphasd",
    dataset  = "01_spn_penicillin_binary",
    baseline = file.path(RESULTS_ROOT, "gwas_spn_penicillin", "inference",
                         "01_spn_penicillin_binary_logistic"),
    metrics  = c(bacc = "balanced accuracy", auc = "AUC", brier = "Brier score"),
    sds      = c(0.075, 1, 1.25, 1.5, 3, 5),
    out      = "logistic_alpha_sweep_summary"
  ),
  list(
    key      = "continuous",
    prefix   = "continuous_inference_alphasd",
    dataset  = "03_spn_penicillin_continuous",
    baseline = file.path(RESULTS_ROOT, "gwas_spn_penicillin", "inference",
                         "03_spn_penicillin_continuous_continuous"),
    metrics  = c(r_squared = "R²", rmse = "RMSE", crps = "CRPS"),
    sds      = c(1.5, 3, 5),
    out      = "continuous_alpha_sweep_summary"
  )
)

# -----------------------------------------------------------------------------
# Readers (record missing files so gaps are reported, not fatal)
# -----------------------------------------------------------------------------

missing_files <- character(0)

# Read named metric columns from a one-row metrics CSV as numerics. A missing
# file or column yields NA (and the file is recorded so we can report gaps).
read_metrics <- function(csv_path, cols) {
  if (!file.exists(csv_path)) {
    missing_files <<- c(missing_files, csv_path)
    return(setNames(rep(NA_real_, length(cols)), cols))
  }
  row <- tryCatch(
    read.csv(csv_path, stringsAsFactors = FALSE, check.names = FALSE)[1, , drop = FALSE],
    error = function(e) NULL)
  if (is.null(row) || nrow(row) == 0) return(setNames(rep(NA_real_, length(cols)), cols))
  vapply(cols, function(c) if (!c %in% names(row)) NA_real_ else suppressWarnings(as.numeric(row[[c]])),
         numeric(1))
}

# Count significant variants in one depruned_variant_effects.csv. These models
# have a single effect per variant (header variant_id,median,signif -- no
# cutpoint), so the count is just the number of signif == "true" rows.
count_signif <- function(csv_path) {
  if (!file.exists(csv_path)) {
    missing_files <<- c(missing_files, csv_path)
    return(NA_real_)
  }
  df <- tryCatch(
    read.csv(csv_path, stringsAsFactors = FALSE,
             colClasses = c(variant_id = "integer", median = "NULL",
                            signif = "character")),
    error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) return(NA_real_)
  sum(tolower(trimws(df$signif)) == "true")
}

# -----------------------------------------------------------------------------
# Shared aesthetics
# -----------------------------------------------------------------------------

colour_cols   <- c(setNames(viridis::viridis(1, option = "D", end = 0.45), SWEEP_LABEL),
                   setNames("#D81B60", REF_LABEL))
colour_levels <- c(SWEEP_LABEL, REF_LABEL)

base_theme <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92", colour = NA),
        legend.position = "bottom")

# -----------------------------------------------------------------------------
# Build + render one family
# -----------------------------------------------------------------------------

build_family <- function(cfg) {
  metric_cols   <- names(cfg$metrics)
  metric_labels <- cfg$metrics
  sds           <- cfg$sds

  # Sweep points: one run per SD.
  rows <- list()
  for (sd in sds) {
    run_dir <- file.path(SWEEP_ROOT, paste0(cfg$prefix, sd_token(sd)), cfg$dataset)
    acc <- read_metrics(file.path(run_dir, "inference_ppc", "prediction_accuracy_metrics.csv"),
                        metric_cols)
    nsig <- count_signif(file.path(run_dir, "fitted_model", "depruned_variant_effects.csv"))
    for (m in metric_cols)
      rows[[length(rows) + 1]] <- data.frame(sd = sd, quantity = m,
                                             value = unname(acc[[m]]), panel = "A",
                                             stringsAsFactors = FALSE)
    rows[[length(rows) + 1]] <- data.frame(sd = sd, quantity = "n_signif_variants",
                                           value = nsig, panel = "B",
                                           stringsAsFactors = FALSE)
  }
  df_long <- do.call(rbind, rows)

  # Baseline (SD = 0.5): accuracy only -> dashed reference line in Panel A.
  base_acc <- read_metrics(file.path(cfg$baseline, "inference_ppc",
                                     "prediction_accuracy_metrics.csv"), metric_cols)
  ref_df <- data.frame(quantity = metric_cols, value = unname(base_acc[metric_cols]),
                       ref = REF_LABEL, stringsAsFactors = FALSE)

  # --- Panel A: accuracy vs SD, one row per metric -------------------------
  dfA <- df_long[df_long$panel == "A", ]
  dfA$quantity <- factor(dfA$quantity, levels = metric_cols)
  ref_df$quantity <- factor(ref_df$quantity, levels = metric_cols)

  panel_a <- ggplot(dfA, aes(x = sd, y = value)) +
    geom_hline(data = ref_df, aes(yintercept = value, colour = ref),
               linetype = "dashed", linewidth = 0.7, inherit.aes = FALSE) +
    geom_line(aes(colour = SWEEP_LABEL), linewidth = 0.7) +
    geom_point(aes(colour = SWEEP_LABEL), size = 2.4) +
    facet_wrap(~ quantity, ncol = 1, scales = "free_y",
               labeller = labeller(quantity = metric_labels)) +
    scale_x_continuous(breaks = sds) +
    scale_colour_manual(values = colour_cols, name = NULL, limits = colour_levels) +
    labs(x = "intercept prior SD", y = "metric value") +
    base_theme

  # --- Panel B: # significant variants vs SD -------------------------------
  dfB <- df_long[df_long$panel == "B", ]
  panel_b <- ggplot(dfB, aes(x = sd, y = value)) +
    geom_line(aes(colour = SWEEP_LABEL), linewidth = 0.7) +
    geom_point(aes(colour = SWEEP_LABEL), size = 2.4) +
    scale_x_continuous(breaks = sds) +
    scale_colour_manual(values = colour_cols, name = NULL, limits = colour_levels) +
    expand_limits(y = 0) +
    labs(x = "intercept prior SD", y = "# significant variants") +
    base_theme

  legend <- cowplot::get_legend(panel_a)
  panels <- cowplot::plot_grid(
    panel_a + theme(legend.position = "none"),
    panel_b + theme(legend.position = "none"),
    ncol = 1, labels = c("A", "B"), label_size = 24, label_fontface = "bold",
    rel_heights = c(2.4, 1))
  figure <- cowplot::plot_grid(panels, legend, ncol = 1, rel_heights = c(1, 0.06))

  png_path <- file.path(OUTPUT_DIR, paste0(cfg$out, ".png"))
  csv_path <- file.path(OUTPUT_DIR, paste0(cfg$out, ".csv"))
  ggsave(png_path, figure, width = 6.5, height = 10, dpi = 300, bg = "white")

  # Saved table: sweep rows + the baseline accuracy rows (sd = 0.5).
  ref_out <- data.frame(sd = BASE_SD, quantity = as.character(ref_df$quantity),
                        value = ref_df$value, panel = "A", stringsAsFactors = FALSE)
  out_df <- rbind(df_long, ref_out)
  write.csv(out_df[order(out_df$panel, out_df$quantity, out_df$sd), ],
            csv_path, row.names = FALSE)

  message("wrote ", png_path)
  message("wrote ", csv_path)
}

# -----------------------------------------------------------------------------
# Run both families
# -----------------------------------------------------------------------------

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
for (cfg in MODELS) build_family(cfg)

if (length(missing_files) > 0) {
  message(sprintf("NOTE: %d input CSV(s) missing (rendered as gaps):", length(missing_files)))
  for (f in missing_files) message("  ", f)
} else {
  message("All input CSVs found.")
}
