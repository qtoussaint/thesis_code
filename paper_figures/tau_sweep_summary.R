#!/usr/bin/env Rscript
# Tau sweep summary figure for the PPOM shrinkage-prior grid.
#
# Summarises what the global shrinkage scale tau buys us: out-of-sample
# prediction accuracy (Panel A) traded off against variant sparsity (Panel B),
# across the fixed-tau grid and two MIC binnings. This is the standard
# diagnostic for choosing a shrinkage level.
#
# The grid is tau in {0.001, 0.01, 0.05, 1, 1.5, 2, 3, 5} crossed with four
# prior configs: slab3/lambda1, slab5/lambda1, slab3/lambda2, slab5/lambda2.
# Coverage is uneven -- the lambda1 series were only run at the low taus
# {0.001, 0.01, 0.05, 1}, slab3/lambda2 spans the full axis, and slab5/lambda2
# was only run at the high taus {1.5, 2, 3, 5}. Runs that don't exist render as
# gaps (the line simply breaks), so each series shows just its own coverage.
# (The standalone estimated-tau variants -- lambda1, lognormaltau02 -- are not
# part of the fixed-tau sweep and are excluded.)
#
# Reads the same per-run files as the prediction pipeline:
#   inference_ppc/prediction_accuracy_metrics.csv   (RPSS, balanced accuracy)
#   fitted_model/depruned_variant_effects.csv       (per-variant signif flags)
#
# Run with:
#   mamba activate gwas_pipeline
#   Rscript paper_figures/tau_sweep_summary.R
#
# Output: <output_dir>/tau_sweep_summary.{png,csv}

suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
  library(viridis)
})

RESULTS_ROOT <- "/nfs/research/jlees/jacqueline/thesis_results"
SWEEP_ROOT   <- file.path(RESULTS_ROOT, "compare_ordinal_models")
OUTPUT_DIR   <- file.path(RESULTS_ROOT, "paper_figures", "parameter_sweeps")
MODEL_PREFIX <- "final_ordered_categorical_PPOM_free_cutpoints_wide_drift"

# -----------------------------------------------------------------------------
# Grid spec
# -----------------------------------------------------------------------------

TAUS <- c(0.001, 0.01, 0.05, 1, 1.5, 2, 3, 5)

# Encode tau -> dir token. Note 0.05 -> "0p05" (NOT "05"): a stray
# ...fixedtau05_slab3 dir exists on disk and must not be matched. The high taus
# use a bare integer token (e.g. 2 -> "2") and 1.5 -> "1p5".
tau_token <- function(tau) {
  switch(as.character(tau),
         "0.001" = "0p001", "0.01" = "0p01",
         "0.05"  = "0p05",  "1"    = "1",
         "1.5"   = "1p5",   "2"    = "2",
         "3"     = "3",     "5"    = "5",
         stop("unmapped tau: ", tau))
}

# Each series -> the dir suffix that follows fixedtau<token>. Coverage differs
# by series (see header); missing runs are rendered as gaps.
SERIES <- list(
  list(label = "slab3/λ1", suffix = "_slab3"),          # slab_scale=3, lambda~C(0,1)
  list(label = "slab5/λ1", suffix = "_slab5"),          # slab_scale=5, lambda~C(0,1)
  list(label = "slab3/λ2", suffix = "_slab3_lambda2"),  # slab_scale=3, lambda~C(0,2)
  list(label = "slab5/λ2", suffix = "_slab5_lambda2")   # slab_scale=5, lambda~C(0,2)
)
series_levels <- vapply(SERIES, `[[`, "", "label")

DATASETS <- list(
  list(id = "02_spn_penicillin_MIC",               label = "doubling (≥5%), K=8", K = 8L),
  list(id = "16_spn_penicillin_MIC_minimabinning", label = "minima, K=5",         K = 5L)
)
dataset_levels <- vapply(DATASETS, `[[`, "", "label")

# Reference baseline PPOM runs: a single estimated-tau fit per dataset. Not part
# of the fixed-tau grid, so they have no tau position and are drawn as horizontal
# reference lines (one per facet) rather than sweep points.
REF_LABEL <- "PPOM (estimated τ)"
REF_RUNS <- list(
  list(dataset = "doubling (≥5%), K=8",
       dir = file.path(RESULTS_ROOT, "gwas_spn_penicillin", "inference",
                       "02_spn_penicillin_MIC_PPOM")),
  list(dataset = "minima, K=5",
       dir = file.path(RESULTS_ROOT, "gwas_spn_penicillin", "inference",
                       "16_spn_penicillin_MIC_minimabinning_PPOM"))
)

# -----------------------------------------------------------------------------
# Path builder
# -----------------------------------------------------------------------------

model_dir <- function(suffix, tau, ds_id)
  file.path(SWEEP_ROOT,
            paste0(MODEL_PREFIX, "_fixedtau", tau_token(tau), suffix),
            ds_id)

metrics_path <- function(suffix, tau, ds_id)
  file.path(model_dir(suffix, tau, ds_id), "inference_ppc", "prediction_accuracy_metrics.csv")

effects_path <- function(suffix, tau, ds_id)
  file.path(model_dir(suffix, tau, ds_id), "fitted_model", "depruned_variant_effects.csv")

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

# Count significant variants in one depruned_variant_effects.csv. signif is the
# lowercase string "true"/"false". colClasses drops the heavy median /
# cutpoint_MIC float columns so the 226k-row files load fast under base R.
# Returns n_unique (distinct variant_ids with any signif), n_pairs (total
# signif (variant,cutpoint) rows).
count_signif <- function(csv_path) {
  na <- c(n_unique = NA_real_, n_pairs = NA_real_)
  if (!file.exists(csv_path)) {
    missing_files <<- c(missing_files, csv_path)
    return(na)
  }
  df <- tryCatch(
    read.csv(csv_path, stringsAsFactors = FALSE,
             colClasses = c(variant_id = "integer", median = "NULL",
                            signif = "character", cutpoint = "integer",
                            cutpoint_MIC = "NULL")),
    error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) return(na)
  is_sig <- tolower(trimws(df$signif)) == "true"
  c(n_unique = length(unique(df$variant_id[is_sig])),
    n_pairs  = sum(is_sig))
}

# -----------------------------------------------------------------------------
# Build tidy long data frame (4 series x 8 tau x 2 datasets, minus the
# series/tau cells with no run on disk, which come through as NA gaps)
# -----------------------------------------------------------------------------

ACC_METRICS <- c("rpss_uniform", "rpss_frequency", "bacc")
rows <- list()

for (s in SERIES) {
  for (tau in TAUS) {
    for (d in DATASETS) {
      acc <- read_metrics(metrics_path(s$suffix, tau, d$id), ACC_METRICS)
      spa <- count_signif(effects_path(s$suffix, tau, d$id))

      common <- data.frame(series = s$label, tau = tau, dataset = d$label,
                           stringsAsFactors = FALSE)
      for (m in ACC_METRICS)
        rows[[length(rows) + 1]] <- cbind(common,
          quantity = m, value = unname(acc[[m]]), panel = "A")
      rows[[length(rows) + 1]] <- cbind(common,
        quantity = "n_unique_variants", value = unname(spa[["n_unique"]]), panel = "B")
      rows[[length(rows) + 1]] <- cbind(common,
        quantity = "n_signif_pairs", value = unname(spa[["n_pairs"]]), panel = "B")
    }
  }
}

df_long <- do.call(rbind, rows)
df_long$series  <- factor(df_long$series,  levels = series_levels)
df_long$dataset <- factor(df_long$dataset, levels = dataset_levels)

# Reference baseline values: one row per (dataset, quantity), tau-independent.
ref_rows <- list()
for (r in REF_RUNS) {
  acc <- read_metrics(file.path(r$dir, "inference_ppc", "prediction_accuracy_metrics.csv"),
                      ACC_METRICS)
  spa <- count_signif(file.path(r$dir, "fitted_model", "depruned_variant_effects.csv"))
  for (m in ACC_METRICS)
    ref_rows[[length(ref_rows) + 1]] <- data.frame(
      dataset = r$dataset, quantity = m, value = unname(acc[[m]]),
      panel = "A", stringsAsFactors = FALSE)
  ref_rows[[length(ref_rows) + 1]] <- data.frame(
    dataset = r$dataset, quantity = "n_unique_variants", value = unname(spa[["n_unique"]]),
    panel = "B", stringsAsFactors = FALSE)
}
ref_df <- do.call(rbind, ref_rows)
ref_df$dataset <- factor(ref_df$dataset, levels = dataset_levels)
ref_df$ref     <- REF_LABEL

if (length(missing_files) > 0) {
  message(sprintf("NOTE: %d input CSV(s) missing (rendered as gaps):", length(missing_files)))
  for (f in missing_files) message("  ", f)
} else {
  message("All input CSVs found.")
}

# -----------------------------------------------------------------------------
# Aesthetics
# -----------------------------------------------------------------------------

series_cols <- viridis::viridis(length(series_levels), option = "D", end = 0.9)
names(series_cols) <- series_levels

# Reference baseline gets a contrasting colour outside the viridis ramp.
colour_cols    <- c(series_cols, setNames("#D81B60", REF_LABEL))
colour_levels  <- c(series_levels, REF_LABEL)

base_theme <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92", colour = NA),
        legend.position = "bottom")

metric_labels <- c(rpss_uniform   = "RPSS (uniform)",
                   rpss_frequency = "RPSS (frequency)",
                   bacc           = "balanced accuracy")

tau_labels <- c("0.001", "0.01", "0.05", "1", "1.5", "2", "3", "5")

# -----------------------------------------------------------------------------
# Panel A: prediction accuracy vs tau
# -----------------------------------------------------------------------------

dfA <- df_long[df_long$panel == "A", ]
dfA$quantity <- factor(dfA$quantity, levels = ACC_METRICS)

refA <- ref_df[ref_df$panel == "A", ]
refA$quantity <- factor(refA$quantity, levels = ACC_METRICS)

panel_a <- ggplot(dfA, aes(x = tau, y = value, colour = series, group = series)) +
  geom_hline(data = refA, aes(yintercept = value, colour = ref),
             linetype = "dashed", linewidth = 0.7, inherit.aes = FALSE) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.4) +
  facet_grid(quantity ~ dataset, scales = "free_y",
             labeller = labeller(quantity = metric_labels)) +
  scale_x_log10(breaks = TAUS, labels = tau_labels) +
  scale_colour_manual(values = colour_cols, name = "prior config",
                      limits = colour_levels) +
  # Zoom each metric in to [0.5, 1] so small differences are legible while the
  # axis still spans up to the theoretical maximum; coord_cartesian clips the
  # view without dropping points.
  coord_cartesian(ylim = c(0.5, 1)) +
  labs(x = expression("global shrinkage scale " * tau), y = "metric value") +
  base_theme

# -----------------------------------------------------------------------------
# Panel B: variant sparsity vs tau (unique significant variants)
# -----------------------------------------------------------------------------

dfB <- df_long[df_long$panel == "B" & df_long$quantity == "n_unique_variants", ]
refB <- ref_df[ref_df$panel == "B", ]

panel_b <- ggplot(dfB, aes(x = tau, y = value, colour = series, group = series)) +
  geom_hline(data = refB, aes(yintercept = value, colour = ref),
             linetype = "dashed", linewidth = 0.7, inherit.aes = FALSE) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.4) +
  facet_wrap(~ dataset, nrow = 1) +
  scale_x_log10(breaks = TAUS, labels = tau_labels) +
  scale_colour_manual(values = colour_cols, name = "prior config",
                      limits = colour_levels) +
  labs(x = expression("global shrinkage scale " * tau),
       y = "# significant variants (unique)") +
  base_theme

# -----------------------------------------------------------------------------
# Assemble + save
# -----------------------------------------------------------------------------

legend <- cowplot::get_legend(panel_a)
panels <- cowplot::plot_grid(
  panel_a + theme(legend.position = "none"),
  panel_b + theme(legend.position = "none"),
  ncol = 1, labels = c("A", "B"), label_size = 24, label_fontface = "bold",
  rel_heights = c(2.4, 1))
figure <- cowplot::plot_grid(panels, legend, ncol = 1, rel_heights = c(1, 0.06))

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
png_path <- file.path(OUTPUT_DIR, "tau_sweep_summary.png")
csv_path <- file.path(OUTPUT_DIR, "tau_sweep_summary.csv")

ggsave(png_path, figure, width = 9, height = 11, dpi = 300, bg = "white")

# Combine sweep + reference rows (reference has no tau) for the saved table.
ref_out <- data.frame(series = REF_LABEL, tau = NA_real_, dataset = ref_df$dataset,
                      quantity = ref_df$quantity, value = ref_df$value,
                      panel = ref_df$panel, stringsAsFactors = FALSE)
out_df <- rbind(df_long, ref_out)
write.csv(out_df[order(out_df$panel, out_df$quantity, out_df$dataset,
                       out_df$series, out_df$tau), ],
          csv_path, row.names = FALSE)

message("wrote ", png_path)
message("wrote ", csv_path)
