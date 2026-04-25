#!/usr/bin/env Rscript
# Cross-dataset synthesis of beta_variant cutpoint inflation across PPOM variants.
#
# Consumes:
#   <results_root>/_compare_PPOM/<dataset>/inflation_ratio_table.csv  (from compare_ppom_variants.R)
#   <results_root>/<variant>/<dataset>/ppc/diagnose/summary_by_k.csv  (from diagnose_ppom_cutpoint_inflation.R)
#   <results_root>/<variant>/<dataset>/ppc/diagnose/metrics_summary.csv
#
# Skips any (variant x dataset) combo whose files are missing -- prints a single
# block to stderr at the end listing skipped combos.
#
# Outputs (to <results_root>/_compare_PPOM/_cross_dataset/):
#   inflation_decision_table.csv  one row per (variant x dataset), median + tail metrics
#   inflation_heatmap.png         tile heatmap, fill = inflation_ratio_max_abs (log)
#   inflation_ratio_lines.png     line plot, x = dataset, colour = variant, faceted by metric
#
# Usage (activate gwas_pipeline first):
#   Rscript analysis/aggregate_compare_PPOM.R \
#     --datasets 02_spn_penicillin_MIC,10_spn_penicillin_MIC_coarse_dilutions,11_spn_penicillin_MIC_large_minbin \
#     --model_subdirs final_ordered_categorical_PPOM_tight_alpha_tau1,final_ordered_categorical_PPOM_tight_slab,final_ordered_categorical_PPOM_poolk,final_ordered_categorical_PPOM_latent_scale,final_ordered_categorical_PPOM_free_cutpoints

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

option_list <- list(
  make_option("--datasets", type = "character", default = NULL,
              help = "Comma-separated dataset nicknames"),
  make_option("--model_subdirs", type = "character", default = NULL,
              help = "Comma-separated model subdir names (under results_root)"),
  make_option("--results_root", type = "character",
              default = "/nfs/research/jlees/jacqueline/thesis_results/compare_ordinal_models",
              help = "Parent dir holding per-model result subdirs"),
  make_option("--out_subdir", type = "character", default = "_compare_PPOM/_cross_dataset",
              help = "Output subdir under results_root [_compare_PPOM/_cross_dataset]"),
  make_option("--heatmap_metric", type = "character", default = "inflation_ratio_max_abs",
              help = "Which ratio drives the heatmap fill [inflation_ratio_max_abs]")
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$datasets) || is.null(opt$model_subdirs))
  stop("Both --datasets and --model_subdirs are required")

datasets <- trimws(strsplit(opt$datasets, ",")[[1]])
models   <- trimws(strsplit(opt$model_subdirs, ",")[[1]])
OUTDIR   <- file.path(opt$results_root, opt$out_subdir)
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

short_label <- function(nm) {
  nm <- sub("^final_ordered_categorical_PPOM_", "", nm)
  nm <- sub("^final_no_horseshoe_", "noHS_", nm)
  nm
}

skipped <- character()
rows <- list()

for (ds in datasets) {
  for (m in models) {
    diag_dir <- file.path(opt$results_root, m, ds, "ppc", "diagnose")
    metrics_f <- file.path(diag_dir, "metrics_summary.csv")
    sumk_f    <- file.path(diag_dir, "summary_by_k.csv")

    if (!file.exists(metrics_f) || !file.exists(sumk_f)) {
      skipped <- c(skipped, sprintf("%s / %s", m, ds))
      next
    }

    metrics <- read.csv(metrics_f)
    sumk    <- read.csv(sumk_f)

    # k=1 and k=K-1 rows from summary_by_k
    K_minus_1 <- max(sumk$k)
    sk1 <- sumk[sumk$k == 1L, ]
    skK <- sumk[sumk$k == K_minus_1, ]
    if (nrow(sk1) != 1 || nrow(skK) != 1) {
      warning("Unexpected summary_by_k row count for ", m, " / ", ds, " -- skipping")
      skipped <- c(skipped, sprintf("%s / %s", m, ds))
      next
    }

    ratio_or_na <- function(num, den) {
      if (length(num) != 1 || length(den) != 1) return(NA_real_)
      if (is.na(num) || is.na(den) || den == 0) return(NA_real_)
      num / den
    }

    rows[[length(rows) + 1L]] <- data.frame(
      dataset                    = ds,
      model_subdir               = m,
      model_short                = short_label(m),
      N                          = metrics$N,
      V                          = metrics$V,
      K                          = metrics$K,
      n_ref_samples              = metrics$n_ref_samples,
      median_abs_beta_k1         = sk1$median_abs_beta,
      median_abs_beta_kKminus1   = skK$median_abs_beta,
      q95_abs_beta_k1            = sk1$q95_abs_beta,
      q95_abs_beta_kKminus1      = skK$q95_abs_beta,
      max_abs_beta_k1            = sk1$max_abs_beta,
      max_abs_beta_kKminus1      = skK$max_abs_beta,
      median_post_sd_k1          = sk1$median_post_sd,
      median_post_sd_kKminus1    = skK$median_post_sd,
      q95_post_sd_k1             = sk1$q95_post_sd,
      q95_post_sd_kKminus1       = skK$q95_post_sd,
      inflation_ratio_median_abs = ratio_or_na(skK$median_abs_beta, sk1$median_abs_beta),
      inflation_ratio_q95_abs    = ratio_or_na(skK$q95_abs_beta,    sk1$q95_abs_beta),
      inflation_ratio_max_abs    = ratio_or_na(skK$max_abs_beta,    sk1$max_abs_beta),
      inflation_ratio_post_sd    = ratio_or_na(skK$median_post_sd,  sk1$median_post_sd),
      frac_above_k1              = metrics$frac_above_k1,
      frac_above_kKminus1        = metrics$frac_above_kKminus1,
      stringsAsFactors           = FALSE
    )
  }
}

if (length(rows) == 0)
  stop("No (variant x dataset) combos had diagnostic outputs -- nothing to aggregate")

decision <- bind_rows(rows) %>%
  arrange(dataset, inflation_ratio_max_abs)

decision_path <- file.path(OUTDIR, "inflation_decision_table.csv")
write.csv(decision, decision_path, row.names = FALSE)

# ---- HEATMAP --------------------------------------------------------------
heat_metric <- opt$heatmap_metric
if (!heat_metric %in% colnames(decision))
  stop("heatmap_metric '", heat_metric, "' not in decision table")

# Order variants by mean rank of the heatmap metric across datasets so the
# best-controlled variants sit on the left.
variant_order <- decision %>%
  group_by(model_short) %>%
  summarise(mean_metric = mean(.data[[heat_metric]], na.rm = TRUE), .groups = "drop") %>%
  arrange(mean_metric) %>%
  pull(model_short)

heat_df <- decision %>%
  mutate(model_short = factor(model_short, levels = variant_order))

p_heat <- ggplot(heat_df, aes(x = model_short, y = dataset, fill = .data[[heat_metric]])) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2g", .data[[heat_metric]])),
            size = 3, colour = "black") +
  scale_fill_gradient(name = heat_metric,
                      low = "#2c7fb8", high = "#e31a1c",
                      trans = "log10", na.value = "grey80") +
  labs(title = sprintf("Cross-dataset cutpoint inflation -- fill = %s", heat_metric),
       subtitle = "Log fill scale. Closer to 1 = less inflation. Cells with NA = combo missing.",
       x = "PPOM variant (sorted by mean metric, best left)",
       y = "dataset") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(file.path(OUTDIR, "inflation_heatmap.png"),
       p_heat, width = 8, height = 4.5, dpi = 150)

# ---- LINE PLOT FACETED BY METRIC ------------------------------------------
metric_long <- decision %>%
  select(dataset, model_short,
         median = inflation_ratio_median_abs,
         q95    = inflation_ratio_q95_abs,
         max    = inflation_ratio_max_abs,
         post_sd = inflation_ratio_post_sd) %>%
  pivot_longer(cols = c(median, q95, max, post_sd),
               names_to = "metric", values_to = "ratio") %>%
  mutate(metric = factor(metric, levels = c("median", "q95", "max", "post_sd")))

p_lines <- ggplot(metric_long,
                  aes(x = dataset, y = ratio, colour = model_short, group = model_short)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.2) +
  scale_y_log10() +
  facet_wrap(~ metric, scales = "free_y") +
  labs(title = "Late-k / early-k inflation ratios by metric and dataset",
       subtitle = "Dashed line = 1 (no inflation). Log y. Lower & flatter = better.",
       x = NULL, y = "inflation ratio (k=K-1 / k=1)",
       colour = "variant") +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(file.path(OUTDIR, "inflation_ratio_lines.png"),
       p_lines, width = 10, height = 6.5, dpi = 150)

# ---- LOG --------------------------------------------------------------------
message("\n=== cross-dataset decision table (sorted by dataset, then inflation_ratio_max_abs) ===")
print(decision %>%
        select(dataset, model_short,
               inflation_ratio_median_abs,
               inflation_ratio_q95_abs,
               inflation_ratio_max_abs,
               inflation_ratio_post_sd),
      row.names = FALSE, digits = 3)

message("\nWritten:")
message("  ", decision_path)
message("  ", file.path(OUTDIR, "inflation_heatmap.png"))
message("  ", file.path(OUTDIR, "inflation_ratio_lines.png"))

if (length(skipped) > 0) {
  message("\nSKIPPED (missing diagnostic outputs):")
  for (s in skipped) message("  - ", s)
}
