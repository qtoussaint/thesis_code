#!/usr/bin/env Rscript
# Cross-dataset synthesis of prediction accuracy metrics across model variants.
#
# Consumes (one row, six columns) per (variant x dataset):
#   <results_root>/<variant>/<dataset>/inference_ppc/prediction_accuracy_metrics.csv
#     columns: rps_mean_scaled, rps_median_scaled, rpss_uniform, rpss_frequency,
#              bacc, ppv
#
# Direction: rps_*_scaled lower is better; rpss_*, bacc, ppv higher is better.
#
# Default behaviour is to auto-discover every available combo. Supply
# --model_subdirs and/or --datasets to restrict.
#
# Outputs (to <results_root>/<out_subdir>/):
#   prediction_accuracy_table.csv         wide, one row per (variant x dataset)
#   prediction_accuracy_long.csv          long form for plotting
#   prediction_accuracy_rank_table.csv    per-variant mean rank across datasets
#   prediction_accuracy_heatmap.png       tile heatmap, faceted by metric
#   prediction_accuracy_lines.png         line plot, x = dataset, faceted metric
#   prediction_accuracy_rank_bars.png     bar plot of mean rank per metric
#
# Usage (activate gwas_pipeline first):
#   Rscript analysis/aggregate_prediction_accuracy.R
#   Rscript analysis/aggregate_prediction_accuracy.R \
#     --model_subdirs final_ordered_categorical_PPOM_poolk,final_ordered_categorical_PPOM_slab50 \
#     --datasets 02_spn_penicillin_MIC,11_spn_penicillin_MIC_large_minbin

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
})

option_list <- list(
  make_option("--results_root", type = "character",
              default = "/nfs/research/jlees/jacqueline/thesis_results/compare_ordinal_models",
              help = "Parent dir holding per-model result subdirs"),
  make_option("--model_subdirs", type = "character", default = NULL,
              help = "Comma-separated model subdir names. Default: auto-discover."),
  make_option("--datasets", type = "character", default = NULL,
              help = "Comma-separated dataset nicknames. Default: auto-discover."),
  make_option("--out_subdir", type = "character",
              default = "_compare_PPOM/_cross_dataset_prediction_accuracy",
              help = "Output subdir under results_root"),
  make_option("--rank_metric", type = "character", default = "rpss_frequency",
              help = "Metric used to order variants in the heatmap [rpss_frequency]")
)
opt <- parse_args(OptionParser(option_list = option_list))

OUTDIR <- file.path(opt$results_root, opt$out_subdir)
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

short_label <- function(nm) {
  nm <- sub("^final_ordered_categorical_PPOM_", "", nm)
  nm <- sub("^final_ordered_categorical_POM_",  "POM_", nm)
  nm <- sub("^final_no_horseshoe_", "noHS_", nm)
  nm
}

# Higher-is-better convention. RPS variants are lower-is-better, so we flip
# their sign for ranking and direction-aware fill.
METRIC_DIRECTION <- c(
  rps_mean_scaled   = -1,
  rps_median_scaled = -1,
  rpss_uniform      = +1,
  rpss_frequency    = +1,
  bacc              = +1,
  ppv               = +1
)
METRIC_ORDER <- names(METRIC_DIRECTION)

# ---- DISCOVER COMBOS ------------------------------------------------------
glob_pat <- file.path(opt$results_root, "*", "*", "inference_ppc",
                      "prediction_accuracy_metrics.csv")
all_paths <- Sys.glob(glob_pat)

if (length(all_paths) == 0)
  stop("No prediction_accuracy_metrics.csv files found under ", opt$results_root)

parse_path <- function(p) {
  rel <- sub(paste0("^", opt$results_root, "/"), "", p)
  parts <- strsplit(rel, "/", fixed = TRUE)[[1]]
  # parts: <variant>/<dataset>/inference_ppc/prediction_accuracy_metrics.csv
  if (length(parts) < 4) return(NULL)
  list(model_subdir = parts[1], dataset = parts[2], path = p)
}
combos <- Filter(Negate(is.null), lapply(all_paths, parse_path))

models_filter <- if (is.null(opt$model_subdirs)) NULL else trimws(strsplit(opt$model_subdirs, ",")[[1]])
ds_filter     <- if (is.null(opt$datasets))      NULL else trimws(strsplit(opt$datasets,      ",")[[1]])

if (!is.null(models_filter))
  combos <- Filter(function(c) c$model_subdir %in% models_filter, combos)
if (!is.null(ds_filter))
  combos <- Filter(function(c) c$dataset %in% ds_filter, combos)

if (length(combos) == 0)
  stop("No combos remain after applying --model_subdirs / --datasets filters")

# ---- LOAD METRICS ---------------------------------------------------------
skipped <- character()
rows <- list()

for (c in combos) {
  d <- tryCatch(read.csv(c$path), error = function(e) NULL)
  if (is.null(d) || nrow(d) != 1) {
    skipped <- c(skipped, sprintf("%s / %s (unreadable or wrong row count)",
                                  c$model_subdir, c$dataset))
    next
  }
  missing_cols <- setdiff(METRIC_ORDER, colnames(d))
  if (length(missing_cols) > 0) {
    skipped <- c(skipped, sprintf("%s / %s (missing columns: %s)",
                                  c$model_subdir, c$dataset,
                                  paste(missing_cols, collapse = ",")))
    next
  }
  rows[[length(rows) + 1L]] <- data.frame(
    model_subdir   = c$model_subdir,
    model_short    = short_label(c$model_subdir),
    dataset        = c$dataset,
    rps_mean_scaled   = d$rps_mean_scaled,
    rps_median_scaled = d$rps_median_scaled,
    rpss_uniform      = d$rpss_uniform,
    rpss_frequency    = d$rpss_frequency,
    bacc              = d$bacc,
    ppv               = d$ppv,
    stringsAsFactors  = FALSE
  )
}

if (length(rows) == 0)
  stop("All discovered files were unreadable; nothing to aggregate")

wide <- bind_rows(rows)

# Report combos that were expected (cross product of all observed
# variants/datasets) but missing.
all_models   <- sort(unique(wide$model_subdir))
all_datasets <- sort(unique(wide$dataset))
expected <- expand.grid(model_subdir = all_models, dataset = all_datasets,
                        stringsAsFactors = FALSE)
present <- paste(wide$model_subdir, wide$dataset, sep = "||")
expected_key <- paste(expected$model_subdir, expected$dataset, sep = "||")
missing_combos <- expected[!expected_key %in% present, , drop = FALSE]
if (nrow(missing_combos) > 0) {
  for (i in seq_len(nrow(missing_combos))) {
    skipped <- c(skipped, sprintf("%s / %s (no metrics file)",
                                  missing_combos$model_subdir[i],
                                  missing_combos$dataset[i]))
  }
}

# ---- WIDE TABLE -----------------------------------------------------------
wide <- wide %>%
  arrange(dataset, desc(rpss_frequency))
wide_path <- file.path(OUTDIR, "prediction_accuracy_table.csv")
write.csv(wide, wide_path, row.names = FALSE)

# ---- LONG TABLE -----------------------------------------------------------
long <- wide %>%
  pivot_longer(cols = all_of(METRIC_ORDER),
               names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric, levels = METRIC_ORDER))
long_path <- file.path(OUTDIR, "prediction_accuracy_long.csv")
write.csv(long, long_path, row.names = FALSE)

# ---- RANK TABLE -----------------------------------------------------------
# Within each (dataset, metric), rank variants such that 1 = best. Then take
# the mean rank per variant per metric across datasets, plus an overall mean
# rank averaged across all metrics.
rank_long <- long %>%
  group_by(dataset, metric) %>%
  mutate(
    signed_value = value * METRIC_DIRECTION[as.character(metric)],
    rank = rank(-signed_value, ties.method = "average", na.last = "keep")
  ) %>%
  ungroup()

rank_per_variant_metric <- rank_long %>%
  group_by(model_short, metric) %>%
  summarise(mean_rank = mean(rank, na.rm = TRUE),
            n_datasets = sum(!is.na(rank)),
            .groups = "drop")

overall_rank <- rank_per_variant_metric %>%
  group_by(model_short) %>%
  summarise(overall_mean_rank = mean(mean_rank, na.rm = TRUE), .groups = "drop")

rank_wide <- rank_per_variant_metric %>%
  select(model_short, metric, mean_rank) %>%
  pivot_wider(names_from = metric, values_from = mean_rank,
              names_prefix = "rank_") %>%
  left_join(overall_rank, by = "model_short") %>%
  arrange(overall_mean_rank)
rank_path <- file.path(OUTDIR, "prediction_accuracy_rank_table.csv")
write.csv(rank_wide, rank_path, row.names = FALSE)

# Order variants for plots by overall_mean_rank (best first).
variant_order <- rank_wide$model_short

# ---- PLOTS ----------------------------------------------------------------
plot_long <- long %>%
  mutate(model_short = factor(model_short, levels = variant_order))

# Heatmap: per-metric rescale of the signed value so that within each metric
# facet, green = best, red = worst.
heat_df <- plot_long %>%
  group_by(metric) %>%
  mutate(
    signed_value = value * METRIC_DIRECTION[as.character(metric)],
    fill_norm = scales::rescale(signed_value, to = c(0, 1))
  ) %>%
  ungroup()

p_heat <- ggplot(heat_df,
                 aes(x = model_short, y = dataset, fill = fill_norm)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.3f", value)), size = 2.7) +
  scale_fill_gradient(name = "rank within metric\n(green = best)",
                      low = "#e31a1c", high = "#1a9850",
                      na.value = "grey80",
                      limits = c(0, 1)) +
  facet_wrap(~ metric, ncol = 2, scales = "free") +
  labs(title = "Prediction accuracy across model variants and datasets",
       subtitle = "Variants ordered left-to-right by overall mean rank (best first). Cell text = raw value.",
       x = "model variant", y = "dataset") +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(OUTDIR, "prediction_accuracy_heatmap.png"),
       p_heat, width = 11, height = 7.5, dpi = 150)

# Lines: x = dataset, y = value, colour = variant, faceted by metric.
p_lines <- ggplot(plot_long,
                  aes(x = dataset, y = value,
                      colour = model_short, group = model_short)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.2) +
  facet_wrap(~ metric, scales = "free_y") +
  labs(title = "Prediction accuracy metrics by dataset and variant",
       subtitle = "rps_*: lower is better. rpss_*, bacc, ppv: higher is better.",
       x = NULL, y = "metric value", colour = "variant") +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(OUTDIR, "prediction_accuracy_lines.png"),
       p_lines, width = 11, height = 7, dpi = 150)

# Rank bars: per-variant mean rank, faceted by metric.
rank_plot_df <- rank_per_variant_metric %>%
  mutate(model_short = factor(model_short, levels = variant_order),
         metric = factor(metric, levels = METRIC_ORDER))
p_ranks <- ggplot(rank_plot_df,
                  aes(x = model_short, y = mean_rank, fill = model_short)) +
  geom_col() +
  geom_text(aes(label = sprintf("%.1f", mean_rank)),
            vjust = -0.3, size = 2.7) +
  facet_wrap(~ metric) +
  scale_y_reverse() +
  labs(title = "Mean rank per variant per metric (1 = best)",
       subtitle = "Averaged across datasets. Y axis reversed so taller bars = better.",
       x = NULL, y = "mean rank (1 = best)") +
  theme_bw(base_size = 10) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(OUTDIR, "prediction_accuracy_rank_bars.png"),
       p_ranks, width = 11, height = 7, dpi = 150)

# ---- LOG ------------------------------------------------------------------
message("\n=== prediction accuracy: wide table (sorted by dataset, rpss_frequency desc) ===")
print(as.data.frame(wide), row.names = FALSE, digits = 3)

message("\n=== rank table (overall_mean_rank: lower = better; 1 = best) ===")
print(as.data.frame(rank_wide), row.names = FALSE, digits = 3)

# Sanity checks on metric ranges.
oddities <- list()
neg_rpss <- wide %>% filter(rpss_uniform < 0 | rpss_frequency < 0)
if (nrow(neg_rpss) > 0)
  oddities[["negative_rpss"]] <- neg_rpss %>%
    select(model_subdir, dataset, rpss_uniform, rpss_frequency)
out_of_range <- wide %>% filter(bacc < 0 | bacc > 1 | ppv < 0 | ppv > 1)
if (nrow(out_of_range) > 0)
  oddities[["bacc_ppv_out_of_range"]] <- out_of_range %>%
    select(model_subdir, dataset, bacc, ppv)
if (length(oddities) > 0) {
  message("\nSANITY-CHECK FLAGS:")
  for (nm in names(oddities)) {
    message("  -- ", nm)
    print(oddities[[nm]], row.names = FALSE, digits = 3)
  }
}

message("\nWritten:")
message("  ", wide_path)
message("  ", long_path)
message("  ", rank_path)
message("  ", file.path(OUTDIR, "prediction_accuracy_heatmap.png"))
message("  ", file.path(OUTDIR, "prediction_accuracy_lines.png"))
message("  ", file.path(OUTDIR, "prediction_accuracy_rank_bars.png"))

if (length(skipped) > 0) {
  message("\nSKIPPED:")
  for (s in skipped) message("  - ", s)
}
