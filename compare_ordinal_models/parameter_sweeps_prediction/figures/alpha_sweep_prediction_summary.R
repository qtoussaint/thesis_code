#!/usr/bin/env Rscript
# OUT-OF-SAMPLE intercept-prior-SD sweep figures for the logistic + continuous
# models.
#
# Prediction twin of paper_figures/alpha_sweep_summary.R: held-out test-split
# accuracy vs alpha_prior_sd, accuracy only (the in-sample figure's Panel B
# variant-sparsity panel is dropped). The production prediction model is drawn as
# a dashed reference line.
#
# Output: <fig_dir>/{logistic,continuous}_alpha_prediction_summary.{png,csv}
#
# Run with:
#   mamba activate gwas_pipeline
#   Rscript figures/alpha_sweep_prediction_summary.R   # after the runs land

suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
  library(viridis)
})

HERE <- dirname(sub("^--file=", "",
                    grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(HERE) || !nzchar(HERE)) HERE <- "."
source(file.path(HERE, "prediction_sweep_helpers.R"))

SWEEP_LABEL <- "swept α-prior SD"
REF_LABEL   <- "production prediction model"

sd_token <- function(sd) gsub("\\.", "p", as.character(sd))

# Per-family config: family key (matches the *_prediction_alphasd<token> run
# names), prediction dataset, accuracy metrics (csv column -> facet label),
# swept SDs, production prediction reference run, output stub.
MODELS <- list(
  list(
    key      = "logistic",
    dataset  = "01_spn_penicillin_binary",
    ref_run  = "01_spn_penicillin_binary_logistic_random",
    metrics  = c(bacc = "balanced accuracy", auc = "AUC", brier = "Brier score"),
    sds      = c(0.075, 1, 1.25, 1.5, 2, 2.5, 3, 5),
    out      = "logistic_alpha_prediction_summary"
  ),
  list(
    key      = "continuous",
    dataset  = "03_spn_penicillin_continuous",
    ref_run  = "03_spn_penicillin_continuous_continuous_random",
    metrics  = c(r_squared = "R²", rmse = "RMSE", crps = "CRPS"),
    sds      = c(1.5, 3, 5),
    out      = "continuous_alpha_prediction_summary"
  )
)

colour_cols   <- c(setNames(viridis::viridis(1, option = "D", end = 0.45), SWEEP_LABEL),
                   setNames("#D81B60", REF_LABEL))
colour_levels <- c(SWEEP_LABEL, REF_LABEL)

base_theme <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92", colour = NA),
        legend.position = "bottom")

build_family <- function(cfg) {
  metric_cols   <- names(cfg$metrics)
  metric_labels <- cfg$metrics

  rows <- list()
  for (sd in cfg$sds) {
    run <- paste0(cfg$key, "_prediction_alphasd", sd_token(sd))
    acc <- read_metrics(pred_metrics_csv(run, cfg$dataset), metric_cols)
    for (m in metric_cols)
      rows[[length(rows) + 1]] <- data.frame(sd = sd, quantity = m,
                                             value = unname(acc[[m]]),
                                             stringsAsFactors = FALSE)
  }
  df_long <- do.call(rbind, rows)
  df_long$quantity <- factor(df_long$quantity, levels = metric_cols)

  # Reference = production prediction run for this family.
  base_acc <- read_metrics(
    file.path(RESULTS_ROOT, "gwas_spn_penicillin", "prediction", cfg$ref_run,
              "prediction_results", "prediction_accuracy_metrics.csv"), metric_cols)
  ref_df <- data.frame(quantity = factor(metric_cols, levels = metric_cols),
                       value = unname(base_acc[metric_cols]), ref = REF_LABEL,
                       stringsAsFactors = FALSE)

  panel <- ggplot(df_long, aes(x = sd, y = value)) +
    geom_hline(data = ref_df, aes(yintercept = value, colour = ref),
               linetype = "dashed", linewidth = 0.7, inherit.aes = FALSE) +
    geom_line(aes(colour = SWEEP_LABEL), linewidth = 0.7) +
    geom_point(aes(colour = SWEEP_LABEL), size = 2.4) +
    facet_wrap(~ quantity, ncol = 1, scales = "free_y",
               labeller = labeller(quantity = metric_labels)) +
    scale_x_continuous(breaks = cfg$sds) +
    scale_colour_manual(values = colour_cols, name = NULL, limits = colour_levels) +
    labs(x = "intercept prior SD", y = "out-of-sample metric value") +
    base_theme

  dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
  png_path <- file.path(FIG_DIR, paste0(cfg$out, ".png"))
  csv_path <- file.path(FIG_DIR, paste0(cfg$out, ".csv"))
  ggsave(png_path, panel, width = 6.5, height = 8, dpi = 300, bg = "white")

  ref_out <- data.frame(sd = NA_real_, quantity = as.character(ref_df$quantity),
                        value = ref_df$value, stringsAsFactors = FALSE)
  out_df <- rbind(df_long, ref_out)
  write.csv(out_df[order(out_df$quantity, out_df$sd), ], csv_path, row.names = FALSE)

  message("wrote ", png_path)
  message("wrote ", csv_path)
}

for (cfg in MODELS) build_family(cfg)
report_missing()
