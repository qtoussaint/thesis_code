#!/usr/bin/env Rscript
# OUT-OF-SAMPLE tau-sweep summary for the PPOM shrinkage-prior grid.
#
# Prediction twin of paper_figures/tau_sweep_summary.R: held-out test-split
# accuracy (RPSS uniform/frequency, balanced accuracy) vs the global shrinkage
# scale tau, across the fixed-tau grid and two MIC binnings. Accuracy only (the
# in-sample figure's Panel B variant-sparsity panel is dropped). The production
# prediction PPOM (estimated tau) is a dashed reference line per facet; bACC
# cells recomputed over present-only categories are starred.
#
# Run with:
#   mamba activate gwas_pipeline
#   Rscript figures/tau_sweep_prediction_summary.R   # after the runs land

suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
  library(viridis)
})

HERE <- dirname(sub("^--file=", "",
                    grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(HERE) || !nzchar(HERE)) HERE <- "."
source(file.path(HERE, "prediction_sweep_helpers.R"))

MODEL_PREFIX <- "final_ordered_categorical_PPOM_free_cutpoints_wide_drift"

TAUS <- c(0.001, 0.01, 0.05, 1, 1.5, 2, 3, 5)

tau_token <- function(tau)
  switch(as.character(tau),
         "0.001" = "0p001", "0.01" = "0p01", "0.05" = "0p05", "1" = "1",
         "1.5" = "1p5", "2" = "2", "3" = "3", "5" = "5",
         stop("unmapped tau: ", tau))

SERIES <- list(
  list(label = "slab3/λ1", suffix = "_slab3"),
  list(label = "slab5/λ1", suffix = "_slab5"),
  list(label = "slab3/λ2", suffix = "_slab3_lambda2"),
  list(label = "slab5/λ2", suffix = "_slab5_lambda2")
)
series_levels <- vapply(SERIES, `[[`, "", "label")

DATASETS <- list(
  list(id = "02_spn_penicillin_MIC",               label = "standard (K=8)", K = 8L),
  list(id = "16_spn_penicillin_MIC_minimabinning", label = "minima (K=5)",   K = 5L)
)
dataset_levels <- vapply(DATASETS, `[[`, "", "label")

REF_LABEL <- "prediction PPOM (estimated τ)"
REF_RUNS  <- list(
  list(dataset = "standard (K=8)", run = "02_spn_penicillin_MIC_PPOM_random"),
  list(dataset = "minima (K=5)",   run = "16_spn_penicillin_MIC_minimabinning_PPOM_random")
)
ACC_METRICS <- c("rpss_uniform", "rpss_frequency", "bacc")

run_name_for <- function(suffix, tau)
  paste0(MODEL_PREFIX, "_fixedtau", tau_token(tau), suffix, "_prediction")

# -----------------------------------------------------------------------------
# Build long data frame
# -----------------------------------------------------------------------------
rows <- list()
for (s in SERIES) {
  for (tau in TAUS) {
    for (d in DATASETS) {
      run <- run_name_for(s$suffix, tau)
      rp  <- read_metrics(pred_metrics_csv(run, d$id), c("rpss_uniform", "rpss_frequency"))
      bc  <- resolve_bacc(run, d$id, d$K)
      common <- data.frame(series = s$label, tau = tau, dataset = d$label,
                           stringsAsFactors = FALSE)
      rows[[length(rows) + 1]] <- cbind(common, quantity = "rpss_uniform",
                                        value = unname(rp[["rpss_uniform"]]), restricted = FALSE)
      rows[[length(rows) + 1]] <- cbind(common, quantity = "rpss_frequency",
                                        value = unname(rp[["rpss_frequency"]]), restricted = FALSE)
      rows[[length(rows) + 1]] <- cbind(common, quantity = "bacc",
                                        value = bc$value, restricted = bc$restricted)
    }
  }
}
df_long <- do.call(rbind, rows)
df_long$series   <- factor(df_long$series,  levels = series_levels)
df_long$dataset  <- factor(df_long$dataset, levels = dataset_levels)
df_long$quantity <- factor(df_long$quantity, levels = ACC_METRICS)

# Reference rows (one per dataset x metric), tau-independent.
ref_rows <- list()
for (r in REF_RUNS) {
  m <- read_metrics(file.path(RESULTS_ROOT, "gwas_spn_penicillin", "prediction",
                              r$run, "prediction_results", "prediction_accuracy_metrics.csv"),
                    ACC_METRICS)
  for (q in ACC_METRICS)
    ref_rows[[length(ref_rows) + 1]] <- data.frame(
      dataset = r$dataset, quantity = q, value = unname(m[[q]]), stringsAsFactors = FALSE)
}
ref_df <- do.call(rbind, ref_rows)
ref_df$dataset  <- factor(ref_df$dataset, levels = dataset_levels)
ref_df$quantity <- factor(ref_df$quantity, levels = ACC_METRICS)
ref_df$ref      <- REF_LABEL

report_missing()

# -----------------------------------------------------------------------------
# Plot (accuracy only)
# -----------------------------------------------------------------------------
series_cols   <- viridis::viridis(length(series_levels), option = "D", end = 0.9)
names(series_cols) <- series_levels
colour_cols   <- c(series_cols, setNames("#D81B60", REF_LABEL))
colour_levels <- c(series_levels, REF_LABEL)

base_theme <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92", colour = NA),
        legend.position = "bottom",
        plot.caption = element_text(hjust = 0, size = 9))

metric_labels <- c(rpss_uniform = "RPSS (uniform)", rpss_frequency = "RPSS (frequency)",
                   bacc = "balanced accuracy")
tau_labels <- c("0.001", "0.01", "0.05", "1", "1.5", "2", "3", "5")

panel <- ggplot(df_long, aes(x = tau, y = value, colour = series, group = series)) +
  geom_hline(data = ref_df, aes(yintercept = value, colour = ref),
             linetype = "dashed", linewidth = 0.7, inherit.aes = FALSE) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.4) +
  geom_text(data = df_long[df_long$quantity == "bacc" & df_long$restricted %in% TRUE, ],
            aes(label = "*"), nudge_y = 0.03, size = 5, show.legend = FALSE) +
  facet_grid(quantity ~ dataset, scales = "free_y",
             labeller = labeller(quantity = metric_labels)) +
  scale_x_log10(breaks = TAUS, labels = tau_labels) +
  scale_colour_manual(values = colour_cols, name = "prior config", limits = colour_levels) +
  coord_cartesian(ylim = c(0.5, 1)) +
  labs(x = expression("global shrinkage scale " * tau), y = "out-of-sample metric value") +
  base_theme

dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
png_path <- file.path(FIG_DIR, "tau_sweep_prediction_summary.png")
csv_path <- file.path(FIG_DIR, "tau_sweep_prediction_summary.csv")

ggsave(png_path, panel, width = 9, height = 9, dpi = 300, bg = "white")

ref_out <- data.frame(series = REF_LABEL, tau = NA_real_, dataset = ref_df$dataset,
                      quantity = ref_df$quantity, value = ref_df$value,
                      restricted = FALSE, stringsAsFactors = FALSE)
out_df <- rbind(df_long, ref_out)
write.csv(out_df[order(out_df$quantity, out_df$dataset, out_df$series, out_df$tau), ],
          csv_path, row.names = FALSE)

message("wrote ", png_path)
message("wrote ", csv_path)
