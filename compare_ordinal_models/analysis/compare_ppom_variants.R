#!/usr/bin/env Rscript
# Side-by-side comparison of PPOM variants on a single dataset.
#
# Consumes the per-model outputs of diagnose_ppom_cutpoint_inflation.R
# (summary_by_k.csv, metrics_summary.csv, ref_sublineage_cdf.csv) and emits
# a stacked visualisation so you can tell at a glance which variant collapses
# the late-k inflation.
#
# Usage (activate gwas_pipeline conda env first):
#   # Run diagnose_ppom_cutpoint_inflation.R on every variant first, then:
#   Rscript compare_ppom_variants.R \
#     --dataset 02_spn_penicillin_MIC \
#     --model_subdirs final_ordered_categorical_PPOM_tight_alpha_tau1,final_ordered_categorical_PPOM_tight_slab,final_ordered_categorical_PPOM_poolk,final_ordered_categorical_PPOM_latent_scale
#
# Outputs (to <results_root>/_compare_PPOM/<dataset>/):
#   abs_beta_by_k.png         median |beta| vs k, line per variant
#   post_sd_by_k.png          posterior SD vs k, line per variant
#   ref_cdf_by_model.png      reference-sublineage model-implied CDF, line per variant
#   inflation_ratio_table.csv inflation ratio (late/early median |beta|) per variant
#
# A well-behaved variant collapses the |beta|-by-k curve toward flat and its
# inflation ratio toward 1. A variant that still climbs means the mechanism
# targeted by that variant was not the dominant one.

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

option_list <- list(
  make_option("--dataset", type = "character", default = NULL,
              help = "Dataset nickname (must match the --dataset used in diagnose_*.R)"),
  make_option("--model_subdirs", type = "character", default = NULL,
              help = "Comma-separated model subdir names (under results_root)"),
  make_option("--results_root", type = "character",
              default = "/nfs/research/jlees/jacqueline/thesis_results/compare_ordinal_models",
              help = "Parent dir holding per-model result subdirs"),
  make_option("--out_subdir", type = "character", default = "_compare_PPOM",
              help = "Comparison output subdir under results_root [_compare_PPOM]")
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$dataset) || is.null(opt$model_subdirs))
  stop("Both --dataset and --model_subdirs are required")

models <- trimws(strsplit(opt$model_subdirs, ",")[[1]])
OUTDIR <- file.path(opt$results_root, opt$out_subdir, opt$dataset)
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

short_label <- function(nm) {
  nm <- sub("^final_ordered_categorical_PPOM_", "", nm)
  nm <- sub("^final_no_horseshoe_", "noHS_", nm)
  nm
}

# ---- LOAD PER-MODEL OUTPUTS ----------------------------------------------
load_one <- function(m) {
  dir <- file.path(opt$results_root, m, opt$dataset, "ppc", "diagnose")
  out <- list()
  sf <- file.path(dir, "summary_by_k.csv")
  if (file.exists(sf)) out$summary <- transform(read.csv(sf), model = m)
  mf <- file.path(dir, "metrics_summary.csv")
  if (file.exists(mf)) out$metrics <- read.csv(mf)
  rf <- file.path(dir, "ref_sublineage_cdf.csv")
  if (file.exists(rf)) out$ref_cdf <- transform(read.csv(rf), model = m)
  out
}

per_model <- lapply(models, load_one)
names(per_model) <- models

missing <- models[sapply(per_model, function(x) is.null(x$summary))]
if (length(missing) > 0) {
  warning("diagnose_ppom_cutpoint_inflation.R outputs not found for: ",
          paste(missing, collapse = ", "),
          " -- skipping these models for dataset ", opt$dataset)
  message("SKIPPED (no diagnostic outputs) for dataset '", opt$dataset, "': ",
          paste(missing, collapse = ", "))
  per_model <- per_model[setdiff(models, missing)]
  models    <- setdiff(models, missing)
  if (length(models) == 0)
    stop("No models with diagnostic outputs found for dataset ", opt$dataset)
}

summary_df <- bind_rows(lapply(per_model, `[[`, "summary"))
metrics_df <- bind_rows(lapply(per_model, `[[`, "metrics"))
ref_df     <- bind_rows(lapply(per_model, `[[`, "ref_cdf"))

summary_df$model_short <- short_label(summary_df$model)
metrics_df$model_short <- short_label(metrics_df$model_subdir)
if (nrow(ref_df) > 0) ref_df$model_short <- short_label(ref_df$model)

# ---- PLOT 1: median |beta| vs k, line per model --------------------------
p_abs <- ggplot(summary_df, aes(x = k, y = median_abs_beta,
                                colour = model_short, group = model_short)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2.2) +
  scale_x_continuous(breaks = sort(unique(summary_df$k))) +
  labs(title = sprintf("median |beta_variant| by cutpoint -- %s", opt$dataset),
       subtitle = "Flat line => no inflation. Rising line => late-k inflation.",
       x = "cutpoint index k", y = "median |beta_variant|",
       colour = "model") +
  theme_bw(base_size = 10) + theme(legend.position = "bottom")
ggsave(file.path(OUTDIR, "abs_beta_by_k.png"),
       p_abs, width = 8, height = 5, dpi = 150)

# ---- PLOT 2: posterior SD vs k, line per model ---------------------------
p_sd <- ggplot(summary_df, aes(x = k, y = median_post_sd,
                               colour = model_short, group = model_short)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2.2) +
  scale_x_continuous(breaks = sort(unique(summary_df$k))) +
  labs(title = sprintf("median posterior SD(beta_variant) by cutpoint -- %s",
                       opt$dataset),
       subtitle = "Rising SD with k indicates weakly-identified late-k effects",
       x = "cutpoint index k", y = "median posterior sd(beta_variant)",
       colour = "model") +
  theme_bw(base_size = 10) + theme(legend.position = "bottom")
ggsave(file.path(OUTDIR, "post_sd_by_k.png"),
       p_sd, width = 8, height = 5, dpi = 150)

# ---- PLOT 3: reference-sublineage CDF by model ---------------------------
if (nrow(ref_df) > 0) {
  # All models share the same empirical curve (dataset fixed) -- take from the
  # first model and plot it once.
  emp_df <- ref_df %>% group_by(k) %>% summarise(P_emp_ref = first(P_emp_ref),
                                                  cutpoint  = first(cutpoint),
                                                  .groups = "drop")
  p_ref <- ggplot() +
    geom_ribbon(data = ref_df,
                aes(x = k, ymin = P_mod_lo, ymax = P_mod_hi,
                    fill = model_short, group = model_short),
                alpha = 0.1) +
    geom_line(data = ref_df,
              aes(x = k, y = P_mod_ref, colour = model_short,
                  group = model_short), linewidth = 0.8) +
    geom_line(data = emp_df, aes(x = k, y = P_emp_ref),
              colour = "black", linewidth = 1.2, linetype = "dashed") +
    geom_point(data = emp_df, aes(x = k, y = P_emp_ref),
               colour = "black", size = 2.5, shape = 17) +
    scale_x_continuous(breaks = sort(unique(ref_df$k))) +
    labs(title = sprintf("Reference-sublineage CDF by model -- %s",
                         opt$dataset),
         subtitle = "Dashed black = empirical; coloured = model-implied (alpha only). Biggest gap = most mechanism-(3) exposure",
         x = "cutpoint index k", y = "P(phen <= k | reference)",
         colour = "model", fill = "model") +
    theme_bw(base_size = 10) + theme(legend.position = "bottom")
  ggsave(file.path(OUTDIR, "ref_cdf_by_model.png"),
         p_ref, width = 8, height = 5, dpi = 150)
}

# ---- TABLE: inflation metric per variant ---------------------------------
inflation_tbl <- metrics_df %>%
  select(model_short, N, V, K, n_ref_samples,
         median_abs_beta_k1, median_abs_beta_kKminus1,
         max_abs_beta_k1, max_abs_beta_kKminus1,
         inflation_ratio_median_abs,
         median_post_sd_k1, median_post_sd_kKminus1,
         frac_above_k1, frac_above_kKminus1) %>%
  arrange(inflation_ratio_median_abs)
write.csv(inflation_tbl, file.path(OUTDIR, "inflation_ratio_table.csv"),
          row.names = FALSE)

message("\n=== inflation ratio (late k / early k median |beta|) ===")
print(inflation_tbl %>%
        select(model_short, median_abs_beta_k1, median_abs_beta_kKminus1,
               inflation_ratio_median_abs),
      row.names = FALSE, digits = 3)
message("\ncomparison outputs written to: ", OUTDIR)
