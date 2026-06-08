#!/usr/bin/env Rscript
# OUT-OF-SAMPLE regularisation-prior comparison for the PPOM free-cutpoints sweep.
#
# Prediction twin of paper_figures/regularization_prior_summary.R: reads the
# held-out test-split metrics produced by the parameter_sweeps_prediction runs
# instead of the in-sample inference_ppc/ metrics, and shows accuracy only (the
# in-sample figure's Panel B variant-sparsity panel is dropped).
#
# Panel: out-of-sample prediction accuracy (RPSS uniform/frequency, balanced
# accuracy) as a Cleveland dot plot, 12 regularisation-prior models x 2 datasets.
# The production prediction PPOM (horseshoe, estimated tau) is drawn as a dashed
# reference line per facet. bACC cells recomputed over present-only categories
# (div-by-zero guard) are starred.
#
# Run with:
#   mamba activate gwas_pipeline
#   Rscript figures/regularization_prior_prediction_summary.R   # after the runs land

suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
})

HERE <- dirname(sub("^--file=", "",
                    grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(HERE) || !nzchar(HERE)) HERE <- "."
source(file.path(HERE, "prediction_sweep_helpers.R"))

MODEL_PREFIX <- "final_ordered_categorical_PPOM_free_cutpoints_wide_drift"

# -----------------------------------------------------------------------------
# Model + dataset spec (same layout as the in-sample figure).
# -----------------------------------------------------------------------------
MODELS <- list(
  list(suffix = "lasso",                               penalty = "lasso", centering = "centred",   label = "lasso"),
  list(suffix = "lasso_no_centering",                  penalty = "lasso", centering = "no-center", label = "lasso · no-center"),
  list(suffix = "lasso_estscale",                      penalty = "lasso", centering = "centred",   label = "lasso · est-σ"),
  list(suffix = "lasso_estscale_no_centering",         penalty = "lasso", centering = "no-center", label = "lasso · est-σ · no-center"),
  list(suffix = "lasso_estscale_mixture",              penalty = "lasso", centering = "centred",   label = "lasso · est-σ · mixture"),
  list(suffix = "lasso_estscale_mixture_no_centering", penalty = "lasso", centering = "no-center", label = "lasso · est-σ · mixture · no-center"),
  list(suffix = "ridge",                               penalty = "ridge", centering = "centred",   label = "ridge"),
  list(suffix = "ridge_no_centering",                  penalty = "ridge", centering = "no-center", label = "ridge · no-center"),
  list(suffix = "ridge_estscale",                      penalty = "ridge", centering = "centred",   label = "ridge · est-σ"),
  list(suffix = "ridge_estscale_no_centering",         penalty = "ridge", centering = "no-center", label = "ridge · est-σ · no-center"),
  list(suffix = "ridge_estscale_ncp",                  penalty = "ridge", centering = "centred",   label = "ridge · est-σ · NCP"),
  list(suffix = "ridge_estscale_ncp_no_centering",     penalty = "ridge", centering = "no-center", label = "ridge · est-σ · NCP · no-center")
)
model_levels   <- rev(vapply(MODELS, `[[`, "", "label"))
penalty_levels <- c("lasso", "ridge")
center_levels  <- c("centred", "no-center")

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

run_name_for <- function(suffix) paste0(MODEL_PREFIX, "_", suffix, "_prediction")

# -----------------------------------------------------------------------------
# Build long data frame (rpss read directly; bacc resolved + restriction-flagged)
# -----------------------------------------------------------------------------
rows <- list()
for (mdl in MODELS) {
  run <- run_name_for(mdl$suffix)
  for (d in DATASETS) {
    rp <- read_metrics(pred_metrics_csv(run, d$id), c("rpss_uniform", "rpss_frequency"))
    bc <- resolve_bacc(run, d$id, d$K)
    common <- data.frame(model = mdl$label, penalty = mdl$penalty,
                         centering = mdl$centering, dataset = d$label,
                         stringsAsFactors = FALSE)
    rows[[length(rows) + 1]] <- cbind(common, quantity = "rpss_uniform",
                                      value = unname(rp[["rpss_uniform"]]), restricted = FALSE)
    rows[[length(rows) + 1]] <- cbind(common, quantity = "rpss_frequency",
                                      value = unname(rp[["rpss_frequency"]]), restricted = FALSE)
    rows[[length(rows) + 1]] <- cbind(common, quantity = "bacc",
                                      value = bc$value, restricted = bc$restricted)
  }
}
df <- do.call(rbind, rows)
df$model     <- factor(df$model,     levels = model_levels)
df$penalty   <- factor(df$penalty,   levels = penalty_levels)
df$centering <- factor(df$centering, levels = center_levels)
df$dataset   <- factor(df$dataset,   levels = dataset_levels)
df$quantity  <- factor(df$quantity,  levels = ACC_METRICS)

# Reference rows (one per dataset x metric); bacc read on-disk (production run).
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
# Plot
# -----------------------------------------------------------------------------
penalty_cols <- c(lasso = "#E69F00", ridge = "#0072B2")
center_shape <- c("centred" = 16, "no-center" = 1)
metric_labels <- c(rpss_uniform = "RPSS (uniform)", rpss_frequency = "RPSS (frequency)",
                   bacc = "balanced accuracy")

base_theme <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(colour = "grey90"),
        strip.background = element_rect(fill = "grey92", colour = NA),
        legend.position = "bottom", legend.box = "vertical",
        legend.margin = margin(0, 0, 0, 0),
        plot.caption = element_text(hjust = 0, size = 9))

panel <- ggplot(df, aes(x = value, y = model)) +
  geom_vline(data = ref_df, aes(xintercept = value, linetype = ref),
             colour = "#555555", linewidth = 0.6, inherit.aes = FALSE) +
  geom_point(aes(colour = penalty, shape = centering), size = 2.4, stroke = 0.9) +
  geom_text(data = df[df$quantity == "bacc" & df$restricted %in% TRUE, ],
            aes(label = "*"), nudge_x = 0.02, size = 5, show.legend = FALSE) +
  facet_grid(dataset ~ quantity, scales = "free_x",
             labeller = labeller(quantity = metric_labels)) +
  scale_colour_manual(values = penalty_cols, breaks = penalty_levels, name = "penalty") +
  scale_shape_manual(values = center_shape, name = "centering") +
  scale_linetype_manual(values = c("dashed"), labels = REF_LABEL, name = NULL) +
  guides(colour = guide_legend(order = 1, override.aes = list(shape = 16)),
         shape  = guide_legend(order = 2),
         linetype = guide_legend(order = 3)) +
  labs(x = "out-of-sample metric value", y = NULL,
       caption = "* bACC computed on a reduced category set (categories absent from the test split excluded)") +
  base_theme

dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
png_path <- file.path(FIG_DIR, "regularization_prior_prediction_summary.png")
csv_path <- file.path(FIG_DIR, "regularization_prior_prediction_summary.csv")

ggsave(png_path, panel, width = 11, height = 7, dpi = 300, bg = "white")

ref_out <- data.frame(model = REF_LABEL, penalty = NA_character_, centering = NA_character_,
                      dataset = ref_df$dataset, quantity = ref_df$quantity,
                      value = ref_df$value, restricted = FALSE, stringsAsFactors = FALSE)
out_df <- rbind(df, ref_out)
write.csv(out_df[order(out_df$quantity, out_df$dataset, out_df$model), ],
          csv_path, row.names = FALSE)

message("wrote ", png_path)
message("wrote ", csv_path)
