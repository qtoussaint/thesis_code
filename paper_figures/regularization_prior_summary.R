#!/usr/bin/env Rscript
# Regularisation-prior comparison figure for the PPOM free-cutpoints model.
#
# Unlike the tau sweep, these runs are NOT a continuous sweep -- they are a set
# of discrete prior variants for the variant-effect coefficients, so this is a
# categorical model comparison (a Cleveland dot plot), not a line sweep.
#
# Each model varies along three axes that are NOT a clean factorial, so we list
# every model explicitly rather than forcing a shared x-position:
#   penalty   : lasso (double-exponential) vs ridge (normal)
#   scale     : fixed global scale vs estimated scale ("est-sigma")
#   extra     : none; plus mixture (lasso only, a spike-and-slab prior) or
#               NCP (ridge only, a non-centred reparametrisation -- same model,
#               different geometry)
#   centering : variant genotypes centred vs not ("no-center")
#
# Panel A: out-of-sample prediction accuracy (RPSS uniform/frequency, balanced
#          accuracy). Panel B: variant sparsity (# unique significant variants).
# The horseshoe PPOM (estimated tau) fit is drawn as a reference line in every
# facet -- the regularisation priors are alternatives to that horseshoe.
#
# Reads the same per-run files as the prediction pipeline:
#   inference_ppc/prediction_accuracy_metrics.csv   (RPSS, balanced accuracy)
#   fitted_model/depruned_variant_effects.csv       (per-variant signif flags)
#
# Run with:
#   mamba activate gwas_pipeline
#   Rscript paper_figures/regularization_prior_summary.R
#
# Output: <output_dir>/regularization_prior_summary.{png,csv}

suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
})

RESULTS_ROOT <- "/nfs/research/jlees/jacqueline/thesis_results"
SWEEP_ROOT   <- file.path(RESULTS_ROOT, "compare_ordinal_models")
OUTPUT_DIR   <- file.path(RESULTS_ROOT, "paper_figures", "parameter_sweeps")
MODEL_PREFIX <- "final_ordered_categorical_PPOM_free_cutpoints_wide_drift"

# -----------------------------------------------------------------------------
# Model spec
# -----------------------------------------------------------------------------

# Each model -> dir suffix (after MODEL_PREFIX_), penalty, centering, and a
# compact label. Listed lasso-first; within each penalty, base -> est-sigma ->
# +extra, each as a centred/no-center pair. The y-axis preserves this order.
MODELS <- list(
  list(suffix = "lasso",                                 penalty = "lasso", centering = "centred",    label = "lasso"),
  list(suffix = "lasso_no_centering",                    penalty = "lasso", centering = "no-center",  label = "lasso · no-center"),
  list(suffix = "lasso_estscale",                        penalty = "lasso", centering = "centred",    label = "lasso · est-σ"),
  list(suffix = "lasso_estscale_no_centering",           penalty = "lasso", centering = "no-center",  label = "lasso · est-σ · no-center"),
  list(suffix = "lasso_estscale_mixture",                penalty = "lasso", centering = "centred",    label = "lasso · est-σ · mixture"),
  list(suffix = "lasso_estscale_mixture_no_centering",   penalty = "lasso", centering = "no-center",  label = "lasso · est-σ · mixture · no-center"),
  list(suffix = "ridge",                                 penalty = "ridge", centering = "centred",    label = "ridge"),
  list(suffix = "ridge_no_centering",                    penalty = "ridge", centering = "no-center",  label = "ridge · no-center"),
  list(suffix = "ridge_estscale",                        penalty = "ridge", centering = "centred",    label = "ridge · est-σ"),
  list(suffix = "ridge_estscale_no_centering",           penalty = "ridge", centering = "no-center",  label = "ridge · est-σ · no-center"),
  list(suffix = "ridge_estscale_ncp",                    penalty = "ridge", centering = "centred",    label = "ridge · est-σ · NCP"),
  list(suffix = "ridge_estscale_ncp_no_centering",       penalty = "ridge", centering = "no-center",  label = "ridge · est-σ · NCP · no-center")
)
# Factor levels: rev() so the first model listed sits at the TOP of the y-axis.
model_levels    <- rev(vapply(MODELS, `[[`, "", "label"))
penalty_levels  <- c("lasso", "ridge")
center_levels   <- c("centred", "no-center")

DATASETS <- list(
  list(id = "02_spn_penicillin_MIC",               label = "doubling (≥5%), K=8", K = 8L),
  list(id = "16_spn_penicillin_MIC_minimabinning", label = "minima, K=5",         K = 5L)
)
dataset_levels <- vapply(DATASETS, `[[`, "", "label")

# Reference baseline: the horseshoe PPOM (estimated tau) fit per dataset, drawn
# as a reference line in every facet.
REF_LABEL <- "horseshoe PPOM (estimated τ)"
REF_RUNS <- list(
  list(dataset = "doubling (≥5%), K=8",
       dir = file.path(RESULTS_ROOT, "gwas_spn_penicillin", "inference",
                       "02_spn_penicillin_MIC_PPOM")),
  list(dataset = "minima, K=5",
       dir = file.path(RESULTS_ROOT, "gwas_spn_penicillin", "inference",
                       "16_spn_penicillin_MIC_minimabinning_PPOM"))
)

# -----------------------------------------------------------------------------
# Path builders
# -----------------------------------------------------------------------------

model_dir <- function(suffix, ds_id)
  file.path(SWEEP_ROOT, paste0(MODEL_PREFIX, "_", suffix), ds_id)

metrics_path <- function(suffix, ds_id)
  file.path(model_dir(suffix, ds_id), "inference_ppc", "prediction_accuracy_metrics.csv")

effects_path <- function(suffix, ds_id)
  file.path(model_dir(suffix, ds_id), "fitted_model", "depruned_variant_effects.csv")

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
# Build tidy long data frame (12 models x 2 datasets)
# -----------------------------------------------------------------------------

ACC_METRICS <- c("rpss_uniform", "rpss_frequency", "bacc")
rows <- list()

for (mdl in MODELS) {
  for (d in DATASETS) {
    acc <- read_metrics(metrics_path(mdl$suffix, d$id), ACC_METRICS)
    spa <- count_signif(effects_path(mdl$suffix, d$id))

    common <- data.frame(model = mdl$label, penalty = mdl$penalty,
                         centering = mdl$centering, dataset = d$label,
                         stringsAsFactors = FALSE)
    for (m in ACC_METRICS)
      rows[[length(rows) + 1]] <- cbind(common,
        quantity = m, value = unname(acc[[m]]), panel = "A")
    rows[[length(rows) + 1]] <- cbind(common,
      quantity = "n_unique_variants", value = unname(spa[["n_unique"]]), panel = "B")
  }
}

df <- do.call(rbind, rows)
df$model     <- factor(df$model,     levels = model_levels)
df$penalty   <- factor(df$penalty,   levels = penalty_levels)
df$centering <- factor(df$centering, levels = center_levels)
df$dataset   <- factor(df$dataset,   levels = dataset_levels)

# Reference baseline values: one row per (dataset, quantity).
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

# Penalty -> colour (Okabe-Ito orange / blue); centering -> point shape.
penalty_cols <- c(lasso = "#E69F00", ridge = "#0072B2")
center_shape <- c("centred" = 16, "no-center" = 1)

base_theme <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(colour = "grey90"),
        strip.background = element_rect(fill = "grey92", colour = NA),
        legend.position = "bottom", legend.box = "vertical",
        legend.margin = margin(0, 0, 0, 0))

metric_labels <- c(rpss_uniform   = "RPSS (uniform)",
                   rpss_frequency = "RPSS (frequency)",
                   bacc           = "balanced accuracy")

# -----------------------------------------------------------------------------
# Panel A: prediction accuracy per model
# -----------------------------------------------------------------------------

dfA <- df[df$panel == "A", ]
dfA$quantity <- factor(dfA$quantity, levels = ACC_METRICS)

refA <- ref_df[ref_df$panel == "A", ]
refA$quantity <- factor(refA$quantity, levels = ACC_METRICS)

# The reference line is mapped to linetype (not colour) so it gets its own
# legend key without polluting the two-entry penalty colour scale.
panel_a <- ggplot(dfA, aes(x = value, y = model)) +
  geom_vline(data = refA, aes(xintercept = value, linetype = ref),
             colour = "#555555", linewidth = 0.6, inherit.aes = FALSE) +
  geom_point(aes(colour = penalty, shape = centering), size = 2.4, stroke = 0.9) +
  facet_grid(dataset ~ quantity, scales = "free_x",
             labeller = labeller(quantity = metric_labels)) +
  scale_colour_manual(values = penalty_cols, breaks = penalty_levels, name = "penalty") +
  scale_shape_manual(values = center_shape, name = "centering") +
  scale_linetype_manual(values = c("dashed"), labels = REF_LABEL, name = NULL) +
  guides(colour = guide_legend(order = 1, override.aes = list(shape = 16)),
         shape  = guide_legend(order = 2),
         linetype = guide_legend(order = 3)) +
  labs(x = "metric value", y = NULL) +
  base_theme

# -----------------------------------------------------------------------------
# Panel B: variant sparsity per model (unique significant variants)
# -----------------------------------------------------------------------------

dfB  <- df[df$panel == "B" & df$quantity == "n_unique_variants", ]
refB <- ref_df[ref_df$panel == "B", ]

panel_b <- ggplot(dfB, aes(x = value, y = model)) +
  geom_vline(data = refB, aes(xintercept = value), colour = "#555555",
             linetype = "dashed", linewidth = 0.6, inherit.aes = FALSE) +
  geom_point(aes(colour = penalty, shape = centering), size = 2.4, stroke = 0.9) +
  facet_wrap(~ dataset, nrow = 1, scales = "free_x") +
  scale_colour_manual(values = penalty_cols, breaks = penalty_levels, name = "penalty") +
  scale_shape_manual(values = center_shape, name = "centering") +
  labs(x = "# significant variants (unique)", y = NULL) +
  base_theme + theme(legend.position = "none")

# -----------------------------------------------------------------------------
# Assemble + save
# -----------------------------------------------------------------------------

legend <- cowplot::get_legend(panel_a)
panels <- cowplot::plot_grid(
  panel_a + theme(legend.position = "none"),
  panel_b,
  ncol = 1, labels = c("A", "B"), label_size = 24, label_fontface = "bold",
  rel_heights = c(2.0, 1.2))
figure <- cowplot::plot_grid(panels, legend, ncol = 1, rel_heights = c(1, 0.14))

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
png_path <- file.path(OUTPUT_DIR, "regularization_prior_summary.png")
csv_path <- file.path(OUTPUT_DIR, "regularization_prior_summary.csv")

ggsave(png_path, figure, width = 11, height = 10, dpi = 300, bg = "white")

# Combine model + reference rows for the saved table.
ref_out <- data.frame(model = REF_LABEL, penalty = NA_character_,
                      centering = NA_character_, dataset = ref_df$dataset,
                      quantity = ref_df$quantity, value = ref_df$value,
                      panel = ref_df$panel, stringsAsFactors = FALSE)
out_df <- rbind(df, ref_out)
write.csv(out_df[order(out_df$panel, out_df$quantity, out_df$dataset, out_df$model), ],
          csv_path, row.names = FALSE)

message("wrote ", png_path)
message("wrote ", csv_path)
