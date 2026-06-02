#!/usr/bin/env Rscript
# Summary figure comparing prediction accuracy of the ordinal MIC models (POM,
# PPOM) against the binary logistic baseline, across all drugs, MIC binnings and
# evaluation splits (random + lineage/LOSO).
#
# Reads the same per-run prediction_accuracy_metrics.csv files that
# gwas_finalruns/collect_results.R turns into LaTeX tables. Only `bacc` is shared
# between the ordinal and binary models, so balanced accuracy is the common axis;
# RPSS (skill over uniform / frequency baselines) is shown for POM vs PPOM only.
#
# Run with:
#   mamba activate gwas_pipeline
#   Rscript paper_figures/prediction_accuracy_summary.R
#
# Output: <output_dir>/prediction_accuracy_summary.{png,csv}

suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
  library(viridis)
})

RESULTS_ROOT <- "/nfs/research/jlees/jacqueline/thesis_results"
OUTPUT_DIR   <- file.path(RESULTS_ROOT, "paper_figures")

# -----------------------------------------------------------------------------
# Path builder + reader (mirrors collect_results.R, but returns numerics)
# -----------------------------------------------------------------------------

pred_csv <- function(species_dir, run_dir, split) file.path(
  RESULTS_ROOT, paste0("gwas_", species_dir), "prediction",
  paste0(run_dir, "_", split), "prediction_results", "prediction_accuracy_metrics.csv")

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

# -----------------------------------------------------------------------------
# Specs (drug / binning / model layout, from collect_results.R)
# -----------------------------------------------------------------------------

# Ordinal MIC runs: one entry per (drug, binning). Each fits both PPOM and POM.
ordinal_specs <- list(
  list(species_dir="spn_penicillin",   drug="SPN PEN", organism="S. pneumoniae",   binning="standard",     K=8, nn="02", run_stub="spn_penicillin_MIC"),
  list(species_dir="spn_penicillin",   drug="SPN PEN", organism="S. pneumoniae",   binning="coarse",       K=5, nn="10", run_stub="spn_penicillin_MIC_coarse_dilutions"),
  list(species_dir="spn_penicillin",   drug="SPN PEN", organism="S. pneumoniae",   binning="large minbin", K=4, nn="11", run_stub="spn_penicillin_MIC_large_minbin"),
  list(species_dir="spn_penicillin",   drug="SPN PEN", organism="S. pneumoniae",   binning="minima",       K=5, nn="16", run_stub="spn_penicillin_MIC_minimabinning"),
  list(species_dir="spn_trimethoprim", drug="SPN TMP", organism="S. pneumoniae",   binning="standard",     K=5, nn="05", run_stub="spn_trimethoprim_MIC"),
  list(species_dir="spn_trimethoprim", drug="SPN TMP", organism="S. pneumoniae",   binning="coarse",       K=3, nn="12", run_stub="spn_trimethoprim_MIC_coarse_dilutions"),
  list(species_dir="spn_trimethoprim", drug="SPN TMP", organism="S. pneumoniae",   binning="large minbin", K=3, nn="13", run_stub="spn_trimethoprim_MIC_large_minbin"),
  list(species_dir="tb_rifampicin",    drug="TB RIF",  organism="M. tuberculosis", binning="standard",     K=5, nn="08", run_stub="tb_rifampicin_MIC"),
  list(species_dir="tb_rifampicin",    drug="TB RIF",  organism="M. tuberculosis", binning="coarse",       K=4, nn="14", run_stub="tb_rifampicin_MIC_coarse_dilutions"),
  list(species_dir="tb_rifampicin",    drug="TB RIF",  organism="M. tuberculosis", binning="large minbin", K=4, nn="15", run_stub="tb_rifampicin_MIC_large_minbin")
)

# Binary logistic baseline: one per drug. run_stub is the on-disk dir prefix.
binary_specs <- list(
  list(species_dir="spn_penicillin",   drug="SPN PEN", run_stub="01_spn_penicillin_binary_logistic"),
  list(species_dir="spn_trimethoprim", drug="SPN TMP", run_stub="04_spn_trimethoprim_binary_logistic"),
  list(species_dir="tb_rifampicin",    drug="TB RIF",  run_stub="07_tb_rifampicin_binary_logistic")
)

SPLITS <- list(random = "Random", loso = "Lineage")
ORD_METRICS <- c("bacc", "rpss_uniform", "rpss_frequency")

# -----------------------------------------------------------------------------
# Build long data frame
# -----------------------------------------------------------------------------

rows <- list()

for (s in ordinal_specs) {
  for (model in c("POM", "PPOM")) {
    run_dir <- paste0(s$nn, "_", s$run_stub, "_", model)
    for (sp in names(SPLITS)) {
      m <- read_metrics(pred_csv(s$species_dir, run_dir, sp), ORD_METRICS)
      for (metric in ORD_METRICS) {
        rows[[length(rows) + 1]] <- data.frame(
          drug = s$drug, organism = s$organism, binning = s$binning, K = s$K,
          model = model, split = SPLITS[[sp]], metric = metric,
          value = unname(m[[metric]]), stringsAsFactors = FALSE)
      }
    }
  }
}

for (s in binary_specs) {
  for (sp in names(SPLITS)) {
    m <- read_metrics(pred_csv(s$species_dir, s$run_stub, sp), "bacc")
    rows[[length(rows) + 1]] <- data.frame(
      drug = s$drug, organism = NA_character_, binning = NA_character_, K = 2L,
      model = "Logistic", split = SPLITS[[sp]], metric = "bacc",
      value = unname(m[["bacc"]]), stringsAsFactors = FALSE)
  }
}

df <- do.call(rbind, rows)

# Factor orders
drug_levels    <- c("SPN PEN", "SPN TMP", "TB RIF")
binning_levels <- c("standard", "coarse", "large minbin", "minima")
df$drug    <- factor(df$drug, levels = drug_levels)
df$binning <- factor(df$binning, levels = binning_levels)
df$model   <- factor(df$model, levels = c("POM", "PPOM", "Logistic"))
df$split   <- factor(df$split, levels = c("Random", "Lineage"))

if (length(missing_files) > 0) {
  message(sprintf("NOTE: %d metrics CSV(s) missing (rendered as gaps):", length(missing_files)))
  for (f in missing_files) message("  ", f)
} else {
  message("All metrics CSVs found.")
}

# -----------------------------------------------------------------------------
# Aesthetics
# -----------------------------------------------------------------------------

# Three-colour discrete viridis: POM, PPOM, Logistic
model_cols <- viridis::viridis(3, option = "D", end = 0.9)
names(model_cols) <- c("POM", "PPOM", "Logistic")

base_theme <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92", colour = NA),
        legend.position = "bottom")

# -----------------------------------------------------------------------------
# Panel A: balanced accuracy (POM/PPOM bars + logistic reference line)
# -----------------------------------------------------------------------------

bacc_bars <- df[df$metric == "bacc" & df$model %in% c("POM", "PPOM"), ]
bacc_bars$model <- droplevels(bacc_bars$model)

# Logistic is constant across binnings within a drug -> one hline per drug x split.
bacc_logit <- df[df$metric == "bacc" & df$model == "Logistic", ]

panel_a <- ggplot(bacc_bars, aes(x = binning, y = value, fill = model)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  geom_hline(data = bacc_logit,
             aes(yintercept = value, linetype = "Logistic baseline"),
             colour = model_cols[["Logistic"]], linewidth = 0.7) +
  facet_grid(split ~ drug, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = model_cols, name = "ordinal model") +
  scale_linetype_manual(values = c("Logistic baseline" = "dashed"), name = NULL) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = NULL, y = "balanced accuracy") +
  base_theme +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

# -----------------------------------------------------------------------------
# Panel B: RPSS skill scores (POM vs PPOM; uniform + frequency baselines)
# -----------------------------------------------------------------------------

rpss <- df[df$metric %in% c("rpss_uniform", "rpss_frequency"), ]
rpss$baseline <- factor(ifelse(rpss$metric == "rpss_uniform", "uniform", "frequency"),
                        levels = c("uniform", "frequency"))
rpss$model <- droplevels(rpss$model)

panel_b <- ggplot(rpss,
                  aes(x = binning, y = value, colour = model, shape = baseline,
                      group = interaction(model, baseline))) +
  geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.4) +
  geom_point(position = position_dodge(width = 0.6), size = 2.8) +
  facet_grid(split ~ drug, scales = "free_x", space = "free_x") +
  scale_colour_manual(values = model_cols, name = "ordinal model") +
  scale_shape_manual(values = c(uniform = 16, frequency = 17), name = "RPSS baseline") +
  labs(x = NULL, y = "RPSS (skill over baseline)") +
  base_theme +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

# -----------------------------------------------------------------------------
# Assemble + save
# -----------------------------------------------------------------------------

figure <- plot_grid(
  panel_a, panel_b,
  ncol = 1,
  labels = c("A", "B"),
  label_size = 24,
  label_fontface = "bold",
  rel_heights = c(1, 1)
) + theme(plot.margin = margin(5, 5, 5, 10))

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
png_path <- file.path(OUTPUT_DIR, "prediction_accuracy_summary.png")
csv_path <- file.path(OUTPUT_DIR, "prediction_accuracy_summary.csv")

ggsave(png_path, figure, width = 14, height = 11, dpi = 300, bg = "white")
write.csv(df[order(df$metric, df$drug, df$binning, df$model, df$split), ],
          csv_path, row.names = FALSE)

message("wrote ", png_path)
message("wrote ", csv_path)
