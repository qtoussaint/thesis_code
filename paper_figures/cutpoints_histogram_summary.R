#!/usr/bin/env Rscript
# Collect the per-run cutpoint posterior histograms into one figure per dataset
# (SPN penicillin, SPN trimethoprim, TB rifampicin).
#
# Every ordinal inference run writes plots/cutpoints/cutpoints_histogram.png --
# overlapping posterior densities of each latent-scale cutpoint, coloured by
# cutpoint with the MIC breakpoint in the legend. The pipeline writes only PNGs
# (no CSV), so we reuse those PNGs directly rather than recomputing from draws.
#
# Layout per dataset: two columns (POM left, PPOM right), one row per MIC binning.
# Each panel gets a caption naming its binning + model, since the source PNGs
# carry no run identity. Panels are drawn as images (never rescaled in aspect)
# via cowplot::draw_image -- magick text annotation is unusable on this host
# (ImageMagick built without fontconfig), so cowplot renders the labels instead.
#
# Usage:
#   mamba activate gwas_pipeline
#   Rscript paper_figures/cutpoints_histogram_summary.R
#
# Output: paper_figures/cutpoint_histograms/cutpoints_histogram_summary_<key>.png

suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
  library(magick)
})

# Null device so any incidental plotting doesn't leave a stray Rplots.pdf behind.
grDevices::pdf(NULL)
on.exit(grDevices::dev.off(), add = TRUE)

RESULTS_ROOT <- "/nfs/research/jlees/jacqueline/thesis_results"
OUTPUT_DIR   <- file.path(RESULTS_ROOT, "paper_figures", "cutpoint_histograms")

MODELS <- c("POM", "PPOM")

# -----------------------------------------------------------------------------
# Dataset specs: same per-dataset binning structure as heritability_summary.R
# (key, species_dir, binning_specs). Baselines are dropped -- logistic and
# continuous models have no cutpoints. Unitigs runs (17/18) are not listed, so
# they are excluded just as in heritability_summary.R.
# -----------------------------------------------------------------------------

datasets <- list(
  list(
    key = "spn_penicillin", species_dir = "spn_penicillin",
    binning_specs = list(
      list(binning = "doubling (≥5%)",        nn = "02", run_stub = "spn_penicillin_MIC"),
      list(binning = "4-fold (≥5%)", nn = "10", run_stub = "spn_penicillin_MIC_coarse_dilutions"),
      list(binning = "doubling (≥10%)",       nn = "11", run_stub = "spn_penicillin_MIC_large_minbin"),
      list(binning = "minima",                 nn = "16", run_stub = "spn_penicillin_MIC_minimabinning"))
  ),
  list(
    key = "spn_trimethoprim", species_dir = "spn_trimethoprim",
    binning_specs = list(
      list(binning = "doubling (≥5%)",        nn = "05", run_stub = "spn_trimethoprim_MIC"),
      list(binning = "4-fold (≥5%)", nn = "12", run_stub = "spn_trimethoprim_MIC_coarse_dilutions"),
      list(binning = "doubling (≥10%)",       nn = "13", run_stub = "spn_trimethoprim_MIC_large_minbin"))
  ),
  list(
    key = "tb_rifampicin", species_dir = "tb_rifampicin",
    binning_specs = list(
      list(binning = "doubling (≥5%)",        nn = "08", run_stub = "tb_rifampicin_MIC"),
      list(binning = "4-fold (≥5%)", nn = "14", run_stub = "tb_rifampicin_MIC_coarse_dilutions"),
      list(binning = "doubling (≥10%)",       nn = "15", run_stub = "tb_rifampicin_MIC_large_minbin"))
  )
)

# -----------------------------------------------------------------------------
# Per-dataset figure builder
# -----------------------------------------------------------------------------

hist_png_path <- function(species_dir, nn, run_stub, model) file.path(
  RESULTS_ROOT, paste0("gwas_", species_dir), "inference",
  paste0(nn, "_", run_stub, "_", model), "plots", "cutpoints",
  "cutpoints_histogram.png")

# One labelled grid cell: the histogram image with a caption above it. A missing
# run renders as an empty captioned cell so the grid layout stays aligned.
build_cell <- function(binning, model, path) {
  label <- sprintf("%s — %s", binning, model)
  cell <- ggdraw() +
    draw_label(label, x = 0.5, y = 0.97, vjust = 1, size = 13,
               fontface = "bold", colour = "grey20")
  if (file.exists(path)) {
    cell + draw_image(path, x = 0, y = 0, width = 1, height = 0.94)
  } else {
    cell +
      draw_label("*", x = 0.5, y = 0.47, size = 16, colour = "grey60")
  }
}

make_figure <- function(ds) {
  cells <- list()
  missing_files <- character(0)
  for (s in ds$binning_specs) {
    for (model in MODELS) {
      path <- hist_png_path(ds$species_dir, s$nn, s$run_stub, model)
      cells[[length(cells) + 1]] <- build_cell(s$binning, model, path)
      if (!file.exists(path)) missing_files <- c(missing_files, path)
    }
  }

  message(sprintf("[%s] ", ds$key), appendLF = FALSE)
  if (length(missing_files) > 0) {
    message(sprintf("%d cutpoint histogram PNG(s) missing (rendered as gaps):",
                    length(missing_files)))
    for (f in missing_files) message("  ", f)
  } else {
    message("all cutpoint histogram PNGs found.")
  }

  n_bin <- length(ds$binning_specs)
  figure <- plot_grid(plotlist = cells, ncol = length(MODELS), align = "hv")

  dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
  out_path <- file.path(OUTPUT_DIR,
                        paste0("cutpoints_histogram_summary_", ds$key, ".png"))
  # Each source histogram is 10x6in; size the grid to keep panels close to native.
  ggsave(out_path, figure,
         width = length(MODELS) * 7, height = n_bin * 4.4,
         dpi = 200, bg = "white", limitsize = FALSE)
  message("wrote ", out_path)
  invisible(NULL)
}

invisible(lapply(datasets, make_figure))
