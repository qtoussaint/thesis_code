#!/usr/bin/env Rscript
# Compose one PPOM paper figure from a manhattan_short directory and a locus-zoom directory.
#
# Usage:
#   Rscript manhattan_locuszoom_composite.R <manhattan_dir> <locuszoom_dir> [--output-dir <dir>]
#
# Output: <output_dir>/figure_<analysis_name>.png
# where <analysis_name> is the parent directory name of <manhattan_dir>.

suppressPackageStartupMessages({
  library(cowplot)
  library(ggplot2)
  library(magick)
})

DEFAULT_OUTPUT_DIR <- "/nfs/research/jlees/jacqueline/thesis_results/paper_figures"

parse_args <- function(argv) {
  output_dir <- DEFAULT_OUTPUT_DIR
  positional <- character(0)
  i <- 1
  while (i <= length(argv)) {
    a <- argv[i]
    if (a == "--output-dir") {
      if (i == length(argv)) stop("--output-dir requires a value")
      output_dir <- argv[i + 1]
      i <- i + 2
    } else {
      positional <- c(positional, a)
      i <- i + 1
    }
  }
  if (length(positional) != 2) {
    stop("Usage: Rscript manhattan_locuszoom_composite.R <manhattan_dir> <locuszoom_dir> [--output-dir <dir>]")
  }
  list(
    manhattan_dir = positional[1],
    locuszoom_dir = positional[2],
    output_dir = output_dir
  )
}

require_file <- function(path) {
  if (!file.exists(path)) stop("Missing input file: ", path)
  path
}

panel_from_png <- function(path) {
  ggdraw() + draw_image(path)
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))

  manhattan_dir <- normalizePath(args$manhattan_dir, mustWork = TRUE)
  locuszoom_dir <- normalizePath(args$locuszoom_dir, mustWork = TRUE)

  panel_a_path <- require_file(file.path(manhattan_dir, "manhattan_all_cutpoints_overlayed_median_effects.png"))
  panel_c_path <- require_file(file.path(manhattan_dir, "manhattan_all_cutpoints_overlayed_exp_abs_median.png"))
  panel_d_path <- require_file(file.path(manhattan_dir, "manhattan_all_cutpoints_overlayed_RATE.png"))

  lz_files <- sort(list.files(locuszoom_dir, pattern = "\\.png$", full.names = TRUE, ignore.case = TRUE))
  if (length(lz_files) == 0) stop("No PNG files found in locus-zoom directory: ", locuszoom_dir)

  analysis_name <- basename(dirname(manhattan_dir))
  dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)
  output_path <- file.path(args$output_dir, paste0("figure_", analysis_name, ".png"))

  message("Panel A: ", panel_a_path)
  message("Panel B: ", length(lz_files), " locus-zoom plots from ", locuszoom_dir)
  message("Panel C: ", panel_c_path)
  message("Panel D: ", panel_d_path)
  message("Output:  ", output_path)

  panel_a <- panel_from_png(panel_a_path)
  panel_c <- panel_from_png(panel_c_path)
  panel_d <- panel_from_png(panel_d_path)
  lz_panels <- lapply(lz_files, panel_from_png)
  row_b <- plot_grid(plotlist = lz_panels, nrow = 1)

  figure <- plot_grid(
    panel_a, row_b, panel_c, panel_d,
    ncol = 1,
    labels = c("A", "B", "C", "D"),
    label_size = 28,
    label_fontface = "bold",
    rel_heights = c(8, 3, 8, 8)
  )

  ggsave(
    output_path,
    figure,
    width = 16,
    height = 27,
    dpi = 300,
    limitsize = FALSE,
    bg = "white"
  )

  message("Wrote ", output_path)
}

main()
