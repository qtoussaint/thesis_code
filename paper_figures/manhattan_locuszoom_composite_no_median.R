#!/usr/bin/env Rscript
# Compose one PPOM paper figure (without the median_effects manhattan row) from a
# manhattan_short directory and a locus-zoom directory.
#
# Usage:
#   Rscript manhattan_locuszoom_composite_no_median.R <manhattan_dir> <locuszoom_dir> [--output-dir <dir>] [--analysis-name <name>]
#
# Output: <output_dir>/figure_<analysis_name>_no_median.png
# where <analysis_name> defaults to the parent directory name of <manhattan_dir>.

suppressPackageStartupMessages({
  library(cowplot)
  library(ggplot2)
  library(magick)
})

DEFAULT_OUTPUT_DIR <- "/nfs/research/jlees/jacqueline/thesis_results/paper_figures"

parse_args <- function(argv) {
  output_dir <- DEFAULT_OUTPUT_DIR
  analysis_name <- NULL
  positional <- character(0)
  i <- 1
  while (i <= length(argv)) {
    a <- argv[i]
    if (a == "--output-dir") {
      if (i == length(argv)) stop("--output-dir requires a value")
      output_dir <- argv[i + 1]
      i <- i + 2
    } else if (a == "--analysis-name") {
      if (i == length(argv)) stop("--analysis-name requires a value")
      analysis_name <- argv[i + 1]
      i <- i + 2
    } else {
      positional <- c(positional, a)
      i <- i + 1
    }
  }
  if (length(positional) != 2) {
    stop("Usage: Rscript manhattan_locuszoom_composite_no_median.R <manhattan_dir> <locuszoom_dir> [--output-dir <dir>] [--analysis-name <name>]")
  }
  list(
    manhattan_dir = positional[1],
    locuszoom_dir = positional[2],
    output_dir = output_dir,
    analysis_name = analysis_name
  )
}

require_file <- function(path) {
  if (!file.exists(path)) stop("Missing input file: ", path)
  path
}

panel_from_png <- function(path) {
  trimmed <- magick::image_trim(magick::image_read(path))
  info <- magick::image_info(trimmed)
  list(
    plot = ggdraw() + draw_image(trimmed),
    aspect = info$height / info$width
  )
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))

  manhattan_dir <- normalizePath(args$manhattan_dir, mustWork = TRUE)
  locuszoom_dir <- normalizePath(args$locuszoom_dir, mustWork = TRUE)

  panel_a_path <- require_file(file.path(manhattan_dir, "manhattan_all_cutpoints_overlayed_exp_abs_median.png"))
  panel_c_path <- require_file(file.path(manhattan_dir, "manhattan_all_cutpoints_overlayed_RATE.png"))

  lz_genes <- c("pbp2X", "pbp1a", "pbp2b")
  lz_files <- file.path(locuszoom_dir, paste0(lz_genes, "_exp_abs_median_nolabels.png"))
  lz_rate_files <- file.path(locuszoom_dir, paste0(lz_genes, "_rate_nolabels.png"))
  missing <- c(lz_files, lz_rate_files)[!file.exists(c(lz_files, lz_rate_files))]
  if (length(missing) > 0) stop("Missing locus-zoom PNG(s): ", paste(missing, collapse = ", "))

  analysis_name <- if (!is.null(args$analysis_name)) args$analysis_name else basename(dirname(manhattan_dir))
  dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)
  output_path <- file.path(args$output_dir, paste0("figure_", analysis_name, "_no_median.png"))

  message("Panel A: ", panel_a_path)
  message("Panel B: ", length(lz_files), " locus-zoom plots from ", locuszoom_dir)
  message("Panel C: ", panel_c_path)
  message("Panel D: ", length(lz_rate_files), " RATE locus-zoom plots from ", locuszoom_dir)
  message("Output:  ", output_path)

  panel_a <- panel_from_png(panel_a_path)
  panel_c <- panel_from_png(panel_c_path)
  lz_panels <- lapply(lz_files, panel_from_png)
  lz_rate_panels <- lapply(lz_rate_files, panel_from_png)
  row_b <- plot_grid(plotlist = lapply(lz_panels, `[[`, "plot"), nrow = 1)
  row_d <- plot_grid(plotlist = lapply(lz_rate_panels, `[[`, "plot"), nrow = 1)

  fig_width <- 16
  h_a <- fig_width * panel_a$aspect
  h_c <- fig_width * panel_c$aspect
  h_b <- (fig_width / length(lz_panels)) *
    max(vapply(lz_panels, `[[`, numeric(1), "aspect"))
  h_d <- (fig_width / length(lz_rate_panels)) *
    max(vapply(lz_rate_panels, `[[`, numeric(1), "aspect"))

  figure <- plot_grid(
    panel_a$plot, row_b, panel_c$plot, row_d,
    ncol = 1,
    labels = c("A", "B", "C", "D"),
    label_size = 28,
    label_fontface = "bold",
    label_x = -0.01,
    hjust = 0,
    rel_heights = c(h_a, h_b, h_c, h_d)
  ) + theme(plot.margin = margin(5, 5, 5, 25))

  ggsave(
    output_path,
    figure,
    width = fig_width,
    height = h_a + h_b + h_c + h_d,
    dpi = 300,
    limitsize = FALSE,
    bg = "white"
  )

  message("Wrote ", output_path)
}

main()
