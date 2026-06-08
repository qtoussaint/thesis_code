#!/usr/bin/env Rscript
# Compose a cross-drug RATE paper figure from two manhattan_short directories
# and two locus-zoom directories:
#   A  Spn trimethoprim (standard PPOM) RATE manhattan (short)
#   B  locus zoom of the same (folP, folA RATE)
#   C  TB rifampicin RATE manhattan (short)
#   D  locus zoom of the same (rpoA, rpoB, rpoC RATE)
#
# The Spn PPOM run overlays all cutpoints, so its manhattan filename is
# manhattan_all_cutpoints_overlayed_RATE.png. The TB binary GWAS run produces a
# single manhattan per metric, so its filename is manhattan_RATE.png.
#
# Usage:
#   Rscript manhattan_locuszoom_composite_rate_tmp_rif.R \
#     <tmp_manhattan_dir> <tmp_locuszoom_dir> <rif_manhattan_dir> <rif_locuszoom_dir> \
#     [--output-dir <dir>] [--analysis-name <name>]
#
# Output: <output_dir>/figure_<analysis_name>_rate_tmp_rif.png
# where <analysis_name> defaults to the parent directory name of <tmp_manhattan_dir>.

suppressPackageStartupMessages({
  library(cowplot)
  library(ggplot2)
  library(magick)
})

DEFAULT_OUTPUT_DIR <- "/nfs/research/jlees/jacqueline/thesis_results/paper_figures/manhattans_with_locus"

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
  if (length(positional) != 4) {
    stop("Usage: Rscript manhattan_locuszoom_composite_rate_tmp_rif.R <tmp_manhattan_dir> <tmp_locuszoom_dir> <rif_manhattan_dir> <rif_locuszoom_dir> [--output-dir <dir>] [--analysis-name <name>]")
  }
  list(
    tmp_manhattan_dir = positional[1],
    tmp_locuszoom_dir = positional[2],
    rif_manhattan_dir = positional[3],
    rif_locuszoom_dir = positional[4],
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

locuszoom_row <- function(locuszoom_dir, genes) {
  files <- file.path(locuszoom_dir, paste0(genes, "_rate_nolabels.png"))
  missing <- files[!file.exists(files)]
  if (length(missing) > 0) stop("Missing locus-zoom PNG(s): ", paste(missing, collapse = ", "))
  panels <- lapply(files, panel_from_png)
  list(
    row = plot_grid(plotlist = lapply(panels, `[[`, "plot"), nrow = 1),
    aspect = max(vapply(panels, `[[`, numeric(1), "aspect")) / length(panels)
  )
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))

  tmp_manhattan_dir <- normalizePath(args$tmp_manhattan_dir, mustWork = TRUE)
  tmp_locuszoom_dir <- normalizePath(args$tmp_locuszoom_dir, mustWork = TRUE)
  rif_manhattan_dir <- normalizePath(args$rif_manhattan_dir, mustWork = TRUE)
  rif_locuszoom_dir <- normalizePath(args$rif_locuszoom_dir, mustWork = TRUE)

  # Panel A: Spn trimethoprim PPOM RATE manhattan (cutpoints overlayed)
  panel_a_path <- require_file(file.path(tmp_manhattan_dir, "manhattan_all_cutpoints_overlayed_RATE.png"))
  # Panel C: TB rifampicin binary-GWAS RATE manhattan
  panel_c_path <- require_file(file.path(rif_manhattan_dir, "manhattan_RATE.png"))

  analysis_name <- if (!is.null(args$analysis_name)) args$analysis_name else basename(dirname(tmp_manhattan_dir))
  dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)
  output_path <- file.path(args$output_dir, paste0("figure_", analysis_name, "_rate_tmp_rif.png"))

  message("Panel A: ", panel_a_path)
  message("Panel B: trimethoprim RATE locus-zoom plots from ", tmp_locuszoom_dir)
  message("Panel C: ", panel_c_path)
  message("Panel D: rifampicin RATE locus-zoom plots from ", rif_locuszoom_dir)
  message("Output:  ", output_path)

  panel_a <- panel_from_png(panel_a_path)
  panel_b <- locuszoom_row(tmp_locuszoom_dir, c("folP", "folA"))
  panel_c <- panel_from_png(panel_c_path)
  panel_d <- locuszoom_row(rif_locuszoom_dir, c("rpoA", "rpoB", "rpoC"))

  fig_width <- 16
  h_a <- fig_width * panel_a$aspect
  h_b <- fig_width * panel_b$aspect
  h_c <- fig_width * panel_c$aspect
  h_d <- fig_width * panel_d$aspect

  figure <- plot_grid(
    panel_a$plot, panel_b$row, panel_c$plot, panel_d$row,
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
