#!/usr/bin/env Rscript
# Recreate just the RATE overlay Manhattan plot at 2/3 height, matching the
# styling of replot_ppom_overlay_manhattans_short.R. Reads RATE straight from the
# cppRATE outputs (raw RATE_values for the y-axis, phandango BP column for genome
# coordinates), so it works for both PPOM runs (per-cutpoint overlay) and binary
# GWAS runs (single series). The plot is built by make_rate_manhattan() in
# manhattan_builders.R, which the locus-zoom composites also call.
#
# Usage:
#   Rscript replot_rate_manhattan_short.R --run-dir <inference output dir> \
#     [--output-dir <dir>] [--positions-file <variant_index.csv>] \
#     [--annotations <snpEff fields .txt>] [--genes-of-interest <goi .txt>] \
#     [--label-mode top_n|gene_list] [--label-genes "<gene1,gene2,...>"] \
#     [--n-labels <int>]
#
# Gene labels are added only when --annotations and --label-mode are supplied.
# gene_list labels the top variant in each named (display-name) gene; top_n labels
# the top --n-labels unique genes by RATE.
#
# Default --output-dir is
#   /nfs/research/jlees/jacqueline/thesis_results/paper_figures/<basename(run-dir)>/manhattan_short
#
# Output is manhattan_all_cutpoints_overlayed_RATE.png for a PPOM run, or
# manhattan_RATE.png for a binary run.

DEFAULT_PAPER_FIGURES_DIR <- "/nfs/research/jlees/jacqueline/thesis_results/paper_figures"

source("/nfs/research/jlees/jacqueline/thesis_code/paper_figures/manhattan_builders.R")

parse_args <- function(argv) {
  out <- list(run_dir = NULL, output_dir = NULL, positions_file = NULL,
              annotations = NULL, genes_of_interest = NULL,
              label_mode = NULL, label_genes = NULL, n_labels = 10L)
  i <- 1
  while (i <= length(argv)) {
    a <- argv[i]
    if (a == "--run-dir") {
      out$run_dir <- argv[i + 1]; i <- i + 2
    } else if (a == "--output-dir") {
      out$output_dir <- argv[i + 1]; i <- i + 2
    } else if (a == "--positions-file") {
      out$positions_file <- argv[i + 1]; i <- i + 2
    } else if (a == "--annotations") {
      out$annotations <- argv[i + 1]; i <- i + 2
    } else if (a == "--genes-of-interest") {
      out$genes_of_interest <- argv[i + 1]; i <- i + 2
    } else if (a == "--label-mode") {
      out$label_mode <- argv[i + 1]; i <- i + 2
    } else if (a == "--label-genes") {
      out$label_genes <- trimws(strsplit(argv[i + 1], ",")[[1]]); i <- i + 2
    } else if (a == "--n-labels") {
      out$n_labels <- as.integer(argv[i + 1]); i <- i + 2
    } else {
      stop("Unknown argument: ", a)
    }
  }
  if (is.null(out$run_dir)) stop("Required arg: --run-dir <dir>")
  if (is.null(out$output_dir)) {
    out$output_dir <- file.path(DEFAULT_PAPER_FIGURES_DIR,
                                basename(normalizePath(out$run_dir, mustWork = TRUE)),
                                "manhattan_short")
  }
  out
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  run_dir <- normalizePath(args$run_dir, mustWork = TRUE)
  dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)

  message("Run dir:    ", run_dir)
  message("Output dir: ", args$output_dir)

  p <- make_rate_manhattan(
    run_dir, positions_file = args$positions_file,
    annotations = args$annotations, goi = args$genes_of_interest,
    label_mode = args$label_mode, label_genes = args$label_genes,
    n_labels = args$n_labels
  )

  path <- file.path(args$output_dir, attr(p, "outfile"))
  ggplot2::ggsave(path, plot = p, width = MB_FIG_WIDTH, height = MB_SHORT_HEIGHT,
                  dpi = 300)
  message("Wrote ", path)
}

# Run main() only when executed directly (Rscript), not when sourced.
if (identical(environment(), globalenv()) && sys.nframe() == 0L) {
  main()
}
