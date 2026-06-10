#!/usr/bin/env Rscript
# Recreate the three PPOM overlay Manhattan plots used in the paper composite
# (median effects, exp(|median|), RATE) at 2/3 height. The beta panels are built
# by make_beta_manhattans() and the RATE panel by make_rate_manhattan() in
# manhattan_builders.R, which the locus-zoom composites also call.
#
# Usage:
#   Rscript replot_ppom_overlay_manhattans_short.R \
#     --run-dir   <inference output dir> \
#     --phandango <variant_index.csv> \
#     [--output-dir <dir>] \
#     [--annotations <snpEff fields .txt>] [--genes-of-interest <goi .txt>] \
#     [--label-genes "<gene1,gene2,...>"]
#
# Gene labels are added to the beta panels only when --annotations and
# --label-genes are supplied (top variant in each named display-name gene, ranked
# on |median|). The RATE panel is left unlabeled here; the labeled RATE panel is
# produced by replot_rate_manhattan_short.R.
#
# Default --output-dir is
#   /nfs/research/jlees/jacqueline/thesis_results/paper_figures/<basename(run-dir)>/manhattan_short

DEFAULT_PAPER_FIGURES_DIR <- "/nfs/research/jlees/jacqueline/thesis_results/paper_figures"

source("/nfs/research/jlees/jacqueline/thesis_code/paper_figures/manhattan_builders.R")

parse_args <- function(argv) {
  out <- list(run_dir = NULL, phandango = NULL, output_dir = NULL,
              annotations = NULL, genes_of_interest = NULL, label_genes = NULL)
  i <- 1
  while (i <= length(argv)) {
    a <- argv[i]
    if (a == "--run-dir") {
      out$run_dir <- argv[i + 1]; i <- i + 2
    } else if (a == "--phandango") {
      out$phandango <- argv[i + 1]; i <- i + 2
    } else if (a == "--output-dir") {
      out$output_dir <- argv[i + 1]; i <- i + 2
    } else if (a == "--annotations") {
      out$annotations <- argv[i + 1]; i <- i + 2
    } else if (a == "--genes-of-interest") {
      out$genes_of_interest <- argv[i + 1]; i <- i + 2
    } else if (a == "--label-genes") {
      out$label_genes <- trimws(strsplit(argv[i + 1], ",")[[1]]); i <- i + 2
    } else {
      stop("Unknown argument: ", a)
    }
  }
  if (is.null(out$run_dir) || is.null(out$phandango)) {
    stop("Required args: --run-dir <dir> --phandango <csv>")
  }
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
  phandango <- normalizePath(args$phandango, mustWork = TRUE)
  dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)

  message("Run dir:    ", run_dir)
  message("Phandango:  ", phandango)
  message("Output dir: ", args$output_dir)

  save_short <- function(plot, name) {
    path <- file.path(args$output_dir, name)
    ggplot2::ggsave(path, plot = plot, width = MB_FIG_WIDTH,
                    height = MB_SHORT_HEIGHT, dpi = 300)
    message("Wrote ", path)
  }

  betas <- make_beta_manhattans(run_dir, phandango,
                                annotations = args$annotations,
                                goi = args$genes_of_interest,
                                label_genes = args$label_genes)
  save_short(betas$median,  "manhattan_all_cutpoints_overlayed_median_effects.png")
  save_short(betas$exp_abs, "manhattan_all_cutpoints_overlayed_exp_abs_median.png")

  p_rate <- make_rate_manhattan(run_dir)
  save_short(p_rate, "manhattan_all_cutpoints_overlayed_RATE.png")
}

# Run main() only when executed directly (Rscript), not when sourced.
if (identical(environment(), globalenv()) && sys.nframe() == 0L) {
  main()
}
