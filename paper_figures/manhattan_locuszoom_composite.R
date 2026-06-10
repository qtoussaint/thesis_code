#!/usr/bin/env Rscript
# Compose the SPN penicillin PPOM paper figure (manhattans + locus zooms):
#   A  median-effect beta manhattan        — pbp1a / pbp2X / pbp2b / folA labeled
#   B  exp(|median|) beta manhattan         — pbp1a / pbp2X / pbp2b / folA labeled
#   C  locus zoom of B (pbp2X, pbp1a, pbp2b), shared legend at right
#   D  RATE manhattan                       — pbp1a / pbp2X / pbp2b labeled
#   E  locus zoom of D (pbp2X, pbp1a, pbp2b), shared legend at right
#
# Manhattans are built natively (manhattan_builders.R) so the magnifier connector
# lines from each locus zoom up to its manhattan land at exact base-pair
# coordinates. Locus-zoom panels come from make_locuszoom_plot.R --composite_style.
# The dihydrofolate-reductase gene is annotated "dyr" in the genome; on this
# penicillin figure it is relabeled "folA".
#
# Usage:
#   Rscript manhattan_locuszoom_composite.R [--output-dir <dir>] [--analysis-name <name>]
#
# Output: <output_dir>/figure_<analysis_name>.png

suppressPackageStartupMessages({
  library(cowplot)
  library(ggplot2)
  library(magick)
})

CODE_DIR <- "/nfs/research/jlees/jacqueline/thesis_code/paper_figures"
source(file.path(CODE_DIR, "manhattan_builders.R"))
source(file.path(CODE_DIR, "composite_connectors.R"))

RESULTS   <- "/nfs/research/jlees/jacqueline/thesis_results"
GOI       <- "/nfs/research/jlees/jacqueline/thesis_code/gwas_genesofinterest/spn_penicillin_genesofinterest.txt"
SPN_ANNOT <- "/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/genotype/fields_filtered_maf05_multiallelic.txt"

PEN_RUN  <- file.path(RESULTS, "gwas_spn_penicillin/inference/02_spn_penicillin_MIC_PPOM")
PEN_POS  <- file.path(RESULTS, "gwas_datasets/inference/02_spn_penicillin_MIC/02_spn_penicillin_MIC_variant_index.csv")
PEN_LZ   <- file.path(RESULTS, "locus_zoom/spneumoniae/plots/02_spn_penicillin_MIC_PPOM_pbp_top5_composite")

LZ_GENES     <- c("pbp2X", "pbp1a", "pbp2b")
BETA_GENES   <- c("pbp1a", "pbp2X", "pbp2b", "folA")  # folA = relabeled dyr
RATE_GENES   <- c("pbp1a", "pbp2X", "pbp2b")
GENE_ALIASES <- c(dyr = "folA")

DEFAULT_OUTPUT_DIR <- file.path(RESULTS, "paper_figures/manhattans_with_locus")

parse_args <- function(argv) {
  out <- list(output_dir = DEFAULT_OUTPUT_DIR,
              analysis_name = "02_spn_penicillin_MIC_PPOM_labeled",
              layout = "full")   # "full" (A-E) or "smaller" (exp|median| + RATE only)
  i <- 1
  while (i <= length(argv)) {
    a <- argv[i]
    if (a == "--output-dir") { out$output_dir <- argv[i + 1]; i <- i + 2 }
    else if (a == "--analysis-name") { out$analysis_name <- argv[i + 1]; i <- i + 2 }
    else if (a == "--layout") { out$layout <- argv[i + 1]; i <- i + 2 }
    else stop("Unknown argument: ", a)
  }
  if (!out$layout %in% c("full", "smaller")) stop("--layout must be full or smaller")
  out
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)
  suffix <- if (args$layout == "smaller") "_smaller" else ""
  output_path <- file.path(args$output_dir,
                           paste0("figure_", args$analysis_name, suffix, ".png"))

  message("Building exp(|median|) beta + RATE manhattans...")
  # penicillin: label the regional peak within +/-10kb of each gene's top SNP and
  # shift each label right until its text box clears every plotted point
  betas <- make_beta_manhattans(PEN_RUN, PEN_POS, annotations = SPN_ANNOT, goi = GOI,
                                label_genes = BETA_GENES, gene_aliases = GENE_ALIASES,
                                label_window_bp = 10000, label_no_overlap = TRUE)
  pB <- betas$exp_abs
  pD <- make_rate_manhattan(PEN_RUN, annotations = SPN_ANNOT, goi = GOI,
                            label_mode = "gene_list", label_genes = RATE_GENES,
                            label_window_bp = 10000, label_no_overlap = TRUE)

  message("Building locus-zoom rows...")
  rowC <- mb_lz_row(PEN_LZ, LZ_GENES, "exp_abs_median")
  rowE <- mb_lz_row(PEN_LZ, LZ_GENES, "rate")

  if (args$layout == "full") {
    # median beta, exp(|median|) beta + locus zoom, RATE + locus zoom
    units <- list(
      list(manhattan = betas$median, lz = NULL),
      list(manhattan = pB, lz = rowC),
      list(manhattan = pD, lz = rowE)
    )
  } else {
    # "smaller": drop the median beta panel
    units <- list(
      list(manhattan = pB, lz = rowC),
      list(manhattan = pD, lz = rowE)
    )
  }

  res <- mb_assemble(units)
  message("Output: ", output_path)
  ggsave(output_path, res$canvas, width = res$width, height = res$height,
         dpi = 300, limitsize = FALSE, bg = "white")
  message("Wrote ", output_path)
}

main()
