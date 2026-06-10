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
              analysis_name = "02_spn_penicillin_MIC_PPOM_labeled")
  i <- 1
  while (i <= length(argv)) {
    a <- argv[i]
    if (a == "--output-dir") { out$output_dir <- argv[i + 1]; i <- i + 2 }
    else if (a == "--analysis-name") { out$analysis_name <- argv[i + 1]; i <- i + 2 }
    else stop("Unknown argument: ", a)
  }
  out
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)
  full_path    <- file.path(args$output_dir, paste0("figure_", args$analysis_name, ".png"))
  smaller_path <- file.path(args$output_dir, paste0("figure_", args$analysis_name, "_smaller.png"))

  message("Building beta manhattans (median + exp|median|)...")
  betas <- make_beta_manhattans(PEN_RUN, PEN_POS, annotations = SPN_ANNOT, goi = GOI,
                                label_genes = BETA_GENES, gene_aliases = GENE_ALIASES)
  pA <- betas$median
  pB <- betas$exp_abs
  message("Building RATE manhattan...")
  pD <- make_rate_manhattan(PEN_RUN, annotations = SPN_ANNOT, goi = GOI,
                            label_mode = "gene_list", label_genes = RATE_GENES)

  message("Building locus-zoom rows...")
  rowC <- mb_lz_row(PEN_LZ, LZ_GENES, "exp_abs_median")
  rowE <- mb_lz_row(PEN_LZ, LZ_GENES, "rate")

  # Full figure: median beta, exp(|median|) beta + locus zoom, RATE + locus zoom
  full <- mb_assemble(list(
    list(manhattan = pA, lz = NULL),
    list(manhattan = pB, lz = rowC),
    list(manhattan = pD, lz = rowE)
  ))
  message("Output: ", full_path)
  ggsave(full_path, full$canvas, width = full$width, height = full$height,
         dpi = 300, limitsize = FALSE, bg = "white")
  message("Wrote ", full_path)

  # "Smaller" figure: drop the median beta panel — exp(|median|) beta + locus zoom,
  # RATE + locus zoom only. Rebuild the locus-zoom rows so each carries a fresh
  # (single-use) magick image handle.
  rowC2 <- mb_lz_row(PEN_LZ, LZ_GENES, "exp_abs_median")
  rowE2 <- mb_lz_row(PEN_LZ, LZ_GENES, "rate")
  smaller <- mb_assemble(list(
    list(manhattan = pB, lz = rowC2),
    list(manhattan = pD, lz = rowE2)
  ))
  message("Output: ", smaller_path)
  ggsave(smaller_path, smaller$canvas, width = smaller$width, height = smaller$height,
         dpi = 300, limitsize = FALSE, bg = "white")
  message("Wrote ", smaller_path)
}

main()
