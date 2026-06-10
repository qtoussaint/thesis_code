#!/usr/bin/env Rscript
# Compose the cross-drug RATE paper figure:
#   A  Spn trimethoprim (PPOM) RATE manhattan  — top variant in folA / folP labeled
#   B  locus zoom of the same (folP, folA), shared legend at right
#   C  TB rifampicin (binary) RATE manhattan   — top 10 genes labeled
#   D  locus zoom of the same (rpoA, rpoB, rpoC), shared legend at right
#
# The manhattan panels are built natively (manhattan_builders.R) so the magnifier
# connector lines from each locus zoom's genomic window up to its manhattan region
# land at exact base-pair coordinates. The locus-zoom panels come from
# make_locuszoom_plot.R --composite_style (no per-plot legend, no axis titles, plus
# a shared <gene>_<metric>_legend.png and a <gene>_<metric>_region.txt sidecar).
#
# Usage:
#   Rscript manhattan_locuszoom_composite_rate_tmp_rif.R \
#     [--output-dir <dir>] [--analysis-name <name>]
#
# Output: <output_dir>/figure_<analysis_name>_rate_tmp_rif.png

suppressPackageStartupMessages({
  library(cowplot)
  library(ggplot2)
  library(magick)
})

CODE_DIR <- "/nfs/research/jlees/jacqueline/thesis_code/paper_figures"
source(file.path(CODE_DIR, "manhattan_builders.R"))
source(file.path(CODE_DIR, "composite_connectors.R"))

# --- inputs -----------------------------------------------------------------
RESULTS   <- "/nfs/research/jlees/jacqueline/thesis_results"
GOI_DIR   <- "/nfs/research/jlees/jacqueline/thesis_code/gwas_genesofinterest"
SPN_ANNOT <- "/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/genotype/fields_filtered_maf05_multiallelic.txt"
TB_ANNOT  <- "/nfs/research/jlees/jacqueline/gwas_data/tuberculosis/cryptic_regeno_snpeff/fields_filtered.txt"

TMP_RUN <- file.path(RESULTS, "gwas_spn_trimethoprim/inference/05_spn_trimethoprim_MIC_PPOM")
TB_RUN  <- file.path(RESULTS, "gwas_tb_rifampicin/inference/07_tb_rifampicin_binary_logistic")
TMP_LZ  <- file.path(RESULTS, "locus_zoom/spneumoniae/plots/05_spn_trimethoprim_MIC_PPOM_genes_top2_composite")
TB_LZ   <- file.path(RESULTS, "locus_zoom/mtuberculosis/plots/07_tb_rifampicin_binary_logistic_genes_top3_composite")

TMP_GOI <- file.path(GOI_DIR, "spn_trimethoprim_genesofinterest.txt")
TB_GOI  <- file.path(GOI_DIR, "tb_rifampicin_genesofinterest.txt")

TMP_LZ_GENES <- c("folP", "folA")  # panel B, left-to-right
TB_LZ_GENES  <- c("rpoB")          # panel D (rpoB only)

DEFAULT_OUTPUT_DIR <- file.path(RESULTS, "paper_figures/manhattans_with_locus")

parse_args <- function(argv) {
  out <- list(output_dir = DEFAULT_OUTPUT_DIR, analysis_name = "spn_tmp_tb_rif_labeled")
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
  output_path <- file.path(args$output_dir,
                           paste0("figure_", args$analysis_name, "_rate_tmp_rif.png"))

  message("Building panel A (trimethoprim RATE manhattan)...")
  pA <- make_rate_manhattan(TMP_RUN, annotations = SPN_ANNOT, goi = TMP_GOI,
                            label_mode = "gene_list",
                            label_genes = c("folA (dhfR)", "folP"))
  message("Building panel C (TB RATE manhattan)...")
  pC <- make_rate_manhattan(TB_RUN, annotations = TB_ANNOT, goi = TB_GOI,
                            label_mode = "top_n", n_labels = 10L)

  message("Building locus-zoom rows...")
  rowB <- mb_lz_row(TMP_LZ, TMP_LZ_GENES, "rate")
  rowD <- mb_lz_row(TB_LZ,  TB_LZ_GENES,  "rate")

  res <- mb_assemble(list(
    list(manhattan = pA, lz = rowB),
    list(manhattan = pC, lz = rowD)
  ))

  message("Output: ", output_path)
  ggsave(output_path, res$canvas, width = res$width, height = res$height,
         dpi = 300, limitsize = FALSE, bg = "white")
  message("Wrote ", output_path)
}

main()
