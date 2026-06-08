#!/usr/bin/env Rscript
# Combine the two SPN penicillin faceted RATE plots into a single A/B figure:
#   A = minima binning (16_spn_penicillin_MIC_minimabinning_PPOM, 4 cutpoints)
#   B = granular        (02_spn_penicillin_MIC_PPOM, 7 cutpoints)
# Panels are stacked vertically with relative heights proportional to their facet
# rows. Output goes to the faceted_cutpoints/ directory.

suppressPackageStartupMessages({
  library(ggplot2)
  library(viridis)
  library(ggrepel)
  library(cowplot)
})

# Reuse the per-run plot builder (sourcing does not trigger its main()).
source(file.path(dirname(sub("^--file=", "",
        grep("^--file=", commandArgs(FALSE), value = TRUE)[1])),
        "replot_rate_faceted_cutpoints.R"))

RESULTS  <- "/nfs/research/jlees/jacqueline/thesis_results"
ANNOT    <- "/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/genotype/fields_filtered_maf05_multiallelic.txt"
GOI      <- "/nfs/research/jlees/jacqueline/thesis_code/gwas_genesofinterest/spn_penicillin_genesofinterest.txt"
OUT_DIR  <- file.path(RESULTS, "paper_figures", "faceted_cutpoints")
NCOL     <- 2L

minima   <- file.path(RESULTS, "gwas_spn_penicillin/inference/16_spn_penicillin_MIC_minimabinning_PPOM")
granular <- file.path(RESULTS, "gwas_spn_penicillin/inference/02_spn_penicillin_MIC_PPOM")

a <- build_faceted_plot(minima,   ANNOT, GOI, NCOL)   # A: minima
b <- build_faceted_plot(granular, ANNOT, GOI, NCOL)   # B: granular

rows_a <- ceiling(a$n / NCOL)
rows_b <- ceiling(b$n / NCOL)

# A/B labels styled to match the other paper_figures composites (cowplot,
# label_size 24, bold).
figure <- cowplot::plot_grid(
  a$plot, b$plot,
  ncol = 1, labels = c("A", "B"),
  label_size = 36, label_fontface = "bold",
  rel_heights = c(rows_a, rows_b)
)

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
path <- file.path(OUT_DIR, "02_AND_16_spn_penicillin_MIC_PPOM_faceted_RATE_labeled.png")
ggplot2::ggsave(path, plot = figure, width = 22, height = 4.5 * (rows_a + rows_b),
                dpi = 300, limitsize = FALSE)
message("Wrote ", path)
