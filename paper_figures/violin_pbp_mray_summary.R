#!/usr/bin/env Rscript
# Paper-figure violin plots: PBP genes + mraY, nucleotide labels, most-disruptive impact.
#
# Re-creates the impact-coloured per-gene effect-size violins that the inference
# pipeline already emits for the SPN penicillin PPOM runs, but tailored for the paper:
#   * only the 6 genes of interest: 5 PBPs (pbp2x/pbp1a/pbp1b/pbp2a/pbp2b) + mraY
#   * only the two impact-coloured plots: abs_beta_impact and abs_delta_beta_impact
#   * outlier labels from the NUCLEOTIDE annotation (ANN[*].HGVS_C, e.g. c.378T>G)
#     instead of the protein annotation (ANN[*].HGVS_P)
#   * for multi-allelic positions with more than one functional impact, the MOST
#     disruptive one is chosen (snpEff IMPACT severity HIGH > MODERATE > LOW), and that
#     row's EFFECT (colour) and HGVS_C (label) are used.
#
# Reuses gwas_workflow's write_ppom_gene_delta_plots() so the figure style matches the
# rest of the thesis; the three modifications above are achieved purely by constructing
# the gene/impact/label input vectors, so the package is left untouched.
#
# Run with:
#   mamba activate gwas_pipeline
#   Rscript paper_figures/violin_pbp_mray_summary.R
#
# Output: <output_dir>/violin_plots/<run>_gene_violin_abs_beta_impact.png
#         <output_dir>/violin_plots/<run>_gene_violin_abs_delta_beta_impact.png
#         <output_dir>/violin_plots/<run>_violin_data.csv

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(viridis)
})

RESULTS_ROOT <- "/nfs/research/jlees/jacqueline/thesis_results"
OUTPUT_DIR   <- file.path(RESULTS_ROOT, "paper_figures", "violin_plots")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Pull in write_ppom_gene_delta_plots() + the internal .italic_gene_labeller it uses.
source("/nfs/research/jlees/jacqueline/gwas_workflow/code/gwas_workflow/R/manhattan_plots.R")

ANNOTATIONS <- "/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/genotype/fields_filtered_maf05_multiallelic.txt"

# col1 = annotation gene name to filter on; col2 = display label used in facet strips.
GENES_OF_INTEREST <- list(
  c("pbpX", "pbp1A", "pbp1B", "pbp2A", "penA", "mraY"),
  c("pbp2x", "pbp1a", "pbp1b", "pbp2a", "pbp2b", "mraY")
)

# Each run: where the fitted effects live and the variant index giving genomic positions
# (column 2, row-aligned to variant_id) — the same --phandango files the pipeline used.
RUNS <- list(
  list(
    name  = "02_spn_penicillin_MIC_PPOM",
    eff   = file.path(RESULTS_ROOT, "gwas_spn_penicillin", "inference",
                      "02_spn_penicillin_MIC_PPOM", "fitted_model",
                      "depruned_variant_effects.csv"),
    index = file.path(RESULTS_ROOT, "gwas_datasets", "inference",
                      "02_spn_penicillin_MIC",
                      "02_spn_penicillin_MIC_variant_index.csv")
  ),
  list(
    name  = "16_spn_penicillin_MIC_minimabinning_PPOM",
    eff   = file.path(RESULTS_ROOT, "gwas_spn_penicillin", "inference",
                      "16_spn_penicillin_MIC_minimabinning_PPOM", "fitted_model",
                      "depruned_variant_effects.csv"),
    index = file.path(RESULTS_ROOT, "gwas_datasets", "inference",
                      "16_spn_penicillin_MIC_minimabinning",
                      "16_spn_penicillin_MIC_minimabinning_variant_index.csv")
  )
)

# -----------------------------------------------------------------------------
# Most-disruptive annotation per genomic position.
#
# The snpEff table has one row per ALT allele / transcript, so a position can carry
# several functional impacts. For each position we keep the row whose IMPACT is most
# severe (HIGH > MODERATE > LOW), and return that row's gene, EFFECT and HGVS_C.
# Returns vectors indexed positionally by variant_id (the i-th element is variant i).
IMPACT_RANK <- c(MODIFIER = 0L, LOW = 1L, MODERATE = 2L, HIGH = 3L)

build_annotation <- function(positions) {
  ann <- read.delim(ANNOTATIONS, stringsAsFactors = FALSE, check.names = FALSE)
  sev <- IMPACT_RANK[ann[["ANN[*].IMPACT"]]]
  sev[is.na(sev)] <- 0L

  # For each position, the index of the most-disruptive row (first wins on ties).
  ord       <- order(ann$POS, -sev)            # within a POS, highest severity first
  ann_ord   <- ann[ord, ]
  first_hit <- !duplicated(ann_ord$POS)        # first row per POS = most disruptive
  best      <- ann_ord[first_hit, ]
  idx       <- match(positions, best$POS)       # one chosen row per variant position

  gene  <- best[["ANN[*].GENE"]][idx]
  gene[is.na(gene)] <- "MODIFIER"
  list(
    gene    = gene,
    impact  = best[["ANN[*].EFFECT"]][idx],     # functional consequence -> colour
    hgvs_c  = best[["ANN[*].HGVS_C"]][idx]       # nucleotide notation -> outlier label
  )
}

# -----------------------------------------------------------------------------
process_run <- function(run) {
  message("Processing ", run$name)
  eff <- read.csv(run$eff, stringsAsFactors = FALSE)
  pos <- read.csv(run$index, stringsAsFactors = FALSE)[, 2]

  n_var <- max(eff$variant_id)
  stopifnot(length(pos) == n_var)

  ann <- build_annotation(pos)
  ids <- as.character(seq_len(n_var))

  # Nucleotide labels; blanks (".", NA, "") fall back to the variant_id.
  labels <- ann$hgvs_c
  labels[is.na(labels) | labels == "." | !nzchar(labels)] <-
    ids[is.na(labels) | labels == "." | !nzchar(labels)]

  variant_labels  <- setNames(labels, ids)
  variant_impacts <- setNames(ann$impact, ids)

  # Render every violin into a scratch dir, then lift out only the two we want.
  tmp <- tempfile("violin_pbp_mray_")
  dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  write_ppom_gene_delta_plots(
    depruned_all_cutpoints = eff,
    variant_positions      = pos,
    gene_names             = ann$gene,
    genes_of_interest      = GENES_OF_INTEREST,
    output_dir             = tmp,
    all_rate_values        = NULL,            # no RATE plots
    variant_labels         = variant_labels,
    variant_impacts        = variant_impacts
  )

  for (f in c("gene_violin_abs_beta_impact.png",
              "gene_violin_abs_delta_beta_impact.png")) {
    src <- file.path(tmp, "plots", "violin", f)
    if (!file.exists(src))
      stop("expected plot not produced: ", src)
    file.copy(src, file.path(OUTPUT_DIR, paste0(run$name, "_", f)),
              overwrite = TRUE)
  }

  # Companion CSV of the filtered underlying data, mirroring the plot inputs.
  display_map <- setNames(GENES_OF_INTEREST[[2]], GENES_OF_INTEREST[[1]])
  df <- eff
  df$position         <- pos[df$variant_id]
  df$gene             <- ann$gene[df$variant_id]
  df$impact           <- variant_impacts[as.character(df$variant_id)]
  df$nucleotide_label <- variant_labels[as.character(df$variant_id)]
  df <- df[df$gene %in% GENES_OF_INTEREST[[1]], , drop = FALSE]
  df$display_name <- display_map[df$gene]
  df <- df[order(df$variant_id, df$cutpoint), ]
  df$delta_abs <- ave(df$median, df$variant_id,
                      FUN = function(x) c(NA_real_, abs(diff(x))))

  out_cols <- c("variant_id", "position", "gene", "display_name",
                "cutpoint", "cutpoint_MIC", "median", "delta_abs",
                "impact", "nucleotide_label")
  write.csv(df[, out_cols],
            file.path(OUTPUT_DIR, paste0(run$name, "_violin_data.csv")),
            row.names = FALSE)
}

invisible(lapply(RUNS, process_run))
message("Done. Output written to ", OUTPUT_DIR)
