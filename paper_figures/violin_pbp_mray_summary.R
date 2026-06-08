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
#   * only the single top SNP (largest value) at each cutpoint within a gene facet is
#     labelled, instead of every statistical outlier.
#
# The two violins reproduce gwas_workflow's style exactly (and borrow its internal
# .italic_gene_labeller) but are drawn here so the label set can be the per-cutpoint
# top SNP; the package itself is left untouched.
#
# Run with:
#   mamba activate gwas_pipeline
#   Rscript paper_figures/violin_pbp_mray_summary.R
#
# Also emits one combined overview: a 2x2 figure (four facets = the four per-run plots)
# overlaying all 6 genes per cutpoint, coloured by gene, with no variant labels.
#
# Output: <output_dir>/violin_plots/<run>_gene_violin_abs_beta_impact.png
#         <output_dir>/violin_plots/<run>_gene_violin_abs_delta_beta_impact.png
#         <output_dir>/violin_plots/<run>_violin_data.csv
#         <output_dir>/violin_plots/overlay_gene_violin_grid.png

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

  # Build the filtered, per-cutpoint data frame the same way the pipeline does:
  # map gene/impact/label by variant_id, keep the 6 genes, then derive |β̃| and |Δβ̃|.
  display_map <- setNames(GENES_OF_INTEREST[[2]], GENES_OF_INTEREST[[1]])
  df <- eff
  df$position        <- pos[df$variant_id]
  df$gene            <- ann$gene[df$variant_id]
  df$impact          <- gsub("_", " ", variant_impacts[as.character(df$variant_id)])
  df$variant_label   <- variant_labels[as.character(df$variant_id)]
  df <- df[df$gene %in% GENES_OF_INTEREST[[1]], , drop = FALSE]
  df$display_name <- factor(display_map[df$gene],
                            levels = c("pbp1a", "pbp1b", "pbp2a",
                                       "pbp2b", "pbp2x", "mraY"))
  df <- df[order(df$variant_id, df$cutpoint), ]
  df$abs_median   <- abs(df$median)
  df$delta_signed <- ave(df$median, df$variant_id,
                         FUN = function(x) c(NA_real_, diff(x)))
  df$delta_abs    <- ave(df$median, df$variant_id,
                         FUN = function(x) c(NA_real_, abs(diff(x))))
  df$x_val <- as.numeric(as.character(df$cutpoint_MIC))   # MIC breakpoint per cutpoint

  x_axis_title <- expression("MIC breakpoint" ~ (mu * "g·mL"^{-1}))

  # Label only the single top SNP (largest value) within each gene facet at each
  # cutpoint, rather than every statistical outlier.
  top_per_group <- function(d, value_col) {
    do.call(rbind, lapply(
      split(d, interaction(d$display_name, d$x_val, drop = TRUE)),
      function(sub) sub[which.max(sub[[value_col]]), , drop = FALSE]
    ))
  }

  # Shared impact-coloured violin (matches gwas_workflow's style exactly).
  impact_violin <- function(d, y_col, ylabel, tops, out_path) {
    p <- ggplot2::ggplot(d,
        ggplot2::aes(x = factor(x_val), y = .data[[y_col]])) +
      ggplot2::geom_violin(fill = "#cfe4f5", colour = "#0b3d91",
                           alpha = 0.7, linewidth = 0.6) +
      ggplot2::geom_jitter(ggplot2::aes(colour = impact),
                           width = 0.15, alpha = 0.6, size = 1.5) +
      ggplot2::scale_colour_viridis_d(na.value = "grey60",
                                      name = "Functional impact") +
      ggrepel::geom_text_repel(
        data = tops,
        ggplot2::aes(label = variant_label),
        size = 4.5, max.overlaps = 10) +
      ggplot2::facet_wrap(~ display_name, labeller = .italic_gene_labeller) +
      ggplot2::labs(x = x_axis_title, y = ylabel) +
      ggplot2::theme_minimal(base_size = 16)
    ggplot2::ggsave(out_path, plot = p, width = 14, height = 10, dpi = 300)
  }

  # |β̃|: all cutpoints.
  impact_violin(
    df, "abs_median", expression("|" ~ tilde(beta) ~ "|"),
    top_per_group(df, "abs_median"),
    file.path(OUTPUT_DIR,
              paste0(run$name, "_gene_violin_abs_beta_impact.png")))

  # |Δβ̃|: transitions only (first cutpoint per variant has no delta).
  df_delta <- df[!is.na(df$delta_signed), ]
  impact_violin(
    df_delta, "delta_abs", expression("|" ~ Delta ~ tilde(beta) ~ "|"),
    top_per_group(df_delta, "delta_abs"),
    file.path(OUTPUT_DIR,
              paste0(run$name, "_gene_violin_abs_delta_beta_impact.png")))

  # Companion CSV of the filtered underlying data, mirroring the plot inputs.
  out_cols <- c("variant_id", "position", "gene", "display_name",
                "cutpoint", "cutpoint_MIC", "median", "delta_abs",
                "impact", "variant_label")
  out <- df[, out_cols]
  names(out)[names(out) == "variant_label"] <- "nucleotide_label"
  write.csv(out,
            file.path(OUTPUT_DIR, paste0(run$name, "_violin_data.csv")),
            row.names = FALSE)

  # Long-format rows for the combined overlay figure (one metric per row).
  rbind(
    data.frame(run = run$name, metric = "beta",
               x_val = df$x_val,       value = df$abs_median,
               gene = df$display_name, stringsAsFactors = FALSE),
    data.frame(run = run$name, metric = "delta",
               x_val = df_delta$x_val, value = df_delta$delta_abs,
               gene = df_delta$display_name, stringsAsFactors = FALSE)
  )
}

combined <- do.call(rbind, lapply(RUNS, process_run))

# -----------------------------------------------------------------------------
# Combined overlay: all 6 genes in one violin per cutpoint, coloured by gene,
# no variant labels. Four facets = the four existing plots (run × metric).
run_label <- c(
  "02_spn_penicillin_MIC_PPOM"               = "MIC PPOM",
  "16_spn_penicillin_MIC_minimabinning_PPOM" = "minima-binning PPOM"
)
facet_expr <- ifelse(
  combined$metric == "beta",
  paste0("plain('", run_label[combined$run], "')~~'|'*tilde(beta)*'|'"),
  paste0("plain('", run_label[combined$run], "')~~'|'*Delta*tilde(beta)*'|'"))
facet_levels <- c(
  "plain('MIC PPOM')~~'|'*tilde(beta)*'|'",
  "plain('MIC PPOM')~~'|'*Delta*tilde(beta)*'|'",
  "plain('minima-binning PPOM')~~'|'*tilde(beta)*'|'",
  "plain('minima-binning PPOM')~~'|'*Delta*tilde(beta)*'|'")
combined$facet <- factor(facet_expr, levels = facet_levels)

p_overlay <- ggplot2::ggplot(combined,
    ggplot2::aes(x = factor(x_val), y = value)) +
  ggplot2::geom_violin(fill = "#cfe4f5", colour = "#0b3d91",
                       alpha = 0.7, linewidth = 0.6) +
  ggplot2::geom_jitter(ggplot2::aes(colour = gene),
                       width = 0.15, alpha = 0.6, size = 1) +
  ggplot2::scale_colour_viridis_d(name = "gene") +
  ggplot2::facet_wrap(~ facet, scales = "free", ncol = 2,
                      labeller = ggplot2::label_parsed) +
  ggplot2::labs(x = expression("MIC breakpoint" ~ (mu * "g·mL"^{-1})), y = NULL) +
  ggplot2::theme_minimal(base_size = 16)
ggplot2::ggsave(
  file.path(OUTPUT_DIR, "overlay_gene_violin_grid.png"),
  plot = p_overlay, width = 14, height = 10, dpi = 300)

message("Done. Output written to ", OUTPUT_DIR)
