############################################################
## compare_grm_vs_subclusters.R
##
## Compares two binary GWAS runs that differ only in how population structure
## is controlled: a lineage/subcluster hierarchy versus a GRM random effect.
##
##   spn_penicillin  01_spn_penicillin_binary_logistic  vs 19_..._grm_logistic
##   tb_rifampicin   07_tb_rifampicin_binary_logistic   vs 20_..._grm_logistic
##
## Both runs in a pair use the same phenotype, the same genotype matrix and the
## same pipeline settings, so their outputs line up row for row on the original
## variant index. This script checks that agreement for variant effects (median
## beta with 89% CI significance) and for RATE values, and writes scatter plots
## plus a summary table.
##
## Usage:
##   Rscript compare_grm_vs_subclusters.R [spn_penicillin|tb_rifampicin]
## Defaults to spn_penicillin.
############################################################

suppressPackageStartupMessages({
  library(ggplot2)
})

## ---------------------------- paths -------------------------------------- ##

RESULTS_BASE <- "/nfs/research/jlees/jacqueline/thesis_results"
DATASETS_DIR <- file.path(RESULTS_BASE, "gwas_datasets", "inference")

COMPARISONS <- list(
  spn_penicillin = list(
    gwas_dir     = "gwas_spn_penicillin",
    run_sub      = "01_spn_penicillin_binary_logistic",
    run_grm      = "19_spn_penicillin_binary_grm_logistic",
    dataset_sub  = "01_spn_penicillin_binary",
    dataset_grm  = "19_spn_penicillin_binary_grm",
    annotations  = "/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/genotype/fields_filtered_maf05_multiallelic.txt",
    genes        = "/nfs/research/jlees/jacqueline/thesis_code/gwas_genesofinterest/spn_penicillin_genesofinterest.txt"
  ),
  tb_rifampicin = list(
    gwas_dir     = "gwas_tb_rifampicin",
    run_sub      = "07_tb_rifampicin_binary_logistic",
    run_grm      = "20_tb_rifampicin_binary_grm_logistic",
    dataset_sub  = "07_tb_rifampicin_binary",
    dataset_grm  = "20_tb_rifampicin_binary_grm",
    annotations  = "/nfs/research/jlees/jacqueline/gwas_data/tuberculosis/cryptic_regeno_snpeff/fields_filtered.txt",
    genes        = "/nfs/research/jlees/jacqueline/thesis_code/gwas_genesofinterest/tb_rifampicin_genesofinterest.txt"
  )
)

cli <- commandArgs(trailingOnly = TRUE)
COMPARISON <- if (length(cli) >= 1) cli[1] else "spn_penicillin"
stopifnot(COMPARISON %in% names(COMPARISONS))
cfg <- COMPARISONS[[COMPARISON]]
message("Comparison: ", COMPARISON)

GWAS_DIR <- file.path(RESULTS_BASE, cfg$gwas_dir, "inference")
RUN_SUB  <- file.path(GWAS_DIR, cfg$run_sub)
RUN_GRM  <- file.path(GWAS_DIR, cfg$run_grm)

VARIANT_INDEX_SUB <- file.path(DATASETS_DIR, cfg$dataset_sub,
                               paste0(cfg$dataset_sub, "_variant_index.csv"))
VARIANT_INDEX_GRM <- file.path(DATASETS_DIR, cfg$dataset_grm,
                               paste0(cfg$dataset_grm, "_variant_index.csv"))

ANNOTATIONS <- cfg$annotations
GENES_OF_INTEREST <- cfg$genes

OUT_DIR <- file.path(RESULTS_BASE, "compare_grm_subclusters", COMPARISON)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

LABEL_SUB <- "subclusters"
LABEL_GRM <- "GRM"

## ---------------------------- readers ------------------------------------ ##

# Median beta and 89% CI significance, one row per variant on the original
# (pre-pruning) variant index.
read_effects <- function(run_dir) {
  path <- file.path(run_dir, "fitted_model", "depruned_variant_effects.csv")
  stopifnot(file.exists(path))
  df <- read.csv(path)
  stopifnot(all(c("variant_id", "median", "signif") %in% names(df)))
  df$signif <- as.logical(df$signif)
  df
}

# RATE values. The cppRATE output carries three comment lines (#ESS, #Delta and
# the column header) before the body; this is the same idiom the pipeline uses
# in gwas_workflow/R/pipeline.R:.load_rates_from_disk().
read_rates <- function(run_dir) {
  path <- file.path(run_dir, "cppRATE_results", "RATE_values_depruned.txt")
  stopifnot(file.exists(path))
  body <- read.delim(path, sep = "", header = FALSE)[-c(1:3), ]
  data.frame(variant_id = as.integer(body[, 1]),
             rate       = as.numeric(body[, 2]),
             kld        = as.numeric(body[, 3]))
}

# Variant IDs kept by LD pruning (the header line of the BacPrune output).
read_kept_variants <- function(run_dir) {
  path <- file.path(run_dir, "cppRATE_matrices", "bacprune_rust_results.csv")
  if (!file.exists(path)) return(NULL)
  as.integer(sub("^X", "", strsplit(readLines(path, n = 1), ",")[[1]]))
}

italic_gene_expr <- function(x) paste0("italic('", gsub("'", "", x), "')")

## ---------------------------- load ---------------------------------------- ##

message("Reading run outputs...")

vi_sub <- read.csv(VARIANT_INDEX_SUB)
vi_grm <- read.csv(VARIANT_INDEX_GRM)
stopifnot(identical(vi_sub, vi_grm))
message("  Variant indices identical: ", nrow(vi_sub), " variants")

eff_sub <- read_effects(RUN_SUB)
eff_grm <- read_effects(RUN_GRM)
rate_sub <- read_rates(RUN_SUB)
rate_grm <- read_rates(RUN_GRM)

stopifnot(nrow(eff_sub) == nrow(vi_sub), nrow(eff_grm) == nrow(vi_sub),
          nrow(rate_sub) == nrow(vi_sub), nrow(rate_grm) == nrow(vi_sub))
stopifnot(identical(eff_sub$variant_id, eff_grm$variant_id),
          identical(rate_sub$variant_id, rate_grm$variant_id),
          identical(eff_sub$variant_id, rate_sub$variant_id))
message("  Effects and RATE files aligned on the same variant order")

# LD pruning should keep the same variants in both runs (same genotype, same
# threshold). De-pruned values sit on the full index either way, so a mismatch
# is reported rather than fatal.
kept_sub <- read_kept_variants(RUN_SUB)
kept_grm <- read_kept_variants(RUN_GRM)
pruning_identical <- NA
if (!is.null(kept_sub) && !is.null(kept_grm)) {
  pruning_identical <- identical(sort(kept_sub), sort(kept_grm))
  message("  LD pruning kept ", length(kept_sub), " (", LABEL_SUB, ") vs ",
          length(kept_grm), " (", LABEL_GRM, ") variants; identical: ",
          pruning_identical)
  if (!isTRUE(pruning_identical)) {
    warning("LD pruning kept different variant sets in the two runs")
  }
}

## ---------------------------- annotate ------------------------------------ ##

ann <- read.delim(ANNOTATIONS, header = TRUE, check.names = TRUE)
pos_col  <- names(ann)[1]
gene_col <- grep("GENE", names(ann), value = TRUE)[1]
hit <- match(vi_sub$position, ann[[pos_col]])
genes <- ann[[gene_col]][hit]
genes[is.na(genes) | genes == ""] <- NA_character_

goi <- read.csv(GENES_OF_INTEREST, header = FALSE,
                strip.white = TRUE, stringsAsFactors = FALSE)
display_map <- stats::setNames(goi[[2]], goi[[1]])
in_goi <- !is.na(genes) & genes %in% names(display_map)
gene_label <- genes
gene_label[in_goi] <- unname(display_map[genes[in_goi]])

## ---------------------------- assemble ------------------------------------ ##

df <- data.frame(
  variant_id   = eff_sub$variant_id,
  variant_name = vi_sub$variant_name,
  position     = vi_sub$position,
  gene         = genes,
  gene_label   = gene_label,
  gene_of_interest = in_goi,
  beta_sub     = eff_sub$median,
  beta_grm     = eff_grm$median,
  signif_sub   = eff_sub$signif,
  signif_grm   = eff_grm$signif,
  rate_sub     = rate_sub$rate,
  rate_grm     = rate_grm$rate,
  stringsAsFactors = FALSE
)

df$agreement <- with(df, ifelse(
  signif_sub & signif_grm, "both",
  ifelse(signif_sub & !signif_grm, paste0(LABEL_SUB, " only"),
         ifelse(!signif_sub & signif_grm, paste0(LABEL_GRM, " only"), "neither"))))
df$agreement <- factor(df$agreement,
                       levels = c("neither", paste0(LABEL_SUB, " only"),
                                  paste0(LABEL_GRM, " only"), "both"))

write.csv(df, file.path(OUT_DIR, "comparison_per_variant.csv"), row.names = FALSE)

## ---------------------------- statistics ---------------------------------- ##

top_overlap <- function(a, b, n) {
  ta <- order(a, decreasing = TRUE)[seq_len(n)]
  tb <- order(b, decreasing = TRUE)[seq_len(n)]
  length(intersect(ta, tb)) / length(union(ta, tb))
}

eps <- min(df$rate_sub[df$rate_sub > 0], df$rate_grm[df$rate_grm > 0]) / 10

stats <- data.frame(
  metric = c("n_variants",
             "beta_pearson", "beta_spearman",
             "abs_beta_pearson", "abs_beta_spearman",
             "beta_pearson_signif_either",
             "rate_spearman", "log10_rate_pearson",
             "rate_top50_jaccard", "rate_top100_jaccard",
             "n_signif_subclusters", "n_signif_grm",
             "n_signif_both", "n_signif_subclusters_only",
             "n_signif_grm_only",
             "ld_pruning_kept_identical"),
  value = c(
    nrow(df),
    cor(df$beta_sub, df$beta_grm, method = "pearson"),
    cor(df$beta_sub, df$beta_grm, method = "spearman"),
    cor(abs(df$beta_sub), abs(df$beta_grm), method = "pearson"),
    cor(abs(df$beta_sub), abs(df$beta_grm), method = "spearman"),
    {
      s <- df$signif_sub | df$signif_grm
      if (sum(s) > 2) cor(df$beta_sub[s], df$beta_grm[s], method = "pearson") else NA_real_
    },
    cor(df$rate_sub, df$rate_grm, method = "spearman"),
    cor(log10(df$rate_sub + eps), log10(df$rate_grm + eps), method = "pearson"),
    top_overlap(df$rate_sub, df$rate_grm, 50),
    top_overlap(df$rate_sub, df$rate_grm, 100),
    sum(df$signif_sub),
    sum(df$signif_grm),
    sum(df$signif_sub & df$signif_grm),
    sum(df$signif_sub & !df$signif_grm),
    sum(!df$signif_sub & df$signif_grm),
    as.numeric(pruning_identical)
  )
)
write.csv(stats, file.path(OUT_DIR, "comparison_summary.csv"), row.names = FALSE)
message("\nSummary:")
print(stats, row.names = FALSE)

## ---------------------------- plots --------------------------------------- ##

# cppRATE returns a distribution summing to 1 over the variants it was run on,
# which is the LD-pruned set, not the full index (de-pruning copies
# representative values onto pruned partners, so the plotted values sum to
# slightly more than 1). 1/V on the pruned count is therefore the uniform
# share: a variant above this line carries more than its even split of the
# total RATE.
n_rate_variants <- if (!is.null(kept_sub)) length(kept_sub) else nrow(df)
RATE_UNIFORM <- 1 / n_rate_variants
message(sprintf("  RATE uniform reference: 1/%d = %.3g",
                n_rate_variants, RATE_UNIFORM))

# Both quantities span many orders of magnitude (RATE 4e-16 to 0.4, |beta|
# 2e-13 to 1.9), so plots use log10 axes; values below these floors are
# squished onto the axis rather than dropped, since that tail is horseshoe
# shrinkage noise rather than signal.
FLOOR_RATE <- 1e-9
FLOOR_BETA <- 1e-10

base_theme <- theme_minimal(base_size = 14)
agreement_colours <- stats::setNames(
  c("grey70", "#E76254", "#376795", "#1E466E"),
  levels(df$agreement)
)

# Label the genes of interest that are significant in either run, one point
# (largest |beta|) per gene.
label_candidates <- df[df$gene_of_interest & (df$signif_sub | df$signif_grm), ]
if (nrow(label_candidates) > 0) {
  label_candidates <- label_candidates[
    order(pmax(abs(label_candidates$beta_sub), abs(label_candidates$beta_grm)),
          decreasing = TRUE), ]
  label_df <- label_candidates[!duplicated(label_candidates$gene_label), ]
  label_df$gene_expr <- italic_gene_expr(label_df$gene_label)
} else {
  label_df <- NULL
}

# Effect magnitudes. Sign is dropped so the comparison is about how much
# weight each model gives a variant, not the direction of the allele effect.
df$abs_beta_sub <- abs(df$beta_sub)
df$abs_beta_grm <- abs(df$beta_grm)
beta_lim <- c(FLOOR_BETA, max(c(df$abs_beta_sub, df$abs_beta_grm)))

p_beta <- ggplot(df, aes(x = abs_beta_sub, y = abs_beta_grm, colour = agreement)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
  geom_point(alpha = 0.5, size = 1) +
  scale_colour_manual(values = agreement_colours) +
  scale_x_log10(limits = beta_lim, oob = scales::squish,
                labels = scales::label_log()) +
  scale_y_log10(limits = beta_lim, oob = scales::squish,
                labels = scales::label_log()) +
  coord_fixed(ratio = 1) +
  labs(x = bquote(group("|", tilde(beta), "|") ~ "(" * .(LABEL_SUB) * ")"),
       y = bquote(group("|", tilde(beta), "|") ~ "(" * .(LABEL_GRM) * ")"),
       colour = "89% CI excludes 0",
       title = "Variant effect magnitudes: subclusters vs GRM",
       subtitle = sprintf(
         "log10 axes, dashed line y = x\n|beta| Pearson r = %.3f, Spearman rho = %.3f",
         cor(df$abs_beta_sub, df$abs_beta_grm),
         cor(df$abs_beta_sub, df$abs_beta_grm, method = "spearman"))) +
  base_theme

if (!is.null(label_df)) {
  p_beta <- p_beta +
    ggrepel::geom_text_repel(
      data = label_df, aes(x = abs(beta_sub), y = abs(beta_grm), label = gene_expr),
      parse = TRUE, size = 3, colour = "black", inherit.aes = FALSE,
      max.overlaps = 30,
      arrow = grid::arrow(length = grid::unit(0.01, "npc"), type = "open"))
}

ggsave(file.path(OUT_DIR, "beta_scatter.png"), p_beta,
       width = 9, height = 8, dpi = 300)

p_rate <- ggplot(df, aes(x = rate_sub + eps, y = rate_grm + eps,
                         colour = agreement)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
  # Uniform-share reference in each model; the upper-right quadrant is the set
  # of variants carrying more than 1/V of the total RATE under both.
  geom_hline(yintercept = RATE_UNIFORM, linetype = "dotted", colour = "grey25") +
  geom_vline(xintercept = RATE_UNIFORM, linetype = "dotted", colour = "grey25") +
  geom_point(alpha = 0.5, size = 1) +
  scale_colour_manual(values = agreement_colours) +
  scale_x_log10() + scale_y_log10() +
  labs(x = paste0("RATE (", LABEL_SUB, ")"),
       y = paste0("RATE (", LABEL_GRM, ")"),
       colour = "89% CI excludes 0",
       title = "RATE values: subcluster vs GRM population-structure control",
       subtitle = sprintf("Spearman rho = %.3f, top-50 Jaccard = %.2f",
                          cor(df$rate_sub, df$rate_grm, method = "spearman"),
                          top_overlap(df$rate_sub, df$rate_grm, 50))) +
  base_theme

if (!is.null(label_df)) {
  p_rate <- p_rate +
    ggrepel::geom_text_repel(
      data = label_df, aes(x = rate_sub + eps, y = rate_grm + eps,
                           label = gene_expr),
      parse = TRUE, size = 3, colour = "black", inherit.aes = FALSE,
      max.overlaps = 30,
      arrow = grid::arrow(length = grid::unit(0.01, "npc"), type = "open"))
}

ggsave(file.path(OUT_DIR, "rate_scatter.png"), p_rate,
       width = 9, height = 8, dpi = 300)

# Overlaid Manhattans as a visual sanity check.
#
long_beta <- rbind(
  data.frame(position = df$position, value = abs(df$beta_sub), model = LABEL_SUB),
  data.frame(position = df$position, value = abs(df$beta_grm), model = LABEL_GRM)
)
long_rate <- rbind(
  data.frame(position = df$position, value = df$rate_sub, model = LABEL_SUB),
  data.frame(position = df$position, value = df$rate_grm, model = LABEL_GRM)
)

# Interleave the two series so neither model systematically overplots the other
set.seed(1)
long_beta <- long_beta[sample.int(nrow(long_beta)), ]
long_rate <- long_rate[sample.int(nrow(long_rate)), ]

model_colours <- stats::setNames(c("#376795", "#E76254"), c(LABEL_SUB, LABEL_GRM))

# The scale of each panel is carried by its y-axis label rather than a
# subtitle, so log panels plot log10(value) directly on a linear axis. Values
# below the floor are clamped onto it; that tail is horseshoe shrinkage noise
# rather than signal.
manhattan_panel <- function(dat, ylab, floor_val = NULL, href = NULL,
                            href_label = NULL) {
  if (!is.null(floor_val)) {
    dat$value <- log10(pmax(dat$value, floor_val))
    if (!is.null(href)) href <- log10(href)
  }
  p <- ggplot(dat, aes(x = position, y = value, colour = model)) +
    geom_point(alpha = 0.35, size = 0.7)
  if (!is.null(href)) {
    # Label sits on a filled box: plain text is unreadable against the cloud
    p <- p +
      geom_hline(yintercept = href, linetype = "dashed", colour = "grey15") +
      ggplot2::annotate("label", x = min(dat$position), y = href,
                        label = href_label, hjust = 0, vjust = -0.25,
                        size = 3.6, colour = "grey15", fill = "white",
                        alpha = 0.85, label.size = 0)
  }
  p +
    scale_colour_manual(values = model_colours) +
    scale_x_continuous(labels = scales::label_number(scale = 1e-6, suffix = " Mb")) +
    guides(colour = ggplot2::guide_legend(override.aes = list(alpha = 1, size = 3))) +
    labs(x = "genome coordinate", y = ylab, colour = NULL) +
    base_theme
}

p_man_beta_lin <- manhattan_panel(
  long_beta, expression(group("|", tilde(beta), "|")))
p_man_beta <- manhattan_panel(
  long_beta, expression(log[10] * group("|", tilde(beta), "|")), FLOOR_BETA)
p_man_rate <- manhattan_panel(
  long_rate, expression(log[10] * "RATE"), FLOOR_RATE,
  href = RATE_UNIFORM,
  href_label = sprintf("1/V = %.2g", RATE_UNIFORM))

ggsave(file.path(OUT_DIR, "manhattan_overlay.png"),
       patchwork::wrap_plots(p_man_beta_lin, p_man_beta, p_man_rate, ncol = 1),
       width = 12, height = 13, dpi = 300)

# Same panels split by model, for when the overlay is still too dense
ggsave(file.path(OUT_DIR, "manhattan_overlay_faceted.png"),
       patchwork::wrap_plots(
         p_man_beta + ggplot2::facet_wrap(~model, ncol = 1),
         p_man_rate + ggplot2::facet_wrap(~model, ncol = 1),
         ncol = 2),
       width = 16, height = 10, dpi = 300)

# Linear-axis version kept for reference
p_lin_beta <- ggplot(long_beta, aes(x = position, y = value, colour = model)) +
  geom_point(alpha = 0.4, size = 0.8) +
  scale_colour_manual(values = model_colours) +
  labs(x = "genome coordinate (bp)", y = expression(group("|", tilde(beta), "|")),
       colour = NULL, title = "Absolute variant effects (linear axis)") +
  base_theme
p_lin_rate <- ggplot(long_rate, aes(x = position, y = value, colour = model)) +
  geom_point(alpha = 0.4, size = 0.8) +
  scale_colour_manual(values = model_colours) +
  labs(x = "genome coordinate (bp)", y = "RATE", colour = NULL,
       title = "RATE (linear axis)") +
  base_theme

ggsave(file.path(OUT_DIR, "manhattan_overlay_linear.png"),
       patchwork::wrap_plots(p_lin_beta, p_lin_rate, ncol = 1),
       width = 12, height = 9, dpi = 300)

## ---------------------------- heritability -------------------------------- ##

h2_path_sub <- file.path(RUN_SUB, "plots", "heritability", "heritability_summary.csv")
h2_path_grm <- file.path(RUN_GRM, "plots", "heritability", "heritability_summary.csv")
if (file.exists(h2_path_sub) && file.exists(h2_path_grm)) {
  h2_sub <- read.csv(h2_path_sub); h2_sub$model <- LABEL_SUB
  h2_grm <- read.csv(h2_path_grm); h2_grm$model <- LABEL_GRM
  h2 <- rbind(h2_sub, h2_grm)
  write.csv(h2, file.path(OUT_DIR, "heritability_comparison.csv"), row.names = FALSE)
  message("\nHeritability:")
  print(h2[, c("model", "metric", "mean", "median", "q_lower", "q_upper")],
        row.names = FALSE)
}

message("\nWrote outputs to: ", OUT_DIR)
