############################################################
## compare_ld_pruning.R
##
## Compares binary GWAS runs that differ only in how variants are LD pruned
## before fitting. Three arms, all the same dataset, the same logistic model and
## the same pipeline settings, so their outputs line up row for row on the
## original variant index:
##
##   none    no pruning                 all 75,272 variants enter the model
##   dprime  |D'| >= 1 pruning          BacPrune --dprime
##   r2      r-squared >= 1 pruning     BacPrune --r, the as-published run
##
## At threshold 1 the two measures are not interchangeable. r2 >= 1 prunes only
## perfectly correlated variants, i.e. exact duplicates. |D'| >= 1 also prunes
## nested variants, where one of the four haplotypes is simply absent, so it
## prunes strictly more and can remove variants that are only weakly correlated
## with their representative.
##
## That difference is what the figures are for. Whatever a pruned arm reports
## for a pruned variant is copied from its representative, not estimated, so
## the more aggressive the pruning, the more of the genome is spoken for by
## proxy. The no-pruning arm is the reference the other two are read against.
##
## Figures, each produced for both quantities:
##   manhattan_faceted_<beta|rate>.png     one row per arm, shared axes
##   scatter_<beta|rate>_vs_nopruning.png  each pruned arm against no pruning
##
## INPUTS. Every arm is read from compare_ld_pruning/depruned/<arm>.csv, written
## by run_deprune.sh. The pipeline's own de-pruned outputs are not used: they
## assume every representative is itself a kept variant, which fails badly on
## the |D'| arm. deprune_from_directions.R has the details.
##
## Usage:
##   Rscript compare_ld_pruning.R [tb_rifampicin]
## Defaults to tb_rifampicin.
## Prerequisite: sbatch run_deprune.sh
############################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
})
grDevices::pdf(NULL)

CODE_BASE <- "/nfs/research/jlees/jacqueline/thesis_code"
source(file.path(CODE_BASE, "paper_figures", "gene_label_helpers.R"))

## ---------------------------- paths -------------------------------------- ##

RESULTS_BASE <- "/nfs/research/jlees/jacqueline/thesis_results"
DATASETS_DIR <- file.path(RESULTS_BASE, "gwas_datasets", "inference")
COMPARE_DIR  <- file.path(RESULTS_BASE, "compare_ld_pruning")

COMPARISONS <- list(
  tb_rifampicin = list(
    dataset     = "07_tb_rifampicin_binary",
    annotations = "/nfs/research/jlees/jacqueline/gwas_data/tuberculosis/cryptic_regeno_snpeff/fields_filtered.txt",
    genes       = file.path(CODE_BASE, "gwas_genesofinterest",
                            "tb_rifampicin_genesofinterest.txt"),
    on_target   = c("rpoB", "rpoC", "rpoA"),
    arms = list(
      none = list(
        label  = "no pruning",
        pruned = FALSE,
        dir    = file.path(COMPARE_DIR,
                           "07_tb_rifampicin_binary_logistic_nopruning")),
      dprime = list(
        label  = "|D'| pruning",
        pruned = TRUE,
        dir    = file.path(COMPARE_DIR,
                           "07_tb_rifampicin_binary_logistic_dprime")),
      r2 = list(
        label  = "r² pruning",
        pruned = TRUE,
        dir    = file.path(RESULTS_BASE, "gwas_tb_rifampicin", "inference",
                           "07_tb_rifampicin_binary_logistic"))
    ),
    # Written by run_deprune.sh, one file per arm key.
    depruned_dir = file.path(COMPARE_DIR, "depruned")
  )
)

cli <- commandArgs(trailingOnly = TRUE)
COMPARISON <- if (length(cli) >= 1) cli[1] else "tb_rifampicin"
stopifnot(COMPARISON %in% names(COMPARISONS))
cfg <- COMPARISONS[[COMPARISON]]
message("Comparison: ", COMPARISON)

VARIANT_INDEX <- file.path(DATASETS_DIR, cfg$dataset,
                           paste0(cfg$dataset, "_variant_index.csv"))
OUT_DIR <- file.path(COMPARE_DIR, COMPARISON)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

ARM_KEYS   <- names(cfg$arms)
REF_KEY    <- "none"                       # the arm the others are read against
ARM_LABELS <- stats::setNames(vapply(cfg$arms, `[[`, character(1), "label"),
                              ARM_KEYS)

# Colour follows the arm. Okabe-Ito blue / vermillion / bluish-green, the same
# trio validated for compare_lineage_subclusters against the dataviz palette
# checks (lightness band, chroma floor, all-pairs CVD separation, normal-vision
# floor, contrast on white). The house Hiroshige picks used elsewhere in
# thesis_code fail the lightness band and chroma floor.
ARM_COLOURS <- stats::setNames(c("#0072B2", "#D55E00", "#009E73"), ARM_KEYS)
COL_NULL    <- "#9aa0a6"   # variants whose 89% CI covers zero

## ---------------------------- readers ------------------------------------ ##

# Median beta, 89% CI significance and cppRATE, one row per variant on the
# original (pre-pruning) variant index. All three arms read the same table
# produced by run_deprune.sh.
#
# None of the pipeline's own de-pruned outputs are used, for two reasons.
# depruned_variant_effects.csv holds the STANDARDIZED effect (cpprate.R matches
# beta_variant_std ahead of beta_variant), so mixing it in would compare beta
# against beta * sd(x_v). More seriously, both it and RATE_values_depruned.txt
# assume every representative is itself a kept variant, which is true under r2
# but not under |D'|: see deprune_from_directions.R.
arm_table_path <- function(key) file.path(cfg$depruned_dir, paste0(key, ".csv"))

read_arm_table <- function(key) {
  path <- arm_table_path(key)
  stopifnot(file.exists(path))
  df <- data.table::fread(path)
  stopifnot(all(c("variant_id", "median", "signif", "rate") %in% names(df)))
  df$signif <- as.logical(df$signif)
  df
}

# How many variants the arm actually fitted, as opposed to inherited from a
# representative. Taken from the de-pruned table's chain_depth, where 0 marks a
# variant that was fitted in its own right.
n_kept_variants <- function(tbl) sum(tbl$chain_depth == 0L)

read_h2 <- function(run_dir) {
  path <- file.path(run_dir, "plots", "heritability", "heritability_summary.csv")
  if (!file.exists(path)) return(NULL)
  read.csv(path)
}

## ---------------------------- load ---------------------------------------- ##

missing <- !file.exists(vapply(ARM_KEYS, arm_table_path, character(1)))
if (any(missing)) {
  stop("Missing de-pruned table for: ",
       paste(ARM_LABELS[missing], collapse = ", "),
       "\n  expected:\n    ",
       paste(vapply(ARM_KEYS[missing], arm_table_path, character(1)),
             collapse = "\n    "),
       "\n  Run compare_ld_pruning/run_deprune.sh once the fits finish.")
}

message("Reading run outputs...")
vi <- read.csv(VARIANT_INDEX)
n_variants <- nrow(vi)
message("  Variant index: ", n_variants, " variants")

eff <- lapply(ARM_KEYS, read_arm_table)
names(eff) <- ARM_KEYS
for (k in ARM_KEYS) stopifnot(nrow(eff[[k]]) == n_variants)
for (k in ARM_KEYS[-1])
  stopifnot(identical(eff[[k]]$variant_id, eff[[ARM_KEYS[1]]]$variant_id))
message("  Arms aligned on the same variant order across all ",
        length(eff), " arms")

has_rate <- vapply(eff, function(e) any(!is.na(e$rate)), logical(1))
if (!all(has_rate))
  message("  No RATE for: ", paste(ARM_LABELS[!has_rate], collapse = ", "),
          " -- RATE figures will be skipped")

## ---------------------------- annotate ------------------------------------ ##

ann_df <- annotate_genes(data.frame(pos = vi$position), cfg$annotations, cfg$genes)

df <- data.frame(
  variant_id   = eff[[1]]$variant_id,
  variant_name = vi$variant_name,
  pos          = vi$position,
  gene         = ann_df$gene,
  stringsAsFactors = FALSE
)
df$on_target <- df$gene %in% cfg$on_target
for (k in ARM_KEYS) {
  df[[paste0("beta_",   k)]] <- eff[[k]]$median
  df[[paste0("signif_", k)]] <- eff[[k]]$signif
  if (has_rate[[k]]) df[[paste0("rate_", k)]] <- eff[[k]]$rate
}
write.csv(df, file.path(OUT_DIR, "comparison_per_variant.csv"), row.names = FALSE)
message("  On-target variants (", paste(cfg$on_target, collapse = "/"), "): ",
        sum(df$on_target), " of ", nrow(df))

## ---------------------------- statistics ---------------------------------- ##

h2 <- do.call(rbind, lapply(ARM_KEYS, function(k) {
  h <- read_h2(cfg$arms[[k]]$dir)
  if (is.null(h)) return(NULL)
  h$model <- unname(ARM_LABELS[k])
  h
}))
if (!is.null(h2))
  write.csv(h2, file.path(OUT_DIR, "heritability_comparison.csv"), row.names = FALSE)

h2_get <- function(label, metric) {
  if (is.null(h2)) return(NA_real_)
  v <- h2$median[h2$model == label & h2$metric == metric]
  if (length(v) == 1) v else NA_real_
}

arm_summary <- function(k) {
  lab  <- unname(ARM_LABELS[k])
  s    <- df[[paste0("signif_", k)]]
  kept <- n_kept_variants(eff[[k]])
  n_on <- sum(s & df$on_target, na.rm = TRUE)
  n_off <- sum(s & !df$on_target, na.rm = TRUE)
  data.frame(model = lab,
             n_variants_fitted = kept,
             n_variants_inherited = n_variants - kept,
             frac_inherited = (n_variants - kept) / n_variants,
             max_chain_depth = max(eff[[k]]$chain_depth),
             n_signif = sum(s, na.rm = TRUE),
             n_signif_on_target = n_on,
             n_signif_off_target = n_off,
             frac_signif_off_target =
               if (sum(s, na.rm = TRUE) > 0) n_off / sum(s, na.rm = TRUE) else NA_real_,
             h2_narrow = h2_get(lab, "h2_narrow"),
             h2_broad  = h2_get(lab, "h2_broad"),
             stringsAsFactors = FALSE)
}

all_summary <- do.call(rbind, lapply(ARM_KEYS, arm_summary))
write.csv(all_summary, file.path(OUT_DIR, "comparison_summary.csv"),
          row.names = FALSE)
message("\nAll-arm summary:")
print(all_summary, row.names = FALSE)

# Each pruned arm against the no-pruning reference. Spearman is reported for
# RATE as well as Pearson because RATE is normalised over the variants an arm
# actually fitted, and the arms fitted different numbers of them, so the two
# scales are not expected to coincide even where the rankings agree.
pairwise <- do.call(rbind, lapply(setdiff(ARM_KEYS, REF_KEY), function(k) {
  ba <- df[[paste0("beta_", REF_KEY)]]; bb <- df[[paste0("beta_", k)]]
  sa <- df[[paste0("signif_", REF_KEY)]]; sb <- df[[paste0("signif_", k)]]
  ra <- df[[paste0("rate_", REF_KEY)]];  rb <- df[[paste0("rate_", k)]]
  ok_r <- !is.null(ra) && !is.null(rb)
  data.frame(
    reference = unname(ARM_LABELS[REF_KEY]),
    model     = unname(ARM_LABELS[k]),
    beta_pearson     = cor(ba, bb, method = "pearson",  use = "complete.obs"),
    beta_spearman    = cor(ba, bb, method = "spearman", use = "complete.obs"),
    abs_beta_pearson = cor(abs(ba), abs(bb), method = "pearson", use = "complete.obs"),
    rate_pearson  = if (ok_r) cor(ra, rb, method = "pearson",  use = "complete.obs") else NA_real_,
    rate_spearman = if (ok_r) cor(ra, rb, method = "spearman", use = "complete.obs") else NA_real_,
    n_signif_both      = sum(sa & sb, na.rm = TRUE),
    n_signif_ref_only  = sum(sa & !sb, na.rm = TRUE),
    n_signif_arm_only  = sum(!sa & sb, na.rm = TRUE),
    stringsAsFactors = FALSE)
}))
write.csv(pairwise, file.path(OUT_DIR, "comparison_pairwise.csv"), row.names = FALSE)
message("\nAgainst the no-pruning reference:")
print(pairwise, row.names = FALSE)

## ---------------------------- figures -------------------------------------- ##

# Everything is +4pt over the paper_figures base of 14. geom_text/geom_text_repel
# take size in mm rather than points, so those go up by 4/.pt (~1.4) instead.
BASE_SIZE <- 18
PT_BUMP   <- 4 / .pt
base_theme <- theme_minimal(base_size = BASE_SIZE)

ARM_ORDER   <- unname(ARM_LABELS[ARM_KEYS])
FILL_VALUES <- stats::setNames(unname(ARM_COLOURS[ARM_KEYS]), ARM_ORDER)

# Long form: one row per variant per arm, for the faceted panels.
long_for <- function(value_prefix) {
  keys <- ARM_KEYS[vapply(ARM_KEYS, function(k)
    !is.null(df[[paste0(value_prefix, k)]]), logical(1))]
  out <- do.call(rbind, lapply(keys, function(k)
    data.frame(pos       = df$pos,
               value     = df[[paste0(value_prefix, k)]],
               signif    = df[[paste0("signif_", k)]],
               gene      = df$gene,
               on_target = df$on_target,
               model     = unname(ARM_LABELS[k]),
               stringsAsFactors = FALSE)))
  out$model <- factor(out$model, levels = ARM_ORDER)
  out
}

# Peak-|value| row per gene, carrying the value through so the label sits on the
# point actually plotted. Split by arm first: select_label_rows() returns one row
# per gene over whatever it is given, so running it on the stacked frame would
# name each gene on a single facet and leave the other two looking empty there.
peak_rows_for <- function(dat, genes) {
  per_arm <- lapply(split(dat, dat$model), function(d) {
    if (nrow(d) == 0L) return(NULL)
    d$score <- abs(d$value)
    select_label_rows(d, gene_col = "gene", score_col = "score",
                      mode = "gene_list", genes = genes)
  })
  do.call(rbind, per_arm[!vapply(per_arm, is.null, logical(1))])
}

# repel_label_layer() hardcodes max.overlaps = 30, which silently drops labels.
# Fine for three genes of interest, not acceptable on the exhaustive panel where
# a missing label would read as "no hit there".
label_layer <- function(label_df, size, nudge_y) {
  ggrepel::geom_text_repel(
    data = label_df,
    mapping = aes(x = .data[["pos"]], y = .data[["value"]],
                  label = .data[["gene_expr"]]),
    parse = TRUE, size = size, colour = "grey25", inherit.aes = FALSE,
    direction = "y", nudge_y = nudge_y, vjust = 0,
    min.segment.length = Inf, box.padding = 0.15, point.padding = 0.2,
    max.overlaps = Inf, seed = 1)
}

# One row per arm. Variants whose 89% CI covers zero are drawn in grey underneath
# so the arm colour marks signal rather than membership; the y scale is shared so
# the rows are read against each other rather than each against itself.
#
# On the RATE panels the colouring still comes from the beta CI, since RATE
# carries no interval of its own. Read it as "the variants this arm's effect
# estimates called significant", not as a statement about RATE.
manhattan_faceted <- function(dat, ylab, log_y = FALSE, label_genes = NULL) {
  null_df <- dat[!dat$signif %in% TRUE, , drop = FALSE]
  sig_df  <- dat[dat$signif %in% TRUE, , drop = FALSE]

  p <- ggplot(mapping = aes(x = pos, y = value)) +
    geom_point(data = null_df, colour = COL_NULL, alpha = 0.25, size = 0.9) +
    geom_point(data = sig_df, aes(colour = model), alpha = 0.75, size = 1.3) +
    scale_colour_manual(values = FILL_VALUES, guide = "none") +
    facet_wrap(~model, ncol = 1) +
    xlab("genome coordinate (bp)") +
    ylab(ylab) +
    base_theme +
    theme(panel.grid.minor = element_blank(),
          strip.text = element_text(hjust = 0))

  if (log_y) {
    n_zero <- sum(dat$value <= 0, na.rm = TRUE)
    if (n_zero > 0)
      message("    dropping ", n_zero, " points with value <= 0 for the log axis")
    p <- p + scale_y_log10()
  }
  if (is.null(label_genes)) return(p)

  label_df <- peak_rows_for(dat, label_genes)
  if (nrow(label_df) == 0L) return(p)
  message("    labelling ", nrow(label_df), " gene-arm peaks")
  span <- diff(range(dat$value, na.rm = TRUE))
  # Headroom for the labels, but only on a linear axis: adding a continuous
  # scale on a log panel would silently replace scale_y_log10().
  if (!log_y) p <- p + scale_y_continuous(expand = expansion(mult = c(0.05, 0.18)))
  p + label_layer(label_df, size = 5 + PT_BUMP, nudge_y = 0.04 * span)
}

# Each pruned arm against the no-pruning reference, one panel per arm. The dashed
# y = x line is the null: a point off it is a variant the pruned arm reports
# differently from the arm that actually estimated every variant.
#
# xlab_val / ylab_val are taken ready-made rather than pasted together from one
# label, because paste0() on an expression coerces it to its deparsed text and
# the axis then reads "tilde(beta), pruned".
scatter_vs_ref <- function(value_prefix, xlab_val, ylab_val, log_scale = FALSE) {
  keys <- setdiff(ARM_KEYS, REF_KEY)
  keys <- keys[vapply(keys, function(k)
    !is.null(df[[paste0(value_prefix, k)]]), logical(1))]
  if (length(keys) == 0L || is.null(df[[paste0(value_prefix, REF_KEY)]])) return(NULL)

  dat <- do.call(rbind, lapply(keys, function(k)
    data.frame(ref       = df[[paste0(value_prefix, REF_KEY)]],
               value     = df[[paste0(value_prefix, k)]],
               signif    = df[[paste0("signif_", k)]] %in% TRUE |
                           df[[paste0("signif_", REF_KEY)]] %in% TRUE,
               on_target = df$on_target,
               model     = unname(ARM_LABELS[k]),
               stringsAsFactors = FALSE)))
  dat$model <- factor(dat$model, levels = ARM_ORDER[ARM_ORDER %in% dat$model])

  if (log_scale) {
    n_drop <- sum(dat$ref <= 0 | dat$value <= 0, na.rm = TRUE)
    if (n_drop > 0)
      message("    dropping ", n_drop, " points with a non-positive value for the log axes")
    dat <- dat[dat$ref > 0 & dat$value > 0, , drop = FALSE]
  }
  dat <- dat[stats::complete.cases(dat$ref, dat$value), , drop = FALSE]

  # Correlations per panel, placed in the corner rather than left to a caption.
  ann <- do.call(rbind, lapply(levels(dat$model), function(m) {
    d <- dat[dat$model == m, ]
    data.frame(model = m,
               lab = sprintf("r == %.3f * ',' ~ rho == %.3f",
                             cor(d$ref, d$value, method = "pearson"),
                             cor(d$ref, d$value, method = "spearman")),
               stringsAsFactors = FALSE)
  }))
  ann$model <- factor(ann$model, levels = levels(dat$model))

  null_df <- dat[!dat$signif, , drop = FALSE]
  sig_df  <- dat[dat$signif, , drop = FALSE]

  p <- ggplot(mapping = aes(x = ref, y = value)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                colour = "grey40", linewidth = 0.5) +
    geom_point(data = null_df, colour = COL_NULL, alpha = 0.2, size = 0.9) +
    geom_point(data = sig_df, aes(colour = model), alpha = 0.7, size = 1.3) +
    scale_colour_manual(values = FILL_VALUES, guide = "none") +
    facet_wrap(~model, nrow = 1) +
    geom_text(data = ann, aes(label = lab), parse = TRUE,
              x = -Inf, y = Inf, hjust = -0.08, vjust = 1.5,
              size = 4.5 + PT_BUMP, colour = "grey25", inherit.aes = FALSE) +
    xlab(xlab_val) +
    ylab(ylab_val) +
    base_theme +
    theme(panel.grid.minor = element_blank(),
          strip.text = element_text(hjust = 0))

  if (log_scale) p + scale_x_log10() + scale_y_log10() else p
}

message("\nBuilding figures...")

beta_long <- long_for("beta_")
message("  manhattan_faceted_beta")
ggsave(file.path(OUT_DIR, "manhattan_faceted_beta.png"),
       manhattan_faceted(beta_long, expression(tilde(beta)),
                         label_genes = cfg$on_target),
       width = 16, height = 12, dpi = 300)
message("  manhattan_faceted_beta_unlabelled")
ggsave(file.path(OUT_DIR, "manhattan_faceted_beta_unlabelled.png"),
       manhattan_faceted(beta_long, expression(tilde(beta))),
       width = 16, height = 12, dpi = 300)

message("  scatter_beta_vs_nopruning")
p <- scatter_vs_ref("beta_",
                    expression(tilde(beta) * ", no pruning"),
                    expression(tilde(beta) * ", pruned"))
if (!is.null(p))
  ggsave(file.path(OUT_DIR, "scatter_beta_vs_nopruning.png"), p,
         width = 14, height = 7.5, dpi = 300)

if (any(has_rate)) {
  rate_long <- long_for("rate_")
  rate_long <- rate_long[rate_long$value > 0 & !is.na(rate_long$value), ]
  message("  manhattan_faceted_rate")
  ggsave(file.path(OUT_DIR, "manhattan_faceted_rate.png"),
         manhattan_faceted(rate_long, "cppRATE relative centrality", log_y = TRUE),
         width = 16, height = 12, dpi = 300)

  message("  scatter_rate_vs_nopruning")
  p <- scatter_vs_ref("rate_", "cppRATE, no pruning", "cppRATE, pruned",
                      log_scale = TRUE)
  if (!is.null(p))
    ggsave(file.path(OUT_DIR, "scatter_rate_vs_nopruning.png"), p,
           width = 14, height = 7.5, dpi = 300)
} else {
  message("  RATE figures skipped: no arm has cppRATE output")
}

message("\nWrote outputs to: ", OUT_DIR)
