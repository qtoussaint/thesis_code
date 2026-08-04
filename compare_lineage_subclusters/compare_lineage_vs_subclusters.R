############################################################
## compare_lineage_vs_subclusters.R
##
## Compares binary GWAS runs that differ only in how population structure is
## controlled. Three arms, all fitted to the same dataset with the same pipeline
## settings, so their outputs line up row for row on the original variant index:
##
##   A  subclusters only        gwas_finalmodels/logistic_inference.stan
##   B  lineage only            models/logistic_lineage_inference.stan
##   C  lineage + subclusters   models/logistic_lineage_subcluster_inference.stan
##
## Arm A is the as-published model. Its within-lineage sum-to-zero centering
## cancels beta_lineage out of beta_sublineage exactly, and its linear predictor
## carries no lineage design matrix, so its population-structure term is purely
## within-lineage deviations and it corrects no between-lineage structure. Arm C
## is arm A with the lineage level restored as a covariate.
##
## Figures are produced as two pairwise sets, each against arm A:
##   subclusters_vs_lineage              A vs B
##   subclusters_vs_lineage_subclusters  A vs C
##
## For rifampicin the causal locus is essentially known, so significant variants
## outside rpoB/rpoC/rpoA are a direct proxy for residual population-structure
## confounding, while the on-target count confirms real signal was not lost.
##
## Usage:
##   Rscript compare_lineage_vs_subclusters.R [tb_rifampicin]
## Defaults to tb_rifampicin.
############################################################

suppressPackageStartupMessages({
  library(ggplot2)
})

CODE_BASE <- "/nfs/research/jlees/jacqueline/thesis_code"
source(file.path(CODE_BASE, "paper_figures", "gene_label_helpers.R"))

## ---------------------------- paths -------------------------------------- ##

RESULTS_BASE <- "/nfs/research/jlees/jacqueline/thesis_results"
DATASETS_DIR <- file.path(RESULTS_BASE, "gwas_datasets", "inference")
COMPARE_DIR  <- file.path(RESULTS_BASE, "compare_lineage_subclusters")

COMPARISONS <- list(
  tb_rifampicin = list(
    dataset     = "07_tb_rifampicin_binary",
    annotations = "/nfs/research/jlees/jacqueline/gwas_data/tuberculosis/cryptic_regeno_snpeff/fields_filtered.txt",
    genes       = file.path(CODE_BASE, "gwas_genesofinterest",
                            "tb_rifampicin_genesofinterest.txt"),
    on_target   = c("rpoB", "rpoC", "rpoA"),
    arms = list(
      subclusters = list(
        label = "lineage subclusters",
        dir = file.path(RESULTS_BASE, "gwas_tb_rifampicin", "inference",
                        "07_tb_rifampicin_binary_logistic")),
      lineage = list(
        label = "lineage clusters",
        dir = file.path(COMPARE_DIR,
                        "07_tb_rifampicin_binary_logistic_lineage")),
      lineage_sub = list(
        label = "lineage clusters + subclusters",
        dir = file.path(COMPARE_DIR,
                        "07_tb_rifampicin_binary_logistic_lineage_subcluster"))
    ),
    sets = list(
      list(key = "subclusters_vs_lineage",
           arms = c("subclusters", "lineage")),
      list(key = "subclusters_vs_lineage_subclusters",
           arms = c("subclusters", "lineage_sub"))
    )
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
ARM_LABELS <- stats::setNames(vapply(cfg$arms, `[[`, character(1), "label"),
                              ARM_KEYS)

# Colour follows the arm, not its position in a set, so an arm keeps its hue
# across both figure sets. Validated with the dataviz palette checks (lightness
# band, chroma floor, all-pairs CVD separation, normal-vision floor, contrast on
# white) as a trio and as each plotted pair: Okabe-Ito vermillion / blue /
# bluish-green. The house Hiroshige picks used elsewhere in thesis_code fail the
# lightness band and chroma floor.
ARM_COLOURS <- stats::setNames(c("#D55E00", "#0072B2", "#009E73"), ARM_KEYS)

## ---------------------------- readers ------------------------------------ ##

# Median beta and 89% CI significance, one row per variant on the original
# (pre-pruning) variant index.
#
# regen_effects_from_draws.R output is preferred over the pipeline's own
# fitted_model/depruned_variant_effects.csv, because the latter holds the
# STANDARDIZED effect: cpprate.R:304 matches beta_variant_std ahead of
# beta_variant (see that script's header). Mixing the two across arms would
# compare beta against beta * sd(x_v). Significance is unaffected either way,
# since the two differ by a positive per-variant factor.
effects_path <- function(run_dir) {
  from_draws <- file.path(run_dir, "fitted_model",
                          "depruned_variant_effects_from_draws.csv")
  canonical  <- file.path(run_dir, "fitted_model", "depruned_variant_effects.csv")
  if (file.exists(from_draws)) from_draws else canonical
}

read_effects <- function(run_dir) {
  path <- effects_path(run_dir)
  stopifnot(file.exists(path))
  if (!grepl("_from_draws[.]csv$", path))
    warning("  ", basename(run_dir), ": falling back to the pipeline CSV, whose ",
            "effects are STANDARDIZED. Run regen_effects_from_draws.R on this ",
            "run before comparing effect sizes across arms.", call. = FALSE)
  df <- read.csv(path)
  stopifnot(all(c("variant_id", "median", "signif") %in% names(df)))
  df$signif <- as.logical(df$signif)
  df
}

read_h2 <- function(run_dir) {
  path <- file.path(run_dir, "plots", "heritability", "heritability_summary.csv")
  if (!file.exists(path)) return(NULL)
  read.csv(path)
}

## ---------------------------- load ---------------------------------------- ##

missing <- vapply(cfg$arms, function(a) !file.exists(effects_path(a$dir)),
                  logical(1))
if (any(missing)) {
  stop("Missing fitted output for: ",
       paste(ARM_LABELS[missing], collapse = ", "),
       "\n  expected fitted_model/depruned_variant_effects[_preview].csv under:\n    ",
       paste(vapply(cfg$arms[missing], `[[`, character(1), "dir"),
             collapse = "\n    "))
}

message("Reading run outputs...")
vi <- read.csv(VARIANT_INDEX)
message("  Variant index: ", nrow(vi), " variants")

eff <- lapply(cfg$arms, function(a) read_effects(a$dir))
for (k in ARM_KEYS) stopifnot(nrow(eff[[k]]) == nrow(vi))
for (k in ARM_KEYS[-1])
  stopifnot(identical(eff[[k]]$variant_id, eff[[ARM_KEYS[1]]]$variant_id))
message("  Effects aligned on the same variant order across all ",
        length(eff), " arms")

## ---------------------------- annotate ------------------------------------ ##

# annotate_genes() maps POS -> ANN....GENE, fills unannotated with "MODIFIER"
# and applies the genes-of-interest display names.
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
  lab <- unname(ARM_LABELS[k])
  s <- df[[paste0("signif_", k)]]
  n_on <- sum(s & df$on_target); n_off <- sum(s & !df$on_target)
  hn <- h2_get(lab, "h2_narrow"); hb <- h2_get(lab, "h2_broad")
  data.frame(model = lab,
             n_signif = sum(s),
             n_signif_on_target = n_on,
             n_signif_off_target = n_off,
             frac_signif_off_target = if (sum(s) > 0) n_off / sum(s) else NA_real_,
             h2_narrow = hn, h2_broad = hb, V_pop_share = hb - hn,
             stringsAsFactors = FALSE)
}

all_summary <- do.call(rbind, lapply(ARM_KEYS, arm_summary))
write.csv(all_summary, file.path(OUT_DIR, "comparison_summary_all_arms.csv"),
          row.names = FALSE)
message("\nAll-arm summary:")
print(all_summary, row.names = FALSE)

## ---------------------------- figures -------------------------------------- ##

# Format follows paper_figures/manhattan_builders.R:make_beta_manhattans() --
# x = genome coordinate, y = tilde(beta), geom_point(alpha = 0.4) -- with colour
# mapped to model instead of cutpoint so the arms overlay in one panel.
#
# Everything is +4pt over the paper_figures base of 14. geom_text/geom_text_repel
# take size in mm rather than points, so those go up by 4/.pt (~1.4) instead.
BASE_SIZE <- 18
PT_BUMP   <- 4 / .pt

# Magnitude floor below which a gene is left unlabelled on the exhaustive panel.
# Set to 0 to label every gene carrying an off-target significant variant; the
# exclusion band below is what keeps those labels legible. Never a top-N cap, so
# a large-effect gene can never be hidden.
LABEL_MIN_ABS_BETA <- 0

# Labels on the exhaustive panel are kept out of +/- this band, where almost all
# variants sit and label text would be unreadable against the point cloud.
LABEL_EXCLUSION_BAND <- 0.025
base_theme <- theme_minimal(base_size = BASE_SIZE)

# Peak-|beta| row per gene, carrying the signed effect through so the label sits
# on the point actually plotted. select_label_rows() ranks on score_col and
# returns whole rows, so `value` stays the signed beta of the chosen row.
peak_rows_for <- function(dat, genes) {
  lab_src <- dat
  lab_src$score <- abs(dat$value)
  select_label_rows(lab_src, gene_col = "gene", score_col = "score",
                    mode = "gene_list", genes = genes)
}

# repel_label_layer() hardcodes max.overlaps = 30, which silently drops labels.
# Fine when labelling three genes of interest, not acceptable on the exhaustive
# panel where a missing label would read as "no hit there".
#
# Two modes. With few labels (genes of interest) each label sits directly above
# its peak SNP, repelled only along y, and needs no connector. With ~100 labels
# they cannot stay above their own points, so the label is free to move in both
# directions and a leader line ties it back to the variant it names; a ring is
# drawn on that variant so the target is visible inside the point cloud.
label_layers_nocap <- function(label_df, size, nudge_y, segments) {
  if (!segments) {
    return(list(ggrepel::geom_text_repel(
      data = label_df,
      mapping = aes(x = .data[["pos"]], y = .data[["value"]],
                    label = .data[["gene_expr"]]),
      parse = TRUE, size = size, colour = "grey25", inherit.aes = FALSE,
      direction = "y", nudge_y = nudge_y, vjust = 0,
      min.segment.length = Inf, box.padding = 0.15, point.padding = 0.2,
      max.overlaps = Inf, seed = 1)))
  }
  # ggrepel only avoids overlapping the points in its OWN data, not the full
  # scatter, so labels otherwise settle straight onto the dense band around
  # zero. Its ylim argument does constrain placement, so labels are split by the
  # sign of their peak effect and each half is confined outside the band: those
  # naming a positive peak go above +LABEL_EXCLUSION_BAND, those naming a
  # negative peak below -LABEL_EXCLUSION_BAND. The leader line still ties each
  # one back to its variant.
  repel_half <- function(dat, ylim) {
    ggrepel::geom_text_repel(
      data = dat,
      mapping = aes(x = .data[["pos"]], y = .data[["value"]],
                    label = .data[["gene_expr"]]),
      parse = TRUE, size = size, colour = "grey20", inherit.aes = FALSE,
      ylim = ylim,
      # leader line to the named variant; 0 forces one on every label
      min.segment.length = 0, segment.colour = "grey55", segment.size = 0.3,
      segment.alpha = 0.9,
      box.padding = 0.45, point.padding = 0.3, max.overlaps = Inf,
      # ~100 labels need more than ggrepel's 0.5 s default to settle
      max.iter = 100000, max.time = 15, seed = 1)
  }
  pos <- label_df[label_df$value >= 0, , drop = FALSE]
  neg <- label_df[label_df$value <  0, , drop = FALSE]
  layers <- list(
    geom_point(data = label_df,
               mapping = aes(x = .data[["pos"]], y = .data[["value"]]),
               inherit.aes = FALSE, shape = 21, fill = NA,
               colour = "grey20", size = 2.4, stroke = 0.6))
  if (nrow(pos) > 0L)
    layers <- c(layers, list(repel_half(pos, c(LABEL_EXCLUSION_BAND, NA))))
  if (nrow(neg) > 0L)
    layers <- c(layers, list(repel_half(neg, c(NA, -LABEL_EXCLUSION_BAND))))
  layers
}

manhattan_beta <- function(dat, colours, labels = c("goi", "offtarget", "none")) {
  labels <- match.arg(labels)
  p <- ggplot(dat, aes(x = pos, y = value, colour = model)) +
    geom_point(alpha = 0.4) +
    scale_colour_manual(values = colours) +
    guides(colour = guide_legend(override.aes = list(alpha = 1, size = 4.5))) +
    xlab("genome coordinate (bp)") +
    ylab(expression(tilde(beta))) +
    labs(colour = NULL) +
    base_theme +
    theme(legend.position = "top")
  if (labels == "none") return(p)

  genes <- if (labels == "goi") cfg$on_target else {
    off <- unique(dat$gene[dat$signif & !dat$on_target & dat$gene != "MODIFIER"])
    message("    off-target significant hits fall in ", length(off), " distinct genes")
    unique(c(cfg$on_target, off))
  }
  label_df <- peak_rows_for(dat, genes)
  # Drop genes whose peak effect is too small to be worth naming. peak_rows_for()
  # already returns the max-|beta| row per gene, so value is that peak.
  if (labels != "goi") {
    n_before <- nrow(label_df)
    label_df <- label_df[abs(label_df$value) >= LABEL_MIN_ABS_BETA, , drop = FALSE]
    message("    dropped ", n_before - nrow(label_df), " genes with peak |median| < ",
            LABEL_MIN_ABS_BETA)
  }
  if (nrow(label_df) == 0L) return(p)
  message("    labelling ", nrow(label_df), " genes on the '", labels, "' panel")
  # the exhaustive panel needs vertical room on both sides for labels to escape
  # the dense band around zero
  expand_mult <- if (labels == "goi") c(0.05, 0.16) else c(0.18, 0.22)
  p +
    scale_y_continuous(expand = expansion(mult = expand_mult)) +
    label_layers_nocap(label_df,
                       size = (if (labels == "goi") 5 else 2.6) + PT_BUMP,
                       nudge_y = 0.04 * diff(range(dat$value, na.rm = TRUE)),
                       segments = labels != "goi")
}

build_set <- function(set) {
  keys    <- set$arms
  labels  <- unname(ARM_LABELS[keys])
  colours <- stats::setNames(unname(ARM_COLOURS[keys]), labels)
  set_dir <- file.path(OUT_DIR, set$key)
  dir.create(set_dir, showWarnings = FALSE, recursive = TRUE)
  message("\n=== set: ", set$key, "  (", paste(labels, collapse = " vs "), ")")

  # summary + pairwise agreement for this pair
  s <- do.call(rbind, lapply(keys, arm_summary))
  write.csv(s, file.path(set_dir, "comparison_summary.csv"), row.names = FALSE)
  print(s, row.names = FALSE)

  ba <- df[[paste0("beta_", keys[1])]]; bb <- df[[paste0("beta_", keys[2])]]
  sa <- df[[paste0("signif_", keys[1])]]; sb <- df[[paste0("signif_", keys[2])]]
  pw <- data.frame(
    model_a = labels[1], model_b = labels[2],
    beta_pearson     = cor(ba, bb, method = "pearson"),
    beta_spearman    = cor(ba, bb, method = "spearman"),
    abs_beta_pearson = cor(abs(ba), abs(bb), method = "pearson"),
    n_signif_both    = sum(sa & sb),
    n_signif_a_only  = sum(sa & !sb),
    n_signif_b_only  = sum(!sa & sb),
    stringsAsFactors = FALSE)
  write.csv(pw, file.path(set_dir, "comparison_pairwise.csv"), row.names = FALSE)
  print(pw, row.names = FALSE)

  long <- do.call(rbind, lapply(keys, function(k)
    data.frame(pos = df$pos,
               value = df[[paste0("beta_", k)]],
               signif = df[[paste0("signif_", k)]],
               gene = df$gene, on_target = df$on_target,
               model = unname(ARM_LABELS[k]), stringsAsFactors = FALSE)))
  long$model <- factor(long$model, levels = labels)
  # Interleave so neither arm systematically overplots the other
  set.seed(1)
  long <- long[sample.int(nrow(long)), ]

  ggsave(file.path(set_dir, "manhattan_overlay_beta.png"),
         manhattan_beta(long, colours, "goi"), width = 16, height = 6, dpi = 300)
  ggsave(file.path(set_dir, "manhattan_overlay_beta_simple.png"),
         manhattan_beta(long, colours, "none"), width = 16, height = 6, dpi = 300)
  ggsave(file.path(set_dir, "manhattan_overlay_beta_offtarget_labelled.png"),
         manhattan_beta(long, colours, "offtarget"), width = 16, height = 10, dpi = 300)

  # Off-target significant hits: the headline panel. Identity comes from the axis
  # text, so the redundant fill legend is dropped; every bar is directly labelled.
  # Panel titles are plotmath so the gene names render italic; label_parsed
  # evaluates the factor levels as expressions.
  goi_expr <- paste(sprintf('italic("%s")', sort(cfg$on_target)),
                    collapse = ' * "/" * ')
  # Leading ellipsis so each strip reads as a continuation of the y-axis label:
  # "significant variants (89% CI excludes 0) ... within rpoA/rpoB/rpoC"
  lab_off <- paste0('"… not within " * ', goi_expr)
  lab_on  <- paste0('"… within " * ', goi_expr)

  bar_df <- rbind(
    data.frame(model = s$model, count = s$n_signif_off_target, panel = lab_off),
    data.frame(model = s$model, count = s$n_signif_on_target,  panel = lab_on))
  bar_df$model <- factor(bar_df$model, levels = labels)
  bar_df$panel <- factor(bar_df$panel, levels = c(lab_off, lab_on))
  p_bars <- ggplot(bar_df, aes(x = model, y = count, fill = model)) +
    geom_col(width = 0.6) +
    geom_text(aes(label = format(count, big.mark = ",")), vjust = -0.4,
              size = 4.5 + PT_BUMP, colour = "grey25") +
    facet_wrap(~panel, scales = "free_y", labeller = label_parsed) +
    scale_fill_manual(values = colours, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(x = NULL, y = "significant variants (89% CI excludes 0)") +
    base_theme +
    theme(axis.text.x = element_text(angle = 15, hjust = 1),
          panel.grid.major.x = element_blank())
  ggsave(file.path(set_dir, "offtarget_hits.png"), p_bars,
         width = 12, height = 6.5, dpi = 300)

  message("  wrote ", set_dir)
}

invisible(lapply(cfg$sets, build_set))

message("\nWrote outputs to: ", OUT_DIR)
