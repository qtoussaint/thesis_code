#!/usr/bin/env Rscript
# Faceted Bayes factor Manhattans for the four unpruned simulated Klebsiella
# fits, one figure per null:
#   * prior null  -- BF against the ROPE derived from the 89% interval of the
#                    prior-null distribution (rope_significance.R)
#   * point null  -- Savage-Dickey BF at beta = 0
#
# Values come from rope_pd_bf/klebsiella_bf_pd_all_variants.R, which uses the
# same definitions as bf_pd_significance.R and validates them against
# bayestestR. Effects here are on the STANDARDIZED scale, because that is the
# scale the ROPE is defined on.
#
# The causal variant carries its probability of direction as a second label.
#
# Usage:
#   mamba activate gwas_pipeline
#   Rscript paper_figures/klebsiella_bf_pointnull_priornull.R
#
# Outputs (to <RES>/gwas_klebsiella_homoplasic/figures/):
#   klebsiella_effectsize_faceted_bf_priornull.png
#   klebsiella_effectsize_faceted_bf_pointnull.png
#   klebsiella_effectsize_bf_pd_summary.csv

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
})
grDevices::pdf(NULL)

RES      <- "/nfs/research/jlees/jacqueline/thesis_results"
DATASETS <- file.path(RES, "gwas_datasets", "inference")
BFDIR    <- file.path(RES, "gwas_klebsiella_homoplasic", "rope_pd_bf")
OUT_DIR  <- file.path(RES, "gwas_klebsiella_homoplasic", "figures")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

COL_NULL     <- "#6b7280"
COL_CAUSAL   <- "#c1121f"
CAUSAL_LABEL <- "simulated causal variant"
EF_LEVELS    <- c("1.5", "2.5", "10", "30")

ef_facet <- function(x) paste0("effect size = ", x)

LOG_BREAKS <- scales::trans_breaks("log10", function(x) 10^x, n = 8)
log10_labels <- function(x) {
  ifelse(is.na(x), "", format(log10(x), trim = TRUE, drop0trailing = TRUE))
}

truth <- read.csv(file.path(DATASETS, "klebsiella_homoplasic_truth_summary.csv"),
                  stringsAsFactors = FALSE)
truth$ef_chr <- sub("\\.0$", "", format(truth$effect_size, trim = TRUE))
causal_name  <- unique(truth$causal_variant)
stopifnot(length(causal_name) == 1L)

rows <- list()
for (i in seq_len(nrow(truth))) {
  f <- file.path(BFDIR, paste0(truth$dataset[i], "_all_variants_bf_pd.csv"))
  if (!file.exists(f)) {
    warning("Missing ", f, " -- run submit_klebsiella_rope.sh first", call. = FALSE)
    next
  }
  d <- data.table::fread(f)
  d[, ef_chr := truth$ef_chr[i]]
  rows[[length(rows) + 1L]] <- d
}
if (length(rows) == 0L) stop("No BF/pd inputs found under ", BFDIR)

df <- data.table::rbindlist(rows)
df[, ef_fac := factor(ef_chr, levels = EF_LEVELS, labels = ef_facet(EF_LEVELS))]
df[, is_causal := variant_name == causal_name]

#' One BF Manhattan. `col` selects which null.
make_bf_plot <- function(df, col) {
  d <- data.table::copy(df)
  d[, bf := get(col)]

  # A log axis cannot show BF <= 0 or non-finite BF; count rather than drop
  # silently. A point-null BF can legitimately come back Inf when the posterior
  # density at 0 underflows, which is information, not an error -- but it cannot
  # be placed on an axis, so it is reported and excluded.
  n_bad <- d[!is.finite(bf) | bf <= 0, .N]
  if (n_bad > 0) {
    message("[", col, "] ", n_bad, " variants have non-finite or non-positive BF (of ",
            nrow(d), ")")
  }
  plot_df <- d[is.finite(bf) & bf > 0]

  # An infinite BF is a real result, not a missing one: it means no posterior
  # draw fell inside the ROPE, so the odds ratio is undefined at this number of
  # draws. Dropping those points would silently delete the causal variant from
  # every panel. Instead they are pinned above the finite cloud in their own
  # panel and marked, so the reader sees a censored value rather than nothing.
  inf_df <- d[is_causal == TRUE & (!is.finite(bf) | bf <= 0)]
  if (nrow(inf_df) > 0) {
    ceil <- plot_df[, .(y_top = max(bf) * 4), by = ef_fac]
    inf_df <- merge(inf_df[, .(ef_fac, position, pd, ef_chr)], ceil, by = "ef_fac")
    message("[", col, "]   causal variant is infinite at EF ",
            paste(inf_df$ef_chr, collapse = ", "),
            " -- drawn at the panel ceiling and marked")
  }
  causal_df <- plot_df[is_causal == TRUE]

  # pd label for the causal variant, plotmath so "pd" renders italic. Applies to
  # the finite and the censored causal points alike.
  pd_lab <- rbind(
    causal_df[, .(ef_fac, position, y = bf, pd)],
    if (nrow(inf_df) > 0) inf_df[, .(ef_fac, position, y = y_top, pd)] else NULL)
  pd_lab[, txt := sprintf('italic(pd) == "%.4f"', pd)]

  p <- ggplot(plot_df[is_causal == FALSE], aes(x = position, y = bf)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 1/3,  ymax = 3,   alpha = 0.10, fill = "#2166ac") +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 3,    ymax = 10,  alpha = 0.10, fill = "#4d9221") +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 1/10, ymax = 1/3, alpha = 0.10, fill = "#4d9221") +
    geom_point(colour = COL_NULL, alpha = 0.35, size = 0.7, stroke = 0) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "black", linewidth = 0.4) +
    geom_point(data = causal_df, colour = "white", size = 3.6, stroke = 0) +
    geom_point(data = causal_df, colour = COL_CAUSAL, size = 2.6, stroke = 0) +
    # Name label keyed off pd_lab, not causal_df, so it still appears when the
    # causal variant's BF is censored and therefore absent from causal_df.
    geom_text(data = pd_lab, aes(x = position, y = y), label = CAUSAL_LABEL,
              colour = COL_CAUSAL, vjust = -2.4, size = 3.6, fontface = "bold",
              inherit.aes = FALSE, show.legend = FALSE) +
    geom_text(data = pd_lab, aes(x = position, y = y, label = txt),
              parse = TRUE, vjust = -0.9, size = 3.4, colour = COL_CAUSAL,
              inherit.aes = FALSE) +
    scale_y_log10(breaks = LOG_BREAKS, labels = log10_labels,
                  expand = expansion(mult = c(0.05, 0.28))) +
    scale_x_continuous(labels = scales::comma) +
    facet_wrap(~ ef_fac, nrow = 1) +
    labs(x = "variants, positioned by rs number",
         y = expression(log[10] * "(Bayes factor)")) +
    theme_minimal(base_size = 13) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(linewidth = 0.3, colour = "grey88"),
          strip.text = element_text(face = "bold"),
          axis.text.x = element_text(size = 8))

  # Censored causal points: open upward triangle at the ceiling, so it reads as
  # "at least this far up" rather than as a measured value.
  if (nrow(inf_df) > 0) {
    p <- p +
      geom_point(data = inf_df, aes(x = position, y = y_top),
                 shape = 24, size = 3.0, stroke = 1.1,
                 colour = COL_CAUSAL, fill = "white", inherit.aes = FALSE) +
      geom_text(data = inf_df, aes(x = position, y = y_top),
                label = "BF == infinity", parse = TRUE, hjust = -0.22,
                size = 3.4, colour = COL_CAUSAL, inherit.aes = FALSE)
  }
  p
}

for (spec in list(list(col = "bf_rope",  tag = "priornull"),
                  list(col = "bf_point", tag = "pointnull"))) {
  p <- make_bf_plot(df, spec$col)
  png_path <- file.path(
    OUT_DIR, paste0("klebsiella_effectsize_faceted_bf_", spec$tag, ".png"))
  ggsave(png_path, p, width = 15, height = 4.6, dpi = 600, bg = "white")
  message("Wrote ", png_path)
}

summary_out <- df[is_causal == TRUE,
  .(ef_chr, pd, bf_rope, bf_point, median_std, ci_lo_std, ci_hi_std,
    rope_lo, rope_hi, signif_median, signif_ci)]
summary_out[, `:=`(
  bf_rope_rank  = sapply(ef_chr, function(e)
    rank(-df[ef_chr == e]$bf_rope,  ties.method = "min")[df[ef_chr == e]$is_causal]),
  n_signif_ci   = sapply(ef_chr, function(e) sum(df[ef_chr == e]$signif_ci)),
  n_signif_med  = sapply(ef_chr, function(e) sum(df[ef_chr == e]$signif_median)))]
summary_out[, ef_num := as.numeric(ef_chr)]
data.table::setorder(summary_out, ef_num)
summary_out[, ef_num := NULL]

csv <- file.path(OUT_DIR, "klebsiella_effectsize_bf_pd_summary.csv")
data.table::fwrite(summary_out, csv)
message("Wrote ", csv)
print(summary_out)
