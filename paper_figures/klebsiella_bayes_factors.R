#!/usr/bin/env Rscript
# Faceted Bayes factor Manhattan across the simulated Klebsiella effect sizes.
#
# The companion to klebsiella_effectsize_faceted_manhattan.R and deliberately
# built to the same design: one panel per effect size, causal variant rs_26645
# in red and directly labelled, no title or subtitle. One figure per LD-pruning
# arm, so the two are never overlaid. Unlike the effect-size figure these panels
# carry no h2 annotation -- h2 still gates which runs are plotted at all, it is
# just not drawn.
#
# Evidence bands follow the original bayesian_bridge.R figure
# (bayes_factors_klebsiella_100.png): neutral inside (1/3, 3), moderate out to
# (1/10, 10), strong beyond.
#
# Bayes factors are computed once by paper_figures/klebsiella_extract_effects.R
# (which owns the prior simulation and the ROPE definition, and documents why
# neither matches the older bayesian_bridge.R); this script only plots them.
#
# Usage:
#   mamba activate gwas_pipeline
#   Rscript paper_figures/klebsiella_extract_effects.R     # writes the BFs
#   Rscript paper_figures/klebsiella_bayes_factors.R       # draws them
#
# Outputs (to <RES>/gwas_klebsiella_homoplasic/figures/):
#   klebsiella_effectsize_faceted_bayes_factors_pruned.png
#   klebsiella_effectsize_faceted_bayes_factors_nopruning.png
#   klebsiella_effectsize_bf_summary.csv

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
})
grDevices::pdf(NULL)

RES       <- "/nfs/research/jlees/jacqueline/thesis_results"
DATASETS  <- file.path(RES, "gwas_datasets", "inference")
EXTRACTED <- file.path(RES, "gwas_klebsiella_homoplasic", "extracted")
OUT_DIR   <- file.path(RES, "gwas_klebsiella_homoplasic", "figures")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

COL_NULL   <- "#6b7280"
COL_CAUSAL <- "#c1121f"
CAUSAL_LABEL <- "simulated causal variant"
EF_LEVELS  <- c("1.5", "2.5", "10", "30")

MODE_TITLES <- c(pruned    = "with LD pruning",
                 nopruning = "without LD pruning")

ef_facet <- function(x) paste0("effect size = ", x)

# Ticks read as exponents so they agree with the log10(BF) axis title. The
# evidence bands are still specified in BF units and land where they should.
LOG_BREAKS <- scales::trans_breaks("log10", function(x) 10^x, n = 8)
log10_labels <- function(x) {
  ifelse(is.na(x), "", format(log10(x), trim = TRUE, drop0trailing = TRUE))
}

truth <- read.csv(file.path(DATASETS, "klebsiella_homoplasic_truth_summary.csv"),
                  stringsAsFactors = FALSE)
truth$ef_chr <- sub("\\.0$", "", format(truth$effect_size, trim = TRUE))
causal_name  <- unique(truth$causal_variant)
stopifnot(length(causal_name) == 1L)

# Same validity gate as the effect-size figure: a fit that recovered essentially
# no variant-explained variance is a bad ADVI optimum, and its Bayes factors are
# not evidence about anything. Such panels are drawn empty and marked, never
# plotted as though they were results.
extract_summary <- {
  f <- file.path(EXTRACTED, "klebsiella_extracted_summary.csv")
  if (file.exists(f)) data.table::fread(f) else NULL
}

# Base-R subsetting: data.table's i-expression would let the arguments shadow
# the identically named columns.
lookup <- function(ds_name, mode_name, column) {
  if (is.null(extract_summary)) return(NA)
  hit <- extract_summary$dataset == ds_name & extract_summary$mode == mode_name
  if (!any(hit)) return(NA)
  extract_summary[[column]][which(hit)[1]]
}
is_valid_run <- function(ds_name, mode_name) {
  if (is.null(extract_summary)) return(TRUE)
  v <- lookup(ds_name, mode_name, "valid")
  if (is.na(v)) return(FALSE)
  isTRUE(as.logical(v))
}

build_arm <- function(mode) {
  rows <- list(); excluded <- character(0)
  for (i in seq_len(nrow(truth))) {
    ds  <- truth$dataset[i]
    csv <- file.path(EXTRACTED, paste0(ds, "_", mode, "_variant_effects.csv"))
    if (!is_valid_run(ds, mode) || !file.exists(csv)) {
      excluded <- c(excluded, truth$ef_chr[i])
      message("[", mode, "] EF ", truth$ef_chr[i],
              " excluded: fit failed or did not pass the h2 validity gate")
      next
    }
    rows[[length(rows) + 1L]] <- data.table::fread(csv)
  }
  if (length(rows) == 0L) {
    warning("No usable extracts for arm '", mode, "' -- skipping", call. = FALSE)
    return(NULL)
  }
  df <- data.table::rbindlist(rows)
  df[, ef_fac := factor(ef_chr, levels = EF_LEVELS, labels = ef_facet(EF_LEVELS))]

  # A log axis cannot show BF <= 0; the continuity correction in
  # klebsiella_extract_effects.R makes that impossible, but guard and count
  # rather than silently dropping points.
  plot_df   <- df[is.finite(bf) & bf > 0]
  n_dropped <- nrow(df) - nrow(plot_df)
  if (n_dropped > 0) {
    message("[", mode, "] dropped ", n_dropped, " non-finite/non-positive BFs")
  }
  causal_df <- plot_df[is_causal == TRUE]
  if (nrow(causal_df) == 0L) {
    warning("[", mode, "] causal variant has no plottable BF", call. = FALSE)
  }

  list(df = df, plot_df = plot_df, causal_df = causal_df, excluded = excluded)
}

make_bf_plot <- function(plot_df, causal_df, excluded = character(0)) {
  # Keep a facet for every effect size, including excluded ones, so the grid is
  # comparable between arms and a failed fit is visibly labelled rather than
  # silently absent.
  excl_df <- NULL
  if (length(excluded) > 0) {
    y_mid <- exp(mean(log(range(plot_df$bf))))
    excl_df <- data.table(
      ef_fac   = factor(ef_facet(excluded), levels = ef_facet(EF_LEVELS)),
      position = mean(range(plot_df$position)),
      bf       = y_mid)
  }

  p <- ggplot(plot_df[is_causal == FALSE], aes(x = position, y = bf)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 1/3,  ymax = 3,   alpha = 0.10, fill = "#2166ac") +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 3,    ymax = 10,  alpha = 0.10, fill = "#4d9221") +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 1/10, ymax = 1/3, alpha = 0.10, fill = "#4d9221") +
    geom_point(colour = COL_NULL, alpha = 0.35, size = 0.7, stroke = 0) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "black", linewidth = 0.4) +
    # Causal variant drawn last so it sits on top, with a surface ring to keep
    # it legible where the null cloud is dense.
    geom_point(data = causal_df, colour = "white", size = 3.6, stroke = 0) +
    geom_point(data = causal_df, colour = COL_CAUSAL, size = 2.6, stroke = 0) +
    geom_text(data = causal_df, label = CAUSAL_LABEL, colour = COL_CAUSAL,
              vjust = -1.1, size = 3.6, fontface = "bold", show.legend = FALSE) +
    scale_y_log10(breaks = LOG_BREAKS, labels = log10_labels,
                  expand = expansion(mult = c(0.05, 0.22))) +
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

  if (!is.null(excl_df)) {
    p <- p +
      geom_blank(data = excl_df, aes(x = position, y = bf)) +
      geom_text(data = excl_df, aes(x = position, y = bf),
                label = "fit failed\n(no usable result)", colour = "grey35",
                size = 3.6, fontface = "italic", lineheight = 0.95)
  }
  p
}

arms <- list()
for (mode in names(MODE_TITLES)) {
  a <- build_arm(mode)
  if (is.null(a)) next
  arms[[mode]] <- a
  p <- make_bf_plot(a$plot_df, a$causal_df, a$excluded)
  png_path <- file.path(
    OUT_DIR, paste0("klebsiella_effectsize_faceted_bayes_factors_", mode, ".png"))
  ggsave(png_path, p, width = 15, height = 4.6, dpi = 600, bg = "white")
  message("Wrote ", png_path)
}
if (length(arms) == 0L) {
  stop("No extracts found in ", EXTRACTED,
       " -- run paper_figures/klebsiella_extract_effects.R first")
}

# `mode` is a column of df as well as the loop variable, so it goes in `by`
# rather than `j` -- assigning mode = mode inside j resolves to the column and
# silently defeats the aggregation.
bf_summary <- data.table::rbindlist(lapply(names(arms), function(arm) {
  arms[[arm]]$df[, .(
    n_variants = .N,
    causal_bf  = bf[is_causal],
    bf_rank    = rank(-bf, ties.method = "min")[is_causal],
    n_bf_gt10  = sum(bf > 10, na.rm = TRUE)
  ), by = .(mode, ef_chr)]
}))
bf_summary[, ef_num := as.numeric(ef_chr)]
data.table::setorder(bf_summary, mode, ef_num)
bf_summary[, ef_num := NULL]

summary_csv <- file.path(OUT_DIR, "klebsiella_effectsize_bf_summary.csv")
data.table::fwrite(bf_summary, summary_csv)
message("Wrote ", summary_csv)
print(bf_summary)
