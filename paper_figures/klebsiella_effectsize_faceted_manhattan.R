#!/usr/bin/env Rscript
# Faceted Manhattan across simulated Klebsiella effect sizes.
#
# One panel per effect size (EF 1.5, 2.5, 10, 30; homoplasic architecture,
# h2 = 0.9, single causal variant rs_26645). The causal variant is drawn in red
# and directly labelled; every other variant is neutral grey.
#
# The y axis is log10 |median beta|, shared across panels. Shared rather than free
# is the whole point of the figure: the four phenotypes are the same simulation
# rescaled, so the causal effect should climb by the same factor as the effect
# size while the null variants stay put. A free y scale would normalise that
# growth away and every panel would look identical.
#
# Colour is doing an identity job (causal vs not) over two categories. The pair
# was checked rather than eyeballed: OKLab dE 21.8 for normal vision, and 13.9 /
# 18.0 / 121.8 under simulated deuteranopia / protanopia / tritanopia, all above
# the >= 8 CVD floor, on a white surface at 4.8:1 and 6.2:1 contrast. The causal
# point is labelled as well as coloured, so identity is never colour-alone.
#
# Outputs (to <RES>/gwas_klebsiella_homoplasic/figures/, alongside the runs and
# extracts they are built from rather than in the shared paper_figures tree):
#   klebsiella_effectsize_faceted_manhattan.png  -- the faceted Manhattan
#   klebsiella_causal_effect_growth.png          -- recovered vs simulated effect
#   klebsiella_effectsize_causal_summary.csv     -- the numbers behind both
#
# Usage:
#   mamba activate gwas_pipeline
#   Rscript paper_figures/klebsiella_effectsize_faceted_manhattan.R

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
})

grDevices::pdf(NULL)

RES         <- "/nfs/research/jlees/jacqueline/thesis_results"
DATASETS    <- file.path(RES, "gwas_datasets", "inference")
RUNS        <- file.path(RES, "gwas_klebsiella_homoplasic", "inference")
OUT_DIR     <- file.path(RES, "gwas_klebsiella_homoplasic", "figures")
TRUTH_CSV   <- file.path(DATASETS, "klebsiella_homoplasic_truth_summary.csv")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Two-colour identity palette (see header for the CVD numbers).
COL_NULL   <- "#6b7280"
COL_CAUSAL <- "#c1121f"

# Callout text for the one labelled point. The variant's rs_ ID is a simulation
# artefact rather than a real locus name, so the label says what the point IS.
CAUSAL_LABEL <- "simulated causal variant"

# Panel order: ascending effect size, not the alphabetical order the
# directory listing would give.
EF_LEVELS <- c("1.5", "2.5", "10", "30")

# Facet strip text. Defined once and used everywhere a factor is built, so the
# real panels and the failed-fit placeholders cannot drift apart -- mismatched
# labels would silently produce duplicate facets.
ef_facet <- function(x) paste0("effect size = ", x)

# The axis is log10-transformed, so the natural breaks are powers of ten. The
# title reads log10|beta|, so the ticks have to read as the exponents (-6, -4,
# ...) and not as 10^-6 -- otherwise the label claims a quantity the ticks do
# not show.
LOG_BREAKS <- scales::trans_breaks("log10", function(x) 10^x, n = 8)
log10_labels <- function(x) {
  ifelse(is.na(x), "", format(log10(x), trim = TRUE, drop0trailing = TRUE))
}

stopifnot(file.exists(TRUTH_CSV))
truth <- read.csv(TRUTH_CSV, stringsAsFactors = FALSE)

# ---------------------------------------------------------------- load runs --
# Effects come from paper_figures/klebsiella_extract_effects.R rather than the
# pipeline: these runs use --ld_pruning false, so there is no
# depruned_variant_effects.csv (that file is only written in the pipeline's
# ld_pruning == "true" branch).
EXTRACTED <- file.path(RES, "gwas_klebsiella_homoplasic", "extracted")

load_run <- function(dataset, effect_size, mode) {
  csv <- file.path(EXTRACTED, paste0(dataset, "_", mode, "_variant_effects.csv"))
  if (!file.exists(csv)) {
    warning("No extract yet for ", dataset, " (", mode,
            ") -- run klebsiella_extract_effects.R first. Expected: ", csv,
            call. = FALSE)
    return(NULL)
  }
  eff <- data.table::fread(csv)
  data.table(
    effect_size  = effect_size,
    mode         = mode,
    variant_name = eff$variant_name,
    position     = eff$position,
    median       = eff$median,
    abs_median   = abs(eff$median),
    q05          = eff$q05,
    q95          = eff$q95,
    # "significant" here means the 90% credible interval excludes zero; the
    # pipeline's own signif column is unavailable without the pruning branch.
    signif       = (eff$q05 > 0) | (eff$q95 < 0)
  )
}

# Effect size as an ordered factor so facets read left-to-right ascending.
ef_label <- function(x) sub("\\.0$", "", format(x, trim = TRUE))
truth$ef_chr <- ef_label(truth$effect_size)
truth$ef_fac <- factor(truth$ef_chr, levels = EF_LEVELS,
                       labels = ef_facet(EF_LEVELS))

# Causal flag, joined on variant name rather than assuming a row index.
causal_name <- unique(truth$causal_variant)
stopifnot(length(causal_name) == 1L)

# Dashed reference at the observed carrier vs non-carrier phenotype difference.
# Note this is a MARGINAL difference whereas beta_variant is estimated
# CONDITIONAL on lineage and sublineage, so it is a reference rather than a
# ground truth and small gaps should not be read as model error. Drawn on the
# Manhattan row only -- the CI zoom row is deliberately left clean so the
# interval is read on its own terms.
truth_lines <- data.table(ef_fac = truth$ef_fac,
                          simulated = abs(truth$observed_diff))

MODE_TITLES <- c(pruned    = "with LD pruning",
                 nopruning = "without LD pruning")

# Which runs are usable. klebsiella_extract_effects.R gates on recovered
# heritability: every phenotype was simulated at h2 = 0.9, so a fit returning
# essentially zero variant-explained variance is a bad ADVI optimum whatever
# its ELBO said. A run that failed outright has no extract at all. Both cases
# are drawn as an explicitly marked empty panel rather than omitted, so a
# missing effect size can never be mistaken for a null result.
extract_summary <- {
  f <- file.path(EXTRACTED, "klebsiella_extracted_summary.csv")
  if (file.exists(f)) data.table::fread(f) else NULL
}

# Base-R subsetting rather than data.table's i-expression: the arguments would
# otherwise shadow the identically named columns.
is_valid_run <- function(ds_name, mode_name) {
  if (is.null(extract_summary)) return(TRUE)
  hit <- extract_summary$dataset == ds_name & extract_summary$mode == mode_name
  if (!any(hit)) return(FALSE)
  isTRUE(as.logical(extract_summary$valid[which(hit)[1]]))
}

# Recovered narrow-sense heritability, annotated per panel so each fit carries
# its own diagnostic. Same base-R subsetting as above, and NA when the run is
# missing rather than a fabricated value.
h2_of_run <- function(ds_name, mode_name) {
  if (is.null(extract_summary)) return(NA_real_)
  hit <- extract_summary$dataset == ds_name & extract_summary$mode == mode_name
  if (!any(hit)) return(NA_real_)
  as.numeric(extract_summary$h2_narrow[which(hit)[1]])
}

# One figure per LD-pruning arm. Each is self-contained: the two arms fit
# different variant sets, so they are never overlaid on shared axes.
build_arm <- function(mode) {
  rows <- list(); excluded <- character(0); h2 <- list()
  for (i in seq_len(nrow(truth))) {
    ok <- is_valid_run(truth$dataset[i], mode)
    r  <- if (ok) load_run(truth$dataset[i], truth$effect_size[i], mode) else NULL
    if (!ok) {
      excluded <- c(excluded, truth$ef_chr[i])
      message("[", mode, "] EF ", truth$ef_chr[i],
              " excluded: fit failed or did not pass the h2 validity gate")
    }
    if (!is.null(r)) {
      rows[[length(rows) + 1L]] <- r
      h2[[length(h2) + 1L]] <- data.table(
        ef_chr = truth$ef_chr[i], h2 = h2_of_run(truth$dataset[i], mode))
    }
  }
  if (length(rows) == 0L) {
    warning("No usable extracts for arm '", mode, "' -- skipping", call. = FALSE)
    return(NULL)
  }
  df <- data.table::rbindlist(rows)
  df[, ef_chr := ef_label(effect_size)]
  df[, ef_fac := factor(ef_chr, levels = EF_LEVELS,
                        labels = ef_facet(EF_LEVELS))]
  df[, is_causal := variant_name == causal_name]

  # A shared log10 y axis cannot show exact zeros. Variants the model pinned to
  # zero (sd_variant == 0) are dropped from the plot and counted, rather than
  # silently floored to an arbitrary epsilon.
  n_zero <- df[abs_median <= 0, .N]
  if (n_zero > 0) {
    message("[", mode, "] dropping ", n_zero,
            " points with |median| == 0 before the log axis")
  }
  plot_df   <- df[abs_median > 0]
  causal_df <- plot_df[is_causal == TRUE]

  # Credible-interval bounds on the |beta| scale the y axis uses. Taking abs()
  # of each end and re-sorting handles a negative effect; an interval straddling
  # zero has no lower bound on a log axis, so it is dropped and reported rather
  # than floored to an arbitrary epsilon.
  causal_df[, ci_lo := pmin(abs(q05), abs(q95))]
  causal_df[, ci_hi := pmax(abs(q05), abs(q95))]
  spans_zero <- causal_df[q05 <= 0 & q95 >= 0]
  if (nrow(spans_zero) > 0) {
    message("[", mode, "] causal CI straddles zero at EF ",
            paste(spans_zero$ef_chr, collapse = ", "),
            " -- no lower bound on a log axis, interval not drawn")
    causal_df[q05 <= 0 & q95 >= 0, c("ci_lo", "ci_hi") := NA_real_]
  }
  if (nrow(causal_df) == 0L) {
    warning("[", mode, "] causal variant has |median| == 0 in every panel",
            call. = FALSE)
  }
  h2_df <- data.table::rbindlist(h2)
  h2_df[, ef_fac := factor(ef_chr, levels = EF_LEVELS, labels = ef_facet(EF_LEVELS))]
  h2_df <- h2_df[!is.na(h2)]

  list(df = df, plot_df = plot_df, causal_df = causal_df, excluded = excluded,
       h2_df = h2_df)
}

make_manhattan <- function(plot_df, causal_df, mode, excluded = character(0),
                           h2_df = NULL, show_ci = FALSE) {
  # Keep a facet for every effect size, including excluded ones, so the panel
  # grid is comparable between arms and a failed fit is visibly labelled rather
  # than silently absent. One off-scale placeholder row creates the facet; the
  # annotation is what the reader actually sees.
  excl_df <- NULL
  if (length(excluded) > 0) {
    y_mid <- exp(mean(log(range(plot_df$abs_median))))
    excl_df <- data.table(
      ef_fac   = factor(ef_facet(excluded), levels = ef_facet(EF_LEVELS)),
      position = mean(range(plot_df$position)),
      abs_median = y_mid)
  }

  p <- ggplot(plot_df[is_causal == FALSE],
            aes(x = position, y = abs_median)) +
  geom_point(colour = COL_NULL, alpha = 0.35, size = 0.7, stroke = 0) +
  geom_hline(data = truth_lines, aes(yintercept = simulated),
             linetype = "dashed", colour = "grey45", linewidth = 0.4) +
  # Causal variant drawn last so it sits on top, with a surface ring to keep it
  # legible where the null cloud is dense.
  geom_point(data = causal_df, colour = "white", size = 3.6, stroke = 0) +
  geom_point(data = causal_df, colour = COL_CAUSAL, size = 2.6, stroke = 0) +
  # The 90% credible interval, drawn to scale. It is only ~0.013 log10 units
  # tall -- narrower than the marker at this axis range -- so this layer is
  # honest but near-invisible; the zoom row below is what makes it readable.
  (if (show_ci && "ci_lo" %in% names(causal_df))
     geom_linerange(data = causal_df[!is.na(ci_lo)],
                    aes(ymin = ci_lo, ymax = ci_hi),
                    colour = COL_CAUSAL, linewidth = 0.9)
   else NULL) +
  # The causal variant is the only labelled point; every other variant stays
  # unlabelled so the label reads as a callout rather than as a gene annotation.
  geom_text(data = causal_df, label = CAUSAL_LABEL,
            colour = COL_CAUSAL, vjust = -1.1, size = 3.6,
            fontface = "bold", show.legend = FALSE) +
  # Headroom on top so the causal label is never clipped by the panel edge.
  scale_y_log10(breaks = LOG_BREAKS, labels = log10_labels,
                expand = expansion(mult = c(0.05, 0.22))) +
  scale_x_continuous(labels = scales::comma) +
  facet_wrap(~ ef_fac, nrow = 1) +
  labs(
    x = "variants, positioned by rs number",
    y = expression(log[10] * group("|", tilde(beta), "|"))
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(linewidth = 0.3, colour = "grey88"),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(size = 8)
  )

  # Recovered h2 in the top-right corner of each panel. The causal variant sits
  # around the 32nd percentile of the x range, so the right corner stays clear
  # of its callout label. parse = TRUE renders a real superscript.
  if (!is.null(h2_df) && nrow(h2_df) > 0) {
    p <- p + geom_label(
      data = h2_df,
      aes(x = Inf, y = Inf,
          # Quoted so plotmath keeps it a literal: unquoted, parse() evaluates
          # 0.900 as a number and prints "0.9", giving ragged decimal places.
          label = sprintf('italic(h)^2 == "%.3f"', h2)),
      parse = TRUE, hjust = 1.04, vjust = 1.15,
      size = 4.6, colour = "grey15", fontface = "bold",
      fill = "white", alpha = 0.9, linewidth = 0.25,
      label.padding = grid::unit(0.3, "lines"), inherit.aes = FALSE)
  }

  if (!is.null(excl_df)) {
    p <- p +
      geom_blank(data = excl_df, aes(x = position, y = abs_median)) +
      geom_text(data = excl_df, aes(x = position, y = abs_median),
                label = "fit failed\n(no usable result)", colour = "grey35",
                size = 3.6, fontface = "italic", lineheight = 0.95)
  }
  p
}

# ------------------------------------------------- causal CI zoom row --------
# The 90% credible interval on the causal effect is ~3% of the estimate, i.e.
# about 0.013 log10 units, which on the shared 14-decade axis above is thinner
# than the point marker. Drawn at that scale the interval is invisible, so the
# CI version pairs the full Manhattan with this row: one free-y panel per effect
# size, showing the interval at a readable scale. Free y is the
# point here -- these panels are deliberately not comparable to each other, only
# each interval to its own truth line.
make_ci_zoom <- function(causal_df, excluded = character(0)) {
  d <- causal_df[!is.na(ci_lo)]
  d[, x := 0]

  # Pinned to the panel top rather than to ci_hi, so it clears the interval at
  # any zoom level.
  lab <- d[, .(ef_fac, txt = sprintf("%.2f [%.2f, %.2f]", abs_median, ci_lo, ci_hi))]

  p <- ggplot(d, aes(x = x, y = abs_median)) +
    geom_linerange(aes(ymin = ci_lo, ymax = ci_hi),
                   colour = COL_CAUSAL, linewidth = 1.1) +
    geom_point(colour = "white", size = 3.8, stroke = 0) +
    geom_point(colour = COL_CAUSAL, size = 2.8, stroke = 0) +
    geom_text(data = lab, aes(x = 0, y = Inf, label = txt),
              vjust = 1.6, size = 3.4, colour = "grey20", inherit.aes = FALSE) +
    scale_x_continuous(limits = c(-1, 1), breaks = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0.35, 0.75))) +
    facet_wrap(~ ef_fac, nrow = 1, scales = "free_y", drop = FALSE) +
    labs(x = NULL, y = expression(tilde(beta) ~ "(90% CrI)")) +
    theme_minimal(base_size = 13) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(linewidth = 0.3, colour = "grey88"),
          strip.text = element_blank(),
          axis.text.x = element_blank())

  if (length(excluded) > 0) {
    excl <- data.table(
      ef_fac = factor(ef_facet(excluded), levels = ef_facet(EF_LEVELS)),
      x = 0, y = 0)
    p <- p + geom_blank(data = excl, aes(x = x, y = y)) +
      geom_text(data = excl, aes(x = x, y = y), label = "no usable fit",
                colour = "grey35", size = 3.4, fontface = "italic",
                inherit.aes = FALSE)
  }
  p
}

arms <- list()
for (mode in names(MODE_TITLES)) {
  a <- build_arm(mode)
  if (is.null(a)) next
  arms[[mode]] <- a

  p <- make_manhattan(a$plot_df, a$causal_df, mode, a$excluded, a$h2_df)
  png_path <- file.path(OUT_DIR,
                        paste0("klebsiella_effectsize_faceted_manhattan_", mode, ".png"))
  ggsave(png_path, p, width = 15, height = 4.6, dpi = 600, bg = "white")
  message("Wrote ", png_path)

  # CI version: same Manhattan with the interval drawn to scale, over the zoom
  # row that makes it legible.
  p_ci <- make_manhattan(a$plot_df, a$causal_df, mode, a$excluded, a$h2_df,
                         show_ci = TRUE)
  zoom <- make_ci_zoom(a$causal_df, a$excluded)
  combined <- patchwork::wrap_plots(p_ci, zoom, ncol = 1, heights = c(3, 1.35))
  ci_path <- file.path(
    OUT_DIR, paste0("klebsiella_effectsize_faceted_manhattan_", mode, "_ci.png"))
  ggsave(ci_path, combined, width = 15, height = 6.6, dpi = 600, bg = "white")
  message("Wrote ", ci_path)
}
if (length(arms) == 0L) {
  stop("No extracts found under ", EXTRACTED,
       " -- run paper_figures/klebsiella_extract_effects.R first")
}

# ------------------------------------------------- recovered vs simulated ----
# Companion panel, both arms on one axis: does the recovered causal effect track
# the simulated one? observed_diff is the carrier-vs-non-carrier phenotype
# difference measured directly from the simulated data, i.e. what the model
# should recover. Colour here encodes the arm, and both arms are direct-labelled
# by effect size, so identity never rests on colour alone.
truth_dt <- data.table(ef_chr = truth$ef_chr,
                       simulated = abs(truth$observed_diff),
                       effect_size = truth$effect_size)

growth <- data.table::rbindlist(lapply(names(arms), function(mode) {
  merge(arms[[mode]]$causal_df[, .(ef_chr, mode, recovered = abs_median)],
        truth_dt, by = "ef_chr")
}))

if (nrow(growth) > 0L) {
  data.table::setorder(growth, effect_size)
  ARM_COLS <- c(pruned = "#b35806", nopruning = "#542788")

  g <- ggplot(growth, aes(x = simulated, y = recovered, colour = mode)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                colour = "grey55", linewidth = 0.5) +
    geom_point(size = 3.2) +
    ggrepel::geom_text_repel(aes(label = paste0("EF ", ef_chr)),
                             size = 3.2, colour = "grey20",
                             seed = 1L, show.legend = FALSE) +
    scale_colour_manual(values = ARM_COLS,
                        labels = MODE_TITLES[names(ARM_COLS)],
                        name = NULL) +
    scale_x_log10() + scale_y_log10() +
    labs(x = "simulated effect (carrier vs non-carrier difference)",
         y = expression(paste("recovered  ", group("|", tilde(beta), "|"))),
         title = "Recovered vs simulated causal effect",
         subtitle = "Dashed line is exact recovery") +
    theme_minimal(base_size = 13) +
    theme(panel.grid.minor = element_blank(),
          legend.position = "top",
          plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(colour = "grey30", size = 10))

  g_path <- file.path(OUT_DIR, "klebsiella_causal_effect_growth.png")
  ggsave(g_path, g, width = 6.5, height = 5.4, dpi = 600, bg = "white")
  message("Wrote ", g_path)
}

# ------------------------------------------------------------ the numbers ----
# Rank of the causal variant within each panel: the honest check on whether the
# red point is actually the top hit, which the eye cannot confirm at this density.
# `mode` is a column of df as well as the loop variable, so it goes in `by`
# rather than `j` -- assigning mode = mode inside j resolves to the column and
# silently defeats the aggregation.
ranks <- data.table::rbindlist(lapply(names(arms), function(arm) {
  arms[[arm]]$df[, .(
    causal_abs_median = abs_median[is_causal],
    causal_rank       = rank(-abs_median, ties.method = "min")[is_causal],
    n_variants        = .N,
    n_larger          = sum(abs_median > abs_median[is_causal]),
    causal_signif     = signif[is_causal]
  ), by = .(mode, ef_chr, effect_size)]
}))

summary_out <- merge(ranks, truth_dt[, .(ef_chr, simulated)], by = "ef_chr")
data.table::setorder(summary_out, mode, effect_size)
summary_csv <- file.path(OUT_DIR, "klebsiella_effectsize_causal_summary.csv")
data.table::fwrite(summary_out, summary_csv)
message("Wrote ", summary_csv)

print(summary_out)
