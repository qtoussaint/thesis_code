#!/usr/bin/env Rscript
# Distribution of posterior median variant effects (beta) across the eight PPOM
# regularisation priors, on the two S. pneumoniae penicillin MIC datasets.
#
# One violin per prior, one panel per dataset, betas pooled over cutpoints. The
# upper half is a kernel density computed from every effect, with the slice
# inside the 50% credible interval shaded; the lower half is a seeded subsample
# of individual variants, drawn as an evenly jittered band. The median with the
# central 50% and 90% of the effects sits on the centre line between them.
#
# The lasso/ridge priors are the centred variants from the regularisation sweep.
# Their prediction-sweep twins were run with --norate and so have no per-variant
# coefficient files; these are the in-sample fits of the same priors.
#
# Reads:
#   <sweep>/final_ordered_categorical_PPOM_free_cutpoints_wide_drift_<suffix>/<dataset>/fitted_model/depruned_variant_effects.csv
#     suffix in {lasso, lasso_estscale, lasso_estscale_mixture,
#                ridge, ridge_estscale, ridge_estscale_ncp}
#   <results>/gwas_spn_penicillin/inference/<dataset>_PPOM/fitted_model/depruned_variant_effects.csv
#     (regularized horseshoe, estimated tau)
#   <sweep>/final_no_horseshoe_with_centering_PPOM/<dataset>/fitted_model/depruned_variant_effects.csv
#     (unregularised normal; the minima-binned table was rebuilt from the saved
#      fit by compare_ordinal_models/analysis/extract_depruned_effects_ppom.R,
#      because that run's pipeline stopped at the inference PPC step)
#
# Run with:
#   mamba activate gwas_pipeline
#   Rscript prior_predictive_checks/variant_effect_distribution_by_prior_summary.R
#
# Output: <results>/prior_predictive_checks/variant_effect_distribution_by_prior_summary.{png,csv}
#         <results>/prior_predictive_checks/variant_effect_distribution_by_prior_summary_pseudolog.png

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
  library(scales)
})

grDevices::pdf(NULL)
on.exit(grDevices::dev.off(), add = TRUE)

# -----------------------------------------------------------------------------
# Constants
# -----------------------------------------------------------------------------
RESULTS_ROOT <- "/nfs/research/jlees/jacqueline/thesis_results"
SWEEP_ROOT   <- file.path(RESULTS_ROOT, "compare_ordinal_models")
OUTPUT_DIR   <- file.path(RESULTS_ROOT, "prior_predictive_checks")
MODEL_PREFIX <- "final_ordered_categorical_PPOM_free_cutpoints_wide_drift"

# Geometry, in row units: the density fills the half above the row line and the
# variant points mirror it in the half below, both starting from the row line
# itself, where the credible-interval bar is drawn.
POINT_SEED   <- 20260728L
N_POINTS_MAX <- 12000L   # plotted points per model x dataset
CLOUD_HEIGHT <- 0.41     # height of the density at its peak
RAIN_DEPTH   <- 0.38     # depth of the point band below the row line
DENS_N       <- 4096L    # density grid points across the display window
X_QUANTILE   <- 0.99     # pooled |beta| quantile that sets the shared x limit
X_LIM_HALF   <- NULL     # set to a number to override the computed limit
PSEUDO_SIGMA <- 1e-5     # pseudo-log linear-region width, companion figure
# Outer bar spans the central 90% of the effects in a row. These are quantiles of
# the spread of point estimates across variants, not a posterior interval, so the
# conventional 5/95 pair is used rather than the 5.5/94.5 of an 89% interval.
OUTER_PROBS  <- c(0.05, 0.95)
INNER_PROBS  <- c(0.25, 0.75)     # interquartile range, the thick bar

set.seed(POINT_SEED)

# -----------------------------------------------------------------------------
# Models and datasets
# -----------------------------------------------------------------------------
MODELS <- list(
  list(key = "lasso",                  family = "lasso",
       label = "lasso",                     kind = "sweep", suffix = "lasso"),
  list(key = "lasso_estscale",         family = "lasso",
       label = "lasso · est-σ",    kind = "sweep", suffix = "lasso_estscale"),
  list(key = "lasso_estscale_mixture", family = "lasso",
       label = "lasso · est-σ · mixture",
       kind = "sweep", suffix = "lasso_estscale_mixture"),
  list(key = "ridge",                  family = "ridge",
       label = "ridge",                     kind = "sweep", suffix = "ridge"),
  list(key = "ridge_estscale",         family = "ridge",
       label = "ridge · est-σ",    kind = "sweep", suffix = "ridge_estscale"),
  list(key = "ridge_estscale_ncp",     family = "ridge",
       label = "ridge · est-σ · NCP",
       kind = "sweep", suffix = "ridge_estscale_ncp"),
  list(key = "horseshoe",              family = "horseshoe",
       label = "horseshoe (est. τ)",    kind = "horseshoe"),
  list(key = "normal",                 family = "normal",
       label = "normal (unregularised)",    kind = "normal")
)
# Row 8 is the top of the panel; y is numeric so the density ribbons and the
# interval bars can be positioned by hand.
for (i in seq_along(MODELS)) MODELS[[i]]$row <- length(MODELS) - i + 1L
row_breaks <- vapply(MODELS, `[[`, 0L, "row")
row_labels <- vapply(MODELS, `[[`, "", "label")

DATASETS <- list(
  list(id = "02_spn_penicillin_MIC",
       label = "doubling dilutions (≥5% min. frequency), K=8", K = 8L),
  list(id = "16_spn_penicillin_MIC_minimabinning",
       label = "minima, K=5",             K = 5L)
)
dataset_levels <- vapply(DATASETS, `[[`, "", "label")

# Okabe-Ito; lasso and ridge keep the colours used in the sweep figures.
family_cols <- c(lasso     = "#E69F00",
                 ridge     = "#0072B2",
                 horseshoe = "#009E73",
                 normal    = "#999999")
family_levels <- c("lasso", "ridge", "horseshoe", "normal")

effects_path <- function(mdl, ds_id)
  switch(mdl$kind,
         sweep     = file.path(SWEEP_ROOT, paste0(MODEL_PREFIX, "_", mdl$suffix),
                               ds_id, "fitted_model", "depruned_variant_effects.csv"),
         horseshoe = file.path(RESULTS_ROOT, "gwas_spn_penicillin", "inference",
                               paste0(ds_id, "_PPOM"), "fitted_model",
                               "depruned_variant_effects.csv"),
         normal    = file.path(SWEEP_ROOT, "final_no_horseshoe_with_centering_PPOM",
                               ds_id, "fitted_model", "depruned_variant_effects.csv"))

# -----------------------------------------------------------------------------
# Reader
# -----------------------------------------------------------------------------
missing_files <- character(0)

# Reads one depruned_variant_effects.csv, taking only the columns we need.
# fread turns the "true"/"false" text of signif into a logical, while read.csv
# leaves it as character, so both representations are handled here.
read_effects <- function(csv_path) {
  if (!file.exists(csv_path)) {
    missing_files <<- c(missing_files, csv_path)
    return(NULL)
  }
  dt <- tryCatch(
    data.table::fread(csv_path, select = c("variant_id", "median", "signif"),
                      showProgress = FALSE),
    error = function(e) NULL)
  if (is.null(dt) || nrow(dt) == 0L) {
    missing_files <<- c(missing_files, csv_path)
    return(NULL)
  }
  s <- dt$signif
  dt[, signif := if (is.logical(s)) s else tolower(trimws(as.character(s))) == "true"]
  dt[is.finite(median)]
}

report_missing <- function() {
  if (length(missing_files) > 0) {
    message(sprintf("NOTE: %d input CSV(s) missing (rendered as blank rows):",
                    length(missing_files)))
    for (f in missing_files) message("  ", f)
  } else {
    message("All input CSVs found.")
  }
}

# -----------------------------------------------------------------------------
# Load every model x dataset cell
# -----------------------------------------------------------------------------
cells <- list()
for (mdl in MODELS) {
  for (ds in DATASETS) {
    path <- effects_path(mdl, ds$id)
    dt   <- read_effects(path)
    cells[[length(cells) + 1]] <- list(
      model = mdl, dataset = ds, path = path,
      beta       = if (is.null(dt)) numeric(0) else dt$median,
      n_signif   = if (is.null(dt)) NA_integer_ else sum(dt$signif, na.rm = TRUE),
      n_variants = if (is.null(dt)) NA_integer_ else data.table::uniqueN(dt$variant_id),
      missing    = is.null(dt))
    rm(dt)
  }
}

pooled <- unlist(lapply(cells, `[[`, "beta"), use.names = FALSE)
if (length(pooled) == 0L) stop("no variant effects were read; nothing to plot")

# -----------------------------------------------------------------------------
# Shared x limits and clipping accounting
# -----------------------------------------------------------------------------
nice_ceiling <- function(x) {
  e <- floor(log10(x))
  m <- x / 10^e
  10^e * c(1, 2, 5, 10)[which(m <= c(1, 2, 5, 10))[1]]
}

xhalf <- if (is.null(X_LIM_HALF)) {
  nice_ceiling(stats::quantile(abs(pooled), X_QUANTILE, names = FALSE))
} else {
  X_LIM_HALF
}
XLIM <- c(-xhalf, xhalf)

pct_clipped_all <- 100 * sum(abs(pooled) > xhalf) / length(pooled)
message(sprintf("shared x limits: [%g, %g] (%.2f%% of all betas outside)",
                XLIM[1], XLIM[2], pct_clipped_all))

# -----------------------------------------------------------------------------
# Density on the display window only, with all values setting the scale
# -----------------------------------------------------------------------------
# The bandwidth and the weights come from every value, so the curve integrates to
# the fraction of mass inside the window: clipped tails lower the curve instead of
# being renormalised away. The bandwidth floor keeps the most heavily shrunk
# priors from drawing as an invisible needle.
window_density <- function(x, lo, hi, n = DENS_N) {
  bw_raw <- stats::bw.nrd0(x)
  bw     <- max(bw_raw, (hi - lo) / 500)
  pad    <- 4 * bw
  keep   <- x >= lo - pad & x <= hi + pad
  if (!any(keep)) return(NULL)
  d <- suppressWarnings(stats::density(
    x[keep], bw = bw, weights = rep(1 / length(x), sum(keep)),
    from = lo - pad, to = hi + pad, n = n, cut = 0))
  list(x = d$x, y = d$y, bw = bw, bw_floored = bw > bw_raw,
       mass_in_window = sum(d$y) * diff(d$x[1:2]))
}

# -----------------------------------------------------------------------------
# Build the plotting frames
# -----------------------------------------------------------------------------
pts_l <- list(); rib_l <- list(); rib_pl_l <- list()
int_l <- list(); blank_l <- list(); summary_l <- list()

ptrans <- scales::pseudo_log_trans(sigma = PSEUDO_SIGMA, base = 10)
t_pooled <- ptrans$transform(pooled)
TLIM <- range(t_pooled)

message(sprintf("%-30s %-22s %10s %10s %8s",
                "model", "dataset", "n", "n_clipped", "pct"))

for (cell in cells) {
  mdl <- cell$model; ds <- cell$dataset; b <- cell$beta

  if (cell$missing) {
    # All 16 tables are present as of the minima rebuild; this keeps the row
    # labelled rather than silently absent should one ever go missing again.
    blank_l[[length(blank_l) + 1]] <- data.frame(
      row = mdl$row, dataset = ds$label, note = "effects table not written",
      stringsAsFactors = FALSE)
    summary_l[[length(summary_l) + 1]] <- data.frame(
      model = mdl$key, model_label = mdl$label, family = mdl$family, row = mdl$row,
      dataset = ds$id, dataset_label = ds$label, K = ds$K, n_cutpoints = ds$K - 1L,
      input_path = cell$path, file_missing = TRUE,
      n = NA_integer_, n_variants = NA_integer_,
      median = NA_real_, mean = NA_real_, sd = NA_real_, min = NA_real_, max = NA_real_,
      q05 = NA_real_, q95 = NA_real_, q25 = NA_real_, q75 = NA_real_,
      abs_q50 = NA_real_, abs_q90 = NA_real_, abs_q95 = NA_real_, abs_q99 = NA_real_,
      n_signif = NA_integer_, prop_signif = NA_real_,
      x_lim_low = XLIM[1], x_lim_high = XLIM[2],
      n_clipped_low = NA_integer_, n_clipped_high = NA_integer_,
      n_clipped = NA_integer_, pct_clipped = NA_real_,
      n_points_plotted = 0L, point_seed = POINT_SEED,
      bw_used = NA_real_, bw_floored = NA, density_mass_in_window = NA_real_,
      pseudo_log_sigma = PSEUDO_SIGMA, stringsAsFactors = FALSE)
    next
  }

  # point cloud: one seeded subsample, jittered evenly through the band below the
  # row line. Both figures share it, since the depth does not depend on x.
  idx <- if (length(b) > N_POINTS_MAX)
           sort(sample.int(length(b), N_POINTS_MAX)) else seq_along(b)
  pts_l[[length(pts_l) + 1]] <- data.frame(
    x = b[idx], row = mdl$row, family = mdl$family, dataset = ds$label,
    y = mdl$row - stats::runif(length(idx), 0, RAIN_DEPTH),
    stringsAsFactors = FALSE)

  # density on the linear display window
  d <- window_density(b, XLIM[1], XLIM[2])
  if (!is.null(d)) {
    keep <- d$x >= XLIM[1] & d$x <= XLIM[2]
    h <- d$y[keep] / max(d$y)
    rib_l[[length(rib_l) + 1]] <- data.frame(
      x = d$x[keep], row = mdl$row, family = mdl$family, dataset = ds$label,
      ymin = mdl$row,
      ymax = mdl$row + CLOUD_HEIGHT * h,
      stringsAsFactors = FALSE)
  }

  # density in pseudo-log space, then mapped back for the companion figure
  d_pl <- window_density(ptrans$transform(b), TLIM[1], TLIM[2])
  if (!is.null(d_pl)) {
    keep <- d_pl$x >= TLIM[1] & d_pl$x <= TLIM[2]
    h_pl <- d_pl$y[keep] / max(d_pl$y)
    rib_pl_l[[length(rib_pl_l) + 1]] <- data.frame(
      x = ptrans$inverse(d_pl$x[keep]), row = mdl$row, family = mdl$family,
      dataset = ds$label,
      ymin = mdl$row,
      ymax = mdl$row + CLOUD_HEIGHT * h_pl,
      stringsAsFactors = FALSE)
  }

  qs  <- stats::quantile(b, OUTER_PROBS, names = FALSE)
  qs50 <- stats::quantile(b, INNER_PROBS, names = FALSE)
  aqs <- stats::quantile(abs(b), c(0.5, 0.9, 0.95, 0.99), names = FALSE)
  int_l[[length(int_l) + 1]] <- data.frame(
    row = mdl$row, family = mdl$family, dataset = ds$label,
    median = stats::median(b), q_lo = qs[1], q_hi = qs[2],
    q50_lo = qs50[1], q50_hi = qs50[2],
    abs_med = aqs[1], stringsAsFactors = FALSE)

  n_lo <- sum(b < XLIM[1]); n_hi <- sum(b > XLIM[2])
  message(sprintf("%-30s %-22s %10d %10d %7.2f%%",
                  mdl$label, ds$label, length(b), n_lo + n_hi,
                  100 * (n_lo + n_hi) / length(b)))

  summary_l[[length(summary_l) + 1]] <- data.frame(
    model = mdl$key, model_label = mdl$label, family = mdl$family, row = mdl$row,
    dataset = ds$id, dataset_label = ds$label, K = ds$K, n_cutpoints = ds$K - 1L,
    input_path = cell$path, file_missing = FALSE,
    n = length(b), n_variants = cell$n_variants,
    median = stats::median(b), mean = mean(b), sd = stats::sd(b),
    min = min(b), max = max(b), q05 = qs[1], q95 = qs[2],
    q25 = qs50[1], q75 = qs50[2],
    abs_q50 = aqs[1], abs_q90 = aqs[2], abs_q95 = aqs[3], abs_q99 = aqs[4],
    n_signif = cell$n_signif, prop_signif = cell$n_signif / length(b),
    x_lim_low = XLIM[1], x_lim_high = XLIM[2],
    n_clipped_low = n_lo, n_clipped_high = n_hi, n_clipped = n_lo + n_hi,
    pct_clipped = 100 * (n_lo + n_hi) / length(b),
    n_points_plotted = length(idx), point_seed = POINT_SEED,
    bw_used = if (is.null(d)) NA_real_ else d$bw,
    bw_floored = if (is.null(d)) NA else d$bw_floored,
    density_mass_in_window = if (is.null(d)) NA_real_ else d$mass_in_window,
    pseudo_log_sigma = PSEUDO_SIGMA, stringsAsFactors = FALSE)
}

pts    <- do.call(rbind, pts_l)
ribbon <- do.call(rbind, rib_l)
rib_pl <- do.call(rbind, rib_pl_l)
ints   <- do.call(rbind, int_l)
blanks <- if (length(blank_l)) do.call(rbind, blank_l) else NULL
summary_df <- do.call(rbind, summary_l)

pts$dataset    <- factor(pts$dataset,    levels = dataset_levels)
ribbon$dataset <- factor(ribbon$dataset, levels = dataset_levels)
rib_pl$dataset <- factor(rib_pl$dataset, levels = dataset_levels)
ints$dataset   <- factor(ints$dataset,   levels = dataset_levels)
if (!is.null(blanks)) blanks$dataset <- factor(blanks$dataset, levels = dataset_levels)

# median |beta| annotation, which is what keeps the most shrunk rows readable
# Three significant figures everywhere, with anything below 0.01 written as a
# mantissa times a power of ten rather than in e-notation. Numbers are quoted so
# plotmath renders them as written; unquoted they are re-parsed as numeric
# literals and the trailing zeros are lost (0.680 prints as 0.68).
format_effect <- function(x) {
  x     <- signif(x, 3)
  small <- abs(x) < 0.01
  out   <- character(length(x))
  out[!small] <- paste0("'", formatC(x[!small], format = "g", digits = 3,
                                     flag = "#"), "'")
  if (any(small)) {
    e <- floor(log10(abs(x[small])))
    m <- formatC(x[small] / 10^e, format = "f", digits = 2)
    out[small] <- paste0("'", m, "' %.% 10^", e)
  }
  out
}
ints$scale_note <- paste0("med.~group('|', tilde(beta), '|') == ",
                          format_effect(ints$abs_med))

# The slice of each violin lying inside the 50% credible interval, shaded darker.
shade_50 <- function(rib) {
  m <- merge(rib, ints[, c("row", "dataset", "q50_lo", "q50_hi")],
             by = c("row", "dataset"))
  m <- m[m$x >= m$q50_lo & m$x <= m$q50_hi, ]
  m[order(m$row, m$dataset, m$x), ]
}
ribbon50 <- shade_50(ribbon)
rib_pl50 <- shade_50(rib_pl)

# The two edges of that slice, so the 50% bounds stay visible through the points.
shade_edges <- function(sh) {
  parts <- split(sh, list(sh$row, sh$dataset), drop = TRUE)
  do.call(rbind, lapply(parts, function(g) {
    g <- g[c(1L, nrow(g)), ]
    data.frame(row = g$row, dataset = g$dataset, x = g$x, ytop = g$ymax,
               stringsAsFactors = FALSE)
  }))
}
edges50    <- shade_edges(ribbon50)
edges_pl50 <- shade_edges(rib_pl50)

# -----------------------------------------------------------------------------
# Figure
# -----------------------------------------------------------------------------
base_theme <- theme_bw(base_size = 15) +
  theme(panel.grid.minor   = element_blank(),
        panel.grid.major.y = element_blank(),
        strip.background   = element_rect(fill = "grey92", colour = NA),
        legend.position    = "bottom",
        plot.subtitle      = element_text(size = 12, colour = "grey30"),
        panel.spacing.x    = grid::unit(1.1, "lines"),
        axis.text.x        = element_text(size = 12),
        axis.text.y        = element_text(hjust = 1))

build_figure <- function(point_df, ribbon_df, shade_df, edge_df, int_df, xscale,
                         xlim, subtitle, note_x, note_hjust) {
  p <- ggplot() +
    geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.3) +
    # violin body, then the 50% slice shaded on top of it
    geom_ribbon(data = ribbon_df,
                aes(x = x, ymin = ymin, ymax = ymax,
                    group = interaction(row, dataset), fill = family),
                alpha = 0.13) +
    geom_ribbon(data = shade_df,
                aes(x = x, ymin = ymin, ymax = ymax,
                    group = interaction(row, dataset), fill = family),
                alpha = 0.5, show.legend = FALSE) +
    # variants, mirrored into the lower half of the same outline
    geom_point(data = point_df, aes(x = x, y = y, colour = family),
               alpha = 0.18, size = 0.5, stroke = 0, shape = 16,
               show.legend = FALSE) +
    # edges of the 50% interval
    geom_segment(data = edge_df,
                 aes(x = x, xend = x, y = row, yend = ytop),
                 colour = "grey20", linewidth = 0.35, linetype = "22") +
    geom_line(data = ribbon_df,
              aes(x = x, y = ymax, group = interaction(row, dataset),
                  colour = family),
              linewidth = 0.5) +
    # central 90% as a thin line, central 50% as a thick bar, median as a dot
    geom_linerange(data = int_df, aes(y = row, xmin = q_lo, xmax = q_hi),
                   colour = "grey15", linewidth = 0.6) +
    geom_linerange(data = int_df, aes(y = row, xmin = q50_lo, xmax = q50_hi),
                   colour = "grey15", linewidth = 2.4) +
    geom_point(data = int_df, aes(x = median, y = row),
               colour = "grey15", fill = "white", shape = 21, size = 2.2,
               stroke = 0.7) +
    geom_text(data = int_df,
              aes(x = note_x, y = row + CLOUD_HEIGHT,
                  label = scale_note),
              hjust = note_hjust, vjust = 1, size = 3.8, colour = "grey35",
              parse = TRUE)
  if (!is.null(blanks))
    p <- p + geom_text(data = blanks, aes(x = 0, y = row, label = note),
                       size = 4.1, colour = "grey45", fontface = "italic")
  p <- p +
    facet_wrap(~ dataset, nrow = 1) +
    # Each row spans row-RAIN_DEPTH .. row+CLOUD_HEIGHT, so only a sliver of
    # padding is needed; a larger expansion just leaves the panel half empty.
    scale_y_continuous(breaks = row_breaks, labels = row_labels,
                       expand = expansion(add = 0.1)) +
    scale_colour_manual(values = family_cols, breaks = family_levels,
                        name = "prior family") +
    scale_fill_manual(values = family_cols, breaks = family_levels,
                      name = "prior family") +
    # Outline each key in its own family colour, so the swatch carries the same
    # darker border the violins do.
    guides(fill = guide_legend(override.aes = list(
             alpha = 0.5, colour = unname(family_cols[family_levels]),
             linewidth = 0.4)),
           colour = "none") +
    labs(x = expression("posterior median variant effect (" * tilde(beta) * ")"),
         y = NULL, subtitle = subtitle) +
    xscale +
    base_theme
  if (!is.null(xlim)) p <- p + coord_cartesian(xlim = xlim)
  p
}

p_lin <- build_figure(
  pts, ribbon, ribbon50, edges50, ints,
  xscale = scale_x_continuous(breaks = scales::pretty_breaks(7),
                              expand = expansion(mult = 0.01)),
  xlim = XLIM,
  subtitle = NULL,   # clipping and subsampling are stated in the caption and CSV
  note_x = XLIM[2], note_hjust = 1)

# The +/-1e-4 breaks collide with 0 at this font size, so they are left out.
pl_breaks <- c(-10000, -100, -1, -0.01, 0, 0.01, 1, 100)
p_log <- build_figure(
  pts, rib_pl, rib_pl50, edges_pl50, ints,
  xscale = scale_x_continuous(transform = ptrans, breaks = pl_breaks,
                              labels = function(v) formatC(v, format = "g"),
                              expand = expansion(mult = 0.01)),
  xlim = NULL,
  subtitle = sprintf(paste0("signed pseudo-log x axis (sigma = %g); no effects ",
                            "clipped. Companion to the linear figure."),
                     PSEUDO_SIGMA),
  note_x = ptrans$inverse(TLIM[2]), note_hjust = 1)

# -----------------------------------------------------------------------------
# Write
# -----------------------------------------------------------------------------
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
png_path   <- file.path(OUTPUT_DIR, "variant_effect_distribution_by_prior_summary.png")
pslog_path <- file.path(OUTPUT_DIR, "variant_effect_distribution_by_prior_summary_pseudolog.png")
csv_path   <- file.path(OUTPUT_DIR, "variant_effect_distribution_by_prior_summary.csv")

ggsave(png_path,   p_lin, width = 12, height = 9.0, dpi = 300, bg = "white")
ggsave(pslog_path, p_log, width = 12, height = 9.0, dpi = 300, bg = "white")
write.csv(summary_df[order(-summary_df$row, summary_df$dataset), ],
          csv_path, row.names = FALSE)

message("wrote ", png_path)
message("wrote ", pslog_path)
message("wrote ", csv_path)
report_missing()
