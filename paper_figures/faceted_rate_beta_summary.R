#!/usr/bin/env Rscript
# Faceted RATE + |median effect| summary figures, one set per drug/species, reusing the
# labeling/faceting scheme from replot_rate_faceted_cutpoints.R (plasma points, top-10
# italic gene labels per facet, theme_minimal base_size 20).
#
# Per species (tb_rifampicin, spn_penicillin, spn_trimethoprim):
#   <species>_POM_faceted_RATE_beta.png  -- |median beta| (top row) and RATE (bottom row),
#                                           one column per binning strategy (POM: no cutpoints)
#   <species>_PPOM_faceted_RATE.png      -- per-binning blocks stacked vertically, each
#                                           faceted by cutpoint, RATE on y
#   <species>_PPOM_faceted_beta.png      -- same blocks, |median beta| on y
# Plus one figure across all species:
#   logistic_all_species_faceted_RATE_beta.png -- logistic |median beta| (top) and RATE
#                                                  (bottom), one column per species
#
# When RATE and beta share a figure, beta is the top row and RATE the bottom row.
# Unitig runs (17/18) are excluded, matching the other paper_figures summaries.
#
# Usage:
#   mamba activate gwas_pipeline
#   Rscript paper_figures/faceted_rate_beta_summary.R
#
# Output: <RES>/paper_figures/faceted_rate_beta/

suppressPackageStartupMessages({
  library(ggplot2)
  library(viridis)
  library(ggrepel)
  library(cowplot)
})

# Null device so any incidental plotting doesn't leave a stray Rplots.pdf behind.
grDevices::pdf(NULL)
on.exit(grDevices::dev.off(), add = TRUE)

# Reuse the per-run helpers (sourcing does not trigger that script's main()):
# read_positions, read_rate, annotate_genes, top_label_idx, .italic_gene_expr,
# build_df, build_faceted_plot.
source(file.path(dirname(sub("^--file=", "",
        grep("^--file=", commandArgs(FALSE), value = TRUE)[1])),
        "replot_rate_faceted_cutpoints.R"))

RES        <- "/nfs/research/jlees/jacqueline/thesis_results"
SPN_ANNOT  <- "/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/genotype/fields_filtered_maf05_multiallelic.txt"
TB_ANNOT   <- "/nfs/research/jlees/jacqueline/gwas_data/tuberculosis/cryptic_regeno_snpeff/fields_filtered.txt"
GOI_DIR    <- "/nfs/research/jlees/jacqueline/thesis_code/gwas_genesofinterest"
OUT_DIR    <- file.path(RES, "paper_figures", "faceted_rate_beta")
NCOL       <- 2L

SINGLE_COLOR <- viridis::viridis(6, option = "plasma")[2]
# plotmath strings (parsed in the facet strips): beta -> |beta-tilde|, rate -> plain text.
METRIC_LABELS <- c(beta = "group('|', tilde(beta), '|')",
                   rate = "'relative centrality (RATE)'")
metric_labeller <- ggplot2::as_labeller(METRIC_LABELS, ggplot2::label_parsed)

# Dataset map (same binning set as cutpoints_histogram_summary.R; unitigs excluded).
species_cfgs <- list(
  list(key = "tb_rifampicin", display = "bolditalic('M. tuberculosis')~bold('(rifampicin)')",
       annot = TB_ANNOT, goi = file.path(GOI_DIR, "tb_rifampicin_genesofinterest.txt"),
       logistic = "07_tb_rifampicin_binary_logistic", n_label = 20L,
       binnings = list(
         list(label = "standard MIC", nn = "08", stub = "tb_rifampicin_MIC"),
         list(label = "coarse",       nn = "14", stub = "tb_rifampicin_MIC_coarse_dilutions"),
         list(label = "large minbin", nn = "15", stub = "tb_rifampicin_MIC_large_minbin"))),
  list(key = "spn_penicillin", display = "bolditalic('S. pneumoniae')~bold('(penicillin)')",
       annot = SPN_ANNOT, goi = file.path(GOI_DIR, "spn_penicillin_genesofinterest.txt"),
       logistic = "01_spn_penicillin_binary_logistic", n_label = 10L,
       binnings = list(
         list(label = "standard MIC", nn = "02", stub = "spn_penicillin_MIC"),
         list(label = "coarse",       nn = "10", stub = "spn_penicillin_MIC_coarse_dilutions"),
         list(label = "large minbin", nn = "11", stub = "spn_penicillin_MIC_large_minbin"),
         list(label = "minima",       nn = "16", stub = "spn_penicillin_MIC_minimabinning")),
       # spn_penicillin's four binnings are split across two figures.
       groups = list(
         list(suffix = "standard_minima",     labels = c("standard MIC", "minima")),
         list(suffix = "coarse_largeminbin",  labels = c("coarse", "large minbin")))),
  list(key = "spn_trimethoprim", display = "bolditalic('S. pneumoniae')~bold('(trimethoprim)')",
       annot = SPN_ANNOT, goi = file.path(GOI_DIR, "spn_trimethoprim_genesofinterest.txt"),
       logistic = "04_spn_trimethoprim_binary_logistic", n_label = 10L,
       binnings = list(
         list(label = "standard MIC", nn = "05", stub = "spn_trimethoprim_MIC"),
         list(label = "coarse",       nn = "12", stub = "spn_trimethoprim_MIC_coarse_dilutions"),
         list(label = "large minbin", nn = "13", stub = "spn_trimethoprim_MIC_large_minbin")))
)

run_dir <- function(key, dir_name) file.path(RES, paste0("gwas_", key), "inference", dir_name)

# -----------------------------------------------------------------------------
# Single-RATE runs (POM, logistic): one RATE + one median per depruned variant.
# Row order of phandango.plot / RATE_values_depruned.txt / depruned_variant_effects.csv
# all match by depruned-variant order, as in build_df().
# -----------------------------------------------------------------------------
# TRUE if the single-RATE outputs (positions + RATE) exist; runs still in progress have
# an empty cppRATE_results and are skipped.
single_has_rate <- function(rdir) {
  file.exists(file.path(rdir, "cppRATE_results", "phandango.plot")) &&
    file.exists(file.path(rdir, "cppRATE_results", "RATE_values_depruned.txt"))
}

# median (effect size) only exists once depruned_variant_effects.csv is regenerated for
# the POM run; absent -> NA, so the beta panel for that binning is simply empty.
read_single_run <- function(rdir) {
  pos  <- read_positions(file.path(rdir, "cppRATE_results", "phandango.plot"))
  rate <- read_rate(file.path(rdir, "cppRATE_results", "RATE_values_depruned.txt"))
  eff_path <- file.path(rdir, "fitted_model", "depruned_variant_effects.csv")
  med <- if (file.exists(eff_path)) read.csv(eff_path)$median else rep(NA_real_, length(rate))
  if (length(pos) != length(rate) || length(pos) != length(med)) {
    stop("Length mismatch in ", rdir, ": pos=", length(pos),
         " rate=", length(rate), " median=", length(med))
  }
  data.frame(pos = pos, rate = rate, median = med)
}

# Long beta/rate frame for a single annotated run-frame (cols pos, rate, median, gene).
# Beta rows with no effect data (NA) are dropped so their facet panel stays empty.
to_long_metric <- function(df) {
  long <- rbind(
    transform(df, metric = "beta", value = abs(median)),
    transform(df, metric = "rate", value = rate))
  long <- long[!is.na(long$value), ]
  long$metric <- factor(long$metric, levels = c("beta", "rate"))
  long
}

# Per-facet top-n labels for a long frame, split by the given grouping columns.
facet_labels <- function(long, group_cols, n = 10L) {
  parts <- split(long, long[group_cols], drop = TRUE)
  label_df <- do.call(rbind, lapply(parts, function(d) {
    d[top_label_idx(d$gene, d$value, n = n), ]
  }))
  label_df$gene_expr <- .italic_gene_expr(label_df$gene)
  label_df
}

repel_layer <- function(label_df) {
  ggrepel::geom_text_repel(
    data = label_df, ggplot2::aes(x = pos, y = value, label = gene_expr),
    parse = TRUE, size = 5.1,
    arrow = grid::arrow(length = grid::unit(0.01, "npc"), type = "open"),
    colour = "black", inherit.aes = FALSE, max.overlaps = 30)
}

# -----------------------------------------------------------------------------
# (A) POM combined: beta row on top, RATE row on bottom, column per binning.
# -----------------------------------------------------------------------------
build_pom_combined <- function(cfg, binnings) {
  frames <- list()
  for (b in binnings) {
    rdir <- run_dir(cfg$key, paste0(b$nn, "_", b$stub, "_POM"))
    if (!single_has_rate(rdir)) {
      message("    skip POM binning '", b$label, "' (no RATE output yet): ", basename(rdir))
      next
    }
    df <- annotate_genes(read_single_run(rdir), cfg$annot, cfg$goi)
    df$binning <- b$label
    frames[[length(frames) + 1]] <- df
  }
  if (length(frames) == 0L) return(NULL)
  bin_labels <- vapply(frames, function(d) d$binning[1], character(1))
  df <- do.call(rbind, frames)
  df$binning <- factor(df$binning, levels = bin_labels)

  long <- to_long_metric(df)
  long$binning <- factor(long$binning, levels = bin_labels)
  label_df <- facet_labels(long, c("binning", "metric"), n = cfg$n_label)

  ggplot2::ggplot(long, ggplot2::aes(x = pos, y = value)) +
    ggplot2::geom_point(alpha = 0.4, colour = SINGLE_COLOR) +
    ggplot2::facet_grid(metric ~ binning, scales = "free_y", switch = "y",
                        labeller = ggplot2::labeller(metric = metric_labeller)) +
    repel_layer(label_df) +
    ggplot2::xlab("genome coordinate (bp)") + ggplot2::ylab(NULL) +
    ggplot2::theme_minimal(base_size = 20) +
    ggplot2::theme(strip.placement = "outside")
}

# -----------------------------------------------------------------------------
# (B/C) PPOM stacked blocks: one titled block per binning, faceted by cutpoint.
# -----------------------------------------------------------------------------
build_ppom_stacked <- function(cfg, metric, binnings) {
  blocks <- list(); rel_h <- numeric(0)
  for (b in binnings) {
    rdir <- run_dir(cfg$key, paste0(b$nn, "_", b$stub, "_PPOM"))
    if (length(Sys.glob(file.path(rdir, "cppRATE_results",
                                  "RATE_values_cutpoint*_depruned.txt"))) == 0L) {
      message("    skip PPOM binning '", b$label, "' (no cutpoint RATE output yet): ",
              basename(rdir))
      next
    }
    built <- build_faceted_plot(rdir, cfg$annot, cfg$goi, NCOL, metric = metric,
                                n_label = cfg$n_label)
    rows <- ceiling(built$n / NCOL)
    title_h <- 0.5; plot_h <- 4.5 * rows
    title <- cowplot::ggdraw() +
      cowplot::draw_label(b$label, fontface = "bold", size = 26, x = 0.01, hjust = 0)
    blocks[[length(blocks) + 1]] <- cowplot::plot_grid(
      title, built$plot, ncol = 1, rel_heights = c(title_h, plot_h))
    rel_h <- c(rel_h, title_h + plot_h)
  }
  if (length(blocks) == 0L) return(NULL)
  list(figure = cowplot::plot_grid(plotlist = blocks, ncol = 1, rel_heights = rel_h),
       height = sum(rel_h))
}

# One logistic block: metric (rows) x species (cols) facet over a long sub-frame.
# Used for both the 2-species SPN block and the single-species TB beta/rate rows.
logi_panel <- function(long, n = 10L) {
  label_df <- facet_labels(long, c("species", "metric"), n = n)
  ggplot2::ggplot(long, ggplot2::aes(x = pos, y = value)) +
    ggplot2::geom_point(alpha = 0.4, colour = SINGLE_COLOR) +
    ggplot2::facet_grid(metric ~ species, scales = "free", switch = "y",
                        labeller = ggplot2::labeller(metric = metric_labeller,
                                                     species = ggplot2::label_parsed)) +
    repel_layer(label_df) +
    ggplot2::xlab("genome coordinate (bp)") + ggplot2::ylab(NULL) +
    ggplot2::theme_minimal(base_size = 20) +
    ggplot2::theme(strip.placement = "outside")
}

# -----------------------------------------------------------------------------
# (D) Logistic figure. The two SPN species share a 2x2 block (beta row on top,
# RATE row below); TB (much wider genome) gets a full-width beta row and a
# full-width RATE row stacked beneath. Returns {figure, height}.
# -----------------------------------------------------------------------------
build_logistic_all <- function(cfgs) {
  frames <- list()
  for (cfg in cfgs) {
    rdir <- run_dir(cfg$key, cfg$logistic)
    if (!single_has_rate(rdir)) {
      message("    skip logistic '", cfg$key, "' (no RATE output yet)")
      next
    }
    df <- annotate_genes(read_single_run(rdir), cfg$annot, cfg$goi)
    df$species <- cfg$display
    df$species_key <- cfg$key
    frames[[length(frames) + 1]] <- df
  }
  if (length(frames) == 0L) return(NULL)
  long <- to_long_metric(do.call(rbind, frames))

  n_label <- stats::setNames(vapply(cfgs, function(c) c$n_label, integer(1)),
                             vapply(cfgs, function(c) c$key, character(1)))
  spn_keys <- c("spn_penicillin", "spn_trimethoprim")
  spn_labels <- vapply(cfgs[vapply(cfgs, function(c) c$key %in% spn_keys, logical(1))],
                       function(c) c$display, character(1))

  blocks <- list(); rel_h <- numeric(0)

  spn_long <- long[long$species_key %in% spn_keys, ]
  if (nrow(spn_long) > 0L) {
    spn_long$species <- factor(spn_long$species, levels = spn_labels)
    blocks[[length(blocks) + 1]] <- logi_panel(spn_long, n = n_label[["spn_penicillin"]])
    rel_h <- c(rel_h, 2)  # two metric rows
  }

  # TB beta then TB RATE, each a full-width single-metric row with TB's label count.
  tb_long <- long[long$species_key == "tb_rifampicin", ]
  for (m in c("beta", "rate")) {
    sub <- tb_long[tb_long$metric == m, ]
    if (nrow(sub) == 0L) next
    sub$species <- factor(sub$species)
    blocks[[length(blocks) + 1]] <- logi_panel(sub, n = n_label[["tb_rifampicin"]])
    rel_h <- c(rel_h, 1)
  }

  if (length(blocks) == 0L) return(NULL)
  list(figure = cowplot::plot_grid(plotlist = blocks, ncol = 1, rel_heights = rel_h),
       height = 5.5 * sum(rel_h))
}

# -----------------------------------------------------------------------------
# Drive
# -----------------------------------------------------------------------------
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

save_or_skip <- function(fig, path, width, height) {
  if (is.null(fig)) { message("    -> nothing to plot, skipped ", basename(path)); return(invisible()) }
  ggplot2::ggsave(path, fig, width = width, height = height,
                  dpi = 300, bg = "white", limitsize = FALSE)
  message("    wrote ", basename(path))
}

for (cfg in species_cfgs) {
  # Each species produces one figure-set per binning group; species without an
  # explicit `groups` get a single set covering all their binnings (no suffix).
  groups <- cfg$groups
  if (is.null(groups)) {
    groups <- list(list(suffix = NULL,
                        labels = vapply(cfg$binnings, function(b) b$label, character(1))))
  }

  for (g in groups) {
    binnings <- Filter(function(b) b$label %in% g$labels, cfg$binnings)
    sfx <- if (is.null(g$suffix)) "" else paste0("_", g$suffix)
    tag <- paste0(cfg$key, sfx)
    n_bin <- length(binnings)

    message("[", tag, "] POM combined")
    pom <- build_pom_combined(cfg, binnings)
    save_or_skip(pom, file.path(OUT_DIR, paste0(tag, "_POM_faceted_RATE_beta.png")),
                 max(20, 6 * n_bin), 11)

    message("[", tag, "] PPOM RATE")
    rate <- build_ppom_stacked(cfg, "rate", binnings)
    save_or_skip(rate$figure, file.path(OUT_DIR, paste0(tag, "_PPOM_faceted_RATE.png")),
                 22, if (is.null(rate)) 0 else rate$height)

    message("[", tag, "] PPOM beta")
    beta <- build_ppom_stacked(cfg, "abs_beta", binnings)
    save_or_skip(beta$figure, file.path(OUT_DIR, paste0(tag, "_PPOM_faceted_beta.png")),
                 22, if (is.null(beta)) 0 else beta$height)
  }
}

message("[all species] logistic RATE + beta")
logi <- build_logistic_all(species_cfgs)
save_or_skip(logi$figure, file.path(OUT_DIR, "logistic_all_species_faceted_RATE_beta.png"),
             22, if (is.null(logi)) 0 else logi$height)

message("Wrote figures to ", OUT_DIR)
