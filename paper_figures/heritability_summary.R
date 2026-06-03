#!/usr/bin/env Rscript
# Summary figures of inferred heritability (h2), one per dataset (SPN penicillin,
# SPN trimethoprim, TB rifampicin). Each figure compares the MIC binning setups
# and the POM vs PPOM ordinal models (panel A), with the logistic and continuous
# non-ordered models shown alongside (panel B).
#
# POM reports a single h2 per model; PPOM reports one h2 per cutpoint. Both the
# narrow-sense (h2_narrow) and broad-sense (h2_broad) estimates are shown.
#
# Reads the per-run heritability_summary.csv files written by the inference
# pipeline (median + 95% CI from q_lower/q_upper).
#
# Run with:
#   mamba activate gwas_pipeline
#   Rscript paper_figures/heritability_summary.R
#
# Output: paper_figures/heritability_plots/heritability_summary_<dataset>.{png,csv}

suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
  library(viridis)
})

# Null device so any incidental plotting doesn't leave a stray Rplots.pdf behind.
grDevices::pdf(NULL)
on.exit(grDevices::dev.off(), add = TRUE)

RESULTS_ROOT <- "/nfs/research/jlees/jacqueline/thesis_results"
OUTPUT_DIR   <- file.path(RESULTS_ROOT, "paper_figures", "heritability_plots")

KEEP <- c("metric", "cutpoint", "cutpoint_mic", "median", "q_lower", "q_upper")

# -----------------------------------------------------------------------------
# Dataset specs: one entry per dataset. binning_specs fits both POM and PPOM;
# baseline_specs are the single-h2 non-ordered models. `title` is a plotmath
# expression for the figure header (full drug name + italicised species).
# -----------------------------------------------------------------------------

datasets <- list(
  list(
    key = "spn_penicillin", species_dir = "spn_penicillin",
    title = "Penicillin~resistance~(italic('S. pneumoniae'))",
    binning_specs = list(
      list(binning = "standard",     nn = "02", run_stub = "spn_penicillin_MIC"),
      list(binning = "coarse",       nn = "10", run_stub = "spn_penicillin_MIC_coarse_dilutions"),
      list(binning = "large minbin", nn = "11", run_stub = "spn_penicillin_MIC_large_minbin"),
      list(binning = "minima",       nn = "16", run_stub = "spn_penicillin_MIC_minimabinning")),
    baseline_specs = list(
      list(model = "logistic",   run_dir = "01_spn_penicillin_binary_logistic"),
      list(model = "continuous", run_dir = "03_spn_penicillin_continuous_continuous"))
  ),
  list(
    key = "spn_trimethoprim", species_dir = "spn_trimethoprim",
    title = "Trimethoprim~resistance~(italic('S. pneumoniae'))",
    binning_specs = list(
      list(binning = "standard",     nn = "05", run_stub = "spn_trimethoprim_MIC"),
      list(binning = "coarse",       nn = "12", run_stub = "spn_trimethoprim_MIC_coarse_dilutions"),
      list(binning = "large minbin", nn = "13", run_stub = "spn_trimethoprim_MIC_large_minbin")),
    baseline_specs = list(
      list(model = "logistic",   run_dir = "04_spn_trimethoprim_binary_logistic"),
      list(model = "continuous", run_dir = "06_spn_trimethoprim_continuous_continuous"))
  ),
  list(
    key = "tb_rifampicin", species_dir = "tb_rifampicin",
    title = "Rifampicin~resistance~(italic('M. tuberculosis'))",
    binning_specs = list(
      list(binning = "standard",     nn = "08", run_stub = "tb_rifampicin_MIC"),
      list(binning = "coarse",       nn = "14", run_stub = "tb_rifampicin_MIC_coarse_dilutions"),
      list(binning = "large minbin", nn = "15", run_stub = "tb_rifampicin_MIC_large_minbin")),
    baseline_specs = list(
      list(model = "logistic",   run_dir = "07_tb_rifampicin_binary_logistic"),
      list(model = "continuous", run_dir = "09_tb_rifampicin_continuous_continuous"))
  )
)

# Binning order/levels covering every dataset (datasets without `minima` just
# drop that level after factoring).
binning_levels <- c("standard", "coarse", "large minbin", "minima")

# -----------------------------------------------------------------------------
# Shared aesthetics
# -----------------------------------------------------------------------------

base_theme <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92", colour = NA),
        legend.position = "bottom")

POM_COLOUR <- "grey20"   # distinct from the PPOM cutpoint-MIC gradient
YLIM       <- c(0, 1)

# -----------------------------------------------------------------------------
# Per-dataset figure builder
# -----------------------------------------------------------------------------

make_figure <- function(ds, style = "default", suffix = "") {
  herit_csv <- function(run_dir) file.path(
    RESULTS_ROOT, paste0("gwas_", ds$species_dir), "inference", run_dir,
    "plots", "heritability", "heritability_summary.csv")

  missing_files <- character(0)
  read_herit <- function(run_dir) {
    path <- herit_csv(run_dir)
    if (!file.exists(path)) { missing_files <<- c(missing_files, path); return(NULL) }
    tryCatch(read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
             error = function(e) NULL)
  }

  # Build long data frame ------------------------------------------------------
  rows <- list()
  for (s in ds$binning_specs) {
    for (model in c("POM", "PPOM")) {
      d <- read_herit(paste0(s$nn, "_", s$run_stub, "_", model))
      if (is.null(d)) next
      d <- d[, KEEP]
      d$group <- "binning"; d$binning <- s$binning; d$model <- model
      rows[[length(rows) + 1]] <- d
    }
  }
  for (s in ds$baseline_specs) {
    d <- read_herit(s$run_dir)
    if (is.null(d)) next
    d <- d[, KEEP]
    d$group <- "baseline"; d$binning <- NA_character_; d$model <- s$model
    rows[[length(rows) + 1]] <- d
  }
  df <- do.call(rbind, rows)

  df$binning <- factor(df$binning, levels = binning_levels)
  df$model   <- factor(df$model, levels = c("POM", "PPOM", "logistic", "continuous"))
  df$metric  <- factor(ifelse(df$metric == "h2_narrow", "narrow-sense", "broad-sense"),
                       levels = c("narrow-sense", "broad-sense"))

  message(sprintf("[%s] ", ds$key), appendLF = FALSE)
  if (length(missing_files) > 0) {
    message(sprintf("%d heritability CSV(s) missing (rendered as gaps):", length(missing_files)))
    for (f in missing_files) message("  ", f)
  } else {
    message("all heritability CSVs found.")
  }

  # Style-dependent aesthetics -------------------------------------------------
  # "workflow" borrows the gwas_workflow PPOM heritability plot look: Hiroshige
  # (MetBrewer) colours, theme_minimal, italic h^2 facet labels, and capless CI
  # segments + median points (gwas_workflow/R/heritability_plots.R).
  is_wf <- identical(style, "workflow")

  plot_theme <- if (is_wf) {
    theme_minimal(base_size = 13) +
      theme(panel.grid.minor = element_blank(),
            panel.grid.major.x = element_blank(),
            strip.text = element_text(size = 13),
            legend.position = "bottom")
  } else {
    base_theme
  }

  metric_labeller <- if (is_wf) {
    labeller(metric = as_labeller(
      c("narrow-sense" = "italic(h)[narrow]^2",
        "broad-sense"  = "italic(h)[broad]^2"), label_parsed))
  } else {
    "label_value"
  }

  colour_scale <- if (is_wf) {
    scale_colour_gradientn(
      colours = MetBrewer::met.brewer("Hiroshige", 256, type = "continuous"),
      trans = "log10", breaks = c(0.03, 0.1, 0.3, 1),
      name = expression("MIC breakpoint ("*mu*g%.%mL^-1*")"))
  } else {
    scale_colour_viridis_c(trans = "log10", option = "D", end = 0.95,
                           breaks = c(0.03, 0.1, 0.3, 1),
                           name = expression("MIC breakpoint ("*mu*g%.%mL^-1*")"))
  }
  # Hiroshige midpoint colour for the PPOM shape-legend key in workflow style.
  ppom_key_col <- if (is_wf) MetBrewer::met.brewer("Hiroshige", 3,
                                                   type = "continuous")[2L]
                  else viridis::viridis(1, end = 0.5)

  # Panel A: binnings, POM vs PPOM, faceted by metric --------------------------
  binA   <- df[df$group == "binning", ]
  pom_A  <- binA[binA$model == "POM", ]
  ppom_A <- binA[binA$model == "PPOM", ]

  # PPOM cutpoints are dodged across the right half of each binning's slot; POM is
  # nudged just left of that cluster so the single POM marker stays attached to its
  # own binning group (close to its cutpoints) without overlapping any of them.
  ppom_dodge <- position_dodge(width = 0.45)
  pom_nudge  <- position_nudge(x = -0.3)
  n_bin      <- nlevels(droplevels(binA$binning))

  # CI layers: capless segments (workflow) vs whisker error bars (default).
  ppom_ci <- if (is_wf) {
    geom_linerange(data = ppom_A,
                   aes(ymin = q_lower, ymax = q_upper, group = cutpoint_mic),
                   position = ppom_dodge, colour = "grey55", linewidth = 0.7)
  } else {
    geom_errorbar(data = ppom_A,
                  aes(ymin = q_lower, ymax = q_upper, group = cutpoint_mic),
                  position = ppom_dodge, width = 0.45, colour = "grey60", linewidth = 0.4)
  }
  pom_ci <- if (is_wf) {
    geom_linerange(data = pom_A, aes(ymin = q_lower, ymax = q_upper),
                   position = pom_nudge, colour = POM_COLOUR, linewidth = 0.9)
  } else {
    geom_errorbar(data = pom_A, aes(ymin = q_lower, ymax = q_upper),
                  position = pom_nudge, width = 0.2, colour = POM_COLOUR, linewidth = 0.5)
  }

  panel_a <- ggplot(mapping = aes(x = binning)) +
    # Light separators so each binning group is unambiguous.
    geom_vline(xintercept = seq_len(max(n_bin - 1, 0)) + 0.5,
               colour = "grey88", linewidth = 0.4) +
    ppom_ci +
    geom_point(data = ppom_A,
               aes(y = median, colour = cutpoint_mic, group = cutpoint_mic, shape = "PPOM"),
               position = ppom_dodge, size = 2.6) +
    pom_ci +
    geom_point(data = pom_A,
               aes(y = median, shape = "POM"),
               position = pom_nudge, colour = POM_COLOUR, size = 3) +
    facet_grid(metric ~ ., labeller = metric_labeller) +
    colour_scale +
    scale_shape_manual(values = c(POM = 17, PPOM = 16), name = "model",
                       limits = c("POM", "PPOM")) +
    scale_x_discrete(drop = TRUE, expand = expansion(add = 0.55)) +
    coord_cartesian(ylim = YLIM) +
    labs(x = NULL, y = expression("heritability ("*italic(h)^2*")"),
         title = "Ordered logistic models with differing MIC category intervals") +
    plot_theme +
    theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
    guides(shape = guide_legend(order = 1,
                                override.aes = list(colour = c(POM_COLOUR, ppom_key_col))),
           colour = guide_colourbar(order = 2, barwidth = grid::unit(4, "cm"),
                                    title.position = "top", title.hjust = 0.5))

  # Panel B: non-ordered baselines, faceted by metric (same y scale) -----------
  baseB <- df[df$group == "baseline", ]
  baseB$model <- droplevels(baseB$model)

  base_ci <- if (is_wf) {
    geom_linerange(aes(ymin = q_lower, ymax = q_upper),
                   colour = "grey40", linewidth = 0.9)
  } else {
    geom_errorbar(aes(ymin = q_lower, ymax = q_upper),
                  width = 0.2, colour = "grey40", linewidth = 0.5)
  }

  panel_b <- ggplot(baseB, aes(x = model, y = median)) +
    base_ci +
    geom_point(size = 3, colour = "grey20") +
    facet_grid(metric ~ ., labeller = metric_labeller) +
    coord_cartesian(ylim = YLIM) +
    labs(x = NULL, y = NULL, title = "Non-ordered models") +
    plot_theme +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

  # Assemble + save ------------------------------------------------------------
  legend <- get_legend(panel_a)
  panels <- plot_grid(
    panel_a + theme(legend.position = "none"), panel_b,
    nrow = 1, rel_widths = c(3, 1.1), align = "h", axis = "tb",
    labels = c("A", "B"), label_size = 20, label_fontface = "bold")
  figure <- plot_grid(panels, legend, ncol = 1, rel_heights = c(1, 0.13)) +
    theme(plot.margin = margin(5, 5, 5, 10))

  dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
  png_path <- file.path(OUTPUT_DIR, paste0("heritability_summary_", ds$key, suffix, ".png"))
  ggsave(png_path, figure, width = 12, height = 8.3, dpi = 300, bg = "white")
  message("wrote ", png_path)

  # The underlying data is style-independent, so write the CSV once (default).
  if (!is_wf) {
    csv_path <- file.path(OUTPUT_DIR, paste0("heritability_summary_", ds$key, ".csv"))
    write.csv(df[order(df$group, df$binning, df$model, df$metric, df$cutpoint), ],
              csv_path, row.names = FALSE)
    message("wrote ", csv_path)
  }
}

# Default (viridis) and workflow-style (Hiroshige / theme_minimal) versions.
invisible(lapply(datasets, make_figure, style = "default", suffix = ""))
invisible(lapply(datasets, make_figure, style = "workflow", suffix = "_workflowstyle"))
