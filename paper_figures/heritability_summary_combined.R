#!/usr/bin/env Rscript
# Combined heritability (h2) summary figure: the three per-dataset heritability
# figures (SPN penicillin, SPN trimethoprim, TB rifampicin) stacked as three
# labelled panels (A, B, C) in a single figure.
#
# This is the multi-panel companion to heritability_summary.R, which renders the
# same content as three separate files. A single shared legend (POM/PPOM shape +
# one MIC colour bar) sits at the bottom of the figure, so the MIC colour scale is
# common to all three panels (limits span every dataset's cutpoints). Within each
# dataset panel the binning setups and POM vs PPOM models sit on the left, the
# non-ordered baselines on the right, both faceted into narrow-sense and
# broad-sense rows.
#
# Reads the per-run heritability_summary.csv files written by the inference
# pipeline (median + 95% CI from q_lower/q_upper).
#
# Run with:
#   mamba activate gwas_pipeline
#   Rscript paper_figures/heritability_summary_combined.R
#
# Output: paper_figures/heritability_plots/heritability_summary_combined.{png,csv}

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
# expression for the panel header (full drug name + italicised species).
# -----------------------------------------------------------------------------

datasets <- list(
  list(
    key = "spn_penicillin", species_dir = "spn_penicillin",
    title = "Penicillin~resistance~(italic('S. pneumoniae'))",
    binning_specs = list(
      list(binning = "doubling (≥5%)",        nn = "02", run_stub = "spn_penicillin_MIC"),
      list(binning = "4-fold (≥5%)", nn = "10", run_stub = "spn_penicillin_MIC_coarse_dilutions"),
      list(binning = "doubling (≥10%)",       nn = "11", run_stub = "spn_penicillin_MIC_large_minbin"),
      list(binning = "minima",                 nn = "16", run_stub = "spn_penicillin_MIC_minimabinning")),
    baseline_specs = list(
      list(model = "logistic",   run_dir = "01_spn_penicillin_binary_logistic"),
      list(model = "continuous", run_dir = "03_spn_penicillin_continuous_continuous"))
  ),
  list(
    key = "spn_trimethoprim", species_dir = "spn_trimethoprim",
    title = "Trimethoprim~resistance~(italic('S. pneumoniae'))",
    binning_specs = list(
      list(binning = "doubling (≥5%)",        nn = "05", run_stub = "spn_trimethoprim_MIC"),
      list(binning = "4-fold (≥5%)", nn = "12", run_stub = "spn_trimethoprim_MIC_coarse_dilutions"),
      list(binning = "doubling (≥10%)",       nn = "13", run_stub = "spn_trimethoprim_MIC_large_minbin")),
    baseline_specs = list(
      list(model = "logistic",   run_dir = "04_spn_trimethoprim_binary_logistic"),
      list(model = "continuous", run_dir = "06_spn_trimethoprim_continuous_continuous"))
  ),
  list(
    key = "tb_rifampicin", species_dir = "tb_rifampicin",
    title = "Rifampicin~resistance~(italic('M. tuberculosis'))",
    binning_specs = list(
      list(binning = "doubling (≥5%)",        nn = "08", run_stub = "tb_rifampicin_MIC"),
      list(binning = "4-fold (≥5%)", nn = "14", run_stub = "tb_rifampicin_MIC_coarse_dilutions"),
      list(binning = "doubling (≥10%)",       nn = "15", run_stub = "tb_rifampicin_MIC_large_minbin")),
    baseline_specs = list(
      list(model = "logistic",   run_dir = "07_tb_rifampicin_binary_logistic"),
      list(model = "continuous", run_dir = "09_tb_rifampicin_continuous_continuous"))
  )
)

# Binning order/levels covering every dataset (datasets without `minima` just
# drop that level after factoring).
binning_levels <- c("doubling (≥5%)", "4-fold (≥5%)", "doubling (≥10%)", "minima")

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
# Per-dataset long data frame
# -----------------------------------------------------------------------------

load_dataset <- function(ds) {
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
  df$dataset <- ds$key

  message(sprintf("[%s] ", ds$key), appendLF = FALSE)
  if (length(missing_files) > 0) {
    message(sprintf("%d heritability CSV(s) missing (rendered as gaps):", length(missing_files)))
    for (f in missing_files) message("  ", f)
  } else {
    message("all heritability CSVs found.")
  }
  df
}

# -----------------------------------------------------------------------------
# Per-dataset row: title + [panel A | panel B] + its own colour legend
# -----------------------------------------------------------------------------

build_row <- function(ds, df, mic_limits) {
  # Panel A: binnings, POM vs PPOM, faceted by metric ------------------------
  binA   <- df[df$group == "binning", ]
  pom_A  <- binA[binA$model == "POM", ]
  ppom_A <- binA[binA$model == "PPOM", ]

  # PPOM cutpoints are dodged across the right half of each binning's slot; POM is
  # nudged just left of that cluster so the single POM marker stays attached to its
  # own binning group (close to its cutpoints) without overlapping any of them.
  ppom_dodge <- position_dodge(width = 0.45)
  pom_nudge  <- position_nudge(x = -0.3)
  n_bin      <- nlevels(droplevels(binA$binning))

  panel_a <- ggplot(mapping = aes(x = binning)) +
    geom_vline(xintercept = seq_len(max(n_bin - 1, 0)) + 0.5,
               colour = "grey88", linewidth = 0.4) +
    geom_errorbar(data = ppom_A,
                  aes(ymin = q_lower, ymax = q_upper, group = cutpoint_mic),
                  position = ppom_dodge, width = 0.45, colour = "grey60", linewidth = 0.4) +
    geom_point(data = ppom_A,
               aes(y = median, colour = cutpoint_mic, group = cutpoint_mic, shape = "PPOM"),
               position = ppom_dodge, size = 2.6) +
    geom_errorbar(data = pom_A,
                  aes(ymin = q_lower, ymax = q_upper),
                  position = pom_nudge, width = 0.2, colour = POM_COLOUR, linewidth = 0.5) +
    geom_point(data = pom_A,
               aes(y = median, shape = "POM"),
               position = pom_nudge, colour = POM_COLOUR, size = 3) +
    facet_grid(metric ~ .) +
    scale_colour_viridis_c(trans = "log10", option = "D", end = 0.95,
                           limits = mic_limits,
                           breaks = c(0.03, 0.1, 0.3, 1),
                           name = expression("MIC breakpoint ("*mu*g%.%mL^-1*")")) +
    scale_shape_manual(values = c(POM = 17, PPOM = 16), name = "model",
                       limits = c("POM", "PPOM")) +
    scale_x_discrete(drop = TRUE, expand = expansion(add = 0.55)) +
    coord_cartesian(ylim = YLIM) +
    labs(x = NULL, y = expression("heritability ("*italic(h)^2*")"),
         title = "Ordered logistic models with differing MIC category intervals") +
    base_theme +
    theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
    guides(shape = guide_legend(order = 1,
                                override.aes = list(colour = c(POM_COLOUR, viridis::viridis(1, end = 0.5)))),
           colour = guide_colourbar(order = 2, barwidth = grid::unit(4, "cm"),
                                    title.position = "top", title.hjust = 0.5))

  # Panel B: non-ordered baselines, faceted by metric (same y scale) ---------
  baseB <- df[df$group == "baseline", ]
  baseB$model <- droplevels(baseB$model)

  panel_b <- ggplot(baseB, aes(x = model, y = median)) +
    geom_errorbar(aes(ymin = q_lower, ymax = q_upper),
                  width = 0.2, colour = "grey40", linewidth = 0.5) +
    geom_point(size = 3, colour = "grey20") +
    facet_grid(metric ~ .) +
    coord_cartesian(ylim = YLIM) +
    labs(x = NULL, y = NULL, title = "Non-ordered models") +
    base_theme +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

  # Assemble this dataset's row (no legend here -- one shared legend is added once
  # at the bottom of the whole figure). Keep panel_a around so its legend can be
  # extracted for that shared block.
  title   <- ggdraw() +
    draw_label(parse(text = ds$title), fontface = "bold", size = 14, x = 0.5, hjust = 0.5)
  panels  <- plot_grid(
    panel_a + theme(legend.position = "none"), panel_b,
    nrow = 1, rel_widths = c(3, 1.1), align = "h", axis = "tb")
  row <- plot_grid(title, panels, ncol = 1, rel_heights = c(0.10, 1))
  list(row = row, legend_src = panel_a)
}

# -----------------------------------------------------------------------------
# Build all three panels and stack into one figure
# -----------------------------------------------------------------------------

dfs  <- lapply(datasets, load_dataset)

# Global MIC range across every dataset's PPOM cutpoints, so the shared colour bar
# covers all panels on a common scale.
all_mic    <- unlist(lapply(dfs, function(d) d$cutpoint_mic[!is.na(d$cutpoint_mic)]))
mic_limits <- range(all_mic)

built <- Map(function(ds, df) build_row(ds, df, mic_limits), datasets, dfs)
rows  <- lapply(built, `[[`, "row")

panels  <- plot_grid(plotlist = rows, ncol = 1,
                     labels = c("A", "B", "C"), label_size = 22, label_fontface = "bold")
legend  <- get_legend(built[[1]]$legend_src)
figure  <- plot_grid(panels, legend, ncol = 1, rel_heights = c(1, 0.05)) +
  theme(plot.margin = margin(5, 5, 5, 10))

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
png_path <- file.path(OUTPUT_DIR, "heritability_summary_combined.png")
csv_path <- file.path(OUTPUT_DIR, "heritability_summary_combined.csv")

ggsave(png_path, figure, width = 12, height = 22, dpi = 300, bg = "white", limitsize = FALSE)

combined <- do.call(rbind, dfs)
write.csv(combined[order(combined$dataset, combined$group, combined$binning,
                         combined$model, combined$metric, combined$cutpoint), ],
          csv_path, row.names = FALSE)

message("wrote ", png_path)
message("wrote ", csv_path)
