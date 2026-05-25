#!/usr/bin/env Rscript
# Recreate the three PPOM overlay Manhattan plots used in the paper composite
# (median effects, exp(|median|), RATE) at 2/3 height, reusing the same data
# and styling as gwas_workflow::write_ppom_combined_manhattan_plots().
#
# Usage:
#   Rscript replot_ppom_overlay_manhattans_short.R \
#     --run-dir   <inference output dir> \
#     --phandango <variant_index.csv> \
#     [--output-dir <dir>]
#
# Default --output-dir is
#   /nfs/research/jlees/jacqueline/thesis_results/paper_figures/<basename(run-dir)>/manhattan_short

suppressPackageStartupMessages({
  library(ggplot2)
  library(viridis)
  library(data.table)
})

DEFAULT_PAPER_FIGURES_DIR <- "/nfs/research/jlees/jacqueline/thesis_results/paper_figures"

parse_args <- function(argv) {
  out <- list(run_dir = NULL, phandango = NULL, output_dir = NULL)
  i <- 1
  while (i <= length(argv)) {
    a <- argv[i]
    if (a == "--run-dir") {
      out$run_dir <- argv[i + 1]; i <- i + 2
    } else if (a == "--phandango") {
      out$phandango <- argv[i + 1]; i <- i + 2
    } else if (a == "--output-dir") {
      out$output_dir <- argv[i + 1]; i <- i + 2
    } else {
      stop("Unknown argument: ", a)
    }
  }
  if (is.null(out$run_dir) || is.null(out$phandango)) {
    stop("Required args: --run-dir <dir> --phandango <csv>")
  }
  if (is.null(out$output_dir)) {
    out$output_dir <- file.path(DEFAULT_PAPER_FIGURES_DIR,
                                basename(normalizePath(out$run_dir, mustWork = TRUE)),
                                "manhattan_short")
  }
  out
}

load_inputs <- function(run_dir, phandango) {
  variant_positions <- read.csv(phandango)[, 2]

  effects_csv <- file.path(run_dir, "fitted_model", "depruned_variant_effects.csv")
  if (!file.exists(effects_csv)) stop("Missing: ", effects_csv)
  depruned_all_cutpoints <- data.table::fread(effects_csv)

  n_cutpoints <- length(unique(depruned_all_cutpoints$cutpoint))
  all_rates <- vector("list", n_cutpoints)
  for (c in seq_len(n_cutpoints)) {
    rate_file <- file.path(run_dir, "cppRATE_results",
                           sprintf("RATE_values_cutpoint%d_depruned.txt", c))
    if (!file.exists(rate_file)) stop("Missing: ", rate_file)
    all_rates[[c]] <- as.numeric(read.table(rate_file)[, 2])
  }

  list(
    variant_positions      = variant_positions,
    depruned_all_cutpoints = depruned_all_cutpoints,
    all_rates              = all_rates,
    n_cutpoints            = n_cutpoints
  )
}

# Mirrors the data-prep block in
# gwas_workflow/R/manhattan_plots.R:369-399
build_overlay_df <- function(variant_positions, depruned_all_cutpoints, n_cutpoints) {
  has_mic <- !all(is.na(depruned_all_cutpoints$cutpoint_MIC))

  df <- data.frame(
    pos      = rep(as.numeric(variant_positions), times = n_cutpoints),
    median   = as.numeric(depruned_all_cutpoints$median),
    cutpoint = depruned_all_cutpoints$cutpoint
  )

  if (has_mic) {
    df$cutpoint_label <- as.character(depruned_all_cutpoints$cutpoint_MIC)
    legend_title <- expression(atop("breakpoint", (mu * "g" %.% "mL"^{-1})))
  } else {
    df$cutpoint_label <- as.character(depruned_all_cutpoints$cutpoint)
    legend_title <- "cutpoint"
  }

  df$cutpoint_label <- factor(df$cutpoint_label,
    levels = unique(df$cutpoint_label[order(df$cutpoint)]))

  list(df = df, legend_title = legend_title, has_mic = has_mic)
}

build_rate_df <- function(variant_positions, all_rates, n_cutpoints, overlay_df) {
  df_rate <- data.frame(
    pos      = rep(as.numeric(variant_positions), times = n_cutpoints),
    rate     = unlist(all_rates),
    cutpoint = rep(seq_len(n_cutpoints), each = length(variant_positions))
  )
  if (overlay_df$has_mic) {
    breakpoints <- levels(overlay_df$df$cutpoint_label)
    df_rate$cutpoint_label <- factor(
      as.character(breakpoints[df_rate$cutpoint]),
      levels = breakpoints
    )
  } else {
    df_rate$cutpoint_label <- factor(df_rate$cutpoint)
  }
  df_rate
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  run_dir <- normalizePath(args$run_dir, mustWork = TRUE)
  phandango <- normalizePath(args$phandango, mustWork = TRUE)
  dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)

  message("Run dir:    ", run_dir)
  message("Phandango:  ", phandango)
  message("Output dir: ", args$output_dir)

  inputs <- load_inputs(run_dir, phandango)
  overlay <- build_overlay_df(inputs$variant_positions,
                              inputs$depruned_all_cutpoints,
                              inputs$n_cutpoints)
  df_rate <- build_rate_df(inputs$variant_positions, inputs$all_rates,
                           inputs$n_cutpoints, overlay)

  plasma_colors <- viridis::viridis(inputs$n_cutpoints + 1L, option = "plasma")[
    2:(inputs$n_cutpoints + 1L)]
  base_theme <- ggplot2::theme_minimal(base_size = 14)

  short_height <- 8 * 2 / 3

  save_short <- function(plot, name) {
    path <- file.path(args$output_dir, name)
    ggplot2::ggsave(path, plot = plot,
                    width = 16, height = short_height, dpi = 300)
    message("Wrote ", path)
  }

  # Panel A: median effects
  p_median <- ggplot2::ggplot(overlay$df,
      ggplot2::aes(x = pos, y = median, colour = cutpoint_label)) +
    ggplot2::geom_point(alpha = 0.4) +
    ggplot2::scale_colour_manual(values = plasma_colors) +
    ggplot2::xlab("genome coordinate (bp)") +
    ggplot2::ylab(expression(tilde(beta))) +
    ggplot2::labs(colour = overlay$legend_title) +
    base_theme
  save_short(p_median, "manhattan_all_cutpoints_overlayed_median_effects.png")

  # Panel C: exp(|median|)
  p_exp_abs <- ggplot2::ggplot(overlay$df,
      ggplot2::aes(x = pos, y = exp(abs(median)), colour = cutpoint_label)) +
    ggplot2::geom_point(alpha = 0.4) +
    ggplot2::scale_colour_manual(values = plasma_colors) +
    ggplot2::xlab("genome coordinate (bp)") +
    ggplot2::ylab(expression("e"^{abs(tilde(beta))})) +
    ggplot2::labs(colour = overlay$legend_title) +
    base_theme
  save_short(p_exp_abs, "manhattan_all_cutpoints_overlayed_exp_abs_median.png")

  # Panel D: RATE
  p_rate <- ggplot2::ggplot(df_rate,
      ggplot2::aes(x = pos, y = rate, colour = cutpoint_label)) +
    ggplot2::geom_point(alpha = 0.4) +
    ggplot2::scale_colour_manual(values = plasma_colors) +
    ggplot2::xlab("genome coordinate (bp)") +
    ggplot2::ylab("relative centrality (RATE)") +
    ggplot2::labs(colour = overlay$legend_title) +
    base_theme
  save_short(p_rate, "manhattan_all_cutpoints_overlayed_RATE.png")
}

main()
