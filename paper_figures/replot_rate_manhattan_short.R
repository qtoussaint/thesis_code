#!/usr/bin/env Rscript
# Recreate just the RATE overlay Manhattan plot at 2/3 height, matching the
# styling of replot_ppom_overlay_manhattans_short.R. Unlike that script this
# reads RATE straight from the cppRATE outputs (raw RATE_values for the y-axis,
# phandango BP column for genome coordinates), so it works for both PPOM runs
# (per-cutpoint overlay) and binary GWAS runs (single series).
#
# Usage:
#   Rscript replot_rate_manhattan_short.R --run-dir <inference output dir> \
#     [--output-dir <dir>] [--positions-file <variant_index.csv>]
#
# Default --output-dir is
#   /nfs/research/jlees/jacqueline/thesis_results/paper_figures/<basename(run-dir)>/manhattan_short
#
# Genome coordinates come from the phandango BP column by default. Some binary
# GWAS runs write a placeholder BP, so pass --positions-file (CSV, column 2 =
# base-pair position, one row per variant in RATE order, same file the locus
# zoom uses) to override.
#
# Output is manhattan_all_cutpoints_overlayed_RATE.png for a PPOM run, or
# manhattan_RATE.png for a binary run.

suppressPackageStartupMessages({
  library(ggplot2)
  library(viridis)
})

DEFAULT_PAPER_FIGURES_DIR <- "/nfs/research/jlees/jacqueline/thesis_results/paper_figures"

parse_args <- function(argv) {
  out <- list(run_dir = NULL, output_dir = NULL, positions_file = NULL)
  i <- 1
  while (i <= length(argv)) {
    a <- argv[i]
    if (a == "--run-dir") {
      out$run_dir <- argv[i + 1]; i <- i + 2
    } else if (a == "--output-dir") {
      out$output_dir <- argv[i + 1]; i <- i + 2
    } else if (a == "--positions-file") {
      out$positions_file <- argv[i + 1]; i <- i + 2
    } else {
      stop("Unknown argument: ", a)
    }
  }
  if (is.null(out$run_dir)) stop("Required arg: --run-dir <dir>")
  if (is.null(out$output_dir)) {
    out$output_dir <- file.path(DEFAULT_PAPER_FIGURES_DIR,
                                basename(normalizePath(out$run_dir, mustWork = TRUE)),
                                "manhattan_short")
  }
  out
}

# phandango .plot is tab-separated with a BP (genome coordinate) column
read_positions <- function(phandango) {
  if (!file.exists(phandango)) stop("Missing: ", phandango)
  as.numeric(read.table(phandango, header = TRUE, sep = "\t")$BP)
}

# RATE_values_*.txt has #-prefixed header lines then "snp_id RATE KLD"; col 2 is RATE
read_rate <- function(rate_file) {
  if (!file.exists(rate_file)) stop("Missing: ", rate_file)
  as.numeric(read.table(rate_file, comment.char = "#")[, 2])
}

# positions CSV: header row, column 2 = base-pair position (one per variant)
read_positions_file <- function(positions_file) {
  if (!file.exists(positions_file)) stop("Missing: ", positions_file)
  suppressWarnings(as.numeric(read.csv(positions_file, header = TRUE)[[2]]))
}

build_df <- function(run_dir, positions_file = NULL) {
  ppom_rate1 <- file.path(run_dir, "cppRATE_results", "RATE_values_cutpoint1_depruned.txt")
  if (file.exists(ppom_rate1)) {
    # PPOM: one phandango + RATE file per cutpoint, overlayed
    cutpoint_files <- sort(Sys.glob(file.path(run_dir, "cppRATE_results",
                                              "RATE_values_cutpoint*_depruned.txt")))
    n_cutpoints <- length(cutpoint_files)

    # cutpoint -> breakpoint MIC label from the fitted effects table
    effects <- read.csv(file.path(run_dir, "fitted_model", "depruned_variant_effects.csv"))
    has_mic <- "cutpoint_MIC" %in% names(effects) && !all(is.na(effects$cutpoint_MIC))
    if (has_mic) {
      map <- unique(effects[, c("cutpoint", "cutpoint_MIC")])
      map <- map[order(map$cutpoint), ]
      labels <- as.character(map$cutpoint_MIC)
      legend_title <- expression(atop("breakpoint", (mu * "g" %.% "mL"^{-1})))
    } else {
      labels <- as.character(seq_len(n_cutpoints))
      legend_title <- "cutpoint"
    }

    parts <- lapply(seq_len(n_cutpoints), function(c) {
      pos  <- read_positions(file.path(run_dir, "cppRATE_results",
                                       sprintf("phandango_cutpoint%d.plot", c)))
      rate <- read_rate(cutpoint_files[c])
      data.frame(pos = pos, rate = rate, cutpoint_label = labels[c])
    })
    df <- do.call(rbind, parts)
    df$cutpoint_label <- factor(df$cutpoint_label, levels = labels)
    list(df = df, n = n_cutpoints, legend_title = legend_title, overlay = TRUE,
         outfile = "manhattan_all_cutpoints_overlayed_RATE.png")
  } else {
    # Binary GWAS: single RATE file (prefer depruned)
    rate_file <- file.path(run_dir, "cppRATE_results", "RATE_values_depruned.txt")
    if (!file.exists(rate_file)) {
      rate_file <- file.path(run_dir, "cppRATE_results", "RATE_values.txt")
    }
    rate <- read_rate(rate_file)
    pos <- if (!is.null(positions_file)) {
      read_positions_file(positions_file)
    } else {
      read_positions(file.path(run_dir, "cppRATE_results", "phandango.plot"))
    }
    if (length(pos) != length(rate)) {
      stop("Position count (", length(pos), ") != RATE count (", length(rate), ")")
    }
    df <- data.frame(pos = pos, rate = rate)
    df <- df[!is.na(df$pos), ]  # drop variants without a mapped coordinate
    list(df = df, n = 1L, legend_title = NULL,
         overlay = FALSE, outfile = "manhattan_RATE.png")
  }
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  run_dir <- normalizePath(args$run_dir, mustWork = TRUE)
  dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)

  message("Run dir:    ", run_dir)
  message("Output dir: ", args$output_dir)

  built <- build_df(run_dir, args$positions_file)
  plasma_colors <- viridis::viridis(built$n + 1L, option = "plasma")[2:(built$n + 1L)]
  base_theme <- ggplot2::theme_minimal(base_size = 14)
  short_height <- 8 * 2 / 3

  if (built$overlay) {
    p <- ggplot2::ggplot(built$df,
        ggplot2::aes(x = pos, y = rate, colour = cutpoint_label)) +
      ggplot2::geom_point(alpha = 0.4) +
      ggplot2::scale_colour_manual(values = plasma_colors) +
      ggplot2::labs(colour = built$legend_title)
  } else {
    # Single-series (binary GWAS): plot in purple, matching the purple end of
    # the plasma palette used for the PPOM overlay panels.
    p <- ggplot2::ggplot(built$df, ggplot2::aes(x = pos, y = rate)) +
      ggplot2::geom_point(alpha = 0.4, colour = "#6A00A8")
  }
  p <- p +
    ggplot2::xlab("genome coordinate (bp)") +
    ggplot2::ylab("relative centrality (RATE)") +
    base_theme

  path <- file.path(args$output_dir, built$outfile)
  ggplot2::ggsave(path, plot = p, width = 16, height = short_height, dpi = 300)
  message("Wrote ", path)
}

main()
