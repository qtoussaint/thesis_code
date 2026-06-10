# manhattan_builders.R
# Native-ggplot builders for the short paper Manhattan panels (RATE + beta),
# shared by the standalone replot_*_short.R scripts and the locus-zoom composites.
#
# Each builder loads from an inference run dir and returns a ggplot object (so the
# composites can introspect panel geometry for the magnifier connector lines and
# the standalone scripts can ggsave it). Gene-name labels are optional and use the
# helpers in gene_label_helpers.R.

suppressPackageStartupMessages({
  library(ggplot2)
  library(viridis)
  library(data.table)
})

.MB_DIR <- "/nfs/research/jlees/jacqueline/thesis_code/paper_figures"
source(file.path(.MB_DIR, "gene_label_helpers.R"))
# mb_panel_cell_frac() (panel data-area geometry) is reused for the penicillin
# right-shift label placement; safe to source (no back-dependency on this file).
source(file.path(.MB_DIR, "composite_connectors.R"))

MB_BASE_SIZE    <- 14
MB_SHORT_HEIGHT <- 8 * 2 / 3
MB_FIG_WIDTH    <- 16

# Manhattan legend text sized to match the shared locus-zoom legend
MB_LEGEND_TITLE <- 12
MB_LEGEND_TEXT  <- 10
mb_legend_theme <- function() ggplot2::theme(
  legend.title = ggplot2::element_text(size = MB_LEGEND_TITLE),
  legend.text  = ggplot2::element_text(size = MB_LEGEND_TEXT))

# ---------------------------------------------------------------------------
# Low-level readers (phandango BP column, RATE col 2, positions CSV col 2)
# ---------------------------------------------------------------------------
mb_read_positions <- function(phandango) {
  if (!file.exists(phandango)) stop("Missing: ", phandango)
  as.numeric(read.table(phandango, header = TRUE, sep = "\t")$BP)
}
mb_read_rate <- function(rate_file) {
  if (!file.exists(rate_file)) stop("Missing: ", rate_file)
  as.numeric(read.table(rate_file, comment.char = "#")[, 2])
}
mb_read_positions_file <- function(positions_file) {
  if (!file.exists(positions_file)) stop("Missing: ", positions_file)
  suppressWarnings(as.numeric(read.csv(positions_file, header = TRUE)[[2]]))
}

# ---------------------------------------------------------------------------
# RATE data frame: PPOM overlay (per-cutpoint) or binary single series.
# Returns list(df, n, overlay, legend_title, outfile).
# ---------------------------------------------------------------------------
mb_build_rate_df <- function(run_dir, positions_file = NULL) {
  ppom_rate1 <- file.path(run_dir, "cppRATE_results",
                          "RATE_values_cutpoint1_depruned.txt")
  if (file.exists(ppom_rate1)) {
    cutpoint_files <- sort(Sys.glob(file.path(run_dir, "cppRATE_results",
                                              "RATE_values_cutpoint*_depruned.txt")))
    n_cutpoints <- length(cutpoint_files)

    effects <- read.csv(file.path(run_dir, "fitted_model",
                                  "depruned_variant_effects.csv"))
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
      pos  <- mb_read_positions(file.path(run_dir, "cppRATE_results",
                                          sprintf("phandango_cutpoint%d.plot", c)))
      rate <- mb_read_rate(cutpoint_files[c])
      data.frame(pos = pos, rate = rate, cutpoint_label = labels[c])
    })
    df <- do.call(rbind, parts)
    df$cutpoint_label <- factor(df$cutpoint_label, levels = labels)
    list(df = df, n = n_cutpoints, overlay = TRUE, legend_title = legend_title,
         outfile = "manhattan_all_cutpoints_overlayed_RATE.png")
  } else {
    rate_file <- file.path(run_dir, "cppRATE_results", "RATE_values_depruned.txt")
    if (!file.exists(rate_file)) {
      rate_file <- file.path(run_dir, "cppRATE_results", "RATE_values.txt")
    }
    rate <- mb_read_rate(rate_file)
    pos <- if (!is.null(positions_file)) {
      mb_read_positions_file(positions_file)
    } else {
      mb_read_positions(file.path(run_dir, "cppRATE_results", "phandango.plot"))
    }
    if (length(pos) != length(rate)) {
      stop("Position count (", length(pos), ") != RATE count (", length(rate), ")")
    }
    df <- data.frame(pos = pos, rate = rate)
    df <- df[!is.na(df$pos), ]
    list(df = df, n = 1L, overlay = FALSE, legend_title = NULL,
         outfile = "manhattan_RATE.png")
  }
}

# Attach gene labels to a manhattan ggplot. `df` must have `pos` and a `value`
# column (the plotted y). `score_col` ranks candidates. Returns the ggplot with a
# repel layer added (or unchanged when no labels are requested / found).
#
# Penicillin-only options:
#   window_bp   when set, re-anchor each gene label at the top SNP within +/-
#               window_bp of the gene's top SNP (label the regional peak, not the
#               gene's own SNP).
#   no_overlap  when TRUE, place each label just above its peak and shift it to the
#               RIGHT until its text box covers no plotted point or other label.
mb_add_labels <- function(p, df, annotations, goi, mode, genes, n_labels,
                          score_col = "value", size = 5, gene_aliases = NULL,
                          window_bp = NULL, no_overlap = FALSE,
                          cell_w_in = MB_FIG_WIDTH, cell_h_in = MB_SHORT_HEIGHT) {
  if (is.null(annotations) || is.null(mode)) return(p)
  ann <- annotate_genes(data.frame(pos = df$pos), annotations, goi)
  df$gene <- ann$gene
  # optional display renames, e.g. c(dyr = "folA") to relabel a gene on one plot
  if (!is.null(gene_aliases)) {
    hit <- df$gene %in% names(gene_aliases)
    df$gene[hit] <- unname(gene_aliases[df$gene[hit]])
  }
  label_df <- select_label_rows(df, gene_col = "gene", score_col = score_col,
                                mode = mode, genes = genes, n = n_labels)
  if (nrow(label_df) == 0L) return(p)

  # re-anchor each label at the regional peak SNP within +/- window_bp
  if (!is.null(window_bp)) {
    for (i in seq_len(nrow(label_df))) {
      w <- which(df$pos >= label_df$pos[i] - window_bp &
                 df$pos <= label_df$pos[i] + window_bp & !is.na(df$value))
      if (length(w) > 0L) {
        j <- w[which.max(df$value[w])]
        label_df$pos[i]   <- df$pos[j]
        label_df$value[i] <- df$value[j]
      }
    }
  }

  # headroom above the tallest peak so labels do not clip the panel edge
  p <- p + ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(mult = c(0.05, 0.16)))

  if (no_overlap) {
    return(p + mb_rightshift_label_layer(p, df, label_df, size = size,
                                         cell_w_in = cell_w_in,
                                         cell_h_in = cell_h_in))
  }
  # default: a small upward nudge so the label sits on top of the peak SNP
  nudge <- 0.04 * diff(range(df$value, na.rm = TRUE))
  p + repel_label_layer(label_df, size = size, nudge_y = nudge)
}

# Place gene labels just above their peak SNP, shifting each label to the RIGHT
# until its text box no longer covers any plotted point (or an already-placed
# label). Used for the penicillin figure. `p` is the scaled-but-unlabeled plot
# (so its built x/y ranges include the headroom); cell_w_in/cell_h_in are the
# rendered size of the manhattan's grid cell. Returns a geom_text layer.
mb_rightshift_label_layer <- function(p, df, label_df, size = 5,
                                      colour = "grey25",
                                      cell_w_in = MB_FIG_WIDTH,
                                      cell_h_in = MB_SHORT_HEIGHT) {
  b  <- ggplot2::ggplot_build(p)
  pr <- b$layout$panel_params[[1]]
  xr <- pr$x.range; yr <- pr$y.range
  fr <- mb_panel_cell_frac(p, cell_w_in, cell_h_in)
  panel_w_in <- (fr$x_right - fr$x_left) * cell_w_in
  panel_h_in <- (fr$y_top  - fr$y_bottom) * cell_h_in
  dpix <- diff(xr) / panel_w_in            # data-x units per rendered inch
  dpiy <- diff(yr) / panel_h_in            # data-y units per rendered inch

  text_h_in <- size / 25.4                 # geom_text size is in mm
  box_h     <- text_h_in * dpiy
  gap_y     <- 0.30 * box_h                # small gap above the peak

  # claim space for the tallest peaks first; never push a label past the panel
  ord    <- order(label_df$value, decreasing = TRUE)
  placed <- list()
  lx <- numeric(nrow(label_df)); ly <- numeric(nrow(label_df))
  for (k in ord) {
    half_w  <- 0.5 * nchar(label_df$gene[k]) * 0.62 * text_h_in * dpix
    step    <- 0.5 * half_w
    box_bot <- label_df$value[k] + gap_y
    box_top <- box_bot + box_h
    x_max   <- xr[2] - half_w
    cx      <- label_df$pos[k]
    for (it in 0:60) {
      L <- cx - half_w; R <- cx + half_w
      hit_pt <- any(df$pos >= L & df$pos <= R &
                    df$value >= box_bot & df$value <= box_top, na.rm = TRUE)
      hit_lab <- FALSE
      for (q in placed) {
        if (!(R < q$L || L > q$R || box_top < q$bot || box_bot > q$top)) {
          hit_lab <- TRUE; break
        }
      }
      if ((!hit_pt && !hit_lab) || cx >= x_max) break
      cx <- min(cx + step, x_max)
    }
    lx[k] <- cx; ly[k] <- box_bot
    placed[[length(placed) + 1L]] <- list(L = cx - half_w, R = cx + half_w,
                                          bot = box_bot, top = box_top)
  }

  placed_df <- label_df
  placed_df$lx <- lx; placed_df$ly <- ly
  ggplot2::geom_text(
    data    = placed_df,
    mapping = ggplot2::aes(x = .data[["lx"]], y = .data[["ly"]],
                           label = .data[["gene_expr"]]),
    parse   = TRUE, size = size, colour = colour, inherit.aes = FALSE,
    vjust   = 0, hjust = 0.5)
}

# Build the short RATE manhattan ggplot for one run, optionally labeled.
#   label_mode  NULL | "top_n" | "gene_list"
#   label_genes display-name vector for gene_list mode
make_rate_manhattan <- function(run_dir, positions_file = NULL,
                                annotations = NULL, goi = NULL,
                                label_mode = NULL, label_genes = NULL,
                                n_labels = 10L, gene_aliases = NULL,
                                label_window_bp = NULL, label_no_overlap = FALSE) {
  built <- mb_build_rate_df(run_dir, positions_file)
  df <- built$df
  plasma_colors <- viridis::viridis(built$n + 1L, option = "plasma")[2:(built$n + 1L)]

  if (built$overlay) {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = pos, y = rate, colour = cutpoint_label)) +
      ggplot2::geom_point(alpha = 0.4) +
      ggplot2::scale_colour_manual(values = plasma_colors) +
      ggplot2::labs(colour = built$legend_title)
  } else {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = pos, y = rate)) +
      ggplot2::geom_point(alpha = 0.4, colour = "#6A00A8")
  }
  p <- p +
    ggplot2::xlab("genome coordinate (bp)") +
    ggplot2::ylab("relative centrality (RATE)") +
    ggplot2::theme_minimal(base_size = MB_BASE_SIZE) +
    mb_legend_theme()

  df$value <- df$rate
  p <- mb_add_labels(p, df, annotations, goi, label_mode, label_genes, n_labels,
                     score_col = "value", gene_aliases = gene_aliases,
                     window_bp = label_window_bp, no_overlap = label_no_overlap)
  attr(p, "outfile") <- built$outfile
  p
}

# ---------------------------------------------------------------------------
# Beta data: PPOM overlay median effects (one row per variant per cutpoint).
# Mirrors replot_ppom_overlay_manhattans_short.R load_inputs/build_overlay_df.
# ---------------------------------------------------------------------------
mb_build_overlay_df <- function(run_dir, phandango) {
  variant_positions <- read.csv(phandango)[, 2]
  effects_csv <- file.path(run_dir, "fitted_model", "depruned_variant_effects.csv")
  if (!file.exists(effects_csv)) stop("Missing: ", effects_csv)
  eff <- data.table::fread(effects_csv)
  n_cutpoints <- length(unique(eff$cutpoint))

  has_mic <- !all(is.na(eff$cutpoint_MIC))
  df <- data.frame(
    pos      = rep(as.numeric(variant_positions), times = n_cutpoints),
    median   = as.numeric(eff$median),
    cutpoint = eff$cutpoint
  )
  if (has_mic) {
    df$cutpoint_label <- as.character(eff$cutpoint_MIC)
    legend_title <- expression(atop("breakpoint", (mu * "g" %.% "mL"^{-1})))
  } else {
    df$cutpoint_label <- as.character(eff$cutpoint)
    legend_title <- "cutpoint"
  }
  df$cutpoint_label <- factor(df$cutpoint_label,
    levels = unique(df$cutpoint_label[order(df$cutpoint)]))
  list(df = df, n = n_cutpoints, legend_title = legend_title)
}

# Build the two beta manhattan ggplots (median and exp(|median|)), optionally
# labeled by gene_list. Ranking is on |median| (the strongest-effect variant per
# gene). Returns list(median = ggplot, exp_abs = ggplot).
make_beta_manhattans <- function(run_dir, phandango,
                                 annotations = NULL, goi = NULL,
                                 label_genes = NULL, gene_aliases = NULL,
                                 label_window_bp = NULL, label_no_overlap = FALSE) {
  ov <- mb_build_overlay_df(run_dir, phandango)
  df <- ov$df
  plasma_colors <- viridis::viridis(ov$n + 1L, option = "plasma")[2:(ov$n + 1L)]
  base_theme <- ggplot2::theme_minimal(base_size = MB_BASE_SIZE) + mb_legend_theme()
  df$absmedian <- abs(df$median)
  mode <- if (!is.null(label_genes)) "gene_list" else NULL

  p_median <- ggplot2::ggplot(df,
      ggplot2::aes(x = pos, y = median, colour = cutpoint_label)) +
    ggplot2::geom_point(alpha = 0.4) +
    ggplot2::scale_colour_manual(values = plasma_colors) +
    ggplot2::xlab("genome coordinate (bp)") +
    ggplot2::ylab(expression(tilde(beta))) +
    ggplot2::labs(colour = ov$legend_title) +
    base_theme
  dfm <- df; dfm$value <- df$median
  p_median <- mb_add_labels(p_median, dfm, annotations, goi, mode, label_genes,
                            n_labels = 10L, score_col = "value",
                            gene_aliases = gene_aliases,
                            window_bp = label_window_bp,
                            no_overlap = label_no_overlap)

  p_exp_abs <- ggplot2::ggplot(df,
      ggplot2::aes(x = pos, y = exp(abs(median)), colour = cutpoint_label)) +
    ggplot2::geom_point(alpha = 0.4) +
    ggplot2::scale_colour_manual(values = plasma_colors) +
    ggplot2::xlab("genome coordinate (bp)") +
    ggplot2::ylab(expression("e"^{abs(tilde(beta))})) +
    ggplot2::labs(colour = ov$legend_title) +
    base_theme
  dfe <- df; dfe$value <- exp(abs(df$median))
  p_exp_abs <- mb_add_labels(p_exp_abs, dfe, annotations, goi, mode, label_genes,
                             n_labels = 10L, score_col = "absmedian",
                             gene_aliases = gene_aliases,
                             window_bp = label_window_bp,
                             no_overlap = label_no_overlap)

  list(median = p_median, exp_abs = p_exp_abs)
}
