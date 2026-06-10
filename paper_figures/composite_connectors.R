# composite_connectors.R
# Helpers for the locus-zoom composites: turn pre-rendered PNGs into cowplot
# panels, and draw "magnifier" connector lines from a genomic window on a native
# Manhattan ggplot down to the top corners of the locus-zoom panel below it.
#
# The Manhattan panels are native ggplots (so bp -> canvas position is exact); the
# locus-zoom panels are rasterized PNGs (so their corner endpoints are layout
# geometry). All canvas coordinates are npc in [0,1] over the whole figure, which
# is the coordinate system cowplot::ggdraw()/draw_line() use.

suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
  library(grid)
  library(magick)
})

# Trimmed PNG -> list(plot = ggdraw image, aspect = height/width)
mb_panel_from_png <- function(path) {
  trimmed <- magick::image_trim(magick::image_read(path))
  info <- magick::image_info(trimmed)
  list(plot = cowplot::ggdraw() + cowplot::draw_image(trimmed),
       aspect = info$height / info$width)
}

# Data-area rectangle of a ggplot's panel, as fractions of its own layout cell.
# cell_w_in / cell_h_in are the cell's rendered size in inches. The panel is the
# only "null"-sized track, so absolute (non-null) widths/heights convert to inches
# and the panel fills the remainder.
mb_panel_cell_frac <- function(p, cell_w_in, cell_h_in) {
  gt <- ggplot2::ggplotGrob(p)
  wn <- grid::convertWidth(gt$widths,  "in", valueOnly = TRUE)
  hn <- grid::convertHeight(gt$heights, "in", valueOnly = TRUE)
  is_panel <- grepl("^panel", gt$layout$name)
  lc <- min(gt$layout$l[is_panel]); rc <- max(gt$layout$r[is_panel])
  tr <- min(gt$layout$t[is_panel]); br <- max(gt$layout$b[is_panel])
  left_in   <- if (lc > 1)           sum(wn[seq_len(lc - 1)])        else 0
  right_in  <- if (rc < length(wn))  sum(wn[(rc + 1):length(wn)])    else 0
  top_in    <- if (tr > 1)           sum(hn[seq_len(tr - 1)])        else 0
  bottom_in <- if (br < length(hn))  sum(hn[(br + 1):length(hn)])    else 0
  list(
    x_left   = left_in / cell_w_in,
    x_right  = 1 - right_in / cell_w_in,
    y_bottom = bottom_in / cell_h_in,
    y_top    = 1 - top_in / cell_h_in
  )
}

mb_xrange <- function(p) ggplot2::ggplot_build(p)$layout$panel_params[[1]]$x.range

# Vertical bands (top/bottom npc, top = 1) for a stack of rel_heights.
mb_bands <- function(rel_heights) {
  H <- sum(rel_heights)
  bottoms <- 1 - cumsum(rel_heights) / H
  tops    <- c(1, head(bottoms, -1))
  lapply(seq_along(rel_heights), function(i) list(top = tops[i], bottom = bottoms[i]))
}

# Add magnifier connectors for one Manhattan -> locus-zoom-row pair.
#   p_manhattan   native Manhattan ggplot (full figure width)
#   band_m        list(top, bottom) of the Manhattan's vertical band (npc)
#   band_lz       list(top, bottom) of the locus-zoom row's band (npc)
#   regions       list of list(start, end) bp windows, one per gene panel, in row
#                 (left-to-right) order
#   panels_frac   fraction of the row width occupied by the gene panels (the rest
#                 is the shared legend on the right)
#   fig_w_in, total_h_in   rendered figure size in inches
mb_add_connectors <- function(canvas, p_manhattan, band_m, band_lz, regions,
                              panel_x, fig_w_in, total_h_in,
                              line_col = "grey45", line_w = 0.4) {
  cell_h_m <- (band_m$top - band_m$bottom) * total_h_in
  fr <- mb_panel_cell_frac(p_manhattan, fig_w_in, cell_h_m)
  xr <- mb_xrange(p_manhattan)
  bp_to_x <- function(bp) fr$x_left + (bp - xr[1]) / diff(xr) * (fr$x_right - fr$x_left)

  # top endpoints sit at the Manhattan panel's data baseline; bottom endpoints at
  # the top corners of each locus-zoom panel cell (panel_x gives each panel's
  # left/right canvas fraction).
  y_top    <- band_m$bottom + fr$y_bottom * (band_m$top - band_m$bottom)
  y_bottom <- band_lz$top
  for (g in seq_along(regions)) {
    reg <- regions[[g]]
    gx_left  <- panel_x[[g]][1]
    gx_right <- panel_x[[g]][2]
    canvas <- canvas +
      cowplot::draw_line(x = c(bp_to_x(reg$start), gx_left),
                         y = c(y_top, y_bottom), colour = line_col, size = line_w) +
      cowplot::draw_line(x = c(bp_to_x(reg$end), gx_right),
                         y = c(y_top, y_bottom), colour = line_col, size = line_w)
  }
  canvas
}

# Read a <gene>_<metric>_region.txt sidecar -> list(start, end)
mb_read_region <- function(path) {
  kv <- read.table(path, sep = "\t", header = FALSE, stringsAsFactors = FALSE,
                   col.names = c("k", "v"))
  list(start = as.numeric(kv$v[kv$k == "region_start"]),
       end   = as.numeric(kv$v[kv$k == "region_end"]))
}

# Build a locus-zoom row: gene panels on the left, a blank slot on the right where
# the shared legend is later overlaid (mb_assemble draws all legends at one common
# width so their text size matches across rows). Returns the panel grob plus the
# geometry the assembler needs.
mb_lz_row <- function(lz_dir, genes, metric, legend_rel = 0.18) {
  panel_files  <- file.path(lz_dir, sprintf("%s_%s_nolabels.png", genes, metric))
  region_files <- file.path(lz_dir, sprintf("%s_%s_region.txt", genes, metric))
  legend_file  <- file.path(lz_dir, sprintf("%s_%s_legend.png", genes[1], metric))
  missing <- c(panel_files, region_files, legend_file)[
    !file.exists(c(panel_files, region_files, legend_file))]
  if (length(missing) > 0) stop("Missing locus-zoom inputs: ",
                                paste(missing, collapse = ", "))

  # Order panels left-to-right by genomic position so the connector lines, which
  # land at exact bp on the manhattan above, never cross.
  regions0 <- lapply(region_files, mb_read_region)
  ord <- order(vapply(regions0, function(r) r$start, numeric(1)))
  panel_files <- panel_files[ord]
  region_files <- region_files[ord]

  panels <- lapply(panel_files, mb_panel_from_png)
  regions <- lapply(region_files, mb_read_region)
  n <- length(panels)
  panel_aspect <- max(vapply(panels, `[[`, numeric(1), "aspect"))
  panels_frac  <- 1 / (1 + legend_rel)   # row width left of the legend
  blank <- ggplot2::ggplot() + ggplot2::theme_void()

  if (n == 1) {
    # a lone panel is centered at half width so it is not stretched full-bleed
    pw   <- panels_frac / 2
    left <- (panels_frac - pw) / 2
    panel_row <- cowplot::plot_grid(
      blank, panels[[1]]$plot, blank, nrow = 1,
      rel_widths = c(left, pw, panels_frac - left - pw))
    panel_x <- list(c(left, left + pw))
    aspect  <- panel_aspect / 2
  } else {
    panel_row <- cowplot::plot_grid(plotlist = lapply(panels, `[[`, "plot"), nrow = 1)
    panel_x <- lapply(seq_len(n),
                      function(g) c((g - 1) / n * panels_frac, g / n * panels_frac))
    aspect  <- panel_aspect / n
  }

  # blank slot at right keeps the panels in [0, panels_frac]; legend overlaid later
  row <- cowplot::plot_grid(panel_row, blank, nrow = 1, rel_widths = c(1, legend_rel))

  legend_img <- magick::image_trim(magick::image_read(legend_file))
  li <- magick::image_info(legend_img)
  list(grob = row, aspect = aspect, regions = regions, panel_x = panel_x,
       panels_frac = panels_frac, legend_img = legend_img,
       legend_aspect = li$height / li$width)
}

# Assemble a stacked composite from a list of units, each
#   list(manhattan = <ggplot>, lz = <mb_lz_row result or NULL>)
# Manhattans get one band (MB_SHORT_HEIGHT); a unit with an lz row also gets a
# spacer + the lz row, and magnifier connectors from the manhattan to that row.
# Shared legends are overlaid last, all at a single common width so their text
# size is identical across rows. Returns list(canvas, width, height).
mb_assemble <- function(units, fig_w = MB_FIG_WIDTH, gap = 0.55,
                        line_col = "grey45") {
  blank <- ggplot2::ggplot() + ggplot2::theme_void()
  plots <- list(); rels <- numeric(0); labels <- character(0)
  conns <- list(); li <- 1L
  for (u in units) {
    plots[[length(plots) + 1L]] <- u$manhattan
    rels <- c(rels, MB_SHORT_HEIGHT)
    labels <- c(labels, LETTERS[li]); li <- li + 1L
    m_idx <- length(plots)
    if (!is.null(u$lz)) {
      plots[[length(plots) + 1L]] <- blank; rels <- c(rels, gap); labels <- c(labels, "")
      plots[[length(plots) + 1L]] <- u$lz$grob; rels <- c(rels, fig_w * u$lz$aspect)
      labels <- c(labels, LETTERS[li]); li <- li + 1L
      conns[[length(conns) + 1L]] <- list(m = m_idx, lz = length(plots), row = u$lz)
    }
  }

  grid <- cowplot::plot_grid(plotlist = plots, ncol = 1, labels = labels,
                             label_size = 28, label_fontface = "bold",
                             label_x = 0.005, hjust = 0, rel_heights = rels)
  total_h <- sum(rels)
  bands <- mb_bands(rels)
  canvas <- cowplot::ggdraw(grid)

  # connectors
  for (cn in conns) {
    canvas <- mb_add_connectors(canvas, plots[[cn$m]], bands[[cn$m]], bands[[cn$lz]],
                                cn$row$regions, cn$row$panel_x, fig_w, total_h,
                                line_col = line_col)
  }

  # one legend width for all rows => identical legend text size. Width is the
  # largest that still fits every row's band height and the legend gutter.
  if (length(conns) > 0) {
    gutter_in <- (1 - conns[[1]]$row$panels_frac) * fig_w
    fits <- vapply(conns, function(cn) {
      band <- bands[[cn$lz]]
      (band$top - band$bottom) * total_h / cn$row$legend_aspect
    }, numeric(1))
    leg_w_in <- min(c(fits, 0.95 * gutter_in, 2.8))
    for (cn in conns) {
      band <- bands[[cn$lz]]
      pf <- cn$row$panels_frac
      w_npc <- leg_w_in / fig_w
      h_npc <- (leg_w_in * cn$row$legend_aspect) / total_h
      x0 <- pf + ((1 - pf) - w_npc) / 2
      y0 <- (band$top + band$bottom) / 2 - h_npc / 2
      canvas <- canvas +
        cowplot::draw_image(cn$row$legend_img, x = x0, y = y0,
                            width = w_npc, height = h_npc)
    }
  }

  list(canvas = canvas, width = fig_w, height = total_h)
}
