#!/usr/bin/env Rscript
# Combined cutpoint-histogram figure: every ordinal inference run's cutpoint
# posterior on one canvas, grouped by drug/species. Each panel reconstructs the
# pipeline's cutpoints_histogram.png (overlapping low-opacity histograms, one
# Hiroshige colour per cutpoint, a count-scaled density curve, and a coloured
# median line) from the small cutpoint_draws.csv caches written by
# extract_cutpoint_draws.R -- so the legends are live ggplot guides, not baked
# pixels.
#
# Per panel the colour scheme and the per-concentration legend are kept, but the
# legend titles are blanked. The long legend title is factored out into one
# shared mini-legend: a white square with a black outline carrying that title
# text, in the same format as the per-panel titles used to.
#
# Layout: two model columns (POM left, PPOM right), one row per MIC binning,
# blocked by drug/species. Runs whose cutpoint_draws.csv is missing (e.g. still
# inferring) render as an empty captioned cell so the grid stays aligned.
#
# Usage:
#   mamba activate gwas_pipeline
#   Rscript paper_figures/cutpoints_histogram_faceted_summary.R   # after extract_cutpoint_draws.R
#
# Output: <results>/paper_figures/cutpoint_histograms/cutpoints_histogram_faceted_summary.png

suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
  library(MetBrewer)
})

grDevices::pdf(NULL)
on.exit(grDevices::dev.off(), add = TRUE)

RESULTS_ROOT <- "/nfs/research/jlees/jacqueline/thesis_results"
OUTPUT_DIR   <- file.path(RESULTS_ROOT, "paper_figures", "cutpoint_histograms")

MODELS <- c("POM", "PPOM")

# Legend title carried by the pipeline's cutpoint histograms; reproduced verbatim
# (same plotmath format) in the shared mini-legend square.
CUTPOINT_TITLE <- expression("Cutpoint, corresponding MIC breakpoint (" *
                               mu * "g" %.% mL^-1 * ")")

# Discrete Hiroshige colours for per-cutpoint fills (matches gwas_workflow's
# .hiroshige_discrete: n=1 is degenerate for continuous interpolation).
hiroshige_discrete <- function(n) {
  n_safe <- max(as.integer(n), 2L)
  MetBrewer::met.brewer("Hiroshige", n = n_safe, type = "continuous")[seq_len(n)]
}

datasets <- list(
  list(species_dir = "spn_penicillin",   display = expression(italic("S. pneumoniae") * " — penicillin"),
       binning_specs = list(
         list(binning = "doubling (≥5%)",        nn = "02", run_stub = "spn_penicillin_MIC"),
         list(binning = "4-fold (≥5%)", nn = "10", run_stub = "spn_penicillin_MIC_coarse_dilutions"),
         list(binning = "doubling (≥10%)",       nn = "11", run_stub = "spn_penicillin_MIC_large_minbin"),
         list(binning = "minima",                 nn = "16", run_stub = "spn_penicillin_MIC_minimabinning"))),
  list(species_dir = "spn_trimethoprim", display = expression(italic("S. pneumoniae") * " — trimethoprim"),
       binning_specs = list(
         list(binning = "doubling (≥5%)",        nn = "05", run_stub = "spn_trimethoprim_MIC"),
         list(binning = "4-fold (≥5%)", nn = "12", run_stub = "spn_trimethoprim_MIC_coarse_dilutions"),
         list(binning = "doubling (≥10%)",       nn = "13", run_stub = "spn_trimethoprim_MIC_large_minbin"))),
  list(species_dir = "tb_rifampicin",    display = expression(italic("M. tuberculosis") * " — rifampicin"),
       binning_specs = list(
         list(binning = "doubling (≥5%)",        nn = "08", run_stub = "tb_rifampicin_MIC"),
         list(binning = "4-fold (≥5%)", nn = "14", run_stub = "tb_rifampicin_MIC_coarse_dilutions"),
         list(binning = "doubling (≥10%)",       nn = "15", run_stub = "tb_rifampicin_MIC_large_minbin")))
)

cache_path <- function(species_dir, nn, run_stub, model) file.path(
  RESULTS_ROOT, paste0("gwas_", species_dir), "inference",
  paste0(nn, "_", run_stub, "_", model), "plots", "cutpoints", "cutpoint_draws.csv")

fmt_mic <- function(x) ifelse(is.na(x), "", format(x, drop0trailing = TRUE, scientific = FALSE))

# Build one run's histogram panel from its cutpoint_draws.csv. Mirrors
# gwas_workflow's .write_cutpoints_histogram_plot, but with blank legend titles.
# A trailing star on a legend entry flags a cutpoint whose 89% CI overlaps any
# other cutpoint's, exactly as the pipeline figures do.
build_panel <- function(binning, model, csv_path) {
  caption <- sprintf("%s — %s", binning, model)

  if (!file.exists(csv_path)) {
    return(ggdraw() +
      draw_label(caption, x = 0.5, y = 0.93, vjust = 1, size = 11,
                 fontface = "bold", colour = "grey20") +
      draw_label("(no draws yet)", x = 0.5, y = 0.5, size = 10, colour = "grey60"))
  }

  d <- read.csv(csv_path)
  ord     <- sort(unique(d$cutpoint))
  n_cp    <- length(ord)
  n_draws <- nrow(d) / n_cp

  # Per-cutpoint median, 89% CI and MIC, then star for overlapping 89% CIs.
  summ <- do.call(rbind, lapply(ord, function(k) {
    x  <- d$value[d$cutpoint == k]
    qs <- stats::quantile(x, probs = c(0.055, 0.945), names = FALSE)
    data.frame(cutpoint = k, mic = d$cutpoint_mic[d$cutpoint == k][1],
               median = stats::median(x), q_lower = qs[1], q_upper = qs[2])
  }))
  star <- vapply(seq_len(nrow(summ)), function(i) {
    others <- setdiff(seq_len(nrow(summ)), i)
    any(summ$q_lower[i] <= summ$q_upper[others] &
          summ$q_lower[others] <= summ$q_upper[i])
  }, logical(1))
  base_lab <- fmt_mic(summ$mic)
  lab      <- ifelse(star, paste0(base_lab, " ★"), base_lab)
  levels_o <- lab[order(summ$cutpoint)]

  d$cutpoint_label    <- factor(lab[match(d$cutpoint, summ$cutpoint)], levels = levels_o)
  summ$cutpoint_label <- factor(lab, levels = levels_o)

  rng          <- range(d$value)
  binwidth_val <- diff(rng) / 52
  pal          <- hiroshige_discrete(n_cp)

  ggplot(d, aes(x = value)) +
    geom_histogram(aes(fill = cutpoint_label), alpha = 0.3,
                   position = "identity", binwidth = binwidth_val, colour = NA) +
    geom_density(aes(colour = cutpoint_label,
                     y = after_stat(density) * binwidth_val * n_draws),
                 linewidth = 0.5, bw = "nrd0") +
    geom_vline(data = summ, aes(xintercept = median, colour = cutpoint_label),
               linewidth = 0.4, show.legend = FALSE) +
    scale_fill_manual(values = pal, name = NULL) +
    scale_colour_manual(values = pal, name = NULL) +
    labs(x = "Cutpoint value (latent scale)", y = "Count", title = caption) +
    theme_minimal(base_size = 11) +
    theme(plot.title       = element_text(size = 11, face = "bold", colour = "grey20"),
          legend.position  = "right",
          legend.key.size  = grid::unit(0.32, "cm"),
          legend.text      = element_text(size = 8),
          legend.margin    = ggplot2::margin(0, 0, 0, 0))
}

# Shared mini-legend: a single white square with a black outline whose label is
# the (blanked) per-panel legend title, in the same plotmath format.
shared_legend <- function() {
  leg_src <- ggplot(data.frame(x = 1, y = 1, g = factor("k"))) +
    geom_point(aes(x = x, y = y, shape = g),
               fill = "white", colour = "black", size = 5, stroke = 0.9) +
    scale_shape_manual(values = 22, labels = CUTPOINT_TITLE, name = NULL) +
    theme_void() +
    theme(legend.position = "bottom",
          legend.text     = element_text(size = 12),
          legend.key.size = grid::unit(0.6, "cm"))
  cowplot::get_legend(leg_src)
}

# One drug/species block: a bold title strip above its POM/PPOM x binning grid.
make_drug_block <- function(ds) {
  panels <- list()
  for (s in ds$binning_specs) {
    for (model in MODELS) {
      panels[[length(panels) + 1]] <-
        build_panel(s$binning, model, cache_path(ds$species_dir, s$nn, s$run_stub, model))
    }
  }
  grid  <- plot_grid(plotlist = panels, ncol = length(MODELS), align = "hv")
  title <- ggdraw() +
    draw_label(ds$display, x = 0.01, hjust = 0, fontface = "bold", size = 15)
  block <- plot_grid(title, grid, ncol = 1,
                     rel_heights = c(0.35, length(ds$binning_specs)))
  list(block = block, weight = 0.35 + length(ds$binning_specs))
}

blocks  <- lapply(datasets, make_drug_block)
body    <- plot_grid(plotlist = lapply(blocks, `[[`, "block"), ncol = 1,
                     rel_heights = vapply(blocks, `[[`, numeric(1), "weight"))
figure  <- plot_grid(body, shared_legend(), ncol = 1, rel_heights = c(1, 0.045))

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
out_path <- file.path(OUTPUT_DIR, "cutpoints_histogram_faceted_summary.png")

total_rows <- sum(vapply(datasets, function(d) length(d$binning_specs), integer(1)))
ggsave(out_path, figure,
       width = length(MODELS) * 7, height = total_rows * 2.9 + 2.5,
       dpi = 200, bg = "white", limitsize = FALSE)
message("wrote ", out_path)
