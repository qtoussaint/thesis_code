#!/usr/bin/env Rscript
# Single summary figure of the posterior predictive checks (PPCs) for every
# ordinal MIC *inference* model (POM + PPOM), across all drugs and MIC binnings.
#
# Two stacked sections, both laid out as rows = {POM, PPOM} x columns = the ten
# (drug, binning) runs:
#   A  Agreement mosaics  -- regenerated here as vector vcd::agreementplot grobs
#      from each run's inference_ppc/y_rep_vs_true.csv (observed y_true vs the
#      posterior-predictive modal category y_rep_mode), in the style of
#      paper_figures/prediction_accuracy_summary.R.
#   B  Category frequencies -- the pipeline's already-rendered
#      inference_ppc/ppc_category_frequencies.png, tiled into the same grid.
#      (The 89% CI band needs the full y_rep_ppc draws, which live only in the
#      21 GB fit RDS per run, so we reuse the per-run PNG rather than recompute.)
#
# Run with:
#   mamba activate gwas_pipeline
#   Rscript paper_figures/inference_ppc_summary.R
#
# Output: <results>/paper_figures/inference_ppc_summary/inference_ppc_summary.png

suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
  library(vcd)     # agreementplot (section A)
  library(grid)
  library(caret)   # confusionMatrix
})

# Null device so grid.grabExpr (section A) doesn't fall back to R's default
# device and leave a stray Rplots.pdf behind when run non-interactively.
grDevices::pdf(NULL)
on.exit(grDevices::dev.off(), add = TRUE)

RESULTS_ROOT <- "/nfs/research/jlees/jacqueline/thesis_results"
OUTPUT_DIR   <- file.path(RESULTS_ROOT, "paper_figures", "inference_ppc_summary")

# -----------------------------------------------------------------------------
# Specs: one row per (drug, binning); each fits both POM and PPOM. nn/run_stub/K
# mirror ordinal_specs in prediction_accuracy_summary.R; `breaks` are the
# mic_breakpoints (length K-1) read from the dataset JSONs.
# -----------------------------------------------------------------------------
specs <- list(
  list(species_dir="spn_penicillin",   drug="SPN PEN", binning="standard",     K=8, nn="02", run_stub="spn_penicillin_MIC",                   breaks=c(0.016,0.03,0.06,0.12,0.25,1,2)),
  list(species_dir="spn_penicillin",   drug="SPN PEN", binning="coarse",       K=5, nn="10", run_stub="spn_penicillin_MIC_coarse_dilutions",  breaks=c(0.016,0.06,0.25,1)),
  list(species_dir="spn_penicillin",   drug="SPN PEN", binning="large minbin", K=4, nn="11", run_stub="spn_penicillin_MIC_large_minbin",      breaks=c(0.03,0.12,1)),
  list(species_dir="spn_penicillin",   drug="SPN PEN", binning="minima",       K=5, nn="16", run_stub="spn_penicillin_MIC_minimabinning",     breaks=c(0.032,0.065,0.2,2)),
  list(species_dir="spn_trimethoprim", drug="SPN TMP", binning="standard",     K=5, nn="05", run_stub="spn_trimethoprim_MIC",                  breaks=c(0.12,0.25,1,2)),
  list(species_dir="spn_trimethoprim", drug="SPN TMP", binning="coarse",       K=3, nn="12", run_stub="spn_trimethoprim_MIC_coarse_dilutions", breaks=c(0.25,1)),
  list(species_dir="spn_trimethoprim", drug="SPN TMP", binning="large minbin", K=3, nn="13", run_stub="spn_trimethoprim_MIC_large_minbin",     breaks=c(0.25,2)),
  list(species_dir="tb_rifampicin",    drug="TB RIF",  binning="standard",     K=5, nn="08", run_stub="tb_rifampicin_MIC",                     breaks=c(0.06,0.12,2,4)),
  list(species_dir="tb_rifampicin",    drug="TB RIF",  binning="coarse",       K=4, nn="14", run_stub="tb_rifampicin_MIC_coarse_dilutions",    breaks=c(0.06,0.25,4)),
  list(species_dir="tb_rifampicin",    drug="TB RIF",  binning="large minbin", K=4, nn="15", run_stub="tb_rifampicin_MIC_large_minbin",        breaks=c(0.06,2,4))
)

# inference_ppc dir for one (spec, model).
ppc_dir <- function(spec, model) file.path(
  RESULTS_ROOT, paste0("gwas_", spec$species_dir), "inference",
  paste0(spec$nn, "_", spec$run_stub, "_", model), "inference_ppc")

missing_files <- character(0)
note_missing  <- function(path) missing_files <<- c(missing_files, path)

# Placeholder for a run whose PPC outputs aren't on disk yet (job still running
# or failed): a centred asterisk rather than an empty cell.
missing_cell <- function() ggdraw() +
  draw_label("*", size = 28, colour = "grey55", fontface = "bold")

# -----------------------------------------------------------------------------
# Aesthetics
# -----------------------------------------------------------------------------
drug_full <- c("SPN PEN" = "Penicillin~(italic('S. pneumoniae'))",
               "SPN TMP" = "Trimethoprim~(italic('S. pneumoniae'))",
               "TB RIF"  = "Rifampicin~(italic('M. tuberculosis'))")

# Bold rotated row label (POM / PPOM), as in prediction_accuracy_summary.R.
row_label <- function(txt) ggdraw() + draw_label(txt, fontface = "bold", angle = 90, size = 13)

# Two-line column header: drug (parsed, italic species) over binning.
col_header <- function(spec) ggdraw() +
  draw_label(parse(text = drug_full[[spec$drug]]), y = 0.70, size = 9) +
  draw_label(spec$binning, y = 0.30, size = 9, fontface = "italic", colour = "grey25")

# -----------------------------------------------------------------------------
# Section A: agreement mosaics
# -----------------------------------------------------------------------------
# Blank out interval labels whose category band centre is within `min_sep` (npc)
# of the previous kept label, so agreementplot's band-centred labels don't pile
# up where a category is sparse. `freqs` are the marginal counts in order.
thin_labels <- function(freqs, labels, min_sep = 0.05) {
  centres <- cumsum(freqs) / sum(freqs) - (freqs / sum(freqs)) / 2
  last <- -Inf
  for (i in seq_along(labels)) {
    if (centres[i] - last < min_sep) labels[i] <- "" else last <- centres[i]
  }
  labels
}

# One agreement-plot grob from a run's y_rep_vs_true.csv (or blank if missing).
# Observed = y_true, predicted = y_rep_mode (the posterior-predictive modal
# category, already saved for both POM and PPOM).
agreement_cell <- function(spec, model) {
  csv <- file.path(ppc_dir(spec, model), "y_rep_vs_true.csv")
  if (!file.exists(csv)) { note_missing(csv); return(missing_cell()) }

  d <- tryCatch(read.csv(csv), error = function(e) NULL)
  if (is.null(d) || !all(c("y_true", "y_rep_mode") %in% names(d))) return(missing_cell())

  K     <- spec$K
  brk   <- spec$breaks
  # MIC interval labels when breaks are the expected length, else plain indices.
  labs  <- if (length(brk) == K - 1)
             c(paste0("≤", brk), paste0(">", brk[K - 1])) else as.character(seq_len(K))
  truth <- factor(d$y_true,     levels = seq_len(K)); levels(truth) <- labs
  pred  <- factor(d$y_rep_mode, levels = seq_len(K)); levels(pred)  <- labs
  tab   <- caret::confusionMatrix(data = pred, reference = truth)$table
  weights <- c(1, 1 - 1 / (ncol(tab) - 1) ^ 2)

  # Thin labels: columns drive the observed (x) axis, rows the predicted (y).
  colnames(tab) <- thin_labels(colSums(tab), colnames(tab))
  rownames(tab) <- thin_labels(rowSums(tab), rownames(tab))

  g <- grid::grid.grabExpr({
    vcd::agreementplot(
      x        = tab,
      weights  = weights,
      fill_col = function(j)
        grDevices::colorRampPalette(c("steelblue", "lightblue"))(length(weights))[j],
      xlab = "", ylab = "", xscale = FALSE, yscale = FALSE,
      xlab_rot = 90, xlab_just = "right", ylab_rot = 0, ylab_just = "right",
      prefix = "ap", pop = FALSE, newpage = FALSE,
      margins = grid::unit(c(3.0, 2.6, 0.6, 0.6), "lines"),
      gp = grid::gpar(fontsize = 8)
    )
    grid::seekViewport("ap agreementplot")
    grid::grid.text("observed (µg⋅mL⁻¹)",
                    y = grid::unit(-0.20, "npc"), gp = grid::gpar(fontsize = 8))
    grid::grid.text("predicted", x = grid::unit(-0.17, "npc"), rot = 90,
                    gp = grid::gpar(fontsize = 8))
    grid::upViewport(0)
  })
  ggdraw() + draw_grob(g)
}

# -----------------------------------------------------------------------------
# Section B: reuse the pipeline's per-run category-frequency PNG.
# -----------------------------------------------------------------------------
catfreq_cell <- function(spec, model) {
  png <- file.path(ppc_dir(spec, model), "ppc_category_frequencies.png")
  if (!file.exists(png)) { note_missing(png); return(missing_cell()) }
  ggdraw() + draw_image(png)
}

# -----------------------------------------------------------------------------
# Assemble a section: header strip + POM row + PPOM row.
# -----------------------------------------------------------------------------
cw <- c(0.05, rep(1, length(specs)))  # left label strip + one column per spec

section <- function(title, cell_fn) {
  banner <- ggdraw() + draw_label(title, fontface = "bold", size = 15,
                                  x = 0.005, hjust = 0)
  header <- plot_grid(NULL, plotlist = lapply(specs, col_header), nrow = 1, rel_widths = cw)
  pom    <- plot_grid(plotlist = c(list(row_label("POM")),
                                   lapply(specs, cell_fn, model = "POM")),
                      nrow = 1, rel_widths = cw)
  ppom   <- plot_grid(plotlist = c(list(row_label("PPOM")),
                                   lapply(specs, cell_fn, model = "PPOM")),
                      nrow = 1, rel_widths = cw)
  plot_grid(banner, header, pom, ppom, ncol = 1, rel_heights = c(0.12, 0.18, 1, 1))
}

# -----------------------------------------------------------------------------
# Build + save
# -----------------------------------------------------------------------------
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

panel_a <- section("A  Observed versus posterior-predictive category", agreement_cell)
panel_b <- section("B  Observed vs. posterior-predictive category frequencies", catfreq_cell)

figure <- plot_grid(
  panel_a, panel_b,
  ncol = 1,
  rel_heights = c(1, 1)
) + theme(plot.margin = margin(10, 8, 6, 8))

if (length(missing_files) > 0) {
  message(sprintf("NOTE: %d inference PPC file(s) missing (rendered as gaps):",
                  length(missing_files)))
  for (f in unique(missing_files)) message("  ", f)
} else {
  message("All inference PPC files found.")
}

png_path <- file.path(OUTPUT_DIR, "inference_ppc_summary.png")
ggsave(png_path, figure, width = 26, height = 13, dpi = 200, bg = "white", limitsize = FALSE)
message("wrote ", png_path)
