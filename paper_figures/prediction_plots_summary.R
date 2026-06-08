#!/usr/bin/env Rscript
# Per-drug composite of the prediction-accuracy plots. One figure per species/drug
# and evaluation split (random + lineage/LOSO). Each figure has one row per MIC
# binning; within a row the columns are:
#
#   POM RPS distribution | POM agreement | PPOM RPS distribution | PPOM agreement |
#   POM-vs-PPOM RPS density overlay
#
# Everything is reconstructed natively (vector) from the pipeline's per-run CSVs
# rather than embedding the saved PNGs, so the panels match in style and stay
# crisp. The RPS distributions mirror gwas_workflow R/prediction_accuracy.R; the
# agreement mosaics reuse the reconstruction in prediction_accuracy_summary.R; the
# overlay is the "PPOM/POM model comparison plot" from
# bayesian_gwas_paper 03_prediction/RPS_and_RPSS.R.
#
# Run with:
#   mamba activate gwas_pipeline
#   Rscript paper_figures/prediction_plots_summary.R
#
# Output: <output_dir>/prediction_plots/prediction_plots_<drug>_<split>.{png,csv}

suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
  library(viridis)
  library(vcd)      # agreementplot
  library(grid)
  library(caret)    # confusionMatrix
})

# Null device so grid.grabExpr / incidental plotting doesn't leave a stray
# Rplots.pdf behind when run non-interactively.
grDevices::pdf(NULL)
on.exit(grDevices::dev.off(), add = TRUE)

RESULTS_ROOT <- "/nfs/research/jlees/jacqueline/thesis_results"
OUTPUT_DIR   <- file.path(RESULTS_ROOT, "paper_figures", "prediction_plots")

# -----------------------------------------------------------------------------
# Specs: one entry per (drug, binning). nn is the on-disk run-dir number prefix,
# run_stub the rest of it. Mirrors ordinal_specs in prediction_accuracy_summary.R.
# -----------------------------------------------------------------------------

drugs <- list(
  list(key = "spn_penicillin", species_dir = "spn_penicillin",
       title = "Penicillin~resistance~(italic('S. pneumoniae'))",
       binnings = list(
         list(binning = "standard",     nn = "02", run_stub = "spn_penicillin_MIC"),
         list(binning = "coarse",       nn = "10", run_stub = "spn_penicillin_MIC_coarse_dilutions"),
         list(binning = "large minbin", nn = "11", run_stub = "spn_penicillin_MIC_large_minbin"),
         list(binning = "minima",       nn = "16", run_stub = "spn_penicillin_MIC_minimabinning"))),
  list(key = "spn_trimethoprim", species_dir = "spn_trimethoprim",
       title = "Trimethoprim~resistance~(italic('S. pneumoniae'))",
       binnings = list(
         list(binning = "standard",     nn = "05", run_stub = "spn_trimethoprim_MIC"),
         list(binning = "coarse",       nn = "12", run_stub = "spn_trimethoprim_MIC_coarse_dilutions"),
         list(binning = "large minbin", nn = "13", run_stub = "spn_trimethoprim_MIC_large_minbin"))),
  list(key = "tb_rifampicin", species_dir = "tb_rifampicin",
       title = "Rifampicin~resistance~(italic('M. tuberculosis'))",
       binnings = list(
         list(binning = "standard",     nn = "08", run_stub = "tb_rifampicin_MIC"),
         list(binning = "coarse",       nn = "14", run_stub = "tb_rifampicin_MIC_coarse_dilutions"),
         list(binning = "large minbin", nn = "15", run_stub = "tb_rifampicin_MIC_large_minbin")))
)

SPLITS <- list(random = "Random", loso = "Lineage")

# Two-colour viridis matching prediction_accuracy_summary.R (POM, PPOM).
model_cols <- viridis::viridis(3, option = "D", end = 0.9)[1:2]
names(model_cols) <- c("POM", "PPOM")

# -----------------------------------------------------------------------------
# Path builders
# -----------------------------------------------------------------------------

# prediction_results dir for one run.
pred_dir_for <- function(species_dir, nn, run_stub, model, split) file.path(
  RESULTS_ROOT, paste0("gwas_", species_dir), "prediction",
  paste0(nn, "_", run_stub, "_", model, "_", split), "prediction_results")

# Dataset stub dir holds the split-specific true phenotypes + data JSON. The loso
# split lives in a parallel "<stub>_loso" directory with a different test set.
dataset_stub <- function(nn, run_stub, split)
  if (split == "loso") paste0(nn, "_", run_stub, "_loso") else paste0(nn, "_", run_stub)

true_csv_for <- function(nn, run_stub, split) {
  stub <- dataset_stub(nn, run_stub, split)
  file.path(RESULTS_ROOT, "gwas_datasets", "prediction", stub,
            paste0(stub, "_test_phenotypes.csv"))
}

json_for <- function(nn, run_stub, split) {
  stub <- dataset_stub(nn, run_stub, split)
  file.path(RESULTS_ROOT, "gwas_datasets", "prediction", stub, paste0(stub, ".json"))
}

missing <- character(0)
note_missing <- function(path) missing <<- c(missing, path)

# K + mic_breakpoints from the dataset JSON (breaks differ per binning). The JSON
# is ~39 MB (it carries the full genotype matrices), so parsing it whole would OOM
# the login node's memory governor. K and mic_breakpoints are the final two keys,
# so we read only the file's last few KB and pull them out by regex.
read_kbreaks <- function(nn, run_stub, split) {
  path <- json_for(nn, run_stub, split)
  if (!file.exists(path)) { note_missing(path); return(NULL) }
  sz  <- file.info(path)$size
  con <- file(path, "rb"); on.exit(close(con))
  if (sz > 4096) seek(con, sz - 4096)
  tail <- readChar(con, 4096, useBytes = TRUE)
  # (?s) = dotall, so .* spans the newlines in the pretty-printed JSON.
  K  <- suppressWarnings(as.integer(sub('(?s).*"K"\\s*:\\s*([0-9]+).*', '\\1', tail, perl = TRUE)))
  mb <- sub('(?s).*"mic_breakpoints"\\s*:\\s*\\[([^]]*)\\].*', '\\1', tail, perl = TRUE)
  breaks <- suppressWarnings(as.numeric(strsplit(mb, ",")[[1]]))
  if (is.na(K) || anyNA(breaks)) return(NULL)
  list(K = K, breaks = breaks)
}

# rpss_uniform / rpss_frequency for the RPS-distribution annotation.
read_rpss <- function(pred_dir) {
  path <- file.path(pred_dir, "prediction_accuracy_metrics.csv")
  if (!file.exists(path)) { note_missing(path); return(c(uniform = NA, frequency = NA)) }
  row <- tryCatch(read.csv(path, check.names = FALSE)[1, , drop = FALSE],
                  error = function(e) NULL)
  if (is.null(row)) return(c(uniform = NA, frequency = NA))
  c(uniform = suppressWarnings(as.numeric(row[["rpss_uniform"]])),
    frequency = suppressWarnings(as.numeric(row[["rpss_frequency"]])))
}

# Per-sample scaled RPS vector (already computed by the pipeline).
read_rps <- function(pred_dir) {
  path <- file.path(pred_dir, "rps_scores.csv")
  if (!file.exists(path)) { note_missing(path); return(NULL) }
  v <- tryCatch(read.csv(path)$rps_scaled, error = function(e) NULL)
  if (is.null(v) || !length(v)) return(NULL)
  as.numeric(v)
}

# -----------------------------------------------------------------------------
# Cell builders
# -----------------------------------------------------------------------------

base_size <- 11

# RPS distribution: histogram + count-scaled KDE + median (dashed) / mean (dotted)
# vlines + an RPSS text annotation. Mirrors prediction_accuracy.R:152-258.
rps_density_cell <- function(rps, rpss, model) {
  if (is.null(rps)) return(ggdraw())
  med <- median(rps); mu <- mean(rps)
  rng <- diff(range(rps)); bw <- if (rng > 0) rng / 52 else 0.02
  n   <- length(rps)
  # Median/mean lines carried by a linetype scale so they generate a legend
  # (dashed = median, dotted = mean) tucked top-right under the RPSS text.
  stat_lines <- data.frame(
    stat = factor(c("median", "mean"), levels = c("median", "mean")),
    x = c(med, mu))
  ggplot(data.frame(RPS = rps), aes(x = RPS)) +
    geom_histogram(binwidth = bw, fill = "grey85", colour = "grey30", linewidth = 0.2) +
    geom_density(aes(y = after_stat(density) * bw * n),
                 colour = model_cols[[model]], linewidth = 0.7, bw = "nrd0", adjust = 0.5) +
    geom_vline(data = stat_lines, aes(xintercept = x, linetype = stat),
               colour = "grey30", linewidth = 0.5, key_glyph = "path") +
    scale_linetype_manual(values = c(median = "dashed", mean = "dotted"), name = NULL) +
    annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.3, size = 3,
             label = sprintf("RPSS_unif = %.3f\nRPSS_freq = %.3f",
                             rpss[["uniform"]], rpss[["frequency"]])) +
    labs(x = "scaled RPS", y = "count") +
    theme_minimal(base_size = base_size) +
    theme(panel.grid.minor = element_blank(),
          legend.position = c(0.99, 0.66), legend.justification = c(1, 1),
          legend.key.height = grid::unit(0.8, "lines"),
          legend.key.width = grid::unit(1.7, "lines"),
          legend.text = element_text(size = 8),
          legend.background = element_rect(fill = scales::alpha("white", 0.6), colour = NA),
          legend.margin = margin(1, 2, 1, 2))
}

# Point predictions, mirroring the pipeline: POM tabulates the integer category
# draws and takes the modal category; PPOM reuses the saved point predictions.
point_predictions_for <- function(pred_dir, model, K) {
  if (model == "POM") {
    draws <- as.matrix(read.csv(file.path(pred_dir, "predicted_draws.csv")))
    # Modal predicted category per test case: tabulate the integer draws in each
    # column. Done column-wise so we never materialise a 2nd big matrix.
    pp <- apply(draws, 2, function(d) which.max(tabulate(d, nbins = K)))
    rm(draws); gc(FALSE)
    pp
  } else {
    read.csv(file.path(pred_dir, "point_predictions.csv"))$point_prediction
  }
}

# Blank labels whose category-band centre is within `min_sep` of the previous kept
# label, so sparse categories' band-centred labels don't pile up.
thin_labels <- function(freqs, labels, min_sep = 0.05) {
  centres <- cumsum(freqs) / sum(freqs) - (freqs / sum(freqs)) / 2
  last <- -Inf
  for (i in seq_along(labels)) {
    if (centres[i] - last < min_sep) labels[i] <- "" else last <- centres[i]
  }
  labels
}

# Agreement mosaic for one (run, model), captured as a grid grob. Reconstruction
# adapted from prediction_accuracy_summary.R:256-305.
agreement_cell <- function(pred_dir, true_csv, model, K, breaks) {
  need <- if (model == "POM") "predicted_draws.csv" else "point_predictions.csv"
  if (!file.exists(file.path(pred_dir, need))) { note_missing(file.path(pred_dir, need)); return(ggdraw()) }
  if (!file.exists(true_csv)) { note_missing(true_csv); return(ggdraw()) }
  if (length(breaks) != K - 1) return(ggdraw())

  true_ph <- read.csv(true_csv)$true_phenotype
  pred    <- point_predictions_for(pred_dir, model, K)

  mic_intervals <- c(paste0("≤", breaks), paste0(">", breaks[K - 1]))
  mic_truth <- factor(true_ph, levels = seq_len(K)); levels(mic_truth) <- mic_intervals
  mic_pred  <- factor(pred,    levels = seq_len(K)); levels(mic_pred)  <- mic_intervals
  tab     <- caret::confusionMatrix(data = mic_pred, reference = mic_truth)$table
  weights <- c(1, 1 - 1 / (ncol(tab) - 1) ^ 2)

  colnames(tab) <- thin_labels(colSums(tab), colnames(tab))
  rownames(tab) <- thin_labels(rowSums(tab), rownames(tab))

  g <- grid::grid.grabExpr({
    vcd::agreementplot(
      x = tab, weights = weights,
      fill_col = function(j)
        grDevices::colorRampPalette(c("steelblue", "lightblue"))(length(weights))[j],
      xlab = "", ylab = "", xscale = FALSE, yscale = FALSE,
      xlab_rot = 90, xlab_just = "right", ylab_rot = 0, ylab_just = "right",
      prefix = "ap", pop = FALSE, newpage = FALSE,
      margins = grid::unit(c(4.0, 1.8, 0.8, 0.8), "lines"),
      gp = grid::gpar(fontsize = 7))
    grid::seekViewport("ap agreementplot")
    # Titles sit just clear of the band-centred MIC tick labels: the observed (x)
    # title tucks below the rotated bottom labels, the predicted (y) title hugs the
    # left labels.
    grid::grid.text("observed interval (µg⋅mL⁻¹)",
                    y = grid::unit(-0.34, "npc"), gp = grid::gpar(fontsize = 10))
    grid::grid.text("predicted interval", x = grid::unit(-0.04, "npc"), rot = 90,
                    gp = grid::gpar(fontsize = 10))
    grid::upViewport(0)
  })
  ggdraw() + draw_grob(g)
}

# POM-vs-PPOM RPS density overlay (RPS_and_RPSS.R:80-88).
overlay_cell <- function(rps_pom, rps_ppom) {
  parts <- list()
  if (!is.null(rps_pom))  parts[[length(parts) + 1]] <- data.frame(model = "POM",  value = rps_pom)
  if (!is.null(rps_ppom)) parts[[length(parts) + 1]] <- data.frame(model = "PPOM", value = rps_ppom)
  if (!length(parts)) return(ggdraw())
  d <- do.call(rbind, parts)
  d$model <- factor(d$model, levels = c("POM", "PPOM"))
  ggplot(d, aes(x = value, colour = model, fill = model)) +
    geom_density(alpha = 0.1, linewidth = 0.7) +
    scale_colour_manual(values = model_cols) +
    scale_fill_manual(values = model_cols, guide = "none") +
    labs(x = "ranked probability score (RPS)", y = "density", colour = "model") +
    theme_minimal(base_size = base_size) +
    theme(panel.grid.minor = element_blank(), legend.position = c(0.8, 0.85),
          legend.background = element_rect(fill = "white", colour = NA))
}

# -----------------------------------------------------------------------------
# Per-figure builder (one drug x split)
# -----------------------------------------------------------------------------

col_titles <- c("POM\nRPS distribution", "POM\nagreement",
                "PPOM\nRPS distribution", "PPOM\nagreement",
                "POM vs PPOM\nRPS overlay")
# Left binning-label strip + five content columns.
cw <- c(0.10, 1, 1, 1, 1, 1)

# All bands render at the same pixel width so magick can stack them; a memory
# governor on the login node (~2.8 GB) means we must NOT hold every row's grobs +
# a full-figure raster at once, so each band is rasterised and freed in turn.
ROW_W <- 18; DPI <- 150
render_band <- function(plot, height) {
  f <- tempfile(fileext = ".png")
  ggsave(f, plot, width = ROW_W, height = height, dpi = DPI, bg = "white",
         device = ragg::agg_png, limitsize = FALSE)
  f
}

make_figure <- function(drug, split) {
  message(sprintf("[%s / %s] ", drug$key, split), appendLF = FALSE)
  csv_rows <- list()

  row_label <- function(txt) ggdraw() + draw_label(txt, fontface = "bold", angle = 90, size = 12)

  # Title + column-header bands, rendered and freed up front.
  title_band <- render_band(ggdraw() + draw_label(
    parse(text = paste0(drug$title, "*' — ", SPLITS[[split]], " split'")),
    fontface = "bold", size = 15), height = 0.55)
  header_band <- render_band(plot_grid(NULL,
    plotlist = lapply(col_titles, function(t) ggdraw() + draw_label(t, size = 11)),
    nrow = 1, rel_widths = cw), height = 0.6)

  band_paths <- vapply(drug$binnings, function(b) {
    kb <- read_kbreaks(b$nn, b$run_stub, split)
    K  <- if (is.null(kb)) NA_integer_ else kb$K
    br <- if (is.null(kb)) numeric(0) else kb$breaks
    tcsv <- true_csv_for(b$nn, b$run_stub, split)

    pd_pom  <- pred_dir_for(drug$species_dir, b$nn, b$run_stub, "POM",  split)
    pd_ppom <- pred_dir_for(drug$species_dir, b$nn, b$run_stub, "PPOM", split)
    rps_pom  <- read_rps(pd_pom);  rpss_pom  <- read_rpss(pd_pom)
    rps_ppom <- read_rps(pd_ppom); rpss_ppom <- read_rpss(pd_ppom)

    for (mdl in c("POM", "PPOM")) {
      v <- if (mdl == "POM") rps_pom else rps_ppom
      r <- if (mdl == "POM") rpss_pom else rpss_ppom
      if (!is.null(v)) csv_rows[[length(csv_rows) + 1]] <<- data.frame(
        drug = drug$key, split = split, binning = b$binning, model = mdl,
        rps_scaled = v, rpss_uniform = r[["uniform"]], rpss_frequency = r[["frequency"]],
        stringsAsFactors = FALSE)
    }

    cells <- list(
      row_label(b$binning),
      rps_density_cell(rps_pom, rpss_pom, "POM"),
      if (is.na(K)) ggdraw() else agreement_cell(pd_pom, tcsv, "POM", K, br),
      rps_density_cell(rps_ppom, rpss_ppom, "PPOM"),
      if (is.na(K)) ggdraw() else agreement_cell(pd_ppom, tcsv, "PPOM", K, br),
      overlay_cell(rps_pom, rps_ppom))
    row <- plot_grid(plotlist = cells, nrow = 1, rel_widths = cw)
    f <- render_band(row, height = 2.8)
    rm(row, cells); gc(FALSE)
    f
  }, character(1))

  # Stack the bands with magick (cheap — each band is a small raster), so we never
  # hold the full set of grobs or a full-figure raster in memory at once.
  dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
  png_path <- file.path(OUTPUT_DIR, sprintf("prediction_plots_%s_%s.png", drug$key, split))
  csv_path <- file.path(OUTPUT_DIR, sprintf("prediction_plots_%s_%s.csv", drug$key, split))
  stacked <- magick::image_append(
    magick::image_read(c(title_band, header_band, band_paths)), stack = TRUE)
  magick::image_write(stacked, png_path)
  unlink(c(title_band, header_band, band_paths))
  if (length(csv_rows)) write.csv(do.call(rbind, csv_rows), csv_path, row.names = FALSE)
  message("wrote ", png_path)
}

# -----------------------------------------------------------------------------
# Drive
# -----------------------------------------------------------------------------

# Optional PRED_DRUG / PRED_SPLIT env vars render just one figure (handy for
# re-rendering a single drug/split without redoing all six).
DRUG_FILTER  <- Sys.getenv("PRED_DRUG",  "")
SPLIT_FILTER <- Sys.getenv("PRED_SPLIT", "")
for (drug in drugs) for (split in names(SPLITS)) {
  if (nzchar(DRUG_FILTER)  && drug$key != DRUG_FILTER)  next
  if (nzchar(SPLIT_FILTER) && split   != SPLIT_FILTER) next
  make_figure(drug, split)
}

if (length(missing) > 0) {
  message(sprintf("NOTE: %d input file(s) missing (rendered as gaps):", length(missing)))
  for (f in unique(missing)) message("  ", f)
} else {
  message("All inputs found.")
}
