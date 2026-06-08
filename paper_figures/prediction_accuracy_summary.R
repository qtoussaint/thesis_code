#!/usr/bin/env Rscript
# Summary figure comparing prediction accuracy of the ordinal MIC models (POM,
# PPOM) against the binary logistic baseline, across all drugs, MIC binnings and
# evaluation splits (random + lineage/LOSO).
#
# Reads the same per-run prediction_accuracy_metrics.csv files that
# gwas_finalruns/collect_results.R turns into LaTeX tables. Only `bacc` is shared
# between the ordinal and binary models, so balanced accuracy is the common axis;
# RPSS (skill over uniform / frequency baselines) is shown for POM vs PPOM only.
#
# Run with:
#   mamba activate gwas_pipeline
#   Rscript paper_figures/prediction_accuracy_summary.R
#
# Output: <output_dir>/prediction_accuracy_summary.{png,csv}
#         <output_dir>/prediction_accuracy_summary_lasso.{png,csv}
# where <output_dir> = <results>/paper_figures/prediction_accuracy_summary
#
# The lasso version reads the parallel `prediction_lasso/` run directories
# (LASSO-penalised models) and is otherwise identical.

suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
  library(viridis)
  library(vcd)     # agreementplot (Panel C)
  library(grid)
  library(caret)   # confusionMatrix
})

# Open a null device so grid.grabExpr (Panel C) doesn't fall back to R's default
# device and leave a stray Rplots.pdf behind when run non-interactively.
grDevices::pdf(NULL)
on.exit(grDevices::dev.off(), add = TRUE)

RESULTS_ROOT <- "/nfs/research/jlees/jacqueline/thesis_results"
OUTPUT_DIR   <- file.path(RESULTS_ROOT, "paper_figures", "prediction_accuracy_summary")

# -----------------------------------------------------------------------------
# Path builder + reader (mirrors collect_results.R, but returns numerics)
# -----------------------------------------------------------------------------
# `pred_subdir` selects the run family: "prediction" (default models) or
# "prediction_lasso" (LASSO-penalised models).

pred_csv <- function(species_dir, run_dir, split, pred_subdir) file.path(
  RESULTS_ROOT, paste0("gwas_", species_dir), pred_subdir,
  paste0(run_dir, "_", split), "prediction_results", "prediction_accuracy_metrics.csv")

missing_files <- character(0)

# Read named metric columns from a one-row metrics CSV as numerics. A missing
# file or column yields NA (and the file is recorded so we can report gaps).
read_metrics <- function(csv_path, cols) {
  if (!file.exists(csv_path)) {
    missing_files <<- c(missing_files, csv_path)
    return(setNames(rep(NA_real_, length(cols)), cols))
  }
  row <- tryCatch(
    read.csv(csv_path, stringsAsFactors = FALSE, check.names = FALSE)[1, , drop = FALSE],
    error = function(e) NULL)
  if (is.null(row) || nrow(row) == 0) return(setNames(rep(NA_real_, length(cols)), cols))
  vapply(cols, function(c) if (!c %in% names(row)) NA_real_ else suppressWarnings(as.numeric(row[[c]])),
         numeric(1))
}

# -----------------------------------------------------------------------------
# Restricted-category bacc (mirrors collect_results.R / gwas_workflow)
# -----------------------------------------------------------------------------
# bacc is the only shared metric that can hit a div-by-zero when a test split is
# missing whole categories. collect_results.R recomputes it over the categories
# actually present (off the K x K confusion table) and flags the cell with a
# superscript dagger in the LaTeX tables; here we mirror that and mark the bar
# with a star. RPSS and the binary logistic bacc are never restricted.

# Per-sample artifact paths. Test phenotypes are shared across run families;
# point predictions / draws live under the run family's pred_subdir.
test_pheno  <- function(stub, split) {
  d <- if (split == "loso") paste0(stub, "_loso") else stub
  file.path(RESULTS_ROOT, "gwas_datasets", "prediction", d, paste0(d, "_test_phenotypes.csv"))
}
pred_points <- function(species_dir, run_dir, split, pred_subdir) file.path(
  RESULTS_ROOT, paste0("gwas_", species_dir), pred_subdir,
  paste0(run_dir, "_", split), "prediction_results", "point_predictions.csv")
pred_draws  <- function(species_dir, run_dir, split, pred_subdir) file.path(
  RESULTS_ROOT, paste0("gwas_", species_dir), pred_subdir,
  paste0(run_dir, "_", split), "prediction_results", "predicted_draws.csv")

# bacc averaged over the categories with a true test sample (needs two or more,
# else NA), off the full K x K confusion table so every sample stays in play.
restricted_bacc <- function(yt, yp, K) {
  M  <- table(factor(yt, levels = seq_len(K)), factor(yp, levels = seq_len(K)))
  N  <- sum(M); rs <- rowSums(M); cs <- colSums(M); tp <- diag(M)
  present <- which(rs > 0)
  sens <- tp / rs; spec <- (N - rs - (cs - tp)) / (N - rs)
  if (length(present) >= 2) mean((sens[present] + spec[present]) / 2) else NA_real_
}

# Recompute bacc over present categories from the saved per-sample artifacts, or
# NA if any needed file is missing. PPOM saves point predictions directly; POM
# rebuilds them as the modal category over the integer draws.
recompute_bacc <- function(species_dir, run_dir, model, split, stub, K, pred_subdir) {
  tf <- test_pheno(stub, split)
  if (!file.exists(tf)) return(NA_real_)
  yt <- read.csv(tf)$true_phenotype
  if (model == "PPOM") {
    pf <- pred_points(species_dir, run_dir, split, pred_subdir)
    if (!file.exists(pf)) return(NA_real_)
    yp <- read.csv(pf)$point_prediction
  } else {
    dr <- pred_draws(species_dir, run_dir, split, pred_subdir)
    if (!file.exists(dr)) return(NA_real_)
    draws <- as.matrix(read.csv(dr))            # rows = draws, cols = samples
    probs <- t(apply(draws, 2, function(col) tabulate(col, nbins = K) / length(col)))
    yp    <- apply(probs, 1, which.max)
  }
  if (is.null(yt) || is.null(yp) || length(yt) != length(yp)) return(NA_real_)
  restricted_bacc(yt, yp, K)
}

# Resolve bacc for one (spec, model, split): keep the on-disk value (starred when
# its bacc_restricted flag is set), or recompute over present categories when the
# cell is NA (a recomputed value is restricted by construction).
resolve_bacc <- function(s, model, sp, pred_subdir) {
  run_dir  <- paste0(s$nn, "_", s$run_stub, "_", model)
  csv_path <- pred_csv(s$species_dir, run_dir, sp, pred_subdir)
  raw <- if (file.exists(csv_path))
           tryCatch(read.csv(csv_path, check.names = FALSE)[1, , drop = FALSE],
                    error = function(e) NULL) else NULL
  onv <- if (!is.null(raw) && "bacc" %in% names(raw))
           suppressWarnings(as.numeric(raw[["bacc"]])) else NA_real_
  if (!is.na(onv)) {
    starred <- !is.null(raw) && "bacc_restricted" %in% names(raw) &&
               isTRUE(as.logical(raw[["bacc_restricted"]]))
    list(value = onv, restricted = starred)
  } else {
    val <- recompute_bacc(s$species_dir, run_dir, model, sp,
                          paste0(s$nn, "_", s$run_stub), s$K, pred_subdir)
    list(value = val, restricted = !is.na(val))
  }
}

# -----------------------------------------------------------------------------
# Specs (drug / binning / model layout, from collect_results.R)
# -----------------------------------------------------------------------------

# Ordinal MIC runs: one entry per (drug, binning). Each fits both PPOM and POM.
ordinal_specs <- list(
  list(species_dir="spn_penicillin",   drug="SPN PEN", organism="S. pneumoniae",   binning="standard",     K=8, nn="02", run_stub="spn_penicillin_MIC"),
  list(species_dir="spn_penicillin",   drug="SPN PEN", organism="S. pneumoniae",   binning="coarse",       K=5, nn="10", run_stub="spn_penicillin_MIC_coarse_dilutions"),
  list(species_dir="spn_penicillin",   drug="SPN PEN", organism="S. pneumoniae",   binning="large minbin", K=4, nn="11", run_stub="spn_penicillin_MIC_large_minbin"),
  list(species_dir="spn_penicillin",   drug="SPN PEN", organism="S. pneumoniae",   binning="minima",       K=5, nn="16", run_stub="spn_penicillin_MIC_minimabinning"),
  list(species_dir="spn_trimethoprim", drug="SPN TMP", organism="S. pneumoniae",   binning="standard",     K=5, nn="05", run_stub="spn_trimethoprim_MIC"),
  list(species_dir="spn_trimethoprim", drug="SPN TMP", organism="S. pneumoniae",   binning="coarse",       K=3, nn="12", run_stub="spn_trimethoprim_MIC_coarse_dilutions"),
  list(species_dir="spn_trimethoprim", drug="SPN TMP", organism="S. pneumoniae",   binning="large minbin", K=3, nn="13", run_stub="spn_trimethoprim_MIC_large_minbin"),
  list(species_dir="tb_rifampicin",    drug="TB RIF",  organism="M. tuberculosis", binning="standard",     K=5, nn="08", run_stub="tb_rifampicin_MIC"),
  list(species_dir="tb_rifampicin",    drug="TB RIF",  organism="M. tuberculosis", binning="coarse",       K=4, nn="14", run_stub="tb_rifampicin_MIC_coarse_dilutions"),
  list(species_dir="tb_rifampicin",    drug="TB RIF",  organism="M. tuberculosis", binning="large minbin", K=4, nn="15", run_stub="tb_rifampicin_MIC_large_minbin")
)

# Binary logistic baseline: one per drug. run_stub is the on-disk dir prefix.
binary_specs <- list(
  list(species_dir="spn_penicillin",   drug="SPN PEN", run_stub="01_spn_penicillin_binary_logistic"),
  list(species_dir="spn_trimethoprim", drug="SPN TMP", run_stub="04_spn_trimethoprim_binary_logistic"),
  list(species_dir="tb_rifampicin",    drug="TB RIF",  run_stub="07_tb_rifampicin_binary_logistic")
)

SPLITS <- list(random = "Random", loso = "Lineage")
ORD_METRICS <- c("bacc", "rpss_uniform", "rpss_frequency")

# -----------------------------------------------------------------------------
# Build long data frame (for one run family / pred_subdir)
# -----------------------------------------------------------------------------

build_df <- function(pred_subdir) {
  rows <- list()

  for (s in ordinal_specs) {
    for (model in c("POM", "PPOM")) {
      run_dir <- paste0(s$nn, "_", s$run_stub, "_", model)
      for (sp in names(SPLITS)) {
        m <- read_metrics(pred_csv(s$species_dir, run_dir, sp, pred_subdir), ORD_METRICS)
        # bacc may be NA on disk (div-by-zero on a missing category); resolve it
        # over present-only categories and flag whether the value is restricted.
        bacc_res <- resolve_bacc(s, model, sp, pred_subdir)
        for (metric in ORD_METRICS) {
          is_bacc <- metric == "bacc"
          rows[[length(rows) + 1]] <- data.frame(
            drug = s$drug, organism = s$organism, binning = s$binning, K = s$K,
            model = model, split = SPLITS[[sp]], metric = metric,
            value = if (is_bacc) bacc_res$value else unname(m[[metric]]),
            restricted = if (is_bacc) bacc_res$restricted else FALSE,
            stringsAsFactors = FALSE)
        }
      }
    }
  }

  for (s in binary_specs) {
    for (sp in names(SPLITS)) {
      m <- read_metrics(pred_csv(s$species_dir, s$run_stub, sp, pred_subdir), "bacc")
      rows[[length(rows) + 1]] <- data.frame(
        drug = s$drug, organism = NA_character_, binning = NA_character_, K = 2L,
        model = "Logistic", split = SPLITS[[sp]], metric = "bacc",
        value = unname(m[["bacc"]]), restricted = FALSE, stringsAsFactors = FALSE)
    }
  }

  df <- do.call(rbind, rows)

  # Factor orders
  drug_levels    <- c("SPN PEN", "SPN TMP", "TB RIF")
  binning_levels <- c("standard", "coarse", "large minbin", "minima")
  df$drug    <- factor(df$drug, levels = drug_levels)
  df$binning <- factor(df$binning, levels = binning_levels)
  df$model   <- factor(df$model, levels = c("POM", "PPOM", "Logistic"))
  df$split   <- factor(df$split, levels = c("Random", "Lineage"))
  df
}

# -----------------------------------------------------------------------------
# Aesthetics
# -----------------------------------------------------------------------------

# Three-colour discrete viridis: POM, PPOM, Logistic
model_cols <- viridis::viridis(3, option = "D", end = 0.9)
names(model_cols) <- c("POM", "PPOM", "Logistic")

# Facet strip labels: full drug name + italicised species (plotmath expressions).
drug_labels <- c(
  "SPN PEN" = "Penicillin~resistance~(italic('S. pneumoniae'))",
  "SPN TMP" = "Trimethoprim~resistance~(italic('S. pneumoniae'))",
  "TB RIF"  = "Rifampicin~resistance~(italic('M. tuberculosis'))")
drug_labeller <- labeller(drug = as_labeller(drug_labels, label_parsed))

base_theme <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92", colour = NA),
        legend.position = "bottom")

# -----------------------------------------------------------------------------
# Panel C helpers: agreement mosaics (standard binning, random split only)
# -----------------------------------------------------------------------------
# Recreated from the prediction pipeline's agreement plot (gwas_workflow
# R/prediction_accuracy.R: vcd::agreementplot on the MIC-labelled confusion
# matrix), captured as a grid grob so it renders as vector in this figure rather
# than embedding the pipeline's low-res PNG. Standard binning = the finest
# (max-K) binning per drug; rows are POM / PPOM, columns are the three drugs.
#
# Inputs per run (same as the pipeline shell scripts in gwas_finalruns/):
#   true phenotypes : gwas_datasets/prediction/<stub>/<stub>_test_phenotypes.csv
#   predictions     : <run>/prediction_results/{predicted_draws,point_predictions}.csv
# K and mic_breakpoints come from the data JSON; they are stable inputs so we
# inline them here (K already matches ordinal_specs above).
c_specs <- list(
  list(species_dir = "spn_penicillin",   stub = "02_spn_penicillin_MIC",
       K = 8, breaks = c(0.016, 0.03, 0.06, 0.12, 0.25, 1, 2)),
  list(species_dir = "spn_trimethoprim", stub = "05_spn_trimethoprim_MIC",
       K = 5, breaks = c(0.12, 0.25, 1, 2)),
  list(species_dir = "tb_rifampicin",    stub = "08_tb_rifampicin_MIC",
       K = 5, breaks = c(0.06, 0.12, 2, 4))
)

# Blank out labels whose category band centre is within `min_sep` (npc) of the
# previous kept label, so agreementplot's band-centred labels don't pile up on
# top of each other where a category is sparse (its band is a near-invisible
# sliver anyway). `freqs` are the marginal counts in category order.
thin_labels <- function(freqs, labels, min_sep = 0.05) {
  centres <- cumsum(freqs) / sum(freqs) - (freqs / sum(freqs)) / 2
  last <- -Inf
  for (i in seq_along(labels)) {
    if (centres[i] - last < min_sep) labels[i] <- "" else last <- centres[i]
  }
  labels
}

# Point predictions, mirroring the pipeline: POM tabulates the integer category
# draws and takes the modal category; PPOM reuses the pipeline's saved file.
point_predictions_for <- function(pred_dir, model, K) {
  if (model == "POM") {
    draws <- read.csv(file.path(pred_dir, "predicted_draws.csv"))
    probs <- t(apply(t(draws), 1, function(d) tabulate(d, nbins = K) / length(d)))
    apply(probs, 1, which.max)
  } else {
    read.csv(file.path(pred_dir, "point_predictions.csv"))$point_prediction
  }
}

# Build one agreement-plot grob (or a blank ggdraw if inputs are missing).
agreement_cell <- function(spec, model, pred_subdir) {
  pred_dir <- file.path(RESULTS_ROOT, paste0("gwas_", spec$species_dir), pred_subdir,
                        paste0(spec$stub, "_", model, "_random"), "prediction_results")
  true_csv <- file.path(RESULTS_ROOT, "gwas_datasets", "prediction", spec$stub,
                        paste0(spec$stub, "_test_phenotypes.csv"))
  need <- if (model == "POM") "predicted_draws.csv" else "point_predictions.csv"
  if (!file.exists(file.path(pred_dir, need)) || !file.exists(true_csv)) return(ggdraw())

  K       <- spec$K
  true_ph <- read.csv(true_csv)$true_phenotype
  pred    <- point_predictions_for(pred_dir, model, K)

  # MIC interval labels: ≤b[1] .. ≤b[K-1], >b[K-1]  (as in the pipeline)
  mic_intervals <- c(paste0("≤", spec$breaks), paste0(">", spec$breaks[K - 1]))
  mic_truth <- factor(true_ph, levels = seq_len(K)); levels(mic_truth) <- mic_intervals
  mic_pred  <- factor(pred,    levels = seq_len(K)); levels(mic_pred)  <- mic_intervals
  cm_mic    <- caret::confusionMatrix(data = mic_pred, reference = mic_truth)
  tab       <- cm_mic$table
  weights   <- c(1, 1 - 1 / (ncol(tab) - 1) ^ 2)

  # Thin the interval labels so they stay legible: columns drive the observed
  # (x) axis, rows the predicted (y) axis.
  colnames(tab) <- thin_labels(colSums(tab), colnames(tab))
  rownames(tab) <- thin_labels(rowSums(tab), rownames(tab))

  g <- grid::grid.grabExpr({
    # Drop the marginal count numbers + tick axes (xscale/yscale = FALSE) to
    # declutter, rotate the observed labels vertical, and use tight margins so
    # the square mosaic fills its cell.
    vcd::agreementplot(
      x        = tab,
      weights  = weights,
      fill_col = function(j)
        grDevices::colorRampPalette(c("steelblue", "lightblue"))(length(weights))[j],
      xlab = "", ylab = "", xscale = FALSE, yscale = FALSE,
      xlab_rot = 90, xlab_just = "right", ylab_rot = 0, ylab_just = "right",
      prefix = "ap", pop = FALSE, newpage = FALSE,
      margins = grid::unit(c(3.8, 3.2, 0.8, 0.8), "lines"),
      gp = grid::gpar(fontsize = 9)
    )
    grid::seekViewport("ap agreementplot")
    grid::grid.text("observed interval (µg⋅mL⁻¹)",
                    y = grid::unit(-0.205, "npc"), gp = grid::gpar(fontsize = 10))
    grid::grid.text("predicted interval", x = grid::unit(-0.16, "npc"), rot = 90,
                    gp = grid::gpar(fontsize = 10))
    grid::upViewport(0)
  })
  ggdraw() + draw_grob(g)
}

row_cells <- function(model, pred_subdir) lapply(c_specs, function(s) agreement_cell(s, model, pred_subdir))

# Column headers (drug names) and bold row labels (model).
col_titles <- lapply(drug_labels, function(e)
  ggdraw() + draw_label(parse(text = e), size = 12))
row_label  <- function(txt) ggdraw() + draw_label(txt, fontface = "bold", angle = 90, size = 13)

cw <- c(0.07, 1, 1, 1)  # left label strip + three drug columns
header_row <- plot_grid(NULL, plotlist = col_titles, nrow = 1, rel_widths = cw)

# -----------------------------------------------------------------------------
# Figure builder: assemble all three panels for one run family and save.
# pred_subdir selects the run dirs ("prediction" / "prediction_lasso"),
# file_stub names the output {png,csv}.
# -----------------------------------------------------------------------------

make_figure <- function(pred_subdir, file_stub) {
  missing_files <<- character(0)
  df <- build_df(pred_subdir)

  if (length(missing_files) > 0) {
    message(sprintf("[%s] NOTE: %d metrics CSV(s) missing (rendered as gaps):",
                    pred_subdir, length(missing_files)))
    for (f in missing_files) message("  ", f)
  } else {
    message(sprintf("[%s] All metrics CSVs found.", pred_subdir))
  }

  # --- Panel A: balanced accuracy (POM/PPOM bars + logistic reference line) ---
  bacc_bars <- df[df$metric == "bacc" & df$model %in% c("POM", "PPOM"), ]
  bacc_bars$model <- droplevels(bacc_bars$model)

  # Logistic is constant across binnings within a drug -> one hline per drug x split.
  bacc_logit <- df[df$metric == "bacc" & df$model == "Logistic", ]

  # One combined "model" legend: POM/PPOM as filled bars, logistic as a line.
  # A shared 3-level scale (same name + breaks) merges the fill and linetype
  # guides; logistic's fill key is transparent and the POM/PPOM lines are blank,
  # so each key shows only the relevant glyph.
  model_levels <- c("POM", "PPOM", "logistic")
  bacc_logit$model_key <- factor("logistic", levels = model_levels)
  fill_vals <- c(POM = model_cols[["POM"]], PPOM = model_cols[["PPOM"]], logistic = NA)
  lty_vals  <- c(POM = "blank", PPOM = "blank", logistic = "dashed")

  panel_a <- ggplot(bacc_bars, aes(x = binning, y = value, fill = model)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.7) +
    # Star bars whose bacc was computed on a reduced category set (cf. the
    # superscript dagger collect_results.R puts on the matching LaTeX cells).
    geom_text(aes(label = ifelse(restricted, "*", ""), group = model),
              position = position_dodge(width = 0.75),
              vjust = -0.2, size = 6, show.legend = FALSE) +
    geom_hline(data = bacc_logit,
               aes(yintercept = value, linetype = model_key),
               colour = model_cols[["Logistic"]], linewidth = 0.7) +
    facet_grid(split ~ drug, scales = "free_x", space = "free_x", labeller = drug_labeller) +
    scale_fill_manual(values = fill_vals, name = "model", limits = model_levels) +
    scale_linetype_manual(values = lty_vals, name = "model", limits = model_levels) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = NULL, y = "balanced accuracy (bACC)",
         caption = "* bACC computed on a reduced category set (categories absent from the test split excluded)") +
    base_theme +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          plot.caption = element_text(hjust = 0, size = 9))

  # --- Panel B: RPSS skill scores (POM vs PPOM; uniform + frequency baselines) ---
  rpss <- df[df$metric %in% c("rpss_uniform", "rpss_frequency"), ]
  rpss$baseline <- factor(ifelse(rpss$metric == "rpss_uniform", "uniform", "frequency"),
                          levels = c("uniform", "frequency"))
  rpss$model <- droplevels(rpss$model)

  panel_b <- ggplot(rpss,
                    aes(x = binning, y = value, colour = model, shape = baseline,
                        group = interaction(model, baseline))) +
    geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.4) +
    geom_point(position = position_dodge(width = 0.6), size = 2.8) +
    facet_grid(split ~ drug, scales = "free_x", space = "free_x", labeller = drug_labeller) +
    scale_colour_manual(values = model_cols, name = "model") +
    scale_shape_manual(values = c(uniform = 16, frequency = 17), name = "RPSS baseline") +
    labs(x = NULL, y = "Ranked Probability Skill Score (RPSS)") +
    base_theme +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

  # --- Panel C: agreement mosaics ---
  pom_row  <- plot_grid(plotlist = c(list(row_label("POM")),  row_cells("POM", pred_subdir)),
                        nrow = 1, rel_widths = cw)
  ppom_row <- plot_grid(plotlist = c(list(row_label("PPOM")), row_cells("PPOM", pred_subdir)),
                        nrow = 1, rel_widths = cw)
  panel_c <- plot_grid(header_row, pom_row, ppom_row, ncol = 1, rel_heights = c(0.12, 1, 1))

  # --- Assemble + save ---
  figure <- plot_grid(
    panel_a, panel_b, panel_c,
    ncol = 1,
    labels = c("A", "B", "C"),
    label_size = 24,
    label_fontface = "bold",
    rel_heights = c(1, 1, 1.9)
  ) + theme(plot.margin = margin(5, 5, 5, 10))

  png_path <- file.path(OUTPUT_DIR, paste0(file_stub, ".png"))
  csv_path <- file.path(OUTPUT_DIR, paste0(file_stub, ".csv"))

  ggsave(png_path, figure, width = 14, height = 20, dpi = 300, bg = "white")
  write.csv(df[order(df$metric, df$drug, df$binning, df$model, df$split), ],
            csv_path, row.names = FALSE)

  message("wrote ", png_path)
  message("wrote ", csv_path)
}

# -----------------------------------------------------------------------------
# Build both run families: default models + LASSO-penalised models.
# -----------------------------------------------------------------------------

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
make_figure("prediction",       "prediction_accuracy_summary")
make_figure("prediction_lasso", "prediction_accuracy_summary_lasso")
