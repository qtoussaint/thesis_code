#!/usr/bin/env Rscript
# Shared readers for the out-of-sample parameter-sweep figures.
#
# These figures read the OUT-OF-SAMPLE prediction metrics produced by the runs
# under
#   thesis_results/compare_ordinal_models/parameter_sweeps_prediction/runs/
#       <run_name>/<dataset>/prediction_results/prediction_accuracy_metrics.csv
# (the prediction twin of the in-sample inference_ppc/ metrics the original
# paper_figures/parameter_sweeps scripts read).
#
# For the ordinal (PPOM) sweeps, balanced accuracy can hit a div-by-zero when the
# held-out test split is missing whole MIC categories. We mirror the recompute
# used in paper_figures/prediction_accuracy_summary.R / collect_results.R: when
# the on-disk bacc is NA, recompute it over the categories actually present and
# flag the value as restricted (so the figure can star it).

RESULTS_ROOT <- "/nfs/research/jlees/jacqueline/thesis_results"
RUNS_BASE    <- file.path(RESULTS_ROOT, "compare_ordinal_models",
                          "parameter_sweeps_prediction", "runs")
FIG_DIR      <- file.path(RESULTS_ROOT, "compare_ordinal_models",
                          "parameter_sweeps_prediction", "figures")
DATASETS_DIR <- file.path(RESULTS_ROOT, "gwas_datasets", "prediction")

missing_files <- character(0)

# Out-of-sample metrics CSV for one (run_name, dataset).
pred_metrics_csv <- function(run_name, dataset)
  file.path(RUNS_BASE, run_name, dataset, "prediction_results",
            "prediction_accuracy_metrics.csv")

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

# bacc averaged over the categories with a true test sample (needs two or more,
# else NA), off the full K x K confusion table so every sample stays in play.
restricted_bacc <- function(yt, yp, K) {
  M  <- table(factor(yt, levels = seq_len(K)), factor(yp, levels = seq_len(K)))
  N  <- sum(M); rs <- rowSums(M); cs <- colSums(M); tp <- diag(M)
  present <- which(rs > 0)
  sens <- tp / rs; spec <- (N - rs - (cs - tp)) / (N - rs)
  if (length(present) >= 2) mean((sens[present] + spec[present]) / 2) else NA_real_
}

# Resolve out-of-sample bacc for one PPOM (run_name, dataset): keep the on-disk
# value (starred when its bacc_restricted flag is set), or recompute over present
# categories from the saved point predictions when the cell is NA (recomputed =>
# restricted by construction). PPOM saves point_predictions.csv directly.
resolve_bacc <- function(run_name, dataset, K) {
  pred_dir <- file.path(RUNS_BASE, run_name, dataset, "prediction_results")
  csv_path <- file.path(pred_dir, "prediction_accuracy_metrics.csv")
  raw <- if (file.exists(csv_path))
           tryCatch(read.csv(csv_path, check.names = FALSE)[1, , drop = FALSE],
                    error = function(e) NULL) else NULL
  onv <- if (!is.null(raw) && "bacc" %in% names(raw))
           suppressWarnings(as.numeric(raw[["bacc"]])) else NA_real_
  if (!is.na(onv)) {
    starred <- !is.null(raw) && "bacc_restricted" %in% names(raw) &&
               isTRUE(as.logical(raw[["bacc_restricted"]]))
    return(list(value = onv, restricted = starred))
  }
  tf <- file.path(DATASETS_DIR, dataset, paste0(dataset, "_test_phenotypes.csv"))
  pf <- file.path(pred_dir, "point_predictions.csv")
  if (!file.exists(tf) || !file.exists(pf)) return(list(value = NA_real_, restricted = FALSE))
  yt <- tryCatch(read.csv(tf)$true_phenotype, error = function(e) NULL)
  yp <- tryCatch(read.csv(pf)$point_prediction, error = function(e) NULL)
  if (is.null(yt) || is.null(yp) || length(yt) != length(yp))
    return(list(value = NA_real_, restricted = FALSE))
  val <- restricted_bacc(yt, yp, K)
  list(value = val, restricted = !is.na(val))
}

report_missing <- function() {
  if (length(missing_files) > 0) {
    message(sprintf("NOTE: %d metrics CSV(s) missing (rendered as gaps):", length(missing_files)))
    for (f in missing_files) message("  ", f)
  } else {
    message("All metrics CSVs found.")
  }
}
