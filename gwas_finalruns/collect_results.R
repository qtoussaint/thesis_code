#!/usr/bin/env Rscript
# Collect inference + prediction accuracy metrics from all runs under
# gwas_finalruns/{inference,prediction}/<species>/ and emit three LaTeX
# table bodies (binary logistic, continuous, ordinal POM/PPOM).
#
# Run with:
#   mamba activate gwas_pipeline
#   Rscript /nfs/research/jlees/jacqueline/thesis_code/gwas_finalruns/collect_results.R

RESULTS_ROOT <- "/nfs/research/jlees/jacqueline/thesis_results"
OUT_DIR      <- "/nfs/research/jlees/jacqueline/thesis_results/paper_figures/prediction_accuracy_tables"
MISSING      <- "--"            # null / no value
RESTRICTED   <- "$^{\\dagger}$"  # metric computed on a reduced category set

# -----------------------------------------------------------------------------
# Readers
# -----------------------------------------------------------------------------

fmt <- function(x) if (is.na(x)) MISSING else sprintf("%.3f", x)

read_metrics <- function(csv_path, cols) {
  if (!file.exists(csv_path)) return(setNames(rep(MISSING, length(cols)), cols))
  row <- tryCatch(
    read.csv(csv_path, stringsAsFactors = FALSE, check.names = FALSE)[1, , drop = FALSE],
    error = function(e) NULL)
  if (is.null(row) || nrow(row) == 0) return(setNames(rep(MISSING, length(cols)), cols))
  vapply(cols, function(c) {
    if (!c %in% names(row)) MISSING else fmt(as.numeric(row[[c]]))
  }, character(1))
}

# Read scalar keys from the *tail* of a JSON file. The pipeline JSONs include
# the full variant_matrix and routinely exceed 100 MB, but the bookkeeping
# keys (N, N_ppc, N_train, N_test, K) sit at the end of the object. We pull
# the last 8 KB and grep them out.
read_json_scalars <- function(json_path, keys) {
  if (!file.exists(json_path)) return(setNames(rep(MISSING, length(keys)), keys))
  size <- file.info(json_path)$size
  n <- as.integer(min(size, 8192))
  con <- file(json_path, "rb")
  on.exit(close(con))
  if (size > n) seek(con, where = size - n, origin = "start")
  raw <- readBin(con, what = "raw", n = n)
  txt <- rawToChar(raw)
  vapply(keys, function(k) {
    pat <- sprintf('"%s"[[:space:]]*:[[:space:]]*[0-9]+', k)
    m <- regmatches(txt, regexpr(pat, txt))
    if (length(m) == 0) MISSING
    else sub(sprintf('"%s"[[:space:]]*:[[:space:]]*', k), "", m)
  }, character(1))
}

# -----------------------------------------------------------------------------
# Restricted-category bacc/ppv (mirrors gwas_workflow .restricted_ordinal_metrics)
# -----------------------------------------------------------------------------

# Recompute bacc/ppv over the categories actually present in the test set, off the
# K x K confusion table so every sample stays in play. bacc averages per-class
# balanced accuracy over classes with a true sample (needs two or more, else NA);
# ppv weight-averages precision over classes that were predicted at least once.
restricted_ordinal_metrics <- function(yt, yp, K) {
  M  <- table(factor(yt, levels = seq_len(K)), factor(yp, levels = seq_len(K)))
  N  <- sum(M); rs <- rowSums(M); cs <- colSums(M); tp <- diag(M)
  present <- which(rs > 0); predicted <- which(cs > 0)
  sens <- tp / rs; spec <- (N - rs - (cs - tp)) / (N - rs)
  bacc <- if (length(present) >= 2) mean((sens[present] + spec[present]) / 2) else NA_real_
  prec <- tp / cs
  ppv  <- if (sum(rs[predicted]) > 0) sum(prec[predicted] * rs[predicted]) / sum(rs[predicted]) else NA_real_
  list(bacc = bacc, ppv = ppv)
}

# -----------------------------------------------------------------------------
# Path builders
# -----------------------------------------------------------------------------

ppc_csv  <- function(species_dir, run_dir) file.path(
  RESULTS_ROOT, paste0("gwas_", species_dir), "inference",
  run_dir, "inference_ppc", "prediction_accuracy_metrics.csv")

pred_csv <- function(species_dir, run_dir, split, pred_subdir = "prediction") file.path(
  RESULTS_ROOT, paste0("gwas_", species_dir), pred_subdir,
  paste0(run_dir, "_", split), "prediction_results", "prediction_accuracy_metrics.csv")

inf_json  <- function(dataset) file.path(
  RESULTS_ROOT, "gwas_datasets", "inference", dataset, paste0(dataset, ".json"))

pred_json <- function(dataset, split) {
  d <- if (split == "loso") paste0(dataset, "_loso") else dataset
  file.path(RESULTS_ROOT, "gwas_datasets", "prediction", d, paste0(d, ".json"))
}

# Per-sample artifacts used to recompute ordinal bacc/ppv over present categories.
ppc_yrep    <- function(species_dir, run_dir) file.path(
  RESULTS_ROOT, paste0("gwas_", species_dir), "inference",
  run_dir, "inference_ppc", "y_rep_vs_true.csv")

pred_points <- function(species_dir, run_dir, split, pred_subdir = "prediction") file.path(
  RESULTS_ROOT, paste0("gwas_", species_dir), pred_subdir,
  paste0(run_dir, "_", split), "prediction_results", "point_predictions.csv")

pred_draws  <- function(species_dir, run_dir, split, pred_subdir = "prediction") file.path(
  RESULTS_ROOT, paste0("gwas_", species_dir), pred_subdir,
  paste0(run_dir, "_", split), "prediction_results", "predicted_draws.csv")

test_pheno  <- function(dataset, split) {
  d <- if (split == "loso") paste0(dataset, "_loso") else dataset
  file.path(RESULTS_ROOT, "gwas_datasets", "prediction", d,
            paste0(d, "_test_phenotypes.csv"))
}

# Returns list(ntr, nte, metrics) for one (species, run, eval) triple.
collect_one <- function(species_dir, dataset_nn, dataset_base, run_dir, eval, metric_cols,
                        pred_subdir = "prediction") {
  ds_full <- paste0(dataset_nn, "_", dataset_base)
  if (eval == "PPC") {
    n <- read_json_scalars(inf_json(ds_full), c("N", "N_ppc"))
    m <- read_metrics(ppc_csv(species_dir, run_dir), metric_cols)
    list(ntr = n[["N"]], nte = n[["N_ppc"]], m = m)
  } else {
    split <- if (eval == "Random") "random" else "loso"
    n <- read_json_scalars(pred_json(ds_full, split), c("N_train", "N_test"))
    m <- read_metrics(pred_csv(species_dir, run_dir, split, pred_subdir), metric_cols)
    list(ntr = n[["N_train"]], nte = n[["N_test"]], m = m)
  }
}

# Load the per-sample artifacts for one ordinal (run, model, eval) and recompute
# bacc/ppv, or NULL if any needed file is missing. PPC reads y_rep_vs_true.csv;
# prediction reads true labels from the dataset and point predictions from the run
# (PPOM saves them directly, POM should rebuild them from the category draws).
recompute_ordinal_stars <- function(species_dir, run_dir, model, eval, dataset, split, K,
                                    pred_subdir = "prediction") {
  if (eval == "PPC") {
    f <- ppc_yrep(species_dir, run_dir)
    if (!file.exists(f)) return(NULL)
    d <- read.csv(f); yt <- d$y_true; yp <- d$y_rep_mode
  } else {
    tf <- test_pheno(dataset, split)
    if (!file.exists(tf)) return(NULL)
    yt <- read.csv(tf)$true_phenotype
    if (model == "PPOM") {
      pf <- pred_points(species_dir, run_dir, split, pred_subdir)
      if (!file.exists(pf)) return(NULL)
      yp <- read.csv(pf)$point_prediction
    } else {                                   # POM: rebuild point predictions from draws
      df <- pred_draws(species_dir, run_dir, split, pred_subdir)
      if (!file.exists(df)) return(NULL)
      draws <- as.matrix(read.csv(df))         # rows = draws, cols = samples
      probs <- t(apply(draws, 2, function(col) tabulate(col, nbins = K) / length(col)))
      yp    <- apply(probs, 1, which.max)
    }
  }
  if (is.null(yt) || is.null(yp) || length(yt) != length(yp)) return(NULL)
  restricted_ordinal_metrics(yt, yp, K)
}

# Star the bacc/ppv cells only. A non-NA on-disk value is kept and starred when the
# CSV's *_restricted flag is set; an NA cell is recomputed from artifacts and starred
# (a recomputed value is restricted by construction), or left XXX if it cannot be.
star_ordinal_cells <- function(m, csv_path, species_dir, run_dir, model, eval, dataset, split, K,
                               pred_subdir = "prediction") {
  raw <- if (file.exists(csv_path))
           tryCatch(read.csv(csv_path, check.names = FALSE)[1, , drop = FALSE],
                    error = function(e) NULL) else NULL
  rc <- NULL  # recompute lazily, only if a cell is NA
  for (col in c("bacc", "ppv")) {
    onv  <- if (!is.null(raw) && col %in% names(raw)) suppressWarnings(as.numeric(raw[[col]])) else NA_real_
    flag <- paste0(col, "_restricted")
    if (!is.na(onv)) {
      starred <- !is.null(raw) && flag %in% names(raw) && isTRUE(as.logical(raw[[flag]]))
      m[[col]] <- paste0(fmt(onv), if (starred) RESTRICTED else "")
    } else {
      if (is.null(rc)) rc <- recompute_ordinal_stars(species_dir, run_dir, model, eval, dataset, split, K, pred_subdir)
      val <- if (!is.null(rc)) rc[[col]] else NA_real_
      m[[col]] <- if (is.na(val)) MISSING else paste0(fmt(val), RESTRICTED)
    }
  }
  m
}

# -----------------------------------------------------------------------------
# Specs (one entry per drug)
# -----------------------------------------------------------------------------

EVALS <- c("PPC", "Random", "Lineage")

# Binary logistic. The upstream prediction_accuracy_metrics.csv computes
# sensitivity/specificity with SUSCEPTIBLE as the positive class, but the paper
# reports them resistance-positive (sensitivity = recall of resistant). Since
# resistance is simply the other class this is exactly a swap, so we read the
# CSV's specificity into the table's Sensitivity column and its sensitivity into
# Specificity. bacc/auc/brier are symmetric and vme/me are already defined
# relative to resistance, so only these two columns move; f1 is already the
# resistant-class F1.
BINARY_COLS <- c("specificity","sensitivity","bacc","auc","brier","f1","vme","me")
binary_specs <- list(
  list(species_dir="spn_penicillin",   org="\\textit{S. pneumoniae}",   amr="PEN", nn="01", base="spn_penicillin_binary",   run="spn_penicillin_binary_logistic"),
  list(species_dir="spn_trimethoprim", org="\\textit{S. pneumoniae}",   amr="TMP", nn="04", base="spn_trimethoprim_binary", run="spn_trimethoprim_binary_logistic"),
  list(species_dir="tb_rifampicin",    org="\\textit{M. tuberculosis}", amr="RIF", nn="07", base="tb_rifampicin_binary",    run="tb_rifampicin_binary_logistic")
)

# Continuous
CONT_COLS <- c("rmse","crps","r_squared","mae","essential_agreement")
continuous_specs <- list(
  list(species_dir="spn_penicillin",   org="\\textit{S. pneumoniae}",   amr="PEN", nn="03", base="spn_penicillin_continuous",   run="spn_penicillin_continuous_continuous"),
  list(species_dir="spn_trimethoprim", org="\\textit{S. pneumoniae}",   amr="TMP", nn="06", base="spn_trimethoprim_continuous", run="spn_trimethoprim_continuous_continuous"),
  list(species_dir="tb_rifampicin",    org="\\textit{M. tuberculosis}", amr="RIF", nn="09", base="tb_rifampicin_continuous",    run="tb_rifampicin_continuous_continuous")
)

# Ordinal (POM + PPOM share the same dataset; we list one entry per binning)
ORD_COLS <- c("bacc","ppv","rps_mean_scaled","rps_median_scaled","rpss_uniform","rpss_frequency")
ordinal_specs <- list(
  # SPN PEN
  list(species_dir="spn_penicillin", org="\\textit{S. pneumoniae}", amr="PEN", binning="standard ($\\geq$5\\%)",  K=8, nn="02", base="spn_penicillin_MIC",                  run_stub="spn_penicillin_MIC"),
  list(species_dir="spn_penicillin", org="\\textit{S. pneumoniae}", amr="PEN", binning="coarse ($\\geq$5\\%)",    K=5, nn="10", base="spn_penicillin_MIC_coarse_dilutions", run_stub="spn_penicillin_MIC_coarse_dilutions"),
  list(species_dir="spn_penicillin", org="\\textit{S. pneumoniae}", amr="PEN", binning="standard ($\\geq$10\\%)", K=4, nn="11", base="spn_penicillin_MIC_large_minbin",     run_stub="spn_penicillin_MIC_large_minbin"),
  list(species_dir="spn_penicillin", org="\\textit{S. pneumoniae}", amr="PEN", binning="minima ($\\geq$5\\%)",    K=5, nn="16", base="spn_penicillin_MIC_minimabinning",    run_stub="spn_penicillin_MIC_minimabinning"),
  # SPN TMP
  list(species_dir="spn_trimethoprim", org="\\textit{S. pneumoniae}", amr="TMP", binning="standard ($\\geq$5\\%)",  K=5, nn="05", base="spn_trimethoprim_MIC",                  run_stub="spn_trimethoprim_MIC"),
  list(species_dir="spn_trimethoprim", org="\\textit{S. pneumoniae}", amr="TMP", binning="coarse ($\\geq$5\\%)",    K=3, nn="12", base="spn_trimethoprim_MIC_coarse_dilutions", run_stub="spn_trimethoprim_MIC_coarse_dilutions"),
  list(species_dir="spn_trimethoprim", org="\\textit{S. pneumoniae}", amr="TMP", binning="standard ($\\geq$10\\%)", K=3, nn="13", base="spn_trimethoprim_MIC_large_minbin",     run_stub="spn_trimethoprim_MIC_large_minbin"),
  # TB RIF
  list(species_dir="tb_rifampicin", org="\\textit{M. tuberculosis}", amr="RIF", binning="standard ($\\geq$5\\%)",  K=5, nn="08", base="tb_rifampicin_MIC",                  run_stub="tb_rifampicin_MIC"),
  list(species_dir="tb_rifampicin", org="\\textit{M. tuberculosis}", amr="RIF", binning="coarse ($\\geq$5\\%)",    K=4, nn="14", base="tb_rifampicin_MIC_coarse_dilutions", run_stub="tb_rifampicin_MIC_coarse_dilutions"),
  list(species_dir="tb_rifampicin", org="\\textit{M. tuberculosis}", amr="RIF", binning="standard ($\\geq$10\\%)", K=4, nn="15", base="tb_rifampicin_MIC_large_minbin",     run_stub="tb_rifampicin_MIC_large_minbin")
)

# -----------------------------------------------------------------------------
# Row emitters
# -----------------------------------------------------------------------------

# Binary / continuous share the same row shape:
#   <Org> & <AMR> & <Eval> & <ntr> & <nte> & <m1> & <m2> & ... \\
flat_row <- function(org, amr, eval, payload) {
  sprintf("%s & %s & %-7s & %s & %s & %s \\\\",
          org, amr, eval, payload$ntr, payload$nte,
          paste(payload$m, collapse = " & "))
}

build_flat_table <- function(specs, metric_cols, evals = EVALS, pred_subdir = "prediction",
                             midrule_before_tb = TRUE) {
  out <- character()
  for (i in seq_along(specs)) {
    s <- specs[[i]]
    prev <- if (i == 1) NULL else specs[[i - 1]]
    # Row separator before this drug
    if (!is.null(prev)) {
      if (prev$species_dir == s$species_dir) {
        out <- c(out, "\\cmidrule(l){2-13}")
      } else {
        out <- c(out, "\\midrule")
      }
    }
    for (j in seq_along(evals)) {
      eval <- evals[j]
      run_dir <- paste0(s$nn, "_", s$run)
      payload <- collect_one(s$species_dir, s$nn, s$base, run_dir, eval, metric_cols, pred_subdir)
      # Only the first eval row of a drug carries Org + AMR. Within a species
      # the second drug also re-prints "Org" as blank — we mimic the template,
      # which leaves Org blank after the first drug of each species block.
      org <- if (j == 1 && (is.null(prev) || prev$species_dir != s$species_dir)) s$org else ""
      amr <- if (j == 1) s$amr else ""
      out <- c(out, flat_row(org, amr, eval, payload))
    }
  }
  out
}

# Ordinal row shape:
#   <Org> & <AMR> & <Binning> & <K> & <Model> & <Eval> & <ntr> & <nte> & <bacc> & <ppv> & <rps_mean> & <rps_med> & <rpss_unif> & <rpss_freq> \\
#
# Cell-blanking by row index within a binning block (PPOM+POM = 6 rows):
#   Row 1 (PPOM PPC):     org/amr/binning/K/model all set
#   Row 2 (PPOM Random):  all blank, just Eval onward
#   Row 3 (PPOM Lineage): all blank
#   --> \cmidrule(l){5-14}
#   Row 4 (POM PPC):      org/amr/binning/K blank, model=POM
#   Row 5 (POM Random):   all blank
#   Row 6 (POM Lineage):  all blank
ordinal_row <- function(org, amr, binning, K, model, eval, payload) {
  sprintf("%s & %s & %s & %s & %s & %-7s & %s & %s & %s \\\\",
          org, amr, binning, K, model, eval, payload$ntr, payload$nte,
          paste(payload$m, collapse = " & "))
}

build_ordinal_table <- function(specs, metric_cols, evals = EVALS, pred_subdir = "prediction") {
  out <- character()
  for (i in seq_along(specs)) {
    s <- specs[[i]]
    prev <- if (i == 1) NULL else specs[[i - 1]]
    if (!is.null(prev)) {
      if (prev$species_dir != s$species_dir) {
        out <- c(out, "\\midrule")
      } else {
        out <- c(out, "\\cmidrule(l){3-14}")
      }
    }
    new_species_block <- is.null(prev) || prev$species_dir != s$species_dir
    for (model in c("PPOM", "POM")) {
      if (model == "POM") out <- c(out, "\\cmidrule(l){5-14}")
      run_dir <- paste0(s$nn, "_", s$run_stub, "_", model)
      for (j in seq_along(evals)) {
        eval <- evals[j]
        payload <- collect_one(s$species_dir, s$nn, s$base, run_dir, eval, metric_cols, pred_subdir)
        split   <- if (eval == "PPC") NA_character_ else if (eval == "Random") "random" else "loso"
        csv_p   <- if (eval == "PPC") ppc_csv(s$species_dir, run_dir)
                   else                pred_csv(s$species_dir, run_dir, split, pred_subdir)
        payload$m <- star_ordinal_cells(payload$m, csv_p, s$species_dir, run_dir,
                                         model, eval, paste0(s$nn, "_", s$base), split, s$K,
                                         pred_subdir)
        first_in_binning <- (model == "PPOM" && j == 1)
        first_in_model   <- (j == 1)
        org     <- if (first_in_binning && new_species_block) s$org     else ""
        amr     <- if (first_in_binning) s$amr     else ""
        binning <- if (first_in_binning) s$binning else ""
        K_cell  <- if (first_in_binning) as.character(s$K) else ""
        model_cell <- if (first_in_model) model else ""
        out <- c(out, ordinal_row(org, amr, binning, K_cell, model_cell, eval, payload))
      }
    }
  }
  out
}

# -----------------------------------------------------------------------------
# Write outputs
# -----------------------------------------------------------------------------

write_table <- function(lines, path) {
  writeLines(lines, path)
  message(sprintf("wrote %d lines -> %s", length(lines), path))
}

binary_lines     <- build_flat_table(binary_specs,     BINARY_COLS)
continuous_lines <- build_flat_table(continuous_specs, CONT_COLS)
ordinal_lines    <- build_ordinal_table(ordinal_specs, ORD_COLS)

write_table(binary_lines,     file.path(OUT_DIR, "binary_table.tex"))
write_table(continuous_lines, file.path(OUT_DIR, "continuous_table.tex"))
write_table(ordinal_lines,    file.path(OUT_DIR, "ordinal_table.tex"))

# Lasso baselines live under prediction_lasso/ and have prediction-only runs
# (no Bayesian inference, so no PPC eval) — emit the same tables with just the
# Random and Lineage rows.
LASSO_EVALS <- c("Random", "Lineage")

binary_lasso     <- build_flat_table(binary_specs,     BINARY_COLS, LASSO_EVALS, "prediction_lasso")
continuous_lasso <- build_flat_table(continuous_specs, CONT_COLS,   LASSO_EVALS, "prediction_lasso")
ordinal_lasso    <- build_ordinal_table(ordinal_specs, ORD_COLS,    LASSO_EVALS, "prediction_lasso")

write_table(binary_lasso,     file.path(OUT_DIR, "binary_table_lasso.tex"))
write_table(continuous_lasso, file.path(OUT_DIR, "continuous_table_lasso.tex"))
write_table(ordinal_lasso,    file.path(OUT_DIR, "ordinal_table_lasso.tex"))
