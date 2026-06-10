#!/usr/bin/env Rscript
# Extract the cutpoint posterior draws from each ordinal inference run's fitted
# RDS and cache them as a small per-run CSV.
#
# The fitted-model RDS files are 13-117 GB each: cutpoints is a tiny ordered[K-1]
# variable buried among the variant betas, and the pipeline saves no draws CSV
# (plots/cutpoints/ holds only PNGs). So we load each RDS once on a big-mem node,
# pull only the `cutpoints` variable, and write
# plots/cutpoints/cutpoint_draws.csv. cutpoints_histogram_faceted_summary.R then
# builds the combined figure from these caches without ever touching the RDS.
#
# Per-cutpoint MIC breakpoint labels come from depruned_variant_effects.csv in
# the same run dir (guaranteed present whenever the RDS exists).
#
# One run per RDS, so this is run as a SLURM array (one task per run, memory
# sized to that run): extract_cutpoint_draws_array.sh maps SLURM_ARRAY_TASK_ID
# to --index here. Resume: a run whose cutpoint_draws.csv already exists is
# skipped; a run with no RDS yet (still inferring) is skipped with a warning.
#
# Modes:
#   Rscript extract_cutpoint_draws.R              # process every run, serially
#   Rscript extract_cutpoint_draws.R --index N    # process only run N (1-based)
#   Rscript extract_cutpoint_draws.R --list       # print "index<TAB>run_dir<TAB>rds_bytes"
#
# Usage (big-mem node):
#   mamba activate gwas_pipeline
#   Rscript paper_figures/extract_cutpoint_draws.R --index 3

RESULTS_ROOT <- "/nfs/research/jlees/jacqueline/thesis_results"

MODELS <- c("POM", "PPOM")

# Same per-dataset binning structure as cutpoints_histogram_summary.R; unitigs
# runs (17/18) are excluded, matching the existing summary figures.
datasets <- list(
  list(species_dir = "spn_penicillin", binning_specs = list(
      list(nn = "02", run_stub = "spn_penicillin_MIC"),
      list(nn = "10", run_stub = "spn_penicillin_MIC_coarse_dilutions"),
      list(nn = "11", run_stub = "spn_penicillin_MIC_large_minbin"),
      list(nn = "16", run_stub = "spn_penicillin_MIC_minimabinning"))),
  list(species_dir = "spn_trimethoprim", binning_specs = list(
      list(nn = "05", run_stub = "spn_trimethoprim_MIC"),
      list(nn = "12", run_stub = "spn_trimethoprim_MIC_coarse_dilutions"),
      list(nn = "13", run_stub = "spn_trimethoprim_MIC_large_minbin"))),
  list(species_dir = "tb_rifampicin", binning_specs = list(
      list(nn = "08", run_stub = "tb_rifampicin_MIC"),
      list(nn = "14", run_stub = "tb_rifampicin_MIC_coarse_dilutions"),
      list(nn = "15", run_stub = "tb_rifampicin_MIC_large_minbin")))
)

run_dir_path <- function(species_dir, nn, run_stub, model) file.path(
  RESULTS_ROOT, paste0("gwas_", species_dir), "inference",
  paste0(nn, "_", run_stub, "_", model))

# Flatten datasets x models into a stable, index-addressable run manifest. The
# ordering is fixed so SLURM_ARRAY_TASK_ID maps to the same run every time.
build_runs <- function() {
  runs <- list()
  for (ds in datasets) {
    for (s in ds$binning_specs) {
      for (model in MODELS) {
        runs[[length(runs) + 1]] <- list(
          tag = paste(s$nn, s$run_stub, model, sep = "_"),
          nn = s$nn, run_stub = s$run_stub,
          run_dir = run_dir_path(ds$species_dir, s$nn, s$run_stub, model))
      }
    }
  }
  runs
}

rds_path <- function(run_dir) {
  f <- list.files(file.path(run_dir, "fitted_model"),
                  pattern = "\\.RDS$", full.names = TRUE)
  if (length(f)) f[1] else NA_character_
}

# Per-cutpoint MIC breakpoints. Prefer the dataset JSON's mic_breakpoints (the
# canonical source, read from the file tail so the multi-GB genotype is never
# loaded); fall back to the run's depruned_variant_effects.csv (PPOM has
# cutpoint_MIC; regenerated POM effects do not). Returns a numeric vector in
# cutpoint order, or NULL.
cutpoint_mic_vec <- function(run_dir, nn, run_stub) {
  json <- file.path(RESULTS_ROOT, "gwas_datasets", "inference",
                    paste0(nn, "_", run_stub), paste0(nn, "_", run_stub, ".json"))
  if (file.exists(json)) {
    sz  <- file.info(json)$size
    n   <- min(sz, 16384)
    con <- file(json, "rb")
    if (sz > n) seek(con, where = sz - n)
    raw <- readChar(con, n, useBytes = TRUE)
    close(con)
    m <- regmatches(raw, regexpr("\"mic_breakpoints\"\\s*:\\s*\\[[^]]*\\]", raw, perl = TRUE))
    if (length(m)) {
      inside <- sub("\\].*$", "", sub("^[^[]*\\[", "", m))
      return(as.numeric(strsplit(inside, "\\s*,\\s*")[[1]]))
    }
  }
  effects_csv <- file.path(run_dir, "fitted_model", "depruned_variant_effects.csv")
  if (!file.exists(effects_csv)) return(NULL)
  eff <- read.csv(effects_csv)
  if (!all(c("cutpoint", "cutpoint_MIC") %in% names(eff))) return(NULL)
  map <- unique(eff[, c("cutpoint", "cutpoint_MIC")])
  map[order(map$cutpoint), ]$cutpoint_MIC
}

extract_one <- function(run_dir, nn, run_stub) {
  out_csv <- file.path(run_dir, "plots", "cutpoints", "cutpoint_draws.csv")
  if (file.exists(out_csv)) {
    message("  cache exists, skipping: ", out_csv)
    return(invisible(NULL))
  }
  rds <- rds_path(run_dir)
  if (is.na(rds)) {
    warning("no RDS (run not finished?), skipping: ", run_dir)
    return(invisible(NULL))
  }

  suppressPackageStartupMessages({
    library(cmdstanr)
    library(posterior)
  })

  message("  reading ", rds)
  fit <- readRDS(rds)
  draws_mat <- tryCatch(
    posterior::as_draws_matrix(fit$draws(variables = "cutpoints")),
    error = function(e) {
      warning("cutpoints not in fit, skipping: ", run_dir, " -- ",
              conditionMessage(e))
      NULL
    })
  rm(fit); gc()
  if (is.null(draws_mat)) return(invisible(NULL))

  n_cp    <- ncol(draws_mat)
  n_draws <- nrow(draws_mat)

  mic <- cutpoint_mic_vec(run_dir, nn, run_stub)
  if (is.null(mic) || length(mic) != n_cp) mic <- rep(NA_real_, n_cp)

  # Long format: each cutpoint's draws stacked in cutpoint order (column-major
  # flatten of the draws matrix matches rep(seq_len(n_cp), each = n_draws)).
  out <- data.frame(
    cutpoint     = rep(seq_len(n_cp), each = n_draws),
    cutpoint_mic = rep(mic, each = n_draws),
    value        = as.numeric(draws_mat),
    stringsAsFactors = FALSE
  )
  dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
  write.csv(out, out_csv, row.names = FALSE)
  message("  wrote ", out_csv, " (", n_cp, " cutpoints x ", n_draws, " draws)")
  rm(draws_mat); gc()
  invisible(NULL)
}

# ---- CLI -------------------------------------------------------------------
argv  <- commandArgs(trailingOnly = TRUE)
runs  <- build_runs()

if ("--list" %in% argv) {
  for (i in seq_along(runs)) {
    rds <- rds_path(runs[[i]]$run_dir)
    bytes <- if (is.na(rds)) 0 else file.info(rds)$size
    cat(i, runs[[i]]$run_dir, format(bytes, scientific = FALSE), sep = "\t")
    cat("\n")
  }
  quit(save = "no")
}

idx_pos <- match("--index", argv)
if (!is.na(idx_pos)) {
  i <- as.integer(argv[idx_pos + 1])
  if (is.na(i) || i < 1 || i > length(runs)) stop("--index out of range: ", argv[idx_pos + 1])
  message(sprintf("[%d/%d %s]", i, length(runs), runs[[i]]$tag))
  extract_one(runs[[i]]$run_dir, runs[[i]]$nn, runs[[i]]$run_stub)
  message("Done.")
  quit(save = "no")
}

# No --index: process all runs serially (fallback / non-array use).
for (i in seq_along(runs)) {
  message(sprintf("[%d/%d %s]", i, length(runs), runs[[i]]$tag))
  tryCatch(extract_one(runs[[i]]$run_dir, runs[[i]]$nn, runs[[i]]$run_stub),
           error = function(e) warning("failed: ", runs[[i]]$run_dir, " -- ",
                                       conditionMessage(e)))
}
message("Done.")
