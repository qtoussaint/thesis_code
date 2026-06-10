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
# Resume: a run whose cutpoint_draws.csv already exists is skipped. Runs with no
# RDS yet (still inferring) are skipped with a warning.
#
# Usage (big-mem node):
#   mamba activate gwas_pipeline
#   Rscript paper_figures/extract_cutpoint_draws.R

suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
})

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

# Per-cutpoint MIC breakpoint, ordered by cutpoint index, from the run's
# depruned_variant_effects.csv (unique cutpoint -> cutpoint_MIC).
cutpoint_mic_map <- function(run_dir) {
  effects_csv <- file.path(run_dir, "fitted_model", "depruned_variant_effects.csv")
  if (!file.exists(effects_csv)) return(NULL)
  eff <- read.csv(effects_csv)
  if (!all(c("cutpoint", "cutpoint_MIC") %in% names(eff))) return(NULL)
  map <- unique(eff[, c("cutpoint", "cutpoint_MIC")])
  map <- map[order(map$cutpoint), ]
  map
}

extract_one <- function(run_dir) {
  out_csv <- file.path(run_dir, "plots", "cutpoints", "cutpoint_draws.csv")
  if (file.exists(out_csv)) {
    message("  cache exists, skipping: ", out_csv)
    return(invisible(NULL))
  }
  rds <- list.files(file.path(run_dir, "fitted_model"),
                    pattern = "\\.RDS$", full.names = TRUE)
  if (length(rds) == 0L) {
    warning("no RDS (run not finished?), skipping: ", run_dir)
    return(invisible(NULL))
  }

  message("  reading ", rds[1])
  fit <- readRDS(rds[1])
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

  map <- cutpoint_mic_map(run_dir)
  mic <- if (!is.null(map) && nrow(map) == n_cp) map$cutpoint_MIC else rep(NA_real_, n_cp)

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

for (ds in datasets) {
  for (s in ds$binning_specs) {
    for (model in MODELS) {
      run_dir <- run_dir_path(ds$species_dir, s$nn, s$run_stub, model)
      message(sprintf("[%s_%s_%s]", s$nn, s$run_stub, model))
      tryCatch(extract_one(run_dir),
               error = function(e) warning("failed: ", run_dir, " -- ",
                                           conditionMessage(e)))
    }
  }
}

message("Done.")
