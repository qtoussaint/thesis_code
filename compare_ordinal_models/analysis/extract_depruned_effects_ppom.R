#!/usr/bin/env Rscript
# Rebuild fitted_model/depruned_variant_effects.csv for a PPOM inference run
# whose pipeline stopped before the effects step. The no-horseshoe Stan programs
# emit no y_rep_ppc, so .run_inference_ppc() errors out immediately after the fit
# is saved, and the effects step that follows it never runs. The fit itself is
# fine and sits in fitted_model/*.RDS, so the table can be rebuilt from it.
#
# Uses the pipeline's own helpers (parse_ld_summary, read_correlation_directions,
# deprun_effects_ppom_from_draws) so the output matches a normal run exactly.
#
# Usage:
#   Rscript extract_depruned_effects_ppom.R <model_subdir> <dataset>
#
# Needs the whole draws object in memory -- run it through SLURM with a large
# --mem, not on a login node.

suppressPackageStartupMessages({
  library(posterior)
  library(jsonlite)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)
model_subdir <- args[1]
dataset      <- args[2]

RESULTS_ROOT <- "/nfs/research/jlees/jacqueline/thesis_results/compare_ordinal_models"
PKG_ROOT     <- "/nfs/research/jlees/jacqueline/gwas_workflow/code/gwas_workflow"
RUNDIR       <- file.path(RESULTS_ROOT, model_subdir, dataset)

invisible(lapply(list.files(file.path(PKG_ROOT, "R"), full.names = TRUE), source))

# -----------------------------------------------------------------------------
# Inputs
# -----------------------------------------------------------------------------
rds_candidates <- list.files(file.path(RUNDIR, "fitted_model"), pattern = "\\.RDS$",
                             ignore.case = TRUE, full.names = TRUE)
stopifnot(length(rds_candidates) >= 1)
rds_path <- rds_candidates[which.max(file.mtime(rds_candidates))]

json_candidates <- list.files(RUNDIR, pattern = "\\.json$", full.names = TRUE)
stopifnot(length(json_candidates) >= 1)
json_path <- json_candidates[which.max(file.mtime(json_candidates))]

ld_summary_csv <- file.path(RUNDIR, "cppRATE_matrices", "ld_pruning_summary.csv")
direction_csv  <- file.path(RUNDIR, "cppRATE_matrices", "direction_of_correlation.csv")
stopifnot(file.exists(ld_summary_csv), file.exists(direction_csv))

effects_csv <- file.path(RUNDIR, "fitted_model", "depruned_variant_effects.csv")
if (file.exists(effects_csv))
  stop("effects CSV already exists, refusing to overwrite: ", effects_csv)

message("RDS:  ", rds_path)
message("JSON: ", json_path)

# One row per variant in the unpruned set, so this gives the full variant count
# without touching the multi-GB dataset JSON.
n_total_variants <- nrow(data.table::fread(direction_csv, select = 1L,
                                           showProgress = FALSE))
message("total variants: ", n_total_variants)

dat         <- jsonlite::fromJSON(json_path, simplifyVector = TRUE)
breakpoints <- dat$mic_breakpoints
n_cutpoints <- length(breakpoints)
message("cutpoints: ", n_cutpoints, " (MIC ",
        paste(breakpoints, collapse = ", "), ")")
rm(dat); invisible(gc())

# -----------------------------------------------------------------------------
# Draws
# -----------------------------------------------------------------------------
message("loading fit ...")
fit <- readRDS(rds_path)
posterior_variants <- posterior::as_draws_matrix(fit$draws("beta_variant"))
rm(fit); invisible(gc())
message("draws: ", nrow(posterior_variants), " x ", ncol(posterior_variants))

# Stan stores matrix[V, K-1] column-major, and the de-pruning helper slices the
# columns as n_kept consecutive per cutpoint, so order them (cutpoint, variant)
# explicitly rather than relying on the order they come back in.
idx <- regmatches(colnames(posterior_variants),
                  regexec("^beta_variant\\[(\\d+),(\\d+)\\]$",
                          colnames(posterior_variants)))
stopifnot(all(lengths(idx) == 3L))
v <- as.integer(vapply(idx, `[`, "", 2L))
k <- as.integer(vapply(idx, `[`, "", 3L))
stopifnot(max(k) == n_cutpoints)
posterior_variants <- posterior_variants[, order(k, v), drop = FALSE]

n_kept <- ncol(posterior_variants) / n_cutpoints
message("kept variants: ", n_kept)

# -----------------------------------------------------------------------------
# De-prune back to the full variant set
# -----------------------------------------------------------------------------
final_rep_variants <- parse_ld_summary(ld_summary_csv)
flip_sign          <- read_correlation_directions(RUNDIR, n_total_variants)

message("de-pruning ...")
depruned <- deprun_effects_ppom_from_draws(
  posterior_variants = posterior_variants,
  final_rep_variants = final_rep_variants,
  n_total_variants   = n_total_variants,
  n_cutpoints        = n_cutpoints,
  flip_sign          = flip_sign
)

depruned$cutpoint_MIC <- factor(depruned$cutpoint)
levels(depruned$cutpoint_MIC) <- breakpoints

data.table::fwrite(depruned, effects_csv, sep = ",", col.names = TRUE,
                   row.names = FALSE, compress = "none")
message("wrote: ", effects_csv, " (", nrow(depruned), " rows)")

print(summary(depruned$median))
print(table(depruned$signif, useNA = "ifany"))
