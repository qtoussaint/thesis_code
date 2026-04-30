#!/usr/bin/env Rscript
# Extract cutpoint posterior summary for the no-horseshoe flexcutpoints PPOM
# variant. Writes <RUNDIR>/ppc/diagnose/cutpoint_posterior.csv with the same
# schema used by the existing diagnose pipeline, plus drift_median (cutpoint
# minus log2(mic_breakpoints) anchor) and prior_sd.
#
# Usage:
#   Rscript extract_cutpoint_posterior_no_hs.R <model_subdir> <dataset>

suppressPackageStartupMessages({
  library(posterior)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)
model_subdir <- args[1]
dataset      <- args[2]

results_root <- "/nfs/research/jlees/jacqueline/thesis_results/compare_ordinal_models"
RUNDIR <- file.path(results_root, model_subdir, dataset)

rds_candidates <- list.files(file.path(RUNDIR, "fitted_model"),
                             pattern = "\\.RDS$", ignore.case = TRUE,
                             full.names = TRUE)
rds_candidates <- rds_candidates[!grepl("depruned", rds_candidates)]
stopifnot(length(rds_candidates) >= 1)
rds_path <- rds_candidates[which.max(file.mtime(rds_candidates))]

json_candidates <- list.files(RUNDIR, pattern = "\\.json$", full.names = TRUE)
stopifnot(length(json_candidates) >= 1)
json_path <- json_candidates[which.max(file.mtime(json_candidates))]

message("RDS:  ", rds_path)
message("JSON: ", json_path)

dat    <- jsonlite::fromJSON(json_path, simplifyVector = TRUE)
anchor <- log2(as.numeric(dat$mic_breakpoints))
K_1    <- length(anchor)
message("K-1 cutpoints: ", K_1)

fit   <- readRDS(rds_path)
draws <- as_draws_df(tryCatch(fit$draws(), error = function(e) fit))
message("draws: ", nrow(draws))

cp_cols <- sprintf("cutpoints[%d]", seq_len(K_1))
missing <- setdiff(cp_cols, names(draws))
if (length(missing) > 0)
  stop("missing cutpoint columns in draws: ", paste(missing, collapse = ", "))
cp_mat <- as.matrix(draws[, cp_cols, drop = FALSE])

out <- data.frame(
  k            = seq_len(K_1),
  anchor       = anchor,
  median       = apply(cp_mat, 2, median),
  sd           = apply(cp_mat, 2, sd),
  q025         = apply(cp_mat, 2, quantile, 0.025, names = FALSE),
  q975         = apply(cp_mat, 2, quantile, 0.975, names = FALSE),
  drift_median = apply(cp_mat, 2, median) - anchor,
  prior_sd     = 0.5
)

OUTDIR <- file.path(RUNDIR, "ppc", "diagnose")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
out_path <- file.path(OUTDIR, "cutpoint_posterior.csv")
write.csv(out, out_path, row.names = FALSE)
message("wrote: ", out_path)
print(out)
