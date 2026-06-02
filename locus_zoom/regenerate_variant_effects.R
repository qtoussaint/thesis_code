#!/usr/bin/env Rscript
# regenerate_variant_effects.R
#
# Regenerate fitted_model/depruned_variant_effects.csv from the fitted-model
# .RDS for non-PPOM inference runs (binary, POM, logistic, continuous).
#
# Used when the inference pipeline finished the Stan fit and wrote the .RDS
# but did not also emit the depruned variant-effects CSV (e.g. the
# 08_tb_rifampicin_MIC_POM run). The CSV is required by locus-zoom plots
# for the abs_median / exp_abs_median metrics.
#
# Output schema (matches the binary CSV produced by the modern pipeline):
#   variant_id, median, signif
#
# Usage:
#   Rscript regenerate_variant_effects.R --run_dir <inference run dir>

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(cmdstanr)
})

# Helpers (parse_ld_summary, read_correlation_directions, read_pruned_variants,
# deprun_effects, deprun_effects_from_draws) live in the gwas_workflow R
# sources. Source them directly rather than installing the package.
GWAS_PIPELINE_R <- "/nfs/research/jlees/jacqueline/gwas_workflow/code/gwas_workflow/R"
source(file.path(GWAS_PIPELINE_R, "ld_pruning.R"))
source(file.path(GWAS_PIPELINE_R, "cpprate.R"))

option_list <- list(
  make_option("--run_dir",
    type = "character", default = NULL,
    help = "Inference run dir; must contain fitted_model/*.RDS and cppRATE_matrices/"
  ),
  make_option("--n_total_variants",
    type = "integer", default = NULL,
    help = "Total variants V. If omitted, inferred from the pruning_input.csv header row."
  ),
  make_option("--force",
    action = "store_true", default = FALSE,
    help = "Overwrite existing CSV instead of skipping."
  )
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$run_dir)) stop("--run_dir is required")
if (!dir.exists(opt$run_dir)) stop("Run dir not found: ", opt$run_dir)

run_dir      <- normalizePath(opt$run_dir)
fitted_dir   <- file.path(run_dir, "fitted_model")
matrices_dir <- file.path(run_dir, "cppRATE_matrices")
out_csv      <- file.path(fitted_dir, "depruned_variant_effects.csv")

if (file.exists(out_csv) && !opt$force) {
  message("Output already exists, skipping (use --force to overwrite): ", out_csv)
  quit(save = "no", status = 0)
}

rds_files <- list.files(fitted_dir, pattern = "\\.RDS$", full.names = TRUE)
if (length(rds_files) != 1) {
  stop("Expected exactly one .RDS in ", fitted_dir, "; found ", length(rds_files))
}
rds_path <- rds_files[1]

if (is.null(opt$n_total_variants)) {
  pruning_input <- file.path(matrices_dir, "pruning_input.csv")
  if (!file.exists(pruning_input)) {
    stop("Cannot infer n_total_variants: ", pruning_input, " is missing. Pass --n_total_variants.")
  }
  header <- readLines(pruning_input, n = 1)
  opt$n_total_variants <- length(strsplit(header, ",", fixed = TRUE)[[1]])
}
message("V (total variants): ", opt$n_total_variants)

message("Parsing LD pruning summary and correlation directions...")
final_rep_variants <- parse_ld_summary(file.path(matrices_dir, "ld_pruning_summary.csv"))
flip_sign          <- read_correlation_directions(run_dir, opt$n_total_variants)
kept_variants      <- read_pruned_variants(run_dir)
message("  kept variants: ", length(kept_variants),
        "  |  representatives with pruned neighbours: ", nrow(final_rep_variants))

message("Loading fitted model: ", rds_path)
fitted_model <- readRDS(rds_path)

message("Computing per-variant posterior medians and 89% CI...")
effects <- deprun_effects(
  fitted_model       = fitted_model,
  final_rep_variants = final_rep_variants,
  kept_variants      = kept_variants,
  n_total_variants   = opt$n_total_variants,
  flip_sign          = flip_sign
)

message("Writing ", out_csv)
fwrite(
  data.table(
    variant_id = seq_along(effects$medians),
    median     = effects$medians,
    signif     = as.logical(effects$signif)
  ),
  file      = out_csv,
  sep       = ",",
  col.names = TRUE,
  row.names = FALSE
)

message("Done. ", length(effects$medians), " variant rows written.")
