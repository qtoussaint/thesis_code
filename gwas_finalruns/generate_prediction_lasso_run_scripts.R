#!/usr/bin/env Rscript
# Emits ONLY the prediction_lasso SLURM scripts -- the lasso variant of every
# prediction run. Reuses build_script() and the dataset/model/species config from
# generate_run_scripts.R (sourced; its emit loop is guarded so sourcing only
# defines helpers). Run this to regenerate prediction_lasso/ on its own, or run
# generate_run_scripts.R to regenerate everything (which includes prediction_lasso).

source("/nfs/research/jlees/jacqueline/thesis_code/gwas_finalruns/generate_run_scripts.R")

written <- 0L
for (i in seq_len(nrow(DATASETS))) {
  row     <- DATASETS[i, ]
  species <- row$species
  dataset <- row$dataset
  models  <- BINNING_MODELS[[row$binning]]

  pred_lasso_dir <- file.path(ROOT, "prediction_lasso", species)
  dir.create(pred_lasso_dir, recursive = TRUE, showWarnings = FALSE)

  for (m in models) {
    for (split in c("random", "loso")) {
      fname <- sprintf("%s_%s_%s.sh", dataset, m$model, split)
      path  <- file.path(pred_lasso_dir, fname)
      writeLines(build_script(species, dataset, m$model, m$model_type,
                              analysis_type = "prediction", split = split, lasso = TRUE),
                 path)
      Sys.chmod(path, mode = "0755")
      written <- written + 1L
    }
  }
}

message(sprintf("Wrote %d prediction_lasso SLURM scripts under %s/prediction_lasso/",
                written, ROOT))
