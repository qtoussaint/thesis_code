#!/usr/bin/env Rscript
# Generate SLURM run scripts (one per model x dataset) that fit each parameter-
# sweep model on a TRAINING split and score it on a held-out TEST split
# (--analysis_type prediction), plus a submit_all.sh launcher.
#
# Mirrors the existing hand-written prediction run script
#   run_PPOM_models/..._tau5_prediction/run_02_spn_penicillin_MIC.sh
# but points --stan_model at the generated prediction models in ./stan_models and
# writes outputs under thesis_results/.../parameter_sweeps_prediction/runs/.
#
#   32 PPOM models    x {02 standard K=8, 16 minima K=5}  = 64 runs
#    8 logistic models x {01 binary}                       =  8 runs
#    3 continuous models x {03 continuous}                 =  3 runs   (75 total)
#
# Run with:
#   mamba activate gwas_pipeline
#   Rscript generate_run_scripts.R
# Then launch on the cluster:
#   bash run_scripts/submit_all.sh

CODE_DIR     <- "/nfs/research/jlees/jacqueline/thesis_code/compare_ordinal_models/parameter_sweeps_prediction"
STAN_DIR     <- file.path(CODE_DIR, "stan_models")
RUN_DIR      <- file.path(CODE_DIR, "run_scripts")
RESULTS_BASE <- "/nfs/research/jlees/jacqueline/thesis_results/compare_ordinal_models/parameter_sweeps_prediction/runs"
DATASETS_DIR <- "/nfs/research/jlees/jacqueline/thesis_results/gwas_datasets/prediction"

RSCRIPT_PATH <- "/nfs/research/jlees/jacqueline/gwas_workflow/code/gwas_workflow/inst/scripts/run_pipeline.R"
PRUNE_SW     <- "/hps/software/users/jlees/jacqueline/manual_installs/bin/BacPrune-Rust/"
ANNOTATIONS  <- "/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/genotype/fields_filtered_maf05_multiallelic.txt"
GENES        <- "/nfs/research/jlees/jacqueline/thesis_code/gwas_genesofinterest/spn_penicillin_genesofinterest.txt"
PFX          <- "final_ordered_categorical_PPOM_free_cutpoints_wide_drift"

dir.create(RUN_DIR, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# Build the model x dataset registry from the generated stan models.
# ---------------------------------------------------------------------------
stans <- list.files(STAN_DIR, pattern = "\\.stan$")

registry <- list()
add <- function(run_name, model_type, datasets) {
  for (ds in datasets)
    registry[[length(registry) + 1]] <<-
      list(run = run_name, model_type = model_type, dataset = ds)
}
for (s in stans) {
  run <- sub("\\.stan$", "", s)
  if (startsWith(run, PFX))                 add(run, "ppom",       c("02_spn_penicillin_MIC", "16_spn_penicillin_MIC_minimabinning"))
  else if (startsWith(run, "logistic_"))    add(run, "binary",     "01_spn_penicillin_binary")
  else if (startsWith(run, "continuous_"))  add(run, "continuous", "03_spn_penicillin_continuous")
  else stop("unrecognised stan model: ", s)
}

# ---------------------------------------------------------------------------
# Script template (mirrors the hand-written tau5_prediction run script).
# ---------------------------------------------------------------------------
short <- function(run) {
  r <- sub(paste0(PFX, "_"), "", run)   # drop the long PPOM prefix
  r <- sub("_prediction$", "", r)
  substr(r, 1, 40)
}

build_script <- function(run, model_type, dataset) {
  out_dir  <- file.path(RESULTS_BASE, run, dataset)
  ds_dir   <- file.path(DATASETS_DIR, dataset)
  json     <- file.path(ds_dir, paste0(dataset, ".json"))
  test_phe <- file.path(ds_dir, paste0(dataset, "_test_phenotypes.csv"))
  vindex   <- file.path(ds_dir, paste0(dataset, "_variant_index.csv"))
  stan     <- file.path(STAN_DIR, paste0(run, ".stan"))
  nn       <- substr(dataset, 1, 2)
  jobname  <- paste0(nn, "_", short(run), "_pred")

  sprintf(r"---(#!/usr/bin/env bash

#SBATCH --job-name=%s
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=600G
#SBATCH --time=06:00:00
#SBATCH --error=%s/logs/%s.err
#SBATCH --output=%s/logs/%s.out

#################################################################################

source /hps/software/users/jlees/jacqueline/etc/profile.d/conda.sh
conda activate gwas_pipeline

RSCRIPT_PATH="%s"

DATA="--data %s"
STAN_MODEL="--stan_model %s"
ANALYSIS_TYPE="--analysis_type prediction"
ANALYSIS_NICKNAME="--analysis_nickname %s"
OUTPUT_DIR="--output_directory %s"
THREADS="--threads 48"

LD_PRUNING="--ld_pruning true"
PRUNING_SOFTWARE="--pruning_software %s"
MAF_CUTOFF="--maf_cutoff 0"
LD_THRESHOLD="--ld_threshold 1"

PHANDANGO="--phandango %s"
ANNOTATIONS="--annotations %s"
MODEL_TYPE="--model_type %s"
TRUE_PHENOTYPES="--true_phenotypes %s"
GENES_OF_INTEREST="--genes_of_interest %s"
NORATE="--norate"
RESUME="--resume"

mkdir -p %s/logs

Rscript $RSCRIPT_PATH \
$DATA \
$STAN_MODEL \
$ANALYSIS_TYPE \
$ANALYSIS_NICKNAME \
$OUTPUT_DIR \
$THREADS \
$LD_PRUNING \
$PRUNING_SOFTWARE \
$MAF_CUTOFF \
$LD_THRESHOLD \
$PHANDANGO \
$ANNOTATIONS \
$MODEL_TYPE \
$TRUE_PHENOTYPES \
$GENES_OF_INTEREST \
$NORATE \
$RESUME \
)---",
    jobname,
    out_dir, dataset, out_dir, dataset,
    RSCRIPT_PATH,
    json, stan,
    paste0(nn, "_", short(run)),
    out_dir,
    PRUNE_SW,
    vindex, ANNOTATIONS, model_type, test_phe, GENES,
    out_dir)
}

# ---------------------------------------------------------------------------
# Emit run scripts + submit_all.sh.
# ---------------------------------------------------------------------------
paths <- character(0)
for (r in registry) {
  d <- file.path(RUN_DIR, r$run)
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  p <- file.path(d, paste0("run_", r$dataset, ".sh"))
  writeLines(build_script(r$run, r$model_type, r$dataset), p)
  paths <- c(paths, p)
}

submit <- c("#!/usr/bin/env bash",
            "# Submit every parameter-sweep prediction run. Generated by generate_run_scripts.R.",
            "set -euo pipefail",
            paste0("sbatch ", shQuote(paths)))
writeLines(submit, file.path(RUN_DIR, "submit_all.sh"))

message(sprintf("Wrote %d run scripts under %s", length(paths), RUN_DIR))
message(sprintf("  PPOM:       %d", sum(vapply(registry, function(x) x$model_type == "ppom", logical(1)))))
message(sprintf("  logistic:   %d", sum(vapply(registry, function(x) x$model_type == "binary", logical(1)))))
message(sprintf("  continuous: %d", sum(vapply(registry, function(x) x$model_type == "continuous", logical(1)))))
message("Launch with: bash ", file.path(RUN_DIR, "submit_all.sh"))
