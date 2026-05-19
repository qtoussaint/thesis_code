#!/usr/bin/env Rscript
# Emits SLURM scripts for every (dataset x model x inference|prediction x split)
# combination defined below. Run once after editing this file; re-run to regenerate.

suppressPackageStartupMessages({
  library(glue)
})

# ----------------------------- static paths -------------------------------- #
ROOT             <- "/nfs/research/jlees/jacqueline/thesis_code/gwas_finalruns"
MODELS_DIR       <- "/nfs/research/jlees/jacqueline/thesis_code/gwas_finalmodels"
DATASETS_DIR     <- "/nfs/research/jlees/jacqueline/thesis_results/gwas_datasets"
RESULTS_BASE     <- "/nfs/research/jlees/jacqueline/thesis_results"
GENES_DIR        <- "/nfs/research/jlees/jacqueline/thesis_code/gwas_genesofinterest"
PIPELINE_RSCRIPT <- "/nfs/research/jlees/jacqueline/gwas_workflow/code/gwas_workflow/inst/scripts/run_pipeline.R"
PRUNING_BIN      <- "/hps/software/users/jlees/jacqueline/manual_installs/bin/BacPrune-Rust/"
CPPRATE_BIN      <- "/hps/software/users/jlees/jacqueline/manual_installs/bin/cpprate-0.2.0/build/bin/cpprate"

SPN_ANNOT <- "/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/genotype/fields_filtered_maf05_multiallelic.txt"
TB_ANNOT  <- "/nfs/research/jlees/jacqueline/gwas_data/tuberculosis/cryptic_regeno_snpeff/cryptic_regeno_fields_filtered.txt"

SPECIES_META <- list(
  spn_penicillin = list(
    results_dir = file.path(RESULTS_BASE, "gwas_spn_penicillin"),
    annotations = SPN_ANNOT,
    genes       = file.path(GENES_DIR, "spn_penicillin_genesofinterest.txt"),
    job_tag     = "spnpen"
  ),
  spn_trimethoprim = list(
    results_dir = file.path(RESULTS_BASE, "gwas_spn_trimethoprim"),
    annotations = SPN_ANNOT,
    genes       = file.path(GENES_DIR, "spn_trimethoprim_genesofinterest.txt"),
    job_tag     = "spntmp"
  ),
  tb_rifampicin = list(
    results_dir = file.path(RESULTS_BASE, "gwas_tb_rifampicin"),
    annotations = TB_ANNOT,
    genes       = file.path(GENES_DIR, "tb_rifampicin_genesofinterest.txt"),
    job_tag     = "tbrif"
  )
)

# ----------------------------- dataset table ------------------------------- #
# Each row: dataset on disk + which species/binning it belongs to.
# binning -> models is expanded below.
DATASETS <- data.frame(
  number      = c("01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16"),
  dataset     = c(
    "01_spn_penicillin_binary",
    "02_spn_penicillin_MIC",
    "03_spn_penicillin_continuous",
    "04_spn_trimethoprim_binary",
    "05_spn_trimethoprim_MIC",
    "06_spn_trimethoprim_continuous",
    "07_tb_rifampicin_binary",
    "08_tb_rifampicin_MIC",
    "09_tb_rifampicin_continuous",
    "10_spn_penicillin_MIC_coarse_dilutions",
    "11_spn_penicillin_MIC_large_minbin",
    "12_spn_trimethoprim_MIC_coarse_dilutions",
    "13_spn_trimethoprim_MIC_large_minbin",
    "14_tb_rifampicin_MIC_coarse_dilutions",
    "15_tb_rifampicin_MIC_large_minbin",
    "16_spn_penicillin_MIC_minimabinning"
  ),
  species     = c(
    "spn_penicillin","spn_penicillin","spn_penicillin",
    "spn_trimethoprim","spn_trimethoprim","spn_trimethoprim",
    "tb_rifampicin","tb_rifampicin","tb_rifampicin",
    "spn_penicillin","spn_penicillin",
    "spn_trimethoprim","spn_trimethoprim",
    "tb_rifampicin","tb_rifampicin",
    "spn_penicillin"
  ),
  binning     = c(
    "binary","ordinal","continuous",
    "binary","ordinal","continuous",
    "binary","ordinal","continuous",
    "ordinal","ordinal",
    "ordinal","ordinal",
    "ordinal","ordinal",
    "ordinal"
  ),
  stringsAsFactors = FALSE
)

# ----------------------------- model dispatch ------------------------------ #
# Per-binning list of (model stem in stan filename, model_type flag).
BINNING_MODELS <- list(
  binary     = list(list(model = "logistic",   model_type = "binary")),
  continuous = list(list(model = "continuous", model_type = "continuous")),
  ordinal    = list(
    list(model = "POM",  model_type = "pom"),
    list(model = "PPOM", model_type = "ppom")
  )
)

# ----------------------------- resources ----------------------------------- #
# (species_kind, model_kind) -> cpus/mem/time.
# species_kind: spn or tb. model_kind: light (logistic/continuous) or heavy (POM/PPOM).
RESOURCES <- list(
  spn_light = list(cpus = 32, mem = 200, time = 6),
  spn_heavy = list(cpus = 48, mem = 650, time = 12),
  tb_light  = list(cpus = 48, mem = 400, time = 12),
  tb_heavy  = list(cpus = 48, mem = 800, time = 24)
)

species_kind <- function(species) if (startsWith(species, "tb_")) "tb" else "spn"
model_kind   <- function(model)   if (model %in% c("POM", "PPOM")) "heavy" else "light"

resources_for <- function(species, model) {
  key <- paste(species_kind(species), model_kind(model), sep = "_")
  RESOURCES[[key]]
}

# ----------------------------- script template ----------------------------- #
build_script <- function(species, dataset, model, model_type, analysis_type, split) {
  meta <- SPECIES_META[[species]]
  res  <- resources_for(species, model)

  # Prediction LOSO uses the _loso dataset variant; random uses the base dataset.
  json_dataset <- if (analysis_type == "prediction" && split == "loso")
    paste0(dataset, "_loso") else dataset

  slug_suffix <- if (analysis_type == "prediction") {
    paste0("_", model, "_", split)   # split is "random" or "loso"
  } else {
    paste0("_", model)
  }
  slug <- paste0(dataset, slug_suffix)

  analysis_subdir <- analysis_type  # "inference" or "prediction"
  results_dir     <- file.path(meta$results_dir, analysis_subdir, slug)
  json_dir        <- file.path(DATASETS_DIR, analysis_subdir, json_dataset)
  json_file       <- file.path(json_dir, paste0(json_dataset, ".json"))
  variant_index   <- file.path(json_dir, paste0(json_dataset, "_variant_index.csv"))

  stan_model <- file.path(MODELS_DIR, paste0(model, "_", analysis_type, ".stan"))

  job_name <- glue("{meta$job_tag}_{substr(dataset,1,2)}_{model}_{substr(analysis_type,1,4)}",
                   .trim = FALSE)
  if (analysis_type == "prediction") {
    job_name <- paste0(job_name, "_", split)
  }

  header <- glue("#!/usr/bin/env bash
#SBATCH --job-name={job_name}
#SBATCH --nodes=1
#SBATCH --cpus-per-task={res$cpus}
#SBATCH --mem={res$mem}G
#SBATCH --time={sprintf('%02d', res$time)}:00:00
#SBATCH --error={results_dir}/logs/{slug}_%j.err
#SBATCH --output={results_dir}/logs/{slug}_%j.out

#################################################################################

source ~/.bashrc
mamba activate gwas_pipeline

mkdir -p {results_dir}/logs

RSCRIPT_PATH=\"{PIPELINE_RSCRIPT}\"

DATA=\"--data {json_file}\"
STAN_MODEL=\"--stan_model {stan_model}\"
ANALYSIS_TYPE=\"--analysis_type {analysis_type}\"
ANALYSIS_NICKNAME=\"--analysis_nickname {slug}\"
OUTPUT_DIR=\"--output_directory {results_dir}\"
THREADS=\"--threads {res$cpus}\"

LD_PRUNING=\"--ld_pruning true\"
PRUNING_SOFTWARE=\"--pruning_software {PRUNING_BIN}\"
MAF_CUTOFF=\"--maf_cutoff 0\"
LD_THRESHOLD=\"--ld_threshold 1\"

PHANDANGO=\"--phandango {variant_index}\"
ANNOTATIONS=\"--annotations {meta$annotations}\"
MODEL_TYPE=\"--model_type {model_type}\"
GENES_OF_INTEREST=\"--genes_of_interest {meta$genes}\"
RESUME=\"--resume\"
CPPRATE=\"--cpprate_bin {CPPRATE_BIN}\"
")

  # Build the Rscript invocation outside glue so the backslash-newline
  # continuations survive (glue eats `\\\n` as line continuation).
  args <- c("$DATA", "$STAN_MODEL", "$ANALYSIS_TYPE", "$ANALYSIS_NICKNAME",
            "$OUTPUT_DIR", "$THREADS", "$LD_PRUNING", "$PRUNING_SOFTWARE",
            "$MAF_CUTOFF", "$LD_THRESHOLD", "$PHANDANGO", "$ANNOTATIONS",
            "$MODEL_TYPE", "$GENES_OF_INTEREST", "$RESUME", "$CPPRATE")
  bs <- "\\"
  invocation <- c(
    paste("Rscript $RSCRIPT_PATH", bs),
    paste(args[-length(args)], bs),
    args[length(args)]
  )

  paste(c(header, invocation, ""), collapse = "\n")
}

# ----------------------------- emit loop ----------------------------------- #
written <- 0L
for (i in seq_len(nrow(DATASETS))) {
  row     <- DATASETS[i, ]
  species <- row$species
  dataset <- row$dataset
  binning <- row$binning
  models  <- BINNING_MODELS[[binning]]

  for (m in models) {
    # Inference (single split)
    inf_dir <- file.path(ROOT, "inference", species)
    dir.create(inf_dir, recursive = TRUE, showWarnings = FALSE)
    inf_path <- file.path(inf_dir, sprintf("%s_%s.sh", dataset, m$model))
    writeLines(build_script(species, dataset, m$model, m$model_type,
                            analysis_type = "inference", split = "base"),
               inf_path)
    Sys.chmod(inf_path, mode = "0755")
    written <- written + 1L

    # Prediction: base + LOSO
    pred_dir <- file.path(ROOT, "prediction", species)
    dir.create(pred_dir, recursive = TRUE, showWarnings = FALSE)

    for (split in c("random", "loso")) {
      fname <- sprintf("%s_%s_%s.sh", dataset, m$model, split)
      path  <- file.path(pred_dir, fname)
      writeLines(build_script(species, dataset, m$model, m$model_type,
                              analysis_type = "prediction", split = split),
                 path)
      Sys.chmod(path, mode = "0755")
      written <- written + 1L
    }
  }
}

message(sprintf("Wrote %d SLURM scripts under %s/", written, ROOT))
