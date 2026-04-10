#!/usr/bin/env bash

#SBATCH --job-name=paperPOM_02_spn_pen_MIC
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=90G
#SBATCH --time=5:00:00
#SBATCH --error=/nfs/research/jlees/jacqueline/thesis_results/compare_ordinal_models/paper_model_POM/02_spn_penicillin_MIC/logs/02_spn_penicillin_MIC.err
#SBATCH --output=/nfs/research/jlees/jacqueline/thesis_results/compare_ordinal_models/paper_model_POM/02_spn_penicillin_MIC/logs/02_spn_penicillin_MIC.out

#################################################################################

#conda activate gwas_pipeline

RSCRIPT_PATH="/nfs/research/jlees/jacqueline/gwas_workflow/code/gwas_workflow/inst/scripts/run_pipeline.R"

# ---------------------------------------------------------------------------
# Required arguments
# ---------------------------------------------------------------------------

DATA="--data /nfs/research/jlees/jacqueline/thesis_results/gwas_datasets/inference/02_spn_penicillin_MIC/02_spn_penicillin_MIC.json"
STAN_MODEL="--stan_model /nfs/research/jlees/jacqueline/thesis_code/compare_ordinal_models/POM_models/paper_model_POM.stan"
ANALYSIS_TYPE="--analysis_type inference"
ANALYSIS_NICKNAME="--analysis_nickname 02_spn_penicillin_MIC_paperPOM"
OUTPUT_DIR="--output_directory /nfs/research/jlees/jacqueline/thesis_results/compare_ordinal_models/paper_model_POM/02_spn_penicillin_MIC"
THREADS="--threads 24"

# ---------------------------------------------------------------------------
# LD pruning arguments
# ---------------------------------------------------------------------------

LD_PRUNING="--ld_pruning true"
PRUNING_SOFTWARE="--pruning_software /hps/software/users/jlees/jacqueline/manual_installs/bin/BacPrune-Rust/"
MAF_CUTOFF="--maf_cutoff 0"
LD_THRESHOLD="--ld_threshold 1"

# ---------------------------------------------------------------------------
# Optional arguments
# ---------------------------------------------------------------------------

PHANDANGO="--phandango /nfs/research/jlees/jacqueline/thesis_results/gwas_datasets/inference/02_spn_penicillin_MIC/02_spn_penicillin_MIC_variant_index.csv"
ANNOTATIONS="--annotations /nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/genotype/fields_filtered_maf05_multiallelic.txt"
MODEL_TYPE="--model_type pom"
NORATE="--norate"

# ---------------------------------------------------------------------------

mkdir -p /nfs/research/jlees/jacqueline/thesis_results/compare_ordinal_models/paper_model_POM/02_spn_penicillin_MIC/logs

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
$NORATE \
