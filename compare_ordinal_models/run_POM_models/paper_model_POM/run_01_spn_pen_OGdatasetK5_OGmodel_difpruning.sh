#!/usr/bin/env bash

#SBATCH --job-name=paperPOM_OGdataK5_OGmodel_difprun
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=50G
#SBATCH --time=5:00:00
#SBATCH --error=/nfs/research/jlees/jacqueline/thesis_results/compare_ordinal_models/paper_model_POM/spn_pen_OGdatasetK5_OGmodel_difpruning/logs/spn_pen_OGdatasetK5_OGmodel_difpruning.err
#SBATCH --output=/nfs/research/jlees/jacqueline/thesis_results/compare_ordinal_models/paper_model_POM/spn_pen_OGdatasetK5_OGmodel_difpruning/logs/spn_pen_OGdatasetK5_OGmodel_difpruning.out

#################################################################################

#conda activate gwas_pipeline

RSCRIPT_PATH="/nfs/research/jlees/jacqueline/gwas_workflow/code/gwas_workflow/inst/scripts/run_pipeline.R"

# ---------------------------------------------------------------------------
# Required arguments
# ---------------------------------------------------------------------------

DATA="--data /nfs/research/jlees/jacqueline/bayesian_gwas_paper/00_data_v1/inference/01_spn_penicillin_subclusters_K5_3f7c4f1ce71de0.json"
STAN_MODEL="--stan_model /nfs/research/jlees/jacqueline/bayesian_gwas_paper/00_stan_models_v1/ordinal-subcluster-standard-association_SPNPENcutpoints.stan"
ANALYSIS_TYPE="--analysis_type inference"
ANALYSIS_NICKNAME="--analysis_nickname spn_pen_OGdatasetK5_OGmodel_difpruning_paperPOM"
OUTPUT_DIR="--output_directory /nfs/research/jlees/jacqueline/thesis_results/compare_ordinal_models/paper_model_POM/spn_pen_OGdatasetK5_OGmodel_difpruning"
THREADS="--threads 48"

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

PHANDANGO="--phandango /nfs/research/jlees/jacqueline/bayesian_gwas_paper/00_data_v1/inference/01_spn_penicillin_subclusters_K5_variant_index.csv"
ANNOTATIONS="--annotations /nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/genotype/fields_filtered_maf05_multiallelic.txt"
MODEL_TYPE="--model_type pom"
NORATE="--norate"

# ---------------------------------------------------------------------------

mkdir -p /nfs/research/jlees/jacqueline/thesis_results/compare_ordinal_models/paper_model_POM/spn_pen_OGdatasetK5_OGmodel_difpruning/logs

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
