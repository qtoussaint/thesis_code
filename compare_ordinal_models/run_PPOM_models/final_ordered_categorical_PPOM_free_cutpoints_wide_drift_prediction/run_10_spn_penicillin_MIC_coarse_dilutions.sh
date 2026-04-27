#!/usr/bin/env bash

#SBATCH --job-name=freeCutsWD_10_spn_pen_coarse_pred
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=400G
#SBATCH --time=5:00:00
#SBATCH --error=/nfs/research/jlees/jacqueline/thesis_results/compare_ordinal_models/final_ordered_categorical_PPOM_free_cutpoints_wide_drift_prediction/10_spn_penicillin_MIC_coarse_dilutions/logs/10_spn_penicillin_MIC_coarse_dilutions.err
#SBATCH --output=/nfs/research/jlees/jacqueline/thesis_results/compare_ordinal_models/final_ordered_categorical_PPOM_free_cutpoints_wide_drift_prediction/10_spn_penicillin_MIC_coarse_dilutions/logs/10_spn_penicillin_MIC_coarse_dilutions.out

#################################################################################

source /hps/software/users/jlees/jacqueline/etc/profile.d/conda.sh
conda activate gwas_pipeline

RSCRIPT_PATH="/nfs/research/jlees/jacqueline/gwas_workflow/code/gwas_workflow/inst/scripts/run_pipeline.R"

DATA="--data /nfs/research/jlees/jacqueline/thesis_results/gwas_datasets/prediction/10_spn_penicillin_MIC_coarse_dilutions/10_spn_penicillin_MIC_coarse_dilutions.json"
STAN_MODEL="--stan_model /nfs/research/jlees/jacqueline/thesis_code/compare_ordinal_models/PPOM_models/final_ordered_categorical_PPOM_free_cutpoints_wide_drift_prediction.stan"
ANALYSIS_TYPE="--analysis_type prediction"
ANALYSIS_NICKNAME="--analysis_nickname 10_spn_penicillin_MIC_coarse_dilutions_freeCutsWDPPOMPred"
OUTPUT_DIR="--output_directory /nfs/research/jlees/jacqueline/thesis_results/compare_ordinal_models/final_ordered_categorical_PPOM_free_cutpoints_wide_drift_prediction/10_spn_penicillin_MIC_coarse_dilutions"
THREADS="--threads 48"

LD_PRUNING="--ld_pruning true"
PRUNING_SOFTWARE="--pruning_software /hps/software/users/jlees/jacqueline/manual_installs/bin/BacPrune-Rust/"
MAF_CUTOFF="--maf_cutoff 0"
LD_THRESHOLD="--ld_threshold 1"

PHANDANGO="--phandango /nfs/research/jlees/jacqueline/thesis_results/gwas_datasets/prediction/10_spn_penicillin_MIC_coarse_dilutions/10_spn_penicillin_MIC_coarse_dilutions_variant_index.csv"
ANNOTATIONS="--annotations /nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/genotype/fields_filtered_maf05_multiallelic.txt"
MODEL_TYPE="--model_type ppom"
GENES_OF_INTEREST="--genes_of_interest /nfs/research/jlees/jacqueline/thesis_code/gwas_genesofinterest/spn_penicillin_genesofinterest.txt"
NORATE="--norate"
RESUME="--resume"

mkdir -p /nfs/research/jlees/jacqueline/thesis_results/compare_ordinal_models/final_ordered_categorical_PPOM_free_cutpoints_wide_drift_prediction/10_spn_penicillin_MIC_coarse_dilutions/logs

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
$GENES_OF_INTEREST \
$NORATE \
$RESUME \
