#!/usr/bin/env bash

#SBATCH --job-name=01_logistic_prediction_alphasd3_pred
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=600G
#SBATCH --time=06:00:00
#SBATCH --error=/nfs/research/jlees/jacqueline/thesis_results/compare_ordinal_models/parameter_sweeps_prediction/runs/logistic_prediction_alphasd3/01_spn_penicillin_binary/logs/01_spn_penicillin_binary.err
#SBATCH --output=/nfs/research/jlees/jacqueline/thesis_results/compare_ordinal_models/parameter_sweeps_prediction/runs/logistic_prediction_alphasd3/01_spn_penicillin_binary/logs/01_spn_penicillin_binary.out

#################################################################################

source /hps/software/users/jlees/jacqueline/etc/profile.d/conda.sh
conda activate gwas_pipeline

RSCRIPT_PATH="/nfs/research/jlees/jacqueline/gwas_workflow/code/gwas_workflow/inst/scripts/run_pipeline.R"

DATA="--data /nfs/research/jlees/jacqueline/thesis_results/gwas_datasets/prediction/01_spn_penicillin_binary/01_spn_penicillin_binary.json"
STAN_MODEL="--stan_model /nfs/research/jlees/jacqueline/thesis_code/compare_ordinal_models/parameter_sweeps_prediction/stan_models/logistic_prediction_alphasd3.stan"
ANALYSIS_TYPE="--analysis_type prediction"
ANALYSIS_NICKNAME="--analysis_nickname 01_logistic_prediction_alphasd3"
OUTPUT_DIR="--output_directory /nfs/research/jlees/jacqueline/thesis_results/compare_ordinal_models/parameter_sweeps_prediction/runs/logistic_prediction_alphasd3/01_spn_penicillin_binary"
THREADS="--threads 48"

LD_PRUNING="--ld_pruning true"
PRUNING_SOFTWARE="--pruning_software /hps/software/users/jlees/jacqueline/manual_installs/bin/BacPrune-Rust/"
MAF_CUTOFF="--maf_cutoff 0"
LD_THRESHOLD="--ld_threshold 1"

PHANDANGO="--phandango /nfs/research/jlees/jacqueline/thesis_results/gwas_datasets/prediction/01_spn_penicillin_binary/01_spn_penicillin_binary_variant_index.csv"
ANNOTATIONS="--annotations /nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/genotype/fields_filtered_maf05_multiallelic.txt"
MODEL_TYPE="--model_type binary"
TRUE_PHENOTYPES="--true_phenotypes /nfs/research/jlees/jacqueline/thesis_results/gwas_datasets/prediction/01_spn_penicillin_binary/01_spn_penicillin_binary_test_phenotypes.csv"
GENES_OF_INTEREST="--genes_of_interest /nfs/research/jlees/jacqueline/thesis_code/gwas_genesofinterest/spn_penicillin_genesofinterest.txt"
NORATE="--norate"
RESUME="--resume"

mkdir -p /nfs/research/jlees/jacqueline/thesis_results/compare_ordinal_models/parameter_sweeps_prediction/runs/logistic_prediction_alphasd3/01_spn_penicillin_binary/logs

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

