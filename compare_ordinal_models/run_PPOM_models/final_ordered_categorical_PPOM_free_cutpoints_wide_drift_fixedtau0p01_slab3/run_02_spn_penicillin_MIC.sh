#!/usr/bin/env bash

#SBATCH --job-name=freeCutsWDFT_02_spn_pen
#SBATCH --nodes=1
#SBATCH --cpus-per-task=80
#SBATCH --mem=650G
#SBATCH --time=12:00:00
#SBATCH --error=/nfs/research/jlees/jacqueline/thesis_results/compare_ordinal_models/final_ordered_categorical_PPOM_free_cutpoints_wide_drift_fixedtau0p01_slab3/02_spn_penicillin_MIC/logs/02_spn_penicillin_MIC.err
#SBATCH --output=/nfs/research/jlees/jacqueline/thesis_results/compare_ordinal_models/final_ordered_categorical_PPOM_free_cutpoints_wide_drift_fixedtau0p01_slab3/02_spn_penicillin_MIC/logs/02_spn_penicillin_MIC.out

#################################################################################

source ~/.bashrc
mamba activate gwas_pipeline

RSCRIPT_PATH="/nfs/research/jlees/jacqueline/gwas_workflow/code/gwas_workflow/inst/scripts/run_pipeline.R"

DATA="--data /nfs/research/jlees/jacqueline/thesis_results/gwas_datasets/inference/02_spn_penicillin_MIC/02_spn_penicillin_MIC.json"
STAN_MODEL="--stan_model /nfs/research/jlees/jacqueline/thesis_code/compare_ordinal_models/PPOM_models/final_ordered_categorical_PPOM_free_cutpoints_wide_drift_fixedtau0p01_slab3.stan"
ANALYSIS_TYPE="--analysis_type inference"
ANALYSIS_NICKNAME="--analysis_nickname 02_spn_penicillin_MIC_freeCutsWDFTPPOM"
OUTPUT_DIR="--output_directory /nfs/research/jlees/jacqueline/thesis_results/compare_ordinal_models/final_ordered_categorical_PPOM_free_cutpoints_wide_drift_fixedtau0p01_slab3/02_spn_penicillin_MIC"
THREADS="--threads 80"

LD_PRUNING="--ld_pruning true"
PRUNING_SOFTWARE="--pruning_software /hps/software/users/jlees/jacqueline/manual_installs/bin/BacPrune-Rust/"
MAF_CUTOFF="--maf_cutoff 0"
LD_THRESHOLD="--ld_threshold 1"

PHANDANGO="--phandango /nfs/research/jlees/jacqueline/thesis_results/gwas_datasets/inference/02_spn_penicillin_MIC/02_spn_penicillin_MIC_variant_index.csv"
ANNOTATIONS="--annotations /nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/genotype/fields_filtered_maf05_multiallelic.txt"
MODEL_TYPE="--model_type ppom"
GENES_OF_INTEREST="--genes_of_interest /nfs/research/jlees/jacqueline/thesis_code/gwas_genesofinterest/spn_penicillin_genesofinterest.txt"
NORATE="--norate"
RESUME="--resume"

mkdir -p /nfs/research/jlees/jacqueline/thesis_results/compare_ordinal_models/final_ordered_categorical_PPOM_free_cutpoints_wide_drift_fixedtau0p01_slab3/02_spn_penicillin_MIC/logs

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
