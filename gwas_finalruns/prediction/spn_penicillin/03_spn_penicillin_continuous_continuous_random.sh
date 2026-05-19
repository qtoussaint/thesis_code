#!/usr/bin/env bash
#SBATCH --job-name=spnpen_03_continuous_pred_random
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=200G
#SBATCH --time=06:00:00
#SBATCH --error=/nfs/research/jlees/jacqueline/thesis_results/gwas_spn_penicillin/prediction/03_spn_penicillin_continuous_continuous_random/logs/03_spn_penicillin_continuous_continuous_random.err
#SBATCH --output=/nfs/research/jlees/jacqueline/thesis_results/gwas_spn_penicillin/prediction/03_spn_penicillin_continuous_continuous_random/logs/03_spn_penicillin_continuous_continuous_random.out

#################################################################################

source ~/.bashrc
mamba activate gwas_pipeline

mkdir -p /nfs/research/jlees/jacqueline/thesis_results/gwas_spn_penicillin/prediction/03_spn_penicillin_continuous_continuous_random/logs

RSCRIPT_PATH="/nfs/research/jlees/jacqueline/gwas_workflow/code/gwas_workflow/inst/scripts/run_pipeline.R"

DATA="--data /nfs/research/jlees/jacqueline/thesis_results/gwas_datasets/prediction/03_spn_penicillin_continuous/03_spn_penicillin_continuous.json"
STAN_MODEL="--stan_model /nfs/research/jlees/jacqueline/thesis_code/gwas_finalmodels/continuous_prediction.stan"
ANALYSIS_TYPE="--analysis_type prediction"
ANALYSIS_NICKNAME="--analysis_nickname 03_spn_penicillin_continuous_continuous_random"
OUTPUT_DIR="--output_directory /nfs/research/jlees/jacqueline/thesis_results/gwas_spn_penicillin/prediction/03_spn_penicillin_continuous_continuous_random"
THREADS="--threads 32"

LD_PRUNING="--ld_pruning true"
PRUNING_SOFTWARE="--pruning_software /hps/software/users/jlees/jacqueline/manual_installs/bin/BacPrune-Rust/"
MAF_CUTOFF="--maf_cutoff 0"
LD_THRESHOLD="--ld_threshold 1"

PHANDANGO="--phandango /nfs/research/jlees/jacqueline/thesis_results/gwas_datasets/prediction/03_spn_penicillin_continuous/03_spn_penicillin_continuous_variant_index.csv"
ANNOTATIONS="--annotations /nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/genotype/fields_filtered_maf05_multiallelic.txt"
MODEL_TYPE="--model_type continuous"
GENES_OF_INTEREST="--genes_of_interest /nfs/research/jlees/jacqueline/thesis_code/gwas_genesofinterest/spn_penicillin_genesofinterest.txt"
RESUME="--resume"
CPPRATE="--cpprate_bin /hps/software/users/jlees/jacqueline/manual_installs/bin/cpprate-0.2.0/build/bin/cpprate"
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
$RESUME \
$CPPRATE

