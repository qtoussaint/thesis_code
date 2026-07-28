#!/usr/bin/env bash
#SBATCH --job-name=spnpen_19_logistic_infe
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=200G
#SBATCH --time=06:00:00
#SBATCH --error=/nfs/research/jlees/jacqueline/thesis_results/gwas_spn_penicillin/inference/19_spn_penicillin_binary_grm_logistic/logs/19_spn_penicillin_binary_grm_logistic_%j.err
#SBATCH --output=/nfs/research/jlees/jacqueline/thesis_results/gwas_spn_penicillin/inference/19_spn_penicillin_binary_grm_logistic/logs/19_spn_penicillin_binary_grm_logistic_%j.out

#################################################################################
# GRM twin of 01_spn_penicillin_binary_logistic.sh. Same dataset, same pipeline
# settings; the only differences are the GRM Stan model and the dataset 19 JSON
# (dataset 01 plus a GRM field). Hand-written rather than emitted by
# generate_run_scripts.R, which would also emit prediction and lasso scripts
# for datasets and models that do not exist here.
#################################################################################

source ~/.bashrc
mamba activate gwas_pipeline

mkdir -p /nfs/research/jlees/jacqueline/thesis_results/gwas_spn_penicillin/inference/19_spn_penicillin_binary_grm_logistic/logs

RSCRIPT_PATH="/nfs/research/jlees/jacqueline/gwas_workflow/code/gwas_workflow/inst/scripts/run_pipeline.R"

DATA="--data /nfs/research/jlees/jacqueline/thesis_results/gwas_datasets/inference/19_spn_penicillin_binary_grm/19_spn_penicillin_binary_grm.json"
STAN_MODEL="--stan_model /nfs/research/jlees/jacqueline/thesis_code/gwas_finalmodels/logistic_grm_inference.stan"
ANALYSIS_TYPE="--analysis_type inference"
ANALYSIS_NICKNAME="--analysis_nickname 19_spn_penicillin_binary_grm_logistic"
OUTPUT_DIR="--output_directory /nfs/research/jlees/jacqueline/thesis_results/gwas_spn_penicillin/inference/19_spn_penicillin_binary_grm_logistic"
THREADS="--threads 32"

LD_PRUNING="--ld_pruning true"
PRUNING_SOFTWARE="--pruning_software /hps/software/users/jlees/jacqueline/manual_installs/bin/BacPrune-Rust/"
MAF_CUTOFF="--maf_cutoff 0"
LD_THRESHOLD="--ld_threshold 1"

PHANDANGO="--phandango /nfs/research/jlees/jacqueline/thesis_results/gwas_datasets/inference/19_spn_penicillin_binary_grm/19_spn_penicillin_binary_grm_variant_index.csv"
ANNOTATIONS="--annotations /nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/genotype/fields_filtered_maf05_multiallelic.txt"
MODEL_TYPE="--model_type binary"
GENES_OF_INTEREST="--genes_of_interest /nfs/research/jlees/jacqueline/thesis_code/gwas_genesofinterest/spn_penicillin_genesofinterest.txt"
RESUME="--resume"
CPPRATE="--cpprate_bin /hps/software/users/jlees/jacqueline/manual_installs/bin/cpprate-0.2.0/build/bin/cpprate"
# Same VI seed as run 01 so the comparison is matched. At this seed the GRM
# model fails with grad_samples = 1 ("dropped evaluations"); averaging the
# ADVI gradient over 10 samples is what lets it fit while holding the seed.
VI_SEED="--vi_seed 987654321"
GRAD_SAMPLES="--grad_samples 10"
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
$CPPRATE \
$VI_SEED \
$GRAD_SAMPLES
