#!/usr/bin/env bash
#SBATCH --job-name=tbrif_14_PPOM_infe
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=800G
#SBATCH --time=24:00:00
#SBATCH --error=/nfs/research/jlees/jacqueline/thesis_results/gwas_tb_rifampicin/inference/14_tb_rifampicin_MIC_coarse_dilutions_PPOM/logs/14_tb_rifampicin_MIC_coarse_dilutions_PPOM_%j.err
#SBATCH --output=/nfs/research/jlees/jacqueline/thesis_results/gwas_tb_rifampicin/inference/14_tb_rifampicin_MIC_coarse_dilutions_PPOM/logs/14_tb_rifampicin_MIC_coarse_dilutions_PPOM_%j.out

#################################################################################

source ~/.bashrc
mamba activate gwas_pipeline

mkdir -p /nfs/research/jlees/jacqueline/thesis_results/gwas_tb_rifampicin/inference/14_tb_rifampicin_MIC_coarse_dilutions_PPOM/logs

RSCRIPT_PATH="/nfs/research/jlees/jacqueline/gwas_workflow/code/gwas_workflow/inst/scripts/run_pipeline.R"

DATA="--data /nfs/research/jlees/jacqueline/thesis_results/gwas_datasets/inference/14_tb_rifampicin_MIC_coarse_dilutions/14_tb_rifampicin_MIC_coarse_dilutions.json"
STAN_MODEL="--stan_model /nfs/research/jlees/jacqueline/thesis_code/gwas_finalmodels/PPOM_inference.stan"
ANALYSIS_TYPE="--analysis_type inference"
ANALYSIS_NICKNAME="--analysis_nickname 14_tb_rifampicin_MIC_coarse_dilutions_PPOM"
OUTPUT_DIR="--output_directory /nfs/research/jlees/jacqueline/thesis_results/gwas_tb_rifampicin/inference/14_tb_rifampicin_MIC_coarse_dilutions_PPOM"
THREADS="--threads 48"

LD_PRUNING="--ld_pruning true"
PRUNING_SOFTWARE="--pruning_software /hps/software/users/jlees/jacqueline/manual_installs/bin/BacPrune-Rust/"
MAF_CUTOFF="--maf_cutoff 0"
LD_THRESHOLD="--ld_threshold 1"

PHANDANGO="--phandango /nfs/research/jlees/jacqueline/thesis_results/gwas_datasets/inference/14_tb_rifampicin_MIC_coarse_dilutions/14_tb_rifampicin_MIC_coarse_dilutions_variant_index.csv"
ANNOTATIONS="--annotations /nfs/research/jlees/jacqueline/gwas_data/tuberculosis/cryptic_regeno_snpeff/cryptic_regeno_fields_filtered.txt"
MODEL_TYPE="--model_type ppom"
GENES_OF_INTEREST="--genes_of_interest /nfs/research/jlees/jacqueline/thesis_code/gwas_genesofinterest/tb_rifampicin_genesofinterest.txt"
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

