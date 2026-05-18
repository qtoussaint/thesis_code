#!/usr/bin/env bash
#SBATCH --job-name=spnpen_11_POM_infe
#SBATCH --nodes=1
#SBATCH --cpus-per-task=80
#SBATCH --mem=650G
#SBATCH --time=12:00:00
#SBATCH --error=/nfs/research/jlees/jacqueline/thesis_results/gwas_spn_penicillin/inference/11_spn_penicillin_MIC_large_minbin_POM/logs/11_spn_penicillin_MIC_large_minbin_POM.err
#SBATCH --output=/nfs/research/jlees/jacqueline/thesis_results/gwas_spn_penicillin/inference/11_spn_penicillin_MIC_large_minbin_POM/logs/11_spn_penicillin_MIC_large_minbin_POM.out

#################################################################################

source ~/.bashrc
mamba activate gwas_pipeline

mkdir -p /nfs/research/jlees/jacqueline/thesis_results/gwas_spn_penicillin/inference/11_spn_penicillin_MIC_large_minbin_POM/logs

RSCRIPT_PATH="/nfs/research/jlees/jacqueline/gwas_workflow/code/gwas_workflow/inst/scripts/run_pipeline.R"

DATA="--data /nfs/research/jlees/jacqueline/thesis_results/gwas_datasets/inference/11_spn_penicillin_MIC_large_minbin/11_spn_penicillin_MIC_large_minbin.json"
STAN_MODEL="--stan_model /nfs/research/jlees/jacqueline/thesis_code/gwas_finalmodels/POM_inference.stan"
ANALYSIS_TYPE="--analysis_type inference"
ANALYSIS_NICKNAME="--analysis_nickname 11_spn_penicillin_MIC_large_minbin_POM"
OUTPUT_DIR="--output_directory /nfs/research/jlees/jacqueline/thesis_results/gwas_spn_penicillin/inference/11_spn_penicillin_MIC_large_minbin_POM"
THREADS="--threads 80"

LD_PRUNING="--ld_pruning true"
PRUNING_SOFTWARE="--pruning_software /hps/software/users/jlees/jacqueline/manual_installs/bin/BacPrune-Rust/"
MAF_CUTOFF="--maf_cutoff 0"
LD_THRESHOLD="--ld_threshold 1"

PHANDANGO="--phandango /nfs/research/jlees/jacqueline/thesis_results/gwas_datasets/inference/11_spn_penicillin_MIC_large_minbin/11_spn_penicillin_MIC_large_minbin_variant_index.csv"
ANNOTATIONS="--annotations /nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/genotype/fields_filtered_maf05_multiallelic.txt"
MODEL_TYPE="--model_type pom"
GENES_OF_INTEREST="--genes_of_interest /nfs/research/jlees/jacqueline/thesis_code/gwas_genesofinterest/spn_penicillin_genesofinterest.txt"
RESUME="--resume"
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
$RESUME

