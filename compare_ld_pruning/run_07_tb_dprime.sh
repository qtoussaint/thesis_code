#!/usr/bin/env bash
#SBATCH --job-name=tbrif_07_dprime_infe
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=400G
#SBATCH --time=36:00:00
#SBATCH --error=/nfs/research/jlees/jacqueline/thesis_results/compare_ld_pruning/07_tb_rifampicin_binary_logistic_dprime/logs/07_tb_rifampicin_binary_logistic_dprime_%j.err
#SBATCH --output=/nfs/research/jlees/jacqueline/thesis_results/compare_ld_pruning/07_tb_rifampicin_binary_logistic_dprime/logs/07_tb_rifampicin_binary_logistic_dprime_%j.out

#################################################################################
#
# Arm B of the LD pruning comparison: |D'| >= 1 pruning.
#
# Identical to the published run
# thesis_results/gwas_tb_rifampicin/inference/07_tb_rifampicin_binary_logistic
# in every argument except --pruning_metric, which switches BacPrune from r2
# (its default, and the published arm) to Lewontin's |D'|.
#
# At threshold 1 the two measures differ in a way that matters here. r2 >= 1
# only prunes variants that are perfectly correlated, i.e. exact duplicates.
# |D'| >= 1 also prunes nested variants, where one of the four haplotypes is
# simply absent, so it prunes strictly more and keeps fewer variants than the
# 40,941 the r2 arm retains.
#
# Fewer variants than the r2 arm, so the published run's resources are ample.
#
#################################################################################

source ~/.bashrc
mamba activate gwas_pipeline

mkdir -p /nfs/research/jlees/jacqueline/thesis_results/compare_ld_pruning/07_tb_rifampicin_binary_logistic_dprime/logs

RSCRIPT_PATH="/nfs/research/jlees/jacqueline/gwas_workflow/code/gwas_workflow/inst/scripts/run_pipeline.R"

DATA="--data /nfs/research/jlees/jacqueline/thesis_results/gwas_datasets/inference/07_tb_rifampicin_binary/07_tb_rifampicin_binary.json"
STAN_MODEL="--stan_model /nfs/research/jlees/jacqueline/thesis_code/gwas_finalmodels/logistic_inference.stan"
ANALYSIS_TYPE="--analysis_type inference"
ANALYSIS_NICKNAME="--analysis_nickname 07_tb_rifampicin_binary_logistic_dprime"
OUTPUT_DIR="--output_directory /nfs/research/jlees/jacqueline/thesis_results/compare_ld_pruning/07_tb_rifampicin_binary_logistic_dprime"
THREADS="--threads 48"

LD_PRUNING="--ld_pruning true"
PRUNING_SOFTWARE="--pruning_software /hps/software/users/jlees/jacqueline/manual_installs/bin/BacPrune-Rust/"
PRUNING_METRIC="--pruning_metric dprime"
MAF_CUTOFF="--maf_cutoff 0"
LD_THRESHOLD="--ld_threshold 1"

PHANDANGO="--phandango /nfs/research/jlees/jacqueline/thesis_results/gwas_datasets/inference/07_tb_rifampicin_binary/07_tb_rifampicin_binary_variant_index.csv"
ANNOTATIONS="--annotations /nfs/research/jlees/jacqueline/gwas_data/tuberculosis/cryptic_regeno_snpeff/fields_filtered.txt"
MODEL_TYPE="--model_type binary"
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
$PRUNING_METRIC \
$MAF_CUTOFF \
$LD_THRESHOLD \
$PHANDANGO \
$ANNOTATIONS \
$MODEL_TYPE \
$GENES_OF_INTEREST \
$RESUME \
$CPPRATE
