#!/usr/bin/env bash
#SBATCH --job-name=tbrif_07_nopruning_infe
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=900G
#SBATCH --time=5-00:00:00
#SBATCH --error=/nfs/research/jlees/jacqueline/thesis_results/compare_ld_pruning/07_tb_rifampicin_binary_logistic_nopruning/logs/07_tb_rifampicin_binary_logistic_nopruning_%j.err
#SBATCH --output=/nfs/research/jlees/jacqueline/thesis_results/compare_ld_pruning/07_tb_rifampicin_binary_logistic_nopruning/logs/07_tb_rifampicin_binary_logistic_nopruning_%j.out

#################################################################################
#
# Arm A of the LD pruning comparison: no pruning at all.
#
# Reuses dataset 07, the same one behind
# thesis_results/gwas_tb_rifampicin/inference/07_tb_rifampicin_binary_logistic,
# so phenotype, genotype and sample order are identical to that run by
# construction. The model is the as-published logistic model; the only thing
# that changes across the three arms is which variants survive pruning.
#
# All 75,272 variants enter the model, against 40,941 under r2 >= 1 pruning, so
# this arm needs roughly twice the memory and time of the published run. The
# scheduler here rejects an explicit --partition=bigmem; asking for the memory
# is what routes the job to a big-memory node.
#
# Because nothing is pruned, this arm needs no de-pruning: cppRATE writes
# RATE_values.txt directly on the original variant index. The no-pruning branch
# of the pipeline writes no variant effects CSV at all, so betas come from
# run_deprune.sh afterwards, which reduces the posterior draws this run leaves
# in cppRATE_matrices/coefficients.csv.
#
#################################################################################

source ~/.bashrc
mamba activate gwas_pipeline

mkdir -p /nfs/research/jlees/jacqueline/thesis_results/compare_ld_pruning/07_tb_rifampicin_binary_logistic_nopruning/logs

RSCRIPT_PATH="/nfs/research/jlees/jacqueline/gwas_workflow/code/gwas_workflow/inst/scripts/run_pipeline.R"

DATA="--data /nfs/research/jlees/jacqueline/thesis_results/gwas_datasets/inference/07_tb_rifampicin_binary/07_tb_rifampicin_binary.json"
STAN_MODEL="--stan_model /nfs/research/jlees/jacqueline/thesis_code/gwas_finalmodels/logistic_inference.stan"
ANALYSIS_TYPE="--analysis_type inference"
ANALYSIS_NICKNAME="--analysis_nickname 07_tb_rifampicin_binary_logistic_nopruning"
OUTPUT_DIR="--output_directory /nfs/research/jlees/jacqueline/thesis_results/compare_ld_pruning/07_tb_rifampicin_binary_logistic_nopruning"
THREADS="--threads 48"

LD_PRUNING="--ld_pruning false"

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
$PHANDANGO \
$ANNOTATIONS \
$MODEL_TYPE \
$GENES_OF_INTEREST \
$RESUME \
$CPPRATE
