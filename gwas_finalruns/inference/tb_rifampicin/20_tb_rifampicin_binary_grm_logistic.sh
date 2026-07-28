#!/usr/bin/env bash
#SBATCH --job-name=tbrif_20_logistic_infe
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=800G
#SBATCH --time=72:00:00
#SBATCH --error=/nfs/research/jlees/jacqueline/thesis_results/gwas_tb_rifampicin/inference/20_tb_rifampicin_binary_grm_logistic/logs/20_tb_rifampicin_binary_grm_logistic_%j.err
#SBATCH --output=/nfs/research/jlees/jacqueline/thesis_results/gwas_tb_rifampicin/inference/20_tb_rifampicin_binary_grm_logistic/logs/20_tb_rifampicin_binary_grm_logistic_%j.out

#################################################################################
# GRM twin of 07_tb_rifampicin_binary_logistic.sh. Same dataset, same pipeline
# settings; the only differences are the GRM Stan model and the dataset 20 JSON
# (dataset 07 plus a GRM field).
#
# Resources are raised over run 07 (400G/36h): the 11,622 x 11,622 GRM adds a
# one-off Cholesky in transformed data and a dense matrix-vector product on the
# autodiff tape at every gradient evaluation, which is the quadratic scaling in
# sample size that the subcluster encoding avoids.
#################################################################################

source ~/.bashrc
mamba activate gwas_pipeline

mkdir -p /nfs/research/jlees/jacqueline/thesis_results/gwas_tb_rifampicin/inference/20_tb_rifampicin_binary_grm_logistic/logs

RSCRIPT_PATH="/nfs/research/jlees/jacqueline/gwas_workflow/code/gwas_workflow/inst/scripts/run_pipeline.R"

DATA="--data /nfs/research/jlees/jacqueline/thesis_results/gwas_datasets/inference/20_tb_rifampicin_binary_grm/20_tb_rifampicin_binary_grm.json"
STAN_MODEL="--stan_model /nfs/research/jlees/jacqueline/thesis_code/gwas_finalmodels/logistic_grm_inference.stan"
ANALYSIS_TYPE="--analysis_type inference"
ANALYSIS_NICKNAME="--analysis_nickname 20_tb_rifampicin_binary_grm_logistic"
OUTPUT_DIR="--output_directory /nfs/research/jlees/jacqueline/thesis_results/gwas_tb_rifampicin/inference/20_tb_rifampicin_binary_grm_logistic"
THREADS="--threads 48"

LD_PRUNING="--ld_pruning true"
PRUNING_SOFTWARE="--pruning_software /hps/software/users/jlees/jacqueline/manual_installs/bin/BacPrune-Rust/"
MAF_CUTOFF="--maf_cutoff 0"
LD_THRESHOLD="--ld_threshold 1"

PHANDANGO="--phandango /nfs/research/jlees/jacqueline/thesis_results/gwas_datasets/inference/20_tb_rifampicin_binary_grm/20_tb_rifampicin_binary_grm_variant_index.csv"
ANNOTATIONS="--annotations /nfs/research/jlees/jacqueline/gwas_data/tuberculosis/cryptic_regeno_snpeff/fields_filtered.txt"
MODEL_TYPE="--model_type binary"
GENES_OF_INTEREST="--genes_of_interest /nfs/research/jlees/jacqueline/thesis_code/gwas_genesofinterest/tb_rifampicin_genesofinterest.txt"
RESUME="--resume"
CPPRATE="--cpprate_bin /hps/software/users/jlees/jacqueline/manual_installs/bin/cpprate-0.2.0/build/bin/cpprate"
# Same VI seed as run 07 so the comparison is matched. The SPN GRM model failed
# at grad_samples = 1 with this seed; 10 gradient samples is carried over here.
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
