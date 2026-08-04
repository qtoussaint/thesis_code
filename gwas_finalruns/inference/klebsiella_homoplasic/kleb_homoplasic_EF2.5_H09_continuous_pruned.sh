#!/usr/bin/env bash
#SBATCH --job-name=kEF2.5_prn
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=250G
#SBATCH --time=24:00:00
#SBATCH --error=/nfs/research/jlees/jacqueline/thesis_results/gwas_klebsiella_homoplasic/inference/kleb_homoplasic_EF2.5_H09_continuous_pruned/logs/kleb_homoplasic_EF2.5_H09_continuous_pruned_%j.err
#SBATCH --output=/nfs/research/jlees/jacqueline/thesis_results/gwas_klebsiella_homoplasic/inference/kleb_homoplasic_EF2.5_H09_continuous_pruned/logs/kleb_homoplasic_EF2.5_H09_continuous_pruned_%j.out

#################################################################################

source ~/.bashrc
mamba activate gwas_pipeline

mkdir -p /nfs/research/jlees/jacqueline/thesis_results/gwas_klebsiella_homoplasic/inference/kleb_homoplasic_EF2.5_H09_continuous_pruned/logs

RSCRIPT_PATH="/nfs/research/jlees/jacqueline/gwas_workflow/code/gwas_workflow/inst/scripts/run_pipeline.R"

DATA="--data /nfs/research/jlees/jacqueline/thesis_results/gwas_datasets/inference/kleb_homoplasic_EF2.5_H09/kleb_homoplasic_EF2.5_H09.json"
STAN_MODEL="--stan_model /nfs/research/jlees/jacqueline/thesis_code/gwas_finalmodels/continuous_inference.stan"
ANALYSIS_TYPE="--analysis_type inference"
ANALYSIS_NICKNAME="--analysis_nickname kleb_homoplasic_EF2.5_H09_continuous_pruned"
OUTPUT_DIR="--output_directory /nfs/research/jlees/jacqueline/thesis_results/gwas_klebsiella_homoplasic/inference/kleb_homoplasic_EF2.5_H09_continuous_pruned"
THREADS="--threads 32"

LD_PRUNING="--ld_pruning true"
PRUNING_SOFTWARE="--pruning_software /hps/software/users/jlees/jacqueline/manual_installs/bin/BacPrune-Rust/"
MAF_CUTOFF="--maf_cutoff 0"
LD_THRESHOLD="--ld_threshold 1"

PHANDANGO="--phandango /nfs/research/jlees/jacqueline/thesis_results/gwas_datasets/inference/kleb_homoplasic_EF2.5_H09/kleb_homoplasic_EF2.5_H09_variant_index.csv"
ANNOTATIONS="--annotations /nfs/research/jlees/jacqueline/gwas_data/klebsiella/genotype/klebsiella_rs_annotations.txt"
MODEL_TYPE="--model_type continuous"
GRAD_SAMPLES="--grad_samples 10"
VI_SEED="--vi_seed 987654321"
GENES_OF_INTEREST="--genes_of_interest /nfs/research/jlees/jacqueline/thesis_code/gwas_genesofinterest/klebsiella_homoplasic_genesofinterest.txt"
RESUME=""
CPPRATE="--cpprate_bin /hps/software/users/jlees/jacqueline/manual_installs/bin/cpprate-0.2.0/build/bin/cpprate"

Rscript $RSCRIPT_PATH $DATA $STAN_MODEL $ANALYSIS_TYPE $ANALYSIS_NICKNAME $OUTPUT_DIR $THREADS $LD_PRUNING $PRUNING_SOFTWARE $MAF_CUTOFF $LD_THRESHOLD $PHANDANGO $ANNOTATIONS $MODEL_TYPE $GRAD_SAMPLES $VI_SEED $GENES_OF_INTEREST $RESUME $CPPRATE
