#!/usr/bin/env bash

#SBATCH --job-name=gwas_grm_dataset
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --error=/nfs/research/jlees/jacqueline/thesis_results/gwas_datasets/logs/grm_dataset.err
#SBATCH --output=/nfs/research/jlees/jacqueline/thesis_results/gwas_datasets/logs/grm_dataset.out

# Builds dataset 19 (SPN penicillin binary + GRM) from dataset 01.
# Memory note: the 588 x 32406 genotype matrix is read as doubles and
# cross-multiplied, so peak usage is a few hundred MB; 64G is generous.

source ~/.bashrc
mamba activate gwas_pipeline

set -euo pipefail

mkdir -p /nfs/research/jlees/jacqueline/thesis_results/gwas_datasets/logs

Rscript /nfs/research/jlees/jacqueline/thesis_code/gwas_datasets/write_grm_dataset.R
