#!/usr/bin/env bash

#SBATCH --job-name=gwas_grm_dataset_tb
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=600G
#SBATCH --time=12:00:00
#SBATCH --error=/nfs/research/jlees/jacqueline/thesis_results/gwas_datasets/logs/grm_dataset_tb.err
#SBATCH --output=/nfs/research/jlees/jacqueline/thesis_results/gwas_datasets/logs/grm_dataset_tb.out

# Builds dataset 20 (TB rifampicin binary + GRM) from dataset 07.
#
# Memory note: the source JSON is 1.75 GB and holds an 11,622 x 75,272 genotype
# matrix; parsing it and holding a double copy for the cross-product is the peak
# cost. The GRM itself is 11,622^2 doubles (~1.1 GB) and the output JSON is
# ~4.5 GB.
# CPUs: OpenBLAS threads the 1e13-flop cross-product.

source ~/.bashrc
mamba activate gwas_pipeline

set -euo pipefail

mkdir -p /nfs/research/jlees/jacqueline/thesis_results/gwas_datasets/logs

# The GRM is rounded to 4 significant figures to shrink the JSON: at full
# precision it is ~2.7 GB on its own. The diagonal jitter is raised to 1e-3
# inside the script to keep the rounded matrix positive definite.
Rscript /nfs/research/jlees/jacqueline/thesis_code/gwas_datasets/write_grm_dataset.R \
  07_tb_rifampicin_binary \
  20_tb_rifampicin_binary_grm \
  4
