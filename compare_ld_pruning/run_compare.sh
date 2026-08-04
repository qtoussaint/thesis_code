#!/usr/bin/env bash

#SBATCH --job-name=compare_ld_pruning
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --error=/nfs/research/jlees/jacqueline/thesis_results/compare_ld_pruning/logs/compare_%j.err
#SBATCH --output=/nfs/research/jlees/jacqueline/thesis_results/compare_ld_pruning/logs/compare_%j.out

# Compares the three LD pruning arms for one dataset.
#   sbatch run_compare.sh tb_rifampicin
# Defaults to tb_rifampicin. Reads 75,272-variant tables for three arms, hence
# the larger memory than the plotting itself needs.
#
# Run after run_deprune.sh, which is what produces the per-arm beta + RATE
# tables in compare_ld_pruning/depruned/.

source ~/.bashrc
mamba activate gwas_pipeline

set -euo pipefail

COMPARISON="${1:-tb_rifampicin}"

mkdir -p /nfs/research/jlees/jacqueline/thesis_results/compare_ld_pruning/logs

Rscript /nfs/research/jlees/jacqueline/thesis_code/compare_ld_pruning/compare_ld_pruning.R \
  "$COMPARISON"
