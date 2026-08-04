#!/usr/bin/env bash

#SBATCH --job-name=compare_lineage_subclusters
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --error=/nfs/research/jlees/jacqueline/thesis_results/compare_lineage_subclusters/logs/compare_%j.err
#SBATCH --output=/nfs/research/jlees/jacqueline/thesis_results/compare_lineage_subclusters/logs/compare_%j.out

# Compares the three population-structure arms for one dataset.
#   sbatch run_compare.sh tb_rifampicin
# Defaults to tb_rifampicin. The TB comparison reads 75,272-variant tables,
# hence the larger memory than the plotting itself needs.

source ~/.bashrc
mamba activate gwas_pipeline

set -euo pipefail

COMPARISON="${1:-tb_rifampicin}"

mkdir -p /nfs/research/jlees/jacqueline/thesis_results/compare_lineage_subclusters/logs

Rscript /nfs/research/jlees/jacqueline/thesis_code/compare_lineage_subclusters/compare_lineage_vs_subclusters.R \
  "$COMPARISON"
