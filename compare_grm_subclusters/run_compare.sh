#!/usr/bin/env bash

#SBATCH --job-name=compare_grm_subclusters
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --error=/nfs/research/jlees/jacqueline/thesis_results/compare_grm_subclusters/logs/compare_%j.err
#SBATCH --output=/nfs/research/jlees/jacqueline/thesis_results/compare_grm_subclusters/logs/compare_%j.out

# Compares the subcluster and GRM runs for one dataset pair.
#   sbatch run_compare.sh spn_penicillin
#   sbatch run_compare.sh tb_rifampicin
# Defaults to spn_penicillin. The TB comparison reads 75,272-variant tables,
# hence the larger memory than the plotting itself needs.

source ~/.bashrc
mamba activate gwas_pipeline

set -euo pipefail

COMPARISON="${1:-spn_penicillin}"

mkdir -p /nfs/research/jlees/jacqueline/thesis_results/compare_grm_subclusters/logs

Rscript /nfs/research/jlees/jacqueline/thesis_code/compare_grm_subclusters/compare_grm_vs_subclusters.R \
  "$COMPARISON"
