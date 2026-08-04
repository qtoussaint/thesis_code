#!/usr/bin/env bash

#SBATCH --job-name=ldprune_deprune
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=500G
#SBATCH --time=12:00:00
#SBATCH --error=/nfs/research/jlees/jacqueline/thesis_results/compare_ld_pruning/logs/deprune_%j.err
#SBATCH --output=/nfs/research/jlees/jacqueline/thesis_results/compare_ld_pruning/logs/deprune_%j.out

# Builds the per-variant beta + RATE table each arm is compared on, by
# de-pruning that arm's posterior draws along BacPrune's representative chains.
# See deprune_from_directions.R for why the pipeline's own de-pruned outputs are
# not used: they assume every representative is itself kept, which holds under
# r2 but not under |D'|, where it silently drops most of the genome.
#
# Writes into compare_ld_pruning/depruned/ rather than into each run directory,
# so the published r2 run is only ever read.
#
# Memory is sized for the no-pruning arm, whose coefficients.csv covers all
# 75,272 variants and so is roughly twice the 5.4 GB of the published run.
#
# Run once both fits have finished:
#   sbatch run_deprune.sh              # all three arms
#   sbatch run_deprune.sh r2 dprime    # only the named arms

source ~/.bashrc
mamba activate gwas_pipeline

set -euo pipefail

RESULTS=/nfs/research/jlees/jacqueline/thesis_results
COMPARE_DIR=$RESULTS/compare_ld_pruning
DEPRUNE=/nfs/research/jlees/jacqueline/thesis_code/compare_ld_pruning/deprune_from_directions.R

mkdir -p "$COMPARE_DIR/logs" "$COMPARE_DIR/depruned"

arm_dir () {
  case "$1" in
    none)   echo "$COMPARE_DIR/07_tb_rifampicin_binary_logistic_nopruning" ;;
    dprime) echo "$COMPARE_DIR/07_tb_rifampicin_binary_logistic_dprime" ;;
    r2)     echo "$RESULTS/gwas_tb_rifampicin/inference/07_tb_rifampicin_binary_logistic" ;;
    *)      echo "Unknown arm: $1" >&2; return 1 ;;
  esac
}

ARMS=("$@")
if [ ${#ARMS[@]} -eq 0 ]; then ARMS=(none dprime r2); fi

for ARM in "${ARMS[@]}"; do
  echo "=== $ARM ==="
  Rscript "$DEPRUNE" --run-dir "$(arm_dir "$ARM")" \
    --out "$COMPARE_DIR/depruned/${ARM}.csv"
done
