#!/usr/bin/env bash
#SBATCH --job-name=regen_preview_effects
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=03:00:00
#SBATCH --error=/nfs/research/jlees/jacqueline/thesis_results/compare_lineage_subclusters/logs/regen_preview_%j.err
#SBATCH --output=/nfs/research/jlees/jacqueline/thesis_results/compare_lineage_subclusters/logs/regen_preview_%j.out

# Builds depruned_variant_effects_preview.csv for the two new arms from the
# posterior draws already written to cppRATE_matrices/coefficients.csv, so the
# comparison figures can be made while cppRATE is still running. Reads ~5.4 GB
# per arm into memory in one pass, hence the 128G request.
#
# The _preview filename cannot collide with the depruned_variant_effects.csv the
# running pipeline writes later; the comparison script prefers the canonical file
# whenever it exists.

source ~/.bashrc
mamba activate gwas_pipeline

set -euo pipefail

BASE=/nfs/research/jlees/jacqueline/thesis_results/compare_lineage_subclusters
mkdir -p $BASE/logs

# All three arms, including the pre-existing subclusters-only run, so every arm's
# effects come from the same quantity (beta_variant) on the same code path.
for RUN in "$BASE/07_tb_rifampicin_binary_logistic_lineage" \
           "$BASE/07_tb_rifampicin_binary_logistic_lineage_subcluster" \
           "/nfs/research/jlees/jacqueline/thesis_results/gwas_tb_rifampicin/inference/07_tb_rifampicin_binary_logistic"; do
  echo "=============== $RUN"
  if [ ! -f "$RUN/cppRATE_matrices/coefficients.csv" ]; then
    echo "  SKIP: no coefficients.csv yet"
    continue
  fi
  Rscript /nfs/research/jlees/jacqueline/thesis_code/compare_lineage_subclusters/regen_effects_from_draws.R \
    --run-dir "$RUN"
done
