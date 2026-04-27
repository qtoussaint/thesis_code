#!/usr/bin/env bash
# Run compare_ppom_variants.R once per dataset, restricted to PPOM variants
# whose diagnose outputs (summary_by_k.csv) actually exist. Mirrors the model
# and dataset arrays in report_diagnose_status.sh.
#
# Usage:
#   bash run_compare_ppom_variants_all.sh

source /hps/software/users/jlees/jacqueline/etc/profile.d/conda.sh
conda activate gwas_pipeline

set -euo pipefail

ANALYSIS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RESULTS_ROOT=/nfs/research/jlees/jacqueline/thesis_results/compare_ordinal_models

MODELS=(
  final_ordered_categorical_PPOM_tight_alpha_tau1
  final_ordered_categorical_PPOM_tight_slab
  final_ordered_categorical_PPOM_poolk
  final_ordered_categorical_PPOM_latent_scale
  final_ordered_categorical_PPOM_free_cutpoints
  final_ordered_categorical_PPOM_free_cutpoints_wide_drift
  final_ordered_categorical_PPOM_free_cutpoints_slab50
  final_ordered_categorical_PPOM_free_cutpoints_wide_drift_slab50
  final_ordered_categorical_PPOM_slab50
)
DATASETS=(
  02_spn_penicillin_MIC
  10_spn_penicillin_MIC_coarse_dilutions
  11_spn_penicillin_MIC_large_minbin
)

for D in "${DATASETS[@]}"; do
  finished=()
  for M in "${MODELS[@]}"; do
    if [[ -f "$RESULTS_ROOT/$M/$D/ppc/diagnose/summary_by_k.csv" ]]; then
      finished+=("$M")
    fi
  done
  if [[ ${#finished[@]} -eq 0 ]]; then
    echo "[skip] $D: no finished variants"
    continue
  fi
  csv=$(IFS=,; echo "${finished[*]}")
  echo "[run]  $D: ${#finished[@]} variants"
  Rscript "$ANALYSIS_DIR/compare_ppom_variants.R" \
    --dataset "$D" \
    --model_subdirs "$csv"
done
