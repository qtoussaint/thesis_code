#!/usr/bin/env bash
set -euo pipefail

# Submits today's full PPOM shrinkage sweep: 14 models x 2 datasets = 28 jobs.
# DRY_RUN=1 ./submit_today_PPOM_sweep.sh   -> print sbatch lines, don't submit.

DRY_RUN="${DRY_RUN:-0}"

RUN_SCRIPTS_DIR="/nfs/research/jlees/jacqueline/thesis_code/compare_ordinal_models/run_PPOM_models"
PREFIX="final_ordered_categorical_PPOM_free_cutpoints_wide_drift"

MODELS=(
  # fixed tau x slab grid (8)
  fixedtau0p001_slab3
  fixedtau0p001_slab5
  fixedtau0p01_slab3
  fixedtau0p01_slab5
  fixedtau0p05_slab3
  fixedtau0p05_slab5
  fixedtau1_slab3
  fixedtau1_slab5
  # fixed tau x slab3 x lambda2 (4)
  fixedtau0p001_slab3_lambda2
  fixedtau0p01_slab3_lambda2
  fixedtau0p05_slab3_lambda2
  fixedtau1_slab3_lambda2
  # standalone variants (2)
  lambda1
  lognormaltau02
)

DATASETS=(
  02_spn_penicillin_MIC
  16_spn_penicillin_MIC_minimabinning
)

for suffix in "${MODELS[@]}"; do
  model="${PREFIX}_${suffix}"
  for ds in "${DATASETS[@]}"; do
    script="${RUN_SCRIPTS_DIR}/${model}/run_${ds}.sh"
    if [[ ! -f "$script" ]]; then
      echo "MISSING: $script" >&2
      continue
    fi
    if [[ "$DRY_RUN" = "1" ]]; then
      echo "DRY: sbatch $script"
    else
      echo "Submitting: $model / $ds"
      sbatch "$script"
    fi
  done
done
