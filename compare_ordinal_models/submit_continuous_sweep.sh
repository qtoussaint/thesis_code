#!/usr/bin/env bash
set -euo pipefail

# Submits the continuous (log2-MIC) alpha_prior_sd sweep: copies of
# gwas_finalmodels/continuous_inference.stan with alpha_prior_sd widened to
# 1.5 / 3 / 5 (baseline model uses 0.5).
# 3 models x 1 dataset = 3 jobs.
# DRY_RUN=1 ./submit_continuous_sweep.sh   -> print sbatch lines, don't submit.

DRY_RUN="${DRY_RUN:-0}"

RUN_SCRIPTS_DIR="/nfs/research/jlees/jacqueline/thesis_code/compare_ordinal_models/run_continuous_models"

MODELS=(
  continuous_inference_alphasd1p5
  continuous_inference_alphasd3
  continuous_inference_alphasd5
)

DATASETS=(
  03_spn_penicillin_continuous
)

for model in "${MODELS[@]}"; do
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
