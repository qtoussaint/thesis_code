#!/usr/bin/env bash
set -euo pipefail

# Submits the logistic alpha_prior_sd sweep: copies of gwas_finalmodels/logistic_inference.stan
# with alpha_prior_sd set to 0.075 / 1 / 1.25 / 1.5 / 2 / 2.5 / 3 / 5 (baseline model uses 0.5).
# 8 models x 1 dataset = 8 jobs.
# DRY_RUN=1 ./submit_logistic_sweep.sh   -> print sbatch lines, don't submit.

DRY_RUN="${DRY_RUN:-0}"

RUN_SCRIPTS_DIR="/nfs/research/jlees/jacqueline/thesis_code/compare_ordinal_models/run_logistic_models"

MODELS=(
  logistic_inference_alphasd0p075
  logistic_inference_alphasd1
  logistic_inference_alphasd1p25
  logistic_inference_alphasd1p5
  logistic_inference_alphasd2
  logistic_inference_alphasd2p5
  logistic_inference_alphasd3
  logistic_inference_alphasd5
)

DATASETS=(
  01_spn_penicillin_binary
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
