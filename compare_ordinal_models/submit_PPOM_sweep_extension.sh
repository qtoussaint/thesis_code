#!/usr/bin/env bash
set -euo pipefail

# Submits the PPOM sweep extension: high-tau horseshoe grid (lambda2) +
# fixed-scale lasso/ridge + estimated-scale (adaptive) lasso/ridge.
# 20 models x 2 datasets = 40 jobs.
# DRY_RUN=1 ./submit_PPOM_sweep_extension.sh   -> print sbatch lines, don't submit.

DRY_RUN="${DRY_RUN:-0}"

RUN_SCRIPTS_DIR="/nfs/research/jlees/jacqueline/thesis_code/compare_ordinal_models/run_PPOM_models"
PREFIX="final_ordered_categorical_PPOM_free_cutpoints_wide_drift"

MODELS=(
  # high-tau x slab grid, lambda~Cauchy(0,2) (8)
  fixedtau1p5_slab3_lambda2
  fixedtau1p5_slab5_lambda2
  fixedtau2_slab3_lambda2
  fixedtau2_slab5_lambda2
  fixedtau3_slab3_lambda2
  fixedtau3_slab5_lambda2
  fixedtau5_slab3_lambda2
  fixedtau5_slab5_lambda2
  # fixed-scale prior comparators (4): standardized (X_std) and raw-genotype (no_centering)
  lasso
  ridge
  lasso_no_centering
  ridge_no_centering
  # estimated-scale (adaptive) prior comparators (8): centered + non-centered/mixture
  # parameterizations of ridge & lasso, each X_std and raw-genotype (no_centering)
  ridge_estscale
  ridge_estscale_no_centering
  ridge_estscale_ncp
  ridge_estscale_ncp_no_centering
  lasso_estscale
  lasso_estscale_no_centering
  lasso_estscale_mixture
  lasso_estscale_mixture_no_centering
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
