#!/usr/bin/env bash

# Walk the 9 PPOM variants x 3 datasets and report which combos produced
# diagnose outputs. Writes a PASS/MISSING table to stdout. Exits non-zero
# if any combo is missing summary_by_k.csv.
#
# Usage:
#   bash report_diagnose_status.sh

set -o pipefail

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

n_pass=0
n_missing=0

printf "%-50s %-45s %-8s %s\n" "model_subdir" "dataset" "status" "summary_by_k.csv"
printf "%-50s %-45s %-8s %s\n" "------------" "-------" "------" "----------------"

for MODEL in "${MODELS[@]}"; do
  for DATASET in "${DATASETS[@]}"; do
    SUMMARY="$RESULTS_ROOT/$MODEL/$DATASET/ppc/diagnose/summary_by_k.csv"
    if [[ -f "$SUMMARY" ]]; then
      printf "%-50s %-45s %-8s %s\n" "$MODEL" "$DATASET" "PASS" "$SUMMARY"
      n_pass=$(( n_pass + 1 ))
    else
      printf "%-50s %-45s %-8s %s\n" "$MODEL" "$DATASET" "MISSING" "(absent)"
      n_missing=$(( n_missing + 1 ))
    fi
  done
done

echo
echo "summary: $n_pass pass, $n_missing missing (of $(( n_pass + n_missing )) combos)"

if [[ $n_missing -gt 0 ]]; then
  exit 1
fi
exit 0
