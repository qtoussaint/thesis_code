#!/usr/bin/env bash
# Delete all plots from every TB rifampicin run.
#
# Runs are enumerated from the submit-script filenames in
#   gwas_finalruns/inference/tb_rifampicin/   ->  thesis_results/.../inference/<run>/
#   gwas_finalruns/prediction/tb_rifampicin/  ->  thesis_results/.../prediction/<run>/
#
# Plots targeted per run:
#   inference/<run>/plots/                  (whole tree)
#   inference/<run>/inference_ppc/*.png
#   prediction/<run>/prediction_results/*.png
#
# Default is dry-run. Pass --execute to actually delete.

set -euo pipefail

CODE_ROOT="/nfs/research/jlees/jacqueline/thesis_code/gwas_finalruns"
RESULTS_ROOT="/nfs/research/jlees/jacqueline/thesis_results/gwas_tb_rifampicin"

EXECUTE=0
if [[ "${1:-}" == "--execute" ]]; then
  EXECUTE=1
fi

if [[ $EXECUTE -eq 0 ]]; then
  echo "DRY RUN — pass --execute to actually delete."
fi

remove_path() {
  local p="$1"
  [[ -e "$p" ]] || return 0
  echo "  rm: $p"
  if [[ $EXECUTE -eq 1 ]]; then
    rm -rf -- "$p"
  fi
}

# Inference runs
for sh in "$CODE_ROOT/inference/tb_rifampicin/"*.sh; do
  run=$(basename "$sh" .sh)
  run_dir="$RESULTS_ROOT/inference/$run"
  if [[ ! -d "$run_dir" ]]; then
    echo "skip (no results dir): $run_dir"
    continue
  fi
  echo "inference run: $run"
  remove_path "$run_dir/plots"
  shopt -s nullglob
  for png in "$run_dir/inference_ppc/"*.png; do
    remove_path "$png"
  done
  shopt -u nullglob
done

# Prediction runs
for sh in "$CODE_ROOT/prediction/tb_rifampicin/"*.sh; do
  run=$(basename "$sh" .sh)
  run_dir="$RESULTS_ROOT/prediction/$run"
  if [[ ! -d "$run_dir" ]]; then
    echo "skip (no results dir): $run_dir"
    continue
  fi
  echo "prediction run: $run"
  shopt -s nullglob
  for png in "$run_dir/prediction_results/"*.png; do
    remove_path "$png"
  done
  shopt -u nullglob
done

if [[ $EXECUTE -eq 0 ]]; then
  echo "DRY RUN complete. Re-run with --execute to delete."
else
  echo "Deletion complete."
fi
