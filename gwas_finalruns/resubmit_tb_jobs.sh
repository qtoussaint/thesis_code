#!/usr/bin/env bash
# Resubmit every TB rifampicin SLURM job:
#   gwas_finalruns/inference/tb_rifampicin/*.sh
#   gwas_finalruns/prediction/tb_rifampicin/*.sh
#
# Default is dry-run. Pass --execute to actually sbatch.

set -euo pipefail

CODE_ROOT="/nfs/research/jlees/jacqueline/thesis_code/gwas_finalruns"

EXECUTE=0
if [[ "${1:-}" == "--execute" ]]; then
  EXECUTE=1
fi

if [[ $EXECUTE -eq 0 ]]; then
  echo "DRY RUN — pass --execute to actually sbatch."
fi

submit() {
  local sh="$1"
  echo "  sbatch: $sh"
  if [[ $EXECUTE -eq 1 ]]; then
    sbatch "$sh"
  fi
}

echo "=== inference ==="
for sh in "$CODE_ROOT/inference/tb_rifampicin/"*.sh; do
  submit "$sh"
done

echo "=== prediction ==="
for sh in "$CODE_ROOT/prediction/tb_rifampicin/"*.sh; do
  submit "$sh"
done

if [[ $EXECUTE -eq 0 ]]; then
  echo "DRY RUN complete. Re-run with --execute to submit."
else
  echo "Submission complete. Check queue with: squeue -u $USER"
fi
