#!/usr/bin/env bash
# Submits every generated SLURM script under gwas_finalruns/{inference,prediction}/<species>/.
# Usage:
#   bash submit_all.sh            # submit every .sh
#   bash submit_all.sh --dry-run  # list what would be submitted
#   bash submit_all.sh inference/spn_penicillin  # limit to one subdir

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

DRY=0
if [[ "${1-}" == "--dry-run" ]]; then
  DRY=1
  shift
fi

if [[ $# -gt 0 ]]; then
  SEARCH_DIRS=()
  for d in "$@"; do
    SEARCH_DIRS+=("$ROOT/$d")
  done
else
  SEARCH_DIRS=(
    "$ROOT/inference/spn_penicillin"
    "$ROOT/inference/spn_trimethoprim"
    "$ROOT/inference/tb_rifampicin"
    "$ROOT/prediction/spn_penicillin"
    "$ROOT/prediction/spn_trimethoprim"
    "$ROOT/prediction/tb_rifampicin"
  )
fi

count=0
for d in "${SEARCH_DIRS[@]}"; do
  if [[ ! -d "$d" ]]; then
    echo "skip: $d (not a directory)" >&2
    continue
  fi
  for f in "$d"/*.sh; do
    [[ -e "$f" ]] || continue
    if (( DRY )); then
      echo "would sbatch: $f"
    else
      sbatch "$f"
    fi
    count=$((count + 1))
  done
done

echo "${count} script(s) processed (dry-run=${DRY})."
