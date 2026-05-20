#!/usr/bin/env bash
# Deletes the cppRATE_results/ subdir from every model output directory
# referenced by the SLURM scripts under gwas_finalruns/{inference,prediction}/<species>/.
# Leaves cppRATE_matrices/ (the expensive intermediate CSVs) untouched.
#
# Usage:
#   bash delete_cpprate_results.sh            # delete everywhere
#   bash delete_cpprate_results.sh --dry-run  # print what would be deleted
#   bash delete_cpprate_results.sh inference/spn_penicillin  # limit to one subdir

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

deleted=0
skipped_missing=0
skipped_no_outdir=0
would_delete=0

for d in "${SEARCH_DIRS[@]}"; do
  if [[ ! -d "$d" ]]; then
    echo "skip: $d (not a directory)" >&2
    continue
  fi
  for f in "$d"/*.sh; do
    [[ -e "$f" ]] || continue

    # Pull `--output_directory <path>` out of the SLURM script. The path is
    # often inside a shell-quoted assignment (OUTPUT_DIR="--output_directory /…/dir"),
    # so strip any leading/trailing single or double quotes after extraction.
    outdir="$(grep -oE -- '--output_directory[[:space:]]+[^[:space:]]+' "$f" | head -n1 | awk '{print $2}')"
    outdir="${outdir%[\"\']}"
    outdir="${outdir#[\"\']}"
    if [[ -z "${outdir}" ]]; then
      echo "warn: no --output_directory found in $f" >&2
      skipped_no_outdir=$((skipped_no_outdir + 1))
      continue
    fi

    target="${outdir%/}/cppRATE_results"

    # Safety guard: refuse anything whose basename isn't exactly cppRATE_results.
    if [[ "$(basename "$target")" != "cppRATE_results" ]]; then
      echo "refuse: $target (basename guard)" >&2
      continue
    fi

    if [[ ! -d "$target" ]]; then
      echo "skip: $target (already gone)"
      skipped_missing=$((skipped_missing + 1))
      continue
    fi

    if (( DRY )); then
      echo "would rm -rf: $target"
      would_delete=$((would_delete + 1))
    else
      rm -rf "$target"
      echo "deleted: $target"
      deleted=$((deleted + 1))
    fi
  done
done

if (( DRY )); then
  echo "dry-run: ${would_delete} dir(s) would be deleted, ${skipped_missing} already gone, ${skipped_no_outdir} script(s) without --output_directory."
else
  echo "${deleted} dir(s) deleted, ${skipped_missing} already gone, ${skipped_no_outdir} script(s) without --output_directory."
fi
