#!/usr/bin/env bash
# Submit every Panaroo job for the split-strain datasets (analysis A0).
#   ./submit_all.sh                    all species, all parts
#   ./submit_all.sh s_aureus           one species, all parts
#   ./submit_all.sh s_aureus quarter_1 one species, one part
#
# 3 species x 7 parts = 21 jobs. The "all" part is the ground-truth graph and is the
# expensive one (~350 isolates); quarters are ~88 isolates each.
set -euo pipefail
cd "$(dirname "$0")"
source ./config.sh

mkdir -p "${RESULTS_ROOT}/logs"

targets_species=("${SPECIES_LIST[@]}")
targets_parts=("${PARTS[@]}")
[[ $# -ge 1 ]] && targets_species=("$1")
[[ $# -ge 2 ]] && targets_parts=("$2")

# Resources scale with part size (~350 / ~175 / ~88 isolates). Sized from the measured
# quarter_1 run: 84 isolates took 50 s and 424 MB peak RSS. These allocations still leave
# 30-100x headroom; the original 300G/150G/100G guesses were ~250x over and would only
# have bought queue time.
mem_for () {
  case $1 in
    all)       echo 64G ;;
    half_*)    echo 32G ;;
    quarter_*) echo 16G ;;
  esac
}
cpus_for () {
  case $1 in
    all)     echo 16 ;;
    *)       echo 8 ;;
  esac
}

for sp in "${targets_species[@]}"; do
  for part in "${targets_parts[@]}"; do
    out=${RESULTS_ROOT}/${sp}/panaroo_${part}
    if [[ -f "${out}/done.txt" ]]; then
      echo "skip (done): $sp $part"
      continue
    fi
    sbatch \
      --job-name="pan_${sp}_${part}" \
      --mem="$(mem_for "$part")" \
      --cpus-per-task="$(cpus_for "$part")" \
      --output="${RESULTS_ROOT}/logs/panaroo_${sp}_${part}.out" \
      --error="${RESULTS_ROOT}/logs/panaroo_${sp}_${part}.err" \
      run_panaroo.sh "$sp" "$part"
  done
done
