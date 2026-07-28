#!/usr/bin/env bash
# Merge + score all six split-strain configurations (3 species x {split2, split4}).
#   ./submit_sweep.sh
# The split4 runs are the first exercise of the iterative multi-graph merge path; split2
# repeats the 2-graph case the original evaluation used. Because quarters nest inside
# halves, the two rows for a species contain identical isolates.
set -euo pipefail
cd "$(dirname "$0")"

OUT=/nfs/research/jlees/jacqueline/thesis_results/pangenomerge_benchmarking/graph_metrics
mkdir -p "${OUT}/logs"

for sp in s_pneumoniae m_tuberculosis s_aureus; do
  for split in split2 split4; do
    cfg=${sp}_${split}
    if [[ -f "${OUT}/metric_validation_${cfg}.csv" ]]; then
      echo "skip (done): $cfg"
      continue
    fi
    sbatch \
      --job-name="metric_${cfg}" \
      --output="${OUT}/logs/${cfg}.out" \
      --error="${OUT}/logs/${cfg}.err" \
      run_config.sh "$sp" "$split"
  done
done
