#!/usr/bin/env bash
set -euo pipefail
RESULTS_ROOT=/nfs/research/jlees/jacqueline/thesis_results/pangenomerge_benchmarking/context_search_iterations
mkdir -p "${RESULTS_ROOT}/logs"
cd "$(dirname "$0")"
for ITER in 1 2 3 5 10 -1; do
    TAG=${ITER/-1/unlimited}
    sbatch \
        --job-name="ctxiter_${TAG}" \
        --output="${RESULTS_ROOT}/logs/iter_${TAG}.out" \
        --error="${RESULTS_ROOT}/logs/iter_${TAG}.err" \
        run_single.sh "$ITER"
done
