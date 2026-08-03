#!/usr/bin/env bash
# Determinism analysis (A3): replicate merges of one configuration.
#
# Two arms:
#   noise floor   n=10 at a fixed thread count -- how much do identical runs differ?
#   thread sweep  n=3 at 1 / 4 / 16 threads -- does the variation vanish at 1 thread?
#
# The thread sweep tests a specific hypothesis. find_mergeable_pairs sorts candidate pairs
# on (ident, sims[0], sims[1], sims[2]) only; equal-scoring pairs keep whatever order
# compute_scores_parallel returned them in, which depends on which worker finished first.
# The greedy seen_nodes walk then breaks ties differently between runs. If that is the
# cause, single-threaded runs should be identical to each other.
set -euo pipefail
cd "$(dirname "$0")"

SPECIES=${SPECIES:-s_pneumoniae}
SPLIT=${SPLIT:-split2}
OUT=/nfs/research/jlees/jacqueline/thesis_results/pangenomerge_benchmarking/determinism
mkdir -p "${OUT}/logs"

submit () {
  local rep=$1 threads=$2
  local tag=${SPECIES}_${SPLIT}_t${threads}_rep${rep}
  if [[ -f "${OUT}/${tag}/done.txt" ]]; then
    echo "skip (done): $tag"; return
  fi
  sbatch --job-name="det_${tag}" --cpus-per-task="$threads" \
    --output="${OUT}/logs/${tag}.out" --error="${OUT}/logs/${tag}.err" \
    run_replicate.sh "$SPECIES" "$SPLIT" "$rep" "$threads"
}

# noise floor: 10 identical runs
for rep in $(seq 1 10); do submit "$rep" 8; done

# thread sweep: 3 runs at each of 1, 4, 16
for threads in 1 4 16; do
  for rep in $(seq 1 3); do submit "$rep" "$threads"; done
done
