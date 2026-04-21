#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=150G
#SBATCH --time=05:00:00

set -euo pipefail

ITER=${1:--1}   # iteration count from $1, or -1 (unlimited) by default
TAG=${ITER/-1/unlimited}
RESULTS_ROOT=/nfs/research/jlees/jacqueline/thesis_results/pangenomerge_benchmarking/context_search_iterations
OUTDIR=${RESULTS_ROOT}/iter_${TAG}
mkdir -p "$OUTDIR"

cd /hps/software/users/jlees/jacqueline/pangenome_merge

/usr/bin/time -v -o "${OUTDIR}/time.log" \
    python -m pangenomerge \
        --mode run \
        --component-graphs /nfs/research/jlees/jacqueline/atb_analyses/neisseria_bf_data/pangenomes.tsv \
        --outdir "$OUTDIR" \
        --threads 48 \
        --context-search-iterations ${ITER}
