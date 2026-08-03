#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --mem=128G
#SBATCH --time=08:00:00
# One replicate merge for the determinism analysis (A3).
#   usage: sbatch run_replicate.sh <species> <split> <rep> <threads>
#
# Inputs, flags and code are identical across replicates -- only the run differs. Any
# variation in the output is therefore the tool's own non-determinism, which is the thing
# being measured.
set -euo pipefail

SS=/nfs/research/jlees/jacqueline/thesis_results/pangenomerge_benchmarking/split_strain_datasets
OUT=/nfs/research/jlees/jacqueline/thesis_results/pangenomerge_benchmarking/determinism
RUNNER=/hps/software/users/jlees/jacqueline/pangenome_merge/pangenomerge-runner.py

SPECIES=$1; SPLIT=$2; REP=$3; THREADS=$4
TAG=${SPECIES}_${SPLIT}_t${THREADS}_rep${REP}
OUTDIR=${OUT}/${TAG}

mkdir -p "$OUTDIR"
source /hps/software/users/jlees/jacqueline/etc/profile.d/conda.sh
conda activate /hps/software/users/jlees/jacqueline/envs/graph_merge

echo "[replicate] $TAG threads=$THREADS"

/usr/bin/time -v -o "${OUTDIR}/time.log" \
  python3 "$RUNNER" \
    --mode test \
    --outdir "$OUTDIR" \
    --component-graphs "${SS}/${SPECIES}/paths/${SPLIT}.txt" \
    --graph-all "${SS}/${SPECIES}/panaroo_all" \
    --threads "$THREADS" \
  > "${OUTDIR}/run.log" 2>&1

test -f "${OUTDIR}/final_graph.gml"
echo "ok" > "${OUTDIR}/done.txt"
echo "[replicate] complete: $OUTDIR"
