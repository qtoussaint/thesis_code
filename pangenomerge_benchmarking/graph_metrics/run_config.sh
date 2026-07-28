#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=08:00:00
# Merge one split-strain configuration and score it (analyses A2d / A7 input).
#   usage: sbatch run_config.sh <species> <split2|split4>
#
# Two steps in one job because the second consumes the first's output:
#   1. pangenomerge --mode test   -> merged_<split>/final_graph.gml (retains seqIDs)
#   2. validate_metrics.py        -> metric_validation_<config>.csv
#
# Memory-led sizing: test mode forces metadata into the graph, and an earlier run of the
# 200-isolate legacy dataset was OOM-killed at 31 GB. These are ~350 isolates.
set -euo pipefail

CODE=/nfs/research/jlees/jacqueline/thesis_code/pangenomerge_benchmarking/graph_metrics
SS=/nfs/research/jlees/jacqueline/thesis_results/pangenomerge_benchmarking/split_strain_datasets
OUT=/nfs/research/jlees/jacqueline/thesis_results/pangenomerge_benchmarking/graph_metrics
RUNNER=/hps/software/users/jlees/jacqueline/pangenome_merge/pangenomerge-runner.py

SPECIES=$1
SPLIT=$2
CONFIG=${SPECIES}_${SPLIT}
THREADS=${SLURM_CPUS_PER_TASK:-8}

PATHS=${SS}/${SPECIES}/paths/${SPLIT}.txt
GRAPH_ALL=${SS}/${SPECIES}/panaroo_all
MERGEDIR=${SS}/${SPECIES}/merged_${SPLIT}

mkdir -p "$MERGEDIR" "$OUT"

source /hps/software/users/jlees/jacqueline/etc/profile.d/conda.sh
conda activate /hps/software/users/jlees/jacqueline/envs/graph_merge

echo "[merge] config=$CONFIG components=$(wc -l < "$PATHS") threads=$THREADS"

/usr/bin/time -v -o "${MERGEDIR}/time.log" \
  python3 "$RUNNER" \
    --mode test \
    --outdir "$MERGEDIR" \
    --component-graphs "$PATHS" \
    --graph-all "$GRAPH_ALL" \
    --threads "$THREADS"

test -f "${MERGEDIR}/final_graph.gml"
echo "[merge] complete: $MERGEDIR"

mapfile -t COMPONENTS < "$PATHS"

echo "[validate] scoring $CONFIG"
/usr/bin/time -v -o "${OUT}/${CONFIG}_time.log" \
  python3 "${CODE}/validate_metrics.py" \
    --graph-all "$GRAPH_ALL" \
    --component-dirs "${COMPONENTS[@]}" \
    --merged-gml "${MERGEDIR}/final_graph.gml" \
    --label "$CONFIG" \
    --out-csv "${OUT}/metric_validation_${CONFIG}.csv"

echo "[validate] complete: ${OUT}/metric_validation_${CONFIG}.csv"
