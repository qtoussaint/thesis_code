#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=04:00:00
# Metric calibration for one dataset configuration (analysis A2d).
#   usage: sbatch run_validate.sh <config>
#
# Configs are <species>_split2 / <species>_split4 for the datasets built by A0, plus
# "legacy_c25" for the original S. aureus cluster-25 half-split, which is the only one with
# a merged graph until A7 runs.
#
# Scoring is memory-hungry rather than CPU-hungry: adjusted_mutual_info_score and the
# contingency tables are built over ~500k genes for each of six scenarios.
set -euo pipefail

CODE=/nfs/research/jlees/jacqueline/thesis_code/pangenomerge_benchmarking/graph_metrics
SS=/nfs/research/jlees/jacqueline/thesis_results/pangenomerge_benchmarking/split_strain_datasets
LEGACY=/nfs/research/jlees/jacqueline/atb_analyses/pangenomerge_split_strain
OUT=/nfs/research/jlees/jacqueline/thesis_results/pangenomerge_benchmarking/graph_metrics

CONFIG=$1
mkdir -p "$OUT"

source /hps/software/users/jlees/jacqueline/etc/profile.d/conda.sh
conda activate /hps/software/users/jlees/jacqueline/envs/graph_merge

case "$CONFIG" in
  legacy_c25)
    GRAPH_ALL=${LEGACY}/panaroo_25all
    COMPONENTS=("${LEGACY}/panaroo_25a" "${LEGACY}/panaroo_25b")
    MERGED=${LEGACY}/method_eval/baseline/final_graph.gml
    ;;
  *_split2|*_split4)
    SPECIES=${CONFIG%_split*}
    SPLIT=${CONFIG##*_}                       # split2 | split4
    GRAPH_ALL=${SS}/${SPECIES}/panaroo_all
    mapfile -t COMPONENTS < "${SS}/${SPECIES}/paths/${SPLIT}.txt"
    MERGED=${SS}/${SPECIES}/merged_${SPLIT}/final_graph.gml
    ;;
  *)
    echo "unknown config: $CONFIG" >&2; exit 1 ;;
esac

if [[ ! -f "$MERGED" ]]; then
  echo "no merged graph yet for ${CONFIG}: ${MERGED}" >&2
  echo "run pangenomerge --mode test for this configuration first" >&2
  exit 1
fi

echo "[validate] config=$CONFIG components=${#COMPONENTS[@]}"

/usr/bin/time -v -o "${OUT}/${CONFIG}_time.log" \
  python3 "${CODE}/validate_metrics.py" \
    --graph-all "$GRAPH_ALL" \
    --component-dirs "${COMPONENTS[@]}" \
    --merged-gml "$MERGED" \
    --label "$CONFIG" \
    --out-csv "${OUT}/metric_validation_${CONFIG}.csv"

echo "[validate] complete: ${OUT}/metric_validation_${CONFIG}.csv"
