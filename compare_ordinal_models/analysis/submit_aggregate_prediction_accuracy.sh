#!/bin/bash
# Submit aggregate_prediction_accuracy.R over all PPOM_models/*_prediction model variants.
# Usage:
#   bash submit_aggregate_prediction_accuracy.sh                        # submit to SLURM
#   SLURM_JOB_ID=local bash submit_aggregate_prediction_accuracy.sh     # run inline (smoke test)

set -euo pipefail

REPO_ROOT=/nfs/research/jlees/jacqueline/thesis_code/compare_ordinal_models
PPOM_DIR="$REPO_ROOT/PPOM_models"
LOG_DIR=/nfs/research/jlees/jacqueline/thesis_results/compare_ordinal_models/_compare_PPOM/aggregate_logs
mkdir -p "$LOG_DIR"

# Discover *_prediction model names from PPOM_models (use .stan sources, strip extension).
mapfile -t MODELS < <(find "$PPOM_DIR" -maxdepth 1 -mindepth 1 -type f -name '*_prediction.stan' -printf '%f\n' | sed 's/\.stan$//' | sort)
if [ "${#MODELS[@]}" -eq 0 ]; then
    echo "ERROR: no *_prediction.stan files under $PPOM_DIR" >&2
    exit 1
fi
MODEL_CSV=$(IFS=,; echo "${MODELS[*]}")

echo "Found ${#MODELS[@]} _prediction model subdirs:"
printf '  %s\n' "${MODELS[@]}"

# If invoked outside SLURM, re-submit self via sbatch; otherwise run the body.
if [ -z "${SLURM_JOB_ID:-}" ]; then
    sbatch \
        --job-name=agg_pred_acc \
        --nodes=1 \
        --cpus-per-task=2 \
        --mem=16G \
        --time=1:00:00 \
        --error="$LOG_DIR/agg_pred_acc_%j.err" \
        --output="$LOG_DIR/agg_pred_acc_%j.out" \
        --export=ALL,MODEL_CSV="$MODEL_CSV" \
        "$0"
    exit 0
fi

# --- inside the SLURM job ---
# conda's activate/deactivate scripts reference unset env vars; relax `-u` for them.
set +u
source /hps/software/users/jlees/jacqueline/etc/profile.d/conda.sh
conda activate gwas_pipeline
set -u

cd "$REPO_ROOT"
echo "started      : $(date -Iseconds)"
echo "model_subdirs: $MODEL_CSV"

Rscript analysis/aggregate_prediction_accuracy.R \
    --model_subdirs "$MODEL_CSV"

echo "finished     : $(date -Iseconds)"
