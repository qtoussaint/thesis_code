#!/usr/bin/env bash

#SBATCH --job-name=diagnose_PPOM
#SBATCH --array=0-26%27
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=350G
#SBATCH --time=02:00:00
#SBATCH --error=/nfs/research/jlees/jacqueline/thesis_results/compare_ordinal_models/diagnose_logs/array_%A_%a.err
#SBATCH --output=/nfs/research/jlees/jacqueline/thesis_results/compare_ordinal_models/diagnose_logs/array_%A_%a.out

# Run diagnose_ppom_cutpoint_inflation.R across all 9 PPOM variants x 3
# datasets (27 tasks, indices 0..26). One SLURM array task per combo;
# task failures are independent.
#
# Usage:
#   bash run_diagnose_all_PPOMs.sh             # self-submits as sbatch array
#   bash run_diagnose_all_PPOMs.sh --print     # print the 15 (model,dataset) pairs and exit
#   SLURM_ARRAY_TASK_ID=0 bash run_diagnose_all_PPOMs.sh
#                                              # run a single task locally (smoke test)

set -o pipefail

MODELS=(
  final_ordered_categorical_PPOM_tight_alpha_tau1
  final_ordered_categorical_PPOM_tight_slab
  final_ordered_categorical_PPOM_poolk
  final_ordered_categorical_PPOM_latent_scale
  final_ordered_categorical_PPOM_free_cutpoints
  final_ordered_categorical_PPOM_free_cutpoints_wide_drift
  final_ordered_categorical_PPOM_free_cutpoints_slab50
  final_ordered_categorical_PPOM_free_cutpoints_wide_drift_slab50
  final_ordered_categorical_PPOM_slab50
)
DATASETS=(
  02_spn_penicillin_MIC
  10_spn_penicillin_MIC_coarse_dilutions
  11_spn_penicillin_MIC_large_minbin
)

N_MODELS=${#MODELS[@]}
N_DATASETS=${#DATASETS[@]}
N_COMBOS=$(( N_MODELS * N_DATASETS ))

LOG_DIR=/nfs/research/jlees/jacqueline/thesis_results/compare_ordinal_models/diagnose_logs
SCRIPT=/nfs/research/jlees/jacqueline/thesis_code/compare_ordinal_models/analysis/diagnose_ppom_cutpoint_inflation.R

if [[ "${1:-}" == "--print" ]]; then
  printf "%-4s %-50s %s\n" "idx" "model_subdir" "dataset"
  for idx in $(seq 0 $(( N_COMBOS - 1 ))); do
    m=$(( idx / N_DATASETS ))
    d=$(( idx % N_DATASETS ))
    printf "%-4d %-50s %s\n" "$idx" "${MODELS[$m]}" "${DATASETS[$d]}"
  done
  exit 0
fi

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  mkdir -p "$LOG_DIR"
  echo "submitting SLURM array (0-$(( N_COMBOS - 1 ))), $N_COMBOS tasks"
  exec sbatch "$0"
fi

source /hps/software/users/jlees/jacqueline/etc/profile.d/conda.sh
conda activate gwas_pipeline

idx=$SLURM_ARRAY_TASK_ID
m_idx=$(( idx / N_DATASETS ))
d_idx=$(( idx % N_DATASETS ))
MODEL=${MODELS[$m_idx]}
DATASET=${DATASETS[$d_idx]}

echo "=========================================================="
echo "task        : $idx / $(( N_COMBOS - 1 ))"
echo "model_subdir: $MODEL"
echo "dataset     : $DATASET"
echo "host        : $(hostname)"
echo "started     : $(date -Iseconds)"
echo "=========================================================="

Rscript "$SCRIPT" --model_subdir "$MODEL" --dataset "$DATASET"
status=$?

echo "=========================================================="
echo "finished    : $(date -Iseconds)"
echo "exit status : $status"
echo "=========================================================="

exit $status
