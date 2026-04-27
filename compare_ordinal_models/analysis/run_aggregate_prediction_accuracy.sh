#!/bin/bash
set -eo pipefail

source /hps/software/users/jlees/jacqueline/etc/profile.d/conda.sh
conda activate gwas_pipeline

cd /nfs/research/jlees/jacqueline/thesis_code/compare_ordinal_models

Rscript analysis/aggregate_prediction_accuracy.R "$@"
