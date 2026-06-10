#!/usr/bin/env bash
#SBATCH --job-name=extract_cpt
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --error=/nfs/research/jlees/jacqueline/thesis_results/paper_figures/cutpoint_histograms/logs/extract_cpt_%A_%a.err
#SBATCH --output=/nfs/research/jlees/jacqueline/thesis_results/paper_figures/cutpoint_histograms/logs/extract_cpt_%A_%a.out

#################################################################################
# One array task per inference run: loads that run's fitted RDS, caches just the
# cutpoint draws to plots/cutpoints/cutpoint_draws.csv. SLURM_ARRAY_TASK_ID is
# the --index into the run manifest (extract_cutpoint_draws.R --list). Resume-safe.
# --array, --mem and --time are set on the sbatch command line so light and heavy
# RDS tiers can be submitted separately.
#################################################################################

source ~/.bashrc
mamba activate gwas_pipeline
set -euo pipefail

Rscript /nfs/research/jlees/jacqueline/thesis_code/paper_figures/extract_cutpoint_draws.R \
  --index "${SLURM_ARRAY_TASK_ID}"
