#!/usr/bin/env bash
#SBATCH --job-name=extract_cutpoint_draws
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=500G
#SBATCH --time=12:00:00
#SBATCH --error=/nfs/research/jlees/jacqueline/thesis_results/paper_figures/cutpoint_histograms/logs/extract_cutpoint_draws_%j.err
#SBATCH --output=/nfs/research/jlees/jacqueline/thesis_results/paper_figures/cutpoint_histograms/logs/extract_cutpoint_draws_%j.out

#################################################################################
# Loads each fitted-model RDS once on a big-mem node and caches just the cutpoint
# posterior draws to plots/cutpoints/cutpoint_draws.csv per run. Resume-safe: a
# run whose cache already exists is skipped, so re-running only fills in newly
# finished runs. ~600 GB of RDS read serially; peak RAM is one RDS at a time.
#################################################################################

source ~/.bashrc
mamba activate gwas_pipeline
set -euo pipefail

Rscript /nfs/research/jlees/jacqueline/thesis_code/paper_figures/extract_cutpoint_draws.R
