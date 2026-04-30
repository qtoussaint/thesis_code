#!/usr/bin/env bash

#SBATCH --job-name=gwas_inference_16_minimabinning
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --error=/nfs/research/jlees/jacqueline/thesis_results/gwas_datasets/logs/inference_16.err
#SBATCH --output=/nfs/research/jlees/jacqueline/thesis_results/gwas_datasets/logs/inference_16.out

# Dataset 16 only (SPN penicillin MIC, breakpoint minima binning).
# Skips TB loading, so memory/time are well below the full pipeline.

mkdir -p /nfs/research/jlees/jacqueline/thesis_results/gwas_datasets/logs

Rscript /nfs/research/jlees/jacqueline/thesis_code/gwas_datasets/write_inference_jsons_16_only.R
