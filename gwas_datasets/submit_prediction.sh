#!/usr/bin/env bash

#SBATCH --job-name=gwas_prediction_datasets
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=300G
#SBATCH --time=06:00:00
#SBATCH --error=/nfs/research/jlees/jacqueline/thesis_results/gwas_datasets/logs/prediction.err
#SBATCH --output=/nfs/research/jlees/jacqueline/thesis_results/gwas_datasets/logs/prediction.out

mkdir -p /nfs/research/jlees/jacqueline/thesis_results/gwas_datasets/logs

Rscript /nfs/research/jlees/jacqueline/thesis_code/gwas_datasets/write_prediction_jsons.R
