#!/usr/bin/env bash

#SBATCH --job-name=gwas_inference_datasets
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=300G
#SBATCH --time=06:00:00
#SBATCH --error=/nfs/research/jlees/jacqueline/thesis_results/gwas_datasets/logs/inference.err
#SBATCH --output=/nfs/research/jlees/jacqueline/thesis_results/gwas_datasets/logs/inference.out

# Memory note: 300G is needed for the TB rifampicin presence/absence matrix (~1.7 GB TSV,
# expands considerably when loaded into R as a dense integer matrix).
# CPUs: data.table uses multiple threads for fread; 8 is sufficient.

mkdir -p /nfs/research/jlees/jacqueline/thesis_results/gwas_datasets/logs

Rscript /nfs/research/jlees/jacqueline/thesis_code/gwas_datasets/write_inference_jsons.R
