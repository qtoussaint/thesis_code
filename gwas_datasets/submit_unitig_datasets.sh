#!/usr/bin/env bash

#SBATCH --job-name=gwas_unitig_datasets
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=400G
#SBATCH --time=06:00:00
#SBATCH --error=/nfs/research/jlees/jacqueline/thesis_results/gwas_datasets/logs/unitig_datasets.err
#SBATCH --output=/nfs/research/jlees/jacqueline/thesis_results/gwas_datasets/logs/unitig_datasets.out

# Builds ONLY the unitig-genotype SPN penicillin datasets (17, 18; inference +
# prediction). 400G is needed because the full unitig rtab (~900 MB, ~740k x 616)
# expands to a dense integer matrix when read into R before MAF filtering.
# CPUs: data.table uses multiple threads for fread; 8 is sufficient.

mkdir -p /nfs/research/jlees/jacqueline/thesis_results/gwas_datasets/logs

Rscript /nfs/research/jlees/jacqueline/thesis_code/gwas_datasets/write_unitig_datasets.R
