#!/usr/bin/env bash

#SBATCH --job-name=kleb_inference_datasets
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=08:00:00
#SBATCH --error=/nfs/research/jlees/jacqueline/thesis_results/gwas_datasets/logs/klebsiella.err
#SBATCH --output=/nfs/research/jlees/jacqueline/thesis_results/gwas_datasets/logs/klebsiella.out

# Memory note: the Klebsiella genotype is 12,410 x 4,368 (~108 MB on disk, ~217 MB
# as a dense R integer matrix). Four datasets are built in one pass off a single
# in-memory copy, and each write_stan_json_streaming call holds one transposed
# sample-major copy at a time. 200G is comfortable; the spn/tb equivalent needs
# more only because of the much larger TB presence/absence matrix.

mkdir -p /nfs/research/jlees/jacqueline/thesis_results/gwas_datasets/logs

source ~/.bashrc
mamba activate gwas_pipeline

Rscript /nfs/research/jlees/jacqueline/thesis_code/gwas_datasets/write_klebsiella_jsons.R
