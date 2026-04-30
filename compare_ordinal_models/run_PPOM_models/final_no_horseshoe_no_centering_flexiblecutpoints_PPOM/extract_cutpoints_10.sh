#!/usr/bin/env bash

#SBATCH --job-name=cutpoint_extract_noHSflex_10
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=96G
#SBATCH --time=1:00:00
#SBATCH --error=/nfs/research/jlees/jacqueline/thesis_results/compare_ordinal_models/final_no_horseshoe_no_centering_flexiblecutpoints_PPOM/10_spn_penicillin_MIC_coarse_dilutions/logs/cutpoint_extract.err
#SBATCH --output=/nfs/research/jlees/jacqueline/thesis_results/compare_ordinal_models/final_no_horseshoe_no_centering_flexiblecutpoints_PPOM/10_spn_penicillin_MIC_coarse_dilutions/logs/cutpoint_extract.out

source ~/.bashrc
mamba activate gwas_pipeline

Rscript /nfs/research/jlees/jacqueline/thesis_code/compare_ordinal_models/analysis/extract_cutpoint_posterior_no_hs.R \
  final_no_horseshoe_no_centering_flexiblecutpoints_PPOM \
  10_spn_penicillin_MIC_coarse_dilutions
