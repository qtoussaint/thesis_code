#!/usr/bin/env bash

#SBATCH --job-name=effects_extract_noHSwCent_16
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=400G
#SBATCH --time=4:00:00
#SBATCH --error=/nfs/research/jlees/jacqueline/thesis_results/compare_ordinal_models/final_no_horseshoe_with_centering_PPOM/16_spn_penicillin_MIC_minimabinning/logs/effects_extract.err
#SBATCH --output=/nfs/research/jlees/jacqueline/thesis_results/compare_ordinal_models/final_no_horseshoe_with_centering_PPOM/16_spn_penicillin_MIC_minimabinning/logs/effects_extract.out

source ~/.bashrc
mamba activate gwas_pipeline

set -euo pipefail

Rscript /nfs/research/jlees/jacqueline/thesis_code/compare_ordinal_models/analysis/extract_depruned_effects_ppom.R \
  final_no_horseshoe_with_centering_PPOM \
  16_spn_penicillin_MIC_minimabinning
