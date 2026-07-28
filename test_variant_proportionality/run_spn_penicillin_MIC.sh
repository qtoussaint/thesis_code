#!/usr/bin/env bash
#SBATCH --job-name=spnpen_proportionality
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=600G
#SBATCH --time=12:00:00
#SBATCH --error=/nfs/research/jlees/jacqueline/thesis_results/test_variant_proportionality/spn_penicillin_MIC/logs/spn_penicillin_MIC_%j.err
#SBATCH --output=/nfs/research/jlees/jacqueline/thesis_results/test_variant_proportionality/spn_penicillin_MIC/logs/spn_penicillin_MIC_%j.out

#################################################################################
# POM vs PPOM proportional-odds comparison for SPN penicillin MIC.
#
# Loads the fitted-model RDS files (POM ~13G, PPOM ~92G), so this is run as a
# high-memory job. Outputs land under $OUT_DIR.

source ~/.bashrc
mamba activate gwas_pipeline

set -euo pipefail

POM_DIR="/nfs/research/jlees/jacqueline/thesis_results/gwas_spn_penicillin/inference/02_spn_penicillin_MIC_POM"
PPOM_DIR="/nfs/research/jlees/jacqueline/thesis_results/gwas_spn_penicillin/inference/02_spn_penicillin_MIC_PPOM"
OUT_DIR="/nfs/research/jlees/jacqueline/thesis_results/test_variant_proportionality/spn_penicillin_MIC"
SCRIPT_DIR="/nfs/research/jlees/jacqueline/thesis_code/test_variant_proportionality"

mkdir -p "$OUT_DIR/logs"

Rscript "$SCRIPT_DIR/test_variant_proportionality.R" \
    --pom_dir     "$POM_DIR" \
    --ppom_dir    "$PPOM_DIR" \
    --out_dir     "$OUT_DIR" \
    --ci_level    0.89 \
    --top_n_plots 20
