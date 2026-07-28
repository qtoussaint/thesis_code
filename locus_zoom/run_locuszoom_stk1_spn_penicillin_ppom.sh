#!/usr/bin/env bash

#SBATCH --job-name=locuszoom_stk1_spn_pen_ppom
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=1:00:00
#SBATCH --error=/nfs/research/jlees/jacqueline/thesis_code/locus_zoom/locuszoom_stk1_ppom.err
#SBATCH --output=/nfs/research/jlees/jacqueline/thesis_code/locus_zoom/locuszoom_stk1_ppom.out

#################################################################################
# LocusZoom plot of the stk1 region in the spn_penicillin
# 02_spn_penicillin_MIC_PPOM inference run, RATE metric.
#
# Same genomic window as the 01_binary_logistic / 02_POM stk1 plots
# (FM211187.1:1,629,135-1,729,135, lead 1,679,135 +/- 50 kb) so the runs are
# directly comparable. PPOM has per-cutpoint effects, so all cutpoints are
# overlaid on one plot (colour = r^2 with lead, shape = cutpoint).
#################################################################################

source ~/.bashrc
mamba activate gwas_pipeline

MAKE_PLOT="/nfs/research/jlees/jacqueline/thesis_code/locus_zoom/make_locuszoom_plot.R"

# ---------------------------------------------------------------------------
# Reference (S. pneumoniae ATCC 700669, GCA_000026665.1_ASM2666v1, seqname FM211187.1)
# ---------------------------------------------------------------------------
SPECIES_OUTPUT_DIR="/nfs/research/jlees/jacqueline/thesis_results/locus_zoom/spneumoniae"
GFF="$SPECIES_OUTPUT_DIR/reference/GCA_000026665.1_ASM2666v1_genomic.gff"

# ---------------------------------------------------------------------------
# Pipeline run
# ---------------------------------------------------------------------------
PIPELINE_OUTPUT_DIR="/nfs/research/jlees/jacqueline/thesis_results/gwas_spn_penicillin/inference/02_spn_penicillin_MIC_PPOM"
DATASET_DIR="/nfs/research/jlees/jacqueline/thesis_results/gwas_datasets/inference/02_spn_penicillin_MIC"

POSITIONS_FILE="$DATASET_DIR/02_spn_penicillin_MIC_variant_index.csv"
GENOTYPE_MATRIX="$PIPELINE_OUTPUT_DIR/cppRATE_matrices/design_matrix.csv"
VARIANT_EFFECTS="$PIPELINE_OUTPUT_DIR/fitted_model/depruned_variant_effects.csv"
ANNOTATIONS="/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/genotype/fields_filtered_maf05_multiallelic.txt"
RATE_DIR="$PIPELINE_OUTPUT_DIR/cppRATE_results"
GENES_OF_INTEREST="/nfs/research/jlees/jacqueline/thesis_code/gwas_genesofinterest/spn_penicillin_genesofinterest.txt"
GENES_OF_INTEREST_GFF="/nfs/research/jlees/jacqueline/thesis_code/gwas_genesofinterest/spn_penicillin_genesofinterest.txt"

# All per-cutpoint RATE files, comma-separated, for the overlay.
RATE_FILES_CSV=$(ls "$RATE_DIR"/RATE_values_cutpoint*_depruned.txt | paste -sd,)
if [[ -z "$RATE_FILES_CSV" ]]; then
  echo "ERROR: No RATE_values_cutpoint*_depruned.txt files in $RATE_DIR"
  exit 1
fi

# stk1 region — matches the 01_binary_logistic stk1 plot window
# (lead 1,679,135 +/- 50 kb).
REGION="FM211187.1:1629135-1729135"

# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------
OUTPUT_DIR="$SPECIES_OUTPUT_DIR/plots/02_spn_penicillin_MIC_PPOM_stk1"
mkdir -p "$OUTPUT_DIR"

Rscript "$MAKE_PLOT" \
  --y_metric         rate \
  --region           "$REGION" \
  --positions_file   "$POSITIONS_FILE" \
  --genotype_matrix  "$GENOTYPE_MATRIX" \
  --gff              "$GFF" \
  --variant_effects  "$VARIANT_EFFECTS" \
  --annotations      "$ANNOTATIONS" \
  --genes_of_interest "$GENES_OF_INTEREST" \
  --genes_of_interest_gff "$GENES_OF_INTEREST_GFF" \
  --rate_files       "$RATE_FILES_CSV" \
  --title            "S. pneumoniae penicillin (PPOM) — stk1 (rate)" \
  --output           "$OUTPUT_DIR/stk1_rate.png" \
  --width 10 --height 7

echo ""
echo "Done. Plots written to: $OUTPUT_DIR"
