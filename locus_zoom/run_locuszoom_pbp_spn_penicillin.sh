#!/usr/bin/env bash

#SBATCH --job-name=locuszoom_pbp_spn_pen
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --error=/nfs/research/jlees/jacqueline/thesis_code/locus_zoom/locuszoom_pbp.err
#SBATCH --output=/nfs/research/jlees/jacqueline/thesis_code/locus_zoom/locuszoom_pbp.out

#################################################################################
# LocusZoom plots for the top 5 pbp genes in the spn_penicillin
# 02_spn_penicillin_MIC_PPOM inference run.
#
# For each gene we pick the variant with the largest |median| effect across all
# cutpoints (depruned_variant_effects.csv); that variant becomes the lead, its
# cutpoint determines which RATE_values_cutpointN_depruned.txt feeds the
# y-axis, and the window extends 50 kb either side.
#################################################################################

source ~/.bashrc
mamba activate gwas_pipeline

MAKE_PLOT="/nfs/research/jlees/jacqueline/thesis_code/locus_zoom/make_locuszoom_plot.R"
LEAD_FINDER="/nfs/research/jlees/jacqueline/thesis_code/locus_zoom/pbp_lead_variants.R"

# ---------------------------------------------------------------------------
# Reference (S. pneumoniae ATCC 700669, seqname NC_011900.1)
# ---------------------------------------------------------------------------
SPECIES_OUTPUT_DIR="/nfs/research/jlees/jacqueline/thesis_results/locus_zoom/spneumoniae"
GFF="$SPECIES_OUTPUT_DIR/reference/genomic.gff"

# ---------------------------------------------------------------------------
# Pipeline run
# ---------------------------------------------------------------------------
PIPELINE_OUTPUT_DIR="/nfs/research/jlees/jacqueline/thesis_results/gwas_spn_penicillin/inference/02_spn_penicillin_MIC_PPOM"
DATASET_DIR="/nfs/research/jlees/jacqueline/thesis_results/gwas_datasets/inference/02_spn_penicillin_MIC"

POSITIONS_FILE="$DATASET_DIR/02_spn_penicillin_MIC_variant_index.csv"
GENOTYPE_MATRIX="$PIPELINE_OUTPUT_DIR/cppRATE_matrices/design_matrix.csv"
VARIANT_EFFECTS="$PIPELINE_OUTPUT_DIR/fitted_model/depruned_variant_effects.csv"
ANNOTATIONS="/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/genotype/fields_filtered_maf05_multiallelic.txt"

# Top 5 pbp genes from spn_penicillin_genesofinterest.txt
# (GFF Name= values; case-insensitive match in pbp_lead_variants.R)
GENES="pbp2X,pbp1a,pbp1b,pbp2a,pbp2b"

# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------
OUTPUT_DIR="$SPECIES_OUTPUT_DIR/plots/02_spn_penicillin_MIC_PPOM_pbp_top5"
mkdir -p "$OUTPUT_DIR"
LEAD_TSV="$OUTPUT_DIR/lead_variants.tsv"

# ---------------------------------------------------------------------------
# 1. Identify the lead variant per gene
# ---------------------------------------------------------------------------
echo "=== Identifying lead variants for: $GENES ==="
Rscript "$LEAD_FINDER" \
  --variant_effects "$VARIANT_EFFECTS" \
  --positions_file  "$POSITIONS_FILE" \
  --gff             "$GFF" \
  --genes           "$GENES" \
  --output          "$LEAD_TSV"

if [[ ! -s "$LEAD_TSV" ]]; then
  echo "ERROR: lead variant TSV is empty: $LEAD_TSV"
  exit 1
fi

# ---------------------------------------------------------------------------
# 2. One locus zoom plot per gene, using its (variant, cutpoint) pair
# ---------------------------------------------------------------------------
# TSV columns: gene  lead_variant  lead_pos  lead_cutpoint  lead_abs_median
{
  read -r _header  # skip header
  while IFS=$'\t' read -r gene vid pos cp absmed; do
    if [[ -z "$vid" || "$vid" == "NA" ]]; then
      echo "WARNING: no lead variant for ${gene}; skipping"
      continue
    fi

    RATE_FILE="$PIPELINE_OUTPUT_DIR/cppRATE_results/RATE_values_cutpoint${cp}_depruned.txt"
    if [[ ! -f "$RATE_FILE" ]]; then
      echo "WARNING: RATE file missing for ${gene} (cutpoint ${cp}): $RATE_FILE"
      continue
    fi

    echo "=== ${gene}: lead variant ${vid} at ${pos} bp (cutpoint ${cp}, |median|=${absmed}) ==="
    Rscript "$MAKE_PLOT" \
      --rate_file       "$RATE_FILE" \
      --positions_file  "$POSITIONS_FILE" \
      --genotype_matrix "$GENOTYPE_MATRIX" \
      --gff             "$GFF" \
      --variant_effects "$VARIANT_EFFECTS" \
      --annotations     "$ANNOTATIONS" \
      --lead_variant    "$vid" \
      --window          50000 \
      --title           "S. pneumoniae penicillin PPOM — ${gene} (cutpoint ${cp})" \
      --output          "$OUTPUT_DIR/${gene}.png" \
      --width           10 \
      --height          7
  done
} < "$LEAD_TSV"

echo ""
echo "Done. Plots written to: $OUTPUT_DIR"
