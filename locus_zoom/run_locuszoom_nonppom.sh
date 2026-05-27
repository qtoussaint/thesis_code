#!/usr/bin/env bash

#SBATCH --job-name=locuszoom_bacterial
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=1:00:00
#SBATCH --error=/nfs/research/jlees/jacqueline/thesis_code/locus_zoom/locuszoom.err
#SBATCH --output=/nfs/research/jlees/jacqueline/thesis_code/locus_zoom/locuszoom.out

#################################################################################
# LocusZoom-style regional association plots for non-PPOM bacterial GWAS
# (binary / continuous phenotypes — single cutpoint).
#
# Three plots are produced (one per y-axis metric: rate, abs_median,
# exp_abs_median). Lead variant is auto-picked as the global peak of the
# chosen metric. Window: 25 kb either side.
#
# Prerequisite (run once):
#   mamba activate gwas_pipeline
#   mamba install -c conda-forge r-patchwork
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
# Pipeline output
# ---------------------------------------------------------------------------
PIPELINE_OUTPUT_DIR="/nfs/research/jlees/jacqueline/thesis_results/test_gwas_workflow"

POSITIONS_FILE="/nfs/research/jlees/jacqueline/bayesian_gwas_paper/00_data/inference/01_spn_penicillin_subclusters_K5_variant_index.csv"
GENOTYPE_MATRIX="$PIPELINE_OUTPUT_DIR/cppRATE_matrices/design_matrix.csv"
VARIANT_EFFECTS="$PIPELINE_OUTPUT_DIR/fitted_model/depruned_variant_effects.csv"
ANNOTATIONS="/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/genotype/fields_filtered_maf05_multiallelic.txt"
RATE_DIR="$PIPELINE_OUTPUT_DIR/cppRATE_results"
GENES_OF_INTEREST="/nfs/research/jlees/jacqueline/thesis_code/gwas_genesofinterest/spn_penicillin_genesofinterest.txt"
GENES_OF_INTEREST_GFF="/nfs/research/jlees/jacqueline/thesis_code/locus_zoom/spn_penicillin_genesofinterest_updatedrefseq.txt"

# Single-cutpoint RATE file (binary / continuous run, no per-cutpoint suffix)
RATE_FILE="$RATE_DIR/RATE_values_depruned.txt"
if [[ ! -f "$RATE_FILE" ]]; then
  # Fall back to RATE_values.txt if the depruned variant doesn't exist
  RATE_FILE="$RATE_DIR/RATE_values.txt"
fi

# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------
OUTPUT_DIR="$SPECIES_OUTPUT_DIR/plots"
mkdir -p "$OUTPUT_DIR"

# ---------------------------------------------------------------------------
# One plot per metric
# ---------------------------------------------------------------------------
for METRIC in rate abs_median exp_abs_median; do
  echo ""
  echo "============================================================"
  echo "=== Metric: ${METRIC}"
  echo "============================================================"

  LEAD_TSV="$OUTPUT_DIR/nonppom_lead_${METRIC}.tsv"

  Rscript "$LEAD_FINDER" \
    --variant_effects "$VARIANT_EFFECTS" \
    --positions_file  "$POSITIONS_FILE" \
    --y_metric        "$METRIC" \
    --rate_dir        "$RATE_DIR" \
    --whole_genome \
    --output          "$LEAD_TSV"

  read -r LEAD_VID LEAD_POS LEAD_CP < <(
    awk -F'\t' 'NR==2 {print $2, $3, $4}' "$LEAD_TSV"
  )
  if [[ -z "$LEAD_VID" || "$LEAD_VID" == "NA" ]]; then
    echo "ERROR: no global lead variant for metric ${METRIC}"
    continue
  fi
  echo "Global lead: variant ${LEAD_VID} at ${LEAD_POS} bp, cutpoint ${LEAD_CP}"

  if [[ "$METRIC" == "rate" ]]; then
    Y_FLAG=(--rate_files "$RATE_FILE")
  else
    Y_FLAG=(--variant_effects "$VARIANT_EFFECTS")
  fi

  Rscript "$MAKE_PLOT" \
    --y_metric         "$METRIC" \
    --positions_file   "$POSITIONS_FILE" \
    --genotype_matrix  "$GENOTYPE_MATRIX" \
    --gff              "$GFF" \
    --variant_effects  "$VARIANT_EFFECTS" \
    --annotations      "$ANNOTATIONS" \
    --genes_of_interest "$GENES_OF_INTEREST" \
    --genes_of_interest_gff "$GENES_OF_INTEREST_GFF" \
    --lead_variant     "$LEAD_VID" \
    --lead_cutpoint    "$LEAD_CP" \
    --window           5000 \
    "${Y_FLAG[@]}" \
    --title            "S. pneumoniae GWAS — penicillin resistance (${METRIC})" \
    --output           "$OUTPUT_DIR/spneu_penicillin_nonppom_${METRIC}.png" \
    --width 10 --height 7
done

echo ""
echo "Done. Plots written to: $OUTPUT_DIR"
