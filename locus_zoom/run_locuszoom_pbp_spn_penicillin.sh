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
# Each gene is plotted three times — once per y-axis metric (rate, abs_median,
# exp_abs_median). Within a plot, all cutpoints are overlaid (color = r² with
# lead, shape = cutpoint). The lead pair (variant, cutpoint) is the peak of the
# chosen metric inside the gene. Window: 25 kb either side of the lead.
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
RATE_DIR="$PIPELINE_OUTPUT_DIR/cppRATE_results"
GENES_OF_INTEREST="/nfs/research/jlees/jacqueline/thesis_code/gwas_genesofinterest/spn_penicillin_genesofinterest.txt"

# Top 5 pbp genes -- names as they appear in the GFF (which already uses display names here).
GENES="pbp2X,pbp1a,pbp1b,pbp2a,pbp2b"

# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------
OUTPUT_DIR="$SPECIES_OUTPUT_DIR/plots/02_spn_penicillin_MIC_PPOM_pbp_top5"
mkdir -p "$OUTPUT_DIR"

# Comma-separated list of per-cutpoint RATE files, for --rate_files
RATE_FILES_CSV=$(ls "$RATE_DIR"/RATE_values_cutpoint*_depruned.txt | paste -sd,)
if [[ -z "$RATE_FILES_CSV" ]]; then
  echo "ERROR: No RATE_values_cutpoint*_depruned.txt files in $RATE_DIR"
  exit 1
fi

# ---------------------------------------------------------------------------
# One pass per metric
# ---------------------------------------------------------------------------
for METRIC in rate abs_median exp_abs_median; do
  echo ""
  echo "============================================================"
  echo "=== Metric: ${METRIC}"
  echo "============================================================"

  LEAD_TSV="$OUTPUT_DIR/lead_variants_${METRIC}.tsv"

  echo "--- Identifying lead variants per gene ---"
  Rscript "$LEAD_FINDER" \
    --variant_effects "$VARIANT_EFFECTS" \
    --positions_file  "$POSITIONS_FILE" \
    --gff             "$GFF" \
    --genes           "$GENES" \
    --y_metric        "$METRIC" \
    --rate_dir        "$RATE_DIR" \
    --output          "$LEAD_TSV"

  if [[ ! -s "$LEAD_TSV" ]]; then
    echo "ERROR: lead variant TSV is empty: $LEAD_TSV"
    exit 1
  fi

  # TSV cols: gene  lead_variant  lead_pos  lead_cutpoint  lead_metric  lead_metric_value
  {
    read -r _header
    while IFS=$'\t' read -r gene vid pos cp _metric metric_value; do
      if [[ -z "$vid" || "$vid" == "NA" ]]; then
        echo "WARNING: no lead variant for ${gene}; skipping"
        continue
      fi

      echo "--- ${gene}: lead variant ${vid} at ${pos} bp (cutpoint ${cp}, ${METRIC}=${metric_value}) ---"

      if [[ "$METRIC" == "rate" ]]; then
        Y_FLAG=(--rate_files "$RATE_FILES_CSV")
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
        --lead_variant     "$vid" \
        --lead_cutpoint    "$cp" \
        --window           5000 \
        "${Y_FLAG[@]}" \
        --title            "S. pneumoniae penicillin PPOM — ${gene} (${METRIC})" \
        --output           "$OUTPUT_DIR/${gene}_${METRIC}.png" \
        --width 10 --height 7
    done
  } < "$LEAD_TSV"
done

echo ""
echo "Done. Plots written to: $OUTPUT_DIR"
