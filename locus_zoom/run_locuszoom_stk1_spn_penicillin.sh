#!/usr/bin/env bash

#SBATCH --job-name=locuszoom_stk1_spn_pen
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=1:00:00
#SBATCH --error=/nfs/research/jlees/jacqueline/thesis_code/locus_zoom/locuszoom_stk1.err
#SBATCH --output=/nfs/research/jlees/jacqueline/thesis_code/locus_zoom/locuszoom_stk1.out

#################################################################################
# LocusZoom plot of stk1 (StkP / SPN23F17350) in the spn_penicillin
# 01_spn_penicillin_binary_logistic inference run.
#
# This is a binary logistic (non-PPOM) run, so there is a single cutpoint and a
# single RATE_values_depruned.txt. The lead variant is the peak of the RATE
# metric within stk1 (+/- a small window). Only RATE is plotted: the posterior
# effect sizes on this run are all ~0, so exp_abs_median just tracks noise.
#################################################################################

source ~/.bashrc
mamba activate gwas_pipeline

MAKE_PLOT="/nfs/research/jlees/jacqueline/thesis_code/locus_zoom/make_locuszoom_plot.R"
LEAD_FINDER="/nfs/research/jlees/jacqueline/thesis_code/locus_zoom/pbp_lead_variants.R"

# ---------------------------------------------------------------------------
# Reference (S. pneumoniae ATCC 700669, GCA_000026665.1_ASM2666v1, seqname FM211187.1)
# ---------------------------------------------------------------------------
SPECIES_OUTPUT_DIR="/nfs/research/jlees/jacqueline/thesis_results/locus_zoom/spneumoniae"
GFF="$SPECIES_OUTPUT_DIR/reference/GCA_000026665.1_ASM2666v1_genomic.gff"

# ---------------------------------------------------------------------------
# Pipeline run
# ---------------------------------------------------------------------------
PIPELINE_OUTPUT_DIR="/nfs/research/jlees/jacqueline/thesis_results/gwas_spn_penicillin/inference/01_spn_penicillin_binary_logistic"
DATASET_DIR="/nfs/research/jlees/jacqueline/thesis_results/gwas_datasets/inference/01_spn_penicillin_binary"

POSITIONS_FILE="$DATASET_DIR/01_spn_penicillin_binary_variant_index.csv"
GENOTYPE_MATRIX="$PIPELINE_OUTPUT_DIR/cppRATE_matrices/design_matrix.csv"
VARIANT_EFFECTS="$PIPELINE_OUTPUT_DIR/fitted_model/depruned_variant_effects.csv"
ANNOTATIONS="/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/genotype/fields_filtered_maf05_multiallelic.txt"
RATE_DIR="$PIPELINE_OUTPUT_DIR/cppRATE_results"
GENES_OF_INTEREST="/nfs/research/jlees/jacqueline/thesis_code/gwas_genesofinterest/spn_penicillin_genesofinterest.txt"
GENES_OF_INTEREST_GFF="/nfs/research/jlees/jacqueline/thesis_code/gwas_genesofinterest/spn_penicillin_genesofinterest.txt"

# Gene name as it appears in the EMBL GFF (Name=stk1).
GENES="stk1"

# Single-cutpoint RATE file (no per-cutpoint suffix on a binary logistic run).
RATE_FILE="$RATE_DIR/RATE_values_depruned.txt"
if [[ ! -f "$RATE_FILE" ]]; then
  RATE_FILE="$RATE_DIR/RATE_values.txt"
fi

# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------
OUTPUT_DIR="$SPECIES_OUTPUT_DIR/plots/01_spn_penicillin_binary_logistic_stk1"
mkdir -p "$OUTPUT_DIR"

# ---------------------------------------------------------------------------
# One pass per metric
# ---------------------------------------------------------------------------
for METRIC in rate; do
  echo ""
  echo "============================================================"
  echo "=== Metric: ${METRIC}"
  echo "============================================================"

  LEAD_TSV="$OUTPUT_DIR/lead_variants_${METRIC}.tsv"

  echo "--- Identifying lead variant in stk1 ---"
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
        --lead_variant     "$vid" \
        --lead_cutpoint    "$cp" \
        --window           50000 \
        "${Y_FLAG[@]}" \
        --title            "S. pneumoniae penicillin (binary logistic) — ${gene} (${METRIC})" \
        --output           "$OUTPUT_DIR/${gene}_${METRIC}.png" \
        --width 10 --height 7
    done
  } < "$LEAD_TSV"
done

echo ""
echo "Done. Plots written to: $OUTPUT_DIR"
