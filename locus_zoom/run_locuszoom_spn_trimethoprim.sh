#!/usr/bin/env bash

#SBATCH --job-name=locuszoom_spn_trim
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --error=/nfs/research/jlees/jacqueline/thesis_code/locus_zoom/locuszoom_spn_trim.err
#SBATCH --output=/nfs/research/jlees/jacqueline/thesis_code/locus_zoom/locuszoom_spn_trim.out

#################################################################################
# LocusZoom plots for the genes of interest in the spn_trimethoprim
# 05_spn_trimethoprim_MIC_PPOM inference run.
#
# Each gene is plotted three times — once per y-axis metric (rate, abs_median,
# exp_abs_median). Within a plot, all cutpoints are overlaid (color = r² with
# lead, shape = cutpoint). The lead pair (variant, cutpoint) is the peak of the
# chosen metric inside the gene. Window: 25 kb either side of the lead.
#
# Gene names passed to the lead finder must match the GFF. After the loop we
# rename per-gene PNGs to the display names in spn_trimethoprim_genesofinterest.txt
# (col2, with text in parens stripped) — e.g. dyr -> folA.
#################################################################################

source ~/.bashrc
mamba activate gwas_pipeline

MAKE_PLOT="/nfs/research/jlees/jacqueline/thesis_code/locus_zoom/make_locuszoom_plot.R"
LEAD_FINDER="/nfs/research/jlees/jacqueline/thesis_code/locus_zoom/pbp_lead_variants.R"

# ---------------------------------------------------------------------------
# Reference (S. pneumoniae ATCC 700669, seqname NC_011900.1)
# ---------------------------------------------------------------------------
SPECIES_OUTPUT_DIR="/nfs/research/jlees/jacqueline/thesis_results/locus_zoom/spneumoniae"
GFF="$SPECIES_OUTPUT_DIR/reference/GCA_000026665.1_ASM2666v1_genomic.gff"

# ---------------------------------------------------------------------------
# Pipeline run
# ---------------------------------------------------------------------------
PIPELINE_OUTPUT_DIR="/nfs/research/jlees/jacqueline/thesis_results/gwas_spn_trimethoprim/inference/05_spn_trimethoprim_MIC_PPOM"
DATASET_DIR="/nfs/research/jlees/jacqueline/thesis_results/gwas_datasets/inference/05_spn_trimethoprim_MIC"

POSITIONS_FILE="$DATASET_DIR/05_spn_trimethoprim_MIC_variant_index.csv"
GENOTYPE_MATRIX="$PIPELINE_OUTPUT_DIR/cppRATE_matrices/design_matrix.csv"
VARIANT_EFFECTS="$PIPELINE_OUTPUT_DIR/fitted_model/depruned_variant_effects.csv"
ANNOTATIONS="/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/genotype/fields_filtered_maf05_multiallelic.txt"
RATE_DIR="$PIPELINE_OUTPUT_DIR/cppRATE_results"
GENES_OF_INTEREST="/nfs/research/jlees/jacqueline/thesis_code/gwas_genesofinterest/spn_trimethoprim_genesofinterest.txt"
GENES_OF_INTEREST_GFF="/nfs/research/jlees/jacqueline/thesis_code/locus_zoom/spn_penicillin_genesofinterest_updatedrefseq.txt"

# Genes — names as they appear in the GFF. Post-rename below maps dyr -> folA.
GENES="folP,dyr"

# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------
OUTPUT_DIR="$SPECIES_OUTPUT_DIR/plots/05_spn_trimethoprim_MIC_PPOM_genes_top2_composite"
mkdir -p "$OUTPUT_DIR"

# Comma-separated list of per-cutpoint RATE files
RATE_FILES_CSV=$(ls "$RATE_DIR"/RATE_values_cutpoint*_depruned.txt | paste -sd,)
if [[ -z "$RATE_FILES_CSV" ]]; then
  echo "ERROR: No RATE_values_cutpoint*_depruned.txt files in $RATE_DIR"
  exit 1
fi

# ---------------------------------------------------------------------------
# One pass per metric
# ---------------------------------------------------------------------------
for METRIC in rate; do
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
        --genes_of_interest_gff "$GENES_OF_INTEREST_GFF" \
        --lead_variant     "$vid" \
        --lead_cutpoint    "$cp" \
        --window           5000 \
        "${Y_FLAG[@]}" \
        --composite_style \
        --title            "S. pneumoniae trimethoprim PPOM — ${gene} (${METRIC})" \
        --output           "$OUTPUT_DIR/${gene}_${METRIC}.png" \
        --width 10 --height 7
    done
  } < "$LEAD_TSV"
done

# ---------------------------------------------------------------------------
# Post-rename: GFF gene names -> display names from genes-of-interest aliases.
# Strip text in parens; skip identity renames.
# ---------------------------------------------------------------------------
echo ""
echo "=== Renaming PNGs to display names from $GENES_OF_INTEREST ==="
while IFS=, read -r gff_name display_name; do
  gff_name=$(echo "$gff_name" | xargs)
  display_clean=$(echo "$display_name" | sed -E 's/[[:space:]]*\([^)]*\)//' | xargs)
  [[ -z "$gff_name" || -z "$display_clean" ]] && continue
  [[ "$gff_name" == "$display_clean" ]] && continue
  for f in "$OUTPUT_DIR"/"${gff_name}"_*; do
    [[ -e "$f" ]] || continue
    new="${f/${gff_name}_/${display_clean}_}"
    echo "  rename $(basename "$f") -> $(basename "$new")"
    mv -f "$f" "$new"
  done
done < "$GENES_OF_INTEREST"

echo ""
echo "Done. Plots written to: $OUTPUT_DIR"
