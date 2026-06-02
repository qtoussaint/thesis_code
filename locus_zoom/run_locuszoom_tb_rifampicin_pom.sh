#!/usr/bin/env bash

#SBATCH --job-name=locuszoom_tb_rif_pom
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=4:00:00
#SBATCH --error=/nfs/research/jlees/jacqueline/thesis_code/locus_zoom/locuszoom_tb_rif_pom.err
#SBATCH --output=/nfs/research/jlees/jacqueline/thesis_code/locus_zoom/locuszoom_tb_rif_pom.out

#################################################################################
# LocusZoom plots for the genes of interest in the tb_rifampicin
# 08_tb_rifampicin_MIC_POM inference run (POM GWAS — one beta per variant
# shared across cutpoints).
#
# Each gene is plotted three times — once per y-axis metric (rate, abs_median,
# exp_abs_median). The lead variant per gene is chosen as the peak of the
# metric within the gene. Window: 5 kb either side of the lead.
#
# POM has a single per-variant beta (proportional-odds assumption), so it
# writes a single RATE_values_depruned.txt — the input handling mirrors the
# binary script, not the PPOM script.
#
# Pre-step: if fitted_model/depruned_variant_effects.csv is absent, regenerate
# it from the .RDS via regenerate_variant_effects.R. The POM pipeline run did
# not emit this CSV (only the PPOM branch did), so the abs_median /
# exp_abs_median metrics need it produced on the fly. High memory is
# requested for this regeneration step (RDS ≈ 21 GB).
#
# Gene names passed to the lead finder must match the GFF. After the loop we
# rename per-gene PNGs to the display names in tb_rifampicin_genesofinterest.txt
# (col2, with text in parens stripped) — no-op for rpoA/B/C but kept for
# symmetry with the trimethoprim script.
#################################################################################

source ~/.bashrc
mamba activate gwas_pipeline

MAKE_PLOT="/nfs/research/jlees/jacqueline/thesis_code/locus_zoom/make_locuszoom_plot.R"
LEAD_FINDER="/nfs/research/jlees/jacqueline/thesis_code/locus_zoom/pbp_lead_variants.R"
REGEN_EFFECTS="/nfs/research/jlees/jacqueline/thesis_code/locus_zoom/regenerate_variant_effects.R"

# ---------------------------------------------------------------------------
# Reference (M. tuberculosis H37Rv, NC_000962.3)
# ---------------------------------------------------------------------------
SPECIES_OUTPUT_DIR="/nfs/research/jlees/jacqueline/thesis_results/locus_zoom/mtuberculosis"
GFF="$SPECIES_OUTPUT_DIR/reference/tb_ref_NC_000962.3.gff3"

# ---------------------------------------------------------------------------
# Pipeline run
# ---------------------------------------------------------------------------
PIPELINE_OUTPUT_DIR="/nfs/research/jlees/jacqueline/thesis_results/gwas_tb_rifampicin/inference/08_tb_rifampicin_MIC_POM"
DATASET_DIR="/nfs/research/jlees/jacqueline/thesis_results/gwas_datasets/inference/08_tb_rifampicin_MIC"

POSITIONS_FILE="$DATASET_DIR/08_tb_rifampicin_MIC_variant_index.csv"
GENOTYPE_MATRIX="$PIPELINE_OUTPUT_DIR/cppRATE_matrices/design_matrix.csv"
VARIANT_EFFECTS="$PIPELINE_OUTPUT_DIR/fitted_model/depruned_variant_effects.csv"
ANNOTATIONS="/nfs/research/jlees/jacqueline/gwas_data/tuberculosis/genotype/snpEff/fields_filtered.txt"
RATE_DIR="$PIPELINE_OUTPUT_DIR/cppRATE_results"
GENES_OF_INTEREST="/nfs/research/jlees/jacqueline/thesis_code/gwas_genesofinterest/tb_rifampicin_genesofinterest.txt"

# Genes — names as they appear in the GFF.
GENES="rpoA,rpoB,rpoC"

# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------
OUTPUT_DIR="$SPECIES_OUTPUT_DIR/plots/08_tb_rifampicin_MIC_POM_genes_top3"
mkdir -p "$OUTPUT_DIR"

# ---------------------------------------------------------------------------
# Pre-step: regenerate depruned_variant_effects.csv if missing
# ---------------------------------------------------------------------------
if [[ ! -f "$VARIANT_EFFECTS" ]]; then
  echo "=== depruned_variant_effects.csv missing; regenerating from .RDS ==="
  Rscript "$REGEN_EFFECTS" --run_dir "$PIPELINE_OUTPUT_DIR"
  if [[ ! -f "$VARIANT_EFFECTS" ]]; then
    echo "ERROR: regenerate_variant_effects.R did not produce $VARIANT_EFFECTS"
    exit 1
  fi
fi

# Single-cutpoint RATE file (POM, like binary)
RATE_FILE="$RATE_DIR/RATE_values_depruned.txt"
if [[ ! -f "$RATE_FILE" ]]; then
  RATE_FILE="$RATE_DIR/RATE_values.txt"
fi
if [[ ! -f "$RATE_FILE" ]]; then
  echo "ERROR: No RATE_values{_depruned,}.txt in $RATE_DIR"
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
        --lead_variant     "$vid" \
        --lead_cutpoint    "$cp" \
        --window           5000 \
        "${Y_FLAG[@]}" \
        --title            "M. tuberculosis rifampicin (POM) — ${gene} (${METRIC})" \
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
  for f in "$OUTPUT_DIR"/"${gff_name}"_*.png; do
    [[ -e "$f" ]] || continue
    new="${f/${gff_name}_/${display_clean}_}"
    echo "  rename $(basename "$f") -> $(basename "$new")"
    mv -f "$f" "$new"
  done
done < "$GENES_OF_INTEREST"

echo ""
echo "Done. Plots written to: $OUTPUT_DIR"
