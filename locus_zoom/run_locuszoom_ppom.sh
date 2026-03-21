#!/usr/bin/env bash

#SBATCH --job-name=locuszoom_ppom             # job name
#SBATCH --nodes=1                              # number of nodes
#SBATCH --cpus-per-task=4                      # CPUs (r² computation is single-threaded)
#SBATCH --mem=32G                              # memory (genotype matrix can be large)
#SBATCH --time=2:00:00                         # time limit (2h to cover all cutpoints)
#SBATCH --error=/nfs/research/jlees/jacqueline/thesis_code/locus_zoom/locuszoom_ppom.err
#SBATCH --output=/nfs/research/jlees/jacqueline/thesis_code/locus_zoom/locuszoom_ppom.out

#################################################################################
# LocusZoom-style regional association plots for ppom bacterial GWAS
#
# Automatically discovers all cutpoints by globbing for:
#   cppRATE_results/RATE_values_cutpoint*_depruned.txt
#
# Produces one plot per cutpoint. Works for any number of cutpoints (K-1).
#
# Prerequisites (run once):
#   mamba activate gwas_pipeline
#   mamba install -c conda-forge r-patchwork r-ggrepel
#################################################################################

MAKE_PLOT="/nfs/research/jlees/jacqueline/thesis_code/locus_zoom/make_locuszoom_plot_v2.R"

# ---------------------------------------------------------------------------
# Species / reference paths
# ---------------------------------------------------------------------------

# S. pneumoniae ATCC 700669 (GCF_000026665.1, assembly ASM2666v1)
# Seqname in GFF: NC_011900.1
SPECIES_OUTPUT_DIR="/nfs/research/jlees/jacqueline/thesis_results/locus_zoom/spneumoniae"
GFF="$SPECIES_OUTPUT_DIR/reference/genomic.gff"

# ---------------------------------------------------------------------------
# GWAS pipeline output paths (S. pneumoniae penicillin resistance analysis)
# Edit PIPELINE_OUTPUT_DIR to point to your specific pipeline output directory
# ---------------------------------------------------------------------------

PIPELINE_OUTPUT_DIR="/nfs/research/jlees/jacqueline/thesis_results/test_gwas_workflow"

POSITIONS_FILE="/nfs/research/jlees/jacqueline/bayesian_gwas_paper/00_data/inference/01_spn_penicillin_subclusters_K5_variant_index.csv"
GENOTYPE_MATRIX="$PIPELINE_OUTPUT_DIR/cppRATE_matrices/design_matrix.csv"
VARIANT_EFFECTS="$PIPELINE_OUTPUT_DIR/fitted_model/depruned_variant_effects.csv"
ANNOTATIONS="/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/genotype/fields_filtered_maf05_multiallelic.txt"

# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

OUTPUT_DIR="$SPECIES_OUTPUT_DIR/plots"
mkdir -p "$OUTPUT_DIR"

# ---------------------------------------------------------------------------
# Discover cutpoints by globbing for existing depruned RATE files
# ---------------------------------------------------------------------------

RATE_FILES=("$PIPELINE_OUTPUT_DIR"/cppRATE_results/RATE_values_cutpoint*_depruned.txt)

if [[ ! -f "${RATE_FILES[0]}" ]]; then
  echo "ERROR: No RATE_values_cutpoint*_depruned.txt files found in:"
  echo "  $PIPELINE_OUTPUT_DIR/cppRATE_results/"
  echo "Is this a ppom pipeline output directory?"
  exit 1
fi

echo "Found ${#RATE_FILES[@]} cutpoint RATE file(s)."

# ---------------------------------------------------------------------------
# Make one LocusZoom plot per cutpoint
#
# Specify either:
#   --lead_variant <integer index>   peak variant index from RATE output
#   --region <seqname:start-end>     explicit genomic window
#
# The seqname in --region MUST match the sequence name in the GFF (NC_011900.1).
# ---------------------------------------------------------------------------

for RATE_FILE in "${RATE_FILES[@]}"; do
  # Extract cutpoint number from filename, e.g. RATE_values_cutpoint3_depruned.txt → 3
  CP=$(basename "$RATE_FILE" | sed 's/RATE_values_cutpoint\([0-9]*\)_depruned\.txt/\1/')

  echo "=== Making plot for cutpoint ${CP} ==="
  echo "    RATE file: $RATE_FILE"

  Rscript "$MAKE_PLOT" \
    --rate_file       "$RATE_FILE" \
    --positions_file  "$POSITIONS_FILE" \
    --genotype_matrix "$GENOTYPE_MATRIX" \
    --gff             "$GFF" \
    --lead_variant    4521 \
    --window          2000 \
    --variant_effects "$VARIANT_EFFECTS" \
    --annotations     "$ANNOTATIONS" \
    --title           "S. pneumoniae GWAS — penicillin resistance (cutpoint ${CP})" \
    --output          "$OUTPUT_DIR/spneu_penicillin_cutpoint${CP}.png" \
    --width           10 \
    --height          7
done

echo ""
echo "Done. Plots written to: $OUTPUT_DIR"
