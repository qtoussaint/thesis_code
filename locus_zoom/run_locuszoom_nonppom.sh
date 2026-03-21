#!/usr/bin/env bash

#SBATCH --job-name=locuszoom_bacterial             # job name
#SBATCH --nodes=1                                  # number of nodes
#SBATCH --cpus-per-task=4                          # CPUs (r² computation is single-threaded)
#SBATCH --mem=32G                                  # memory (genotype matrix can be large)
#SBATCH --time=1:00:00                             # time limit
#SBATCH --error=/nfs/research/jlees/jacqueline/thesis_code/locus_zoom/locuszoom.err
#SBATCH --output=/nfs/research/jlees/jacqueline/thesis_code/locus_zoom/locuszoom.out

#################################################################################
# LocusZoom-style regional association plots for bacterial GWAS
#
# make_locuszoom_plot_v2.R: pure ggplot2 + patchwork implementation.
# Reads the reference GFF3 directly — no SQLite database build step needed.
# No Bioconductor packages required.
#
# Prerequisite (run once):
#   mamba activate gwas_pipeline
#   mamba install -c conda-forge r-patchwork
#
# conda activate gwas_pipeline   # (activate your environment if needed)
#################################################################################

MAKE_PLOT="/nfs/research/jlees/jacqueline/thesis_code/locus_zoom/make_locuszoom_plot_v2.R"

# ---------------------------------------------------------------------------
# Species / reference paths
# ---------------------------------------------------------------------------

# S. pneumoniae ATCC 700669 (GCF_000026665.1, assembly ASM2666v1)
# Seqname in GFF: NC_011900.1
SPECIES_OUTPUT_DIR="/nfs/research/jlees/jacqueline/thesis_results/locus_zoom/spneumoniae"
GFF="$SPECIES_OUTPUT_DIR/reference/genomic.gff"
SEQNAME="NC_011900.1"   # chromosome seqname from GFF header; used in --region below

# ---------------------------------------------------------------------------
# GWAS pipeline output paths (S. pneumoniae penicillin resistance analysis)
# Edit these to point to your specific pipeline output directory
# ---------------------------------------------------------------------------

PIPELINE_OUTPUT_DIR="/nfs/research/jlees/jacqueline/thesis_results/test_gwas_workflow"

RATE_FILE="$PIPELINE_OUTPUT_DIR/cppRATE_results/RATE_values_depruned.txt"
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
# Make LocusZoom plot(s)
#
# Specify either:
#   --lead_variant <integer index>   peak variant index from RATE output
#   --region <seqname:start-end>     explicit genomic window
#
# The seqname in --region MUST match the sequence name in the GFF (NC_011900.1).
#
# Example A — by lead variant index (e.g. top RATE hit):
#   Rscript "$MAKE_PLOT" --lead_variant 4521 ...
#
# Example B — by genomic region (e.g. pbp2x locus in S. pneumoniae):
#   --region NC_011900.1:1900000-2000000
# ---------------------------------------------------------------------------

echo "=== Making LocusZoom plot ==="
Rscript "$MAKE_PLOT" \
  --rate_file        "$RATE_FILE" \
  --positions_file   "$POSITIONS_FILE" \
  --genotype_matrix  "$GENOTYPE_MATRIX" \
  --gff              "$GFF" \
  --lead_variant     4521 \
  --window           250000 \
  --variant_effects  "$VARIANT_EFFECTS" \
  --annotations      "$ANNOTATIONS" \
  --title            "S. pneumoniae GWAS — penicillin resistance" \
  --output           "$OUTPUT_DIR/spneu_penicillin_locus.png" \
  --width            10 \
  --height           7

# ---------------------------------------------------------------------------
# To plot a second locus, repeat with a new --lead_variant / --region and --output
# ---------------------------------------------------------------------------

# Rscript "$MAKE_PLOT" \
#   --rate_file        "$RATE_FILE" \
#   --positions_file   "$POSITIONS_FILE" \
#   --genotype_matrix  "$GENOTYPE_MATRIX" \
#   --gff              "$GFF" \
#   --region           "${SEQNAME}:800000-1100000" \
#   --variant_effects  "$VARIANT_EFFECTS" \
#   --annotations      "$ANNOTATIONS" \
#   --title            "S. pneumoniae GWAS — second locus" \
#   --output           "$OUTPUT_DIR/spneu_locus2.png"

echo ""
echo "Done. Plots written to: $OUTPUT_DIR"
