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
# Two steps:
#   STEP 1 (run once per species): build_gff_db.R
#     Converts the reference GFF3 to an ensembldb SQLite for gene track annotation.
#     Comment this out after the first run and reuse the .sqlite file.
#
#   STEP 2 (run per locus): make_locuszoom_plot.R
#     Produces the LocusZoom plot for a specific region or lead variant.
#
# NOTE on GFF format: This GFF is NCBI RefSeq GFF3, where CDS features link
# directly to gene (no mRNA intermediate). If ensDbFromGff() fails, run the
# AGAT preprocessing command below to add synthetic mRNA features before Step 1.
#
# conda activate gwas_pipeline   # (activate your environment if needed)
#################################################################################

BUILD_GFF_DB="/nfs/research/jlees/jacqueline/thesis_code/locus_zoom/build_gff_db.R"
MAKE_PLOT="/nfs/research/jlees/jacqueline/thesis_code/locus_zoom/make_locuszoom_plot.R"

# ---------------------------------------------------------------------------
# Species / reference paths
# ---------------------------------------------------------------------------

# S. pneumoniae ATCC 700669 (GCF_000026665.1, assembly ASM2666v1)
# Seqname in GFF: NC_011900.1
SPECIES_OUTPUT_DIR="/nfs/research/jlees/jacqueline/thesis_results/locus_zoom/spneumoniae"
GFF="$SPECIES_OUTPUT_DIR/reference/genomic.gff"
ENSDB_SQLITE="$SPECIES_OUTPUT_DIR/reference/spneu_ATCC700669_ensdb.sqlite"
ORGANISM="Streptococcus pneumoniae"
GENOME_VERSION="ASM2666v1"
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
mkdir -p "$SPECIES_OUTPUT_DIR/reference"
mkdir -p "$SPECIES_OUTPUT_DIR/plots"

# ---------------------------------------------------------------------------
# Copy reference GFF into thesis_results (run once, then comment out)
# ---------------------------------------------------------------------------

echo "=== Copying reference GFF to thesis_results ==="
GFF_SOURCE="/nfs/research/jlees/jacqueline/bayesian_gwas_paper/annotationhub/genomic.gff"
cp "$GFF_SOURCE" "$SPECIES_OUTPUT_DIR/reference/genomic.gff"

# ---------------------------------------------------------------------------
# STEP 0 (optional): Preprocess GFF if ensDbFromGff() fails on NCBI GFF3
#
# NCBI RefSeq GFFs sometimes omit mRNA features for protein-coding genes,
# which ensDbFromGff() requires. AGAT adds them back:
#

# mamba activate agat
# (had to create new environment to install agat)

   agat_convert_sp_gxf2gxf.pl \
     --gxf "$GFF" \
     --output "$SPECIES_OUTPUT_DIR/reference/genomic_agat.gff"
   GFF="$SPECIES_OUTPUT_DIR/reference/genomic_agat.gff"
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# STEP 1: Build ensembldb from GFF (run once per species, then comment out)
# ---------------------------------------------------------------------------

echo "=== STEP 1: Building ensembldb from GFF ==="
Rscript "$BUILD_GFF_DB" \
  --gff            "$GFF" \
  --output         "$ENSDB_SQLITE" \
  --organism       "$ORGANISM" \
  --genome_version "$GENOME_VERSION"

# ---------------------------------------------------------------------------
# STEP 2: Make LocusZoom plot(s)
#
# Specify either:
#   --lead_variant <integer index>   peak variant index from RATE output
#   --region <seqname:start-end>     explicit genomic window
#
# The seqname MUST match the sequence name in the GFF (here: NC_011900.1).
#
# Example A — by lead variant index (e.g. top RATE hit):
#   Rscript "$MAKE_PLOT" --lead_variant 4521 ...
#
# Example B — by genomic region (e.g. pbp2x locus in S. pneumoniae):
#   --region NC_011900.1:1900000-2000000
# ---------------------------------------------------------------------------

echo ""
echo "=== STEP 2: Making LocusZoom plot ==="
Rscript "$MAKE_PLOT" \
  --rate_file        "$RATE_FILE" \
  --positions_file   "$POSITIONS_FILE" \
  --genotype_matrix  "$GENOTYPE_MATRIX" \
  --ensdb            "$ENSDB_SQLITE" \
  --lead_variant     4521 \
  --window           250000 \
  --variant_effects  "$VARIANT_EFFECTS" \
  --annotations      "$ANNOTATIONS" \
  --title            "S. pneumoniae GWAS — penicillin resistance" \
  --output           "$OUTPUT_DIR/spneu_penicillin_locus.png" \
  --width            10 \
  --height           7

# ---------------------------------------------------------------------------
# To plot a second locus (e.g. a different peak), repeat Step 2 with a new
# --lead_variant / --region and --output, reusing the same --ensdb sqlite.
# ---------------------------------------------------------------------------

# Rscript "$MAKE_PLOT" \
#   --rate_file        "$RATE_FILE" \
#   --positions_file   "$POSITIONS_FILE" \
#   --genotype_matrix  "$GENOTYPE_MATRIX" \
#   --ensdb            "$ENSDB_SQLITE" \
#   --region           "${SEQNAME}:800000-1100000" \
#   --variant_effects  "$VARIANT_EFFECTS" \
#   --annotations      "$ANNOTATIONS" \
#   --title            "S. pneumoniae GWAS — second locus" \
#   --output           "$OUTPUT_DIR/spneu_locus2.png"

echo ""
echo "Done. Plots written to: $OUTPUT_DIR"
