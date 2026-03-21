#!/usr/bin/env Rscript
# build_gff_db.R
# One-time helper: converts a reference GFF3/GTF file to a TxDb SQLite
# database for use with locuszoomr.
#
# Uses GenomicFeatures::makeTxDbFromGFF(), which correctly handles NCBI
# RefSeq GFF3 files (unlike ensDbFromGff(), which requires Ensembl-style
# attributes not present in NCBI annotations).
#
# Run once per species, then pass the output .sqlite to make_locuszoom_plot.R
#
# Usage:
#   Rscript build_gff_db.R \
#     --gff /path/to/reference.gff3 \
#     --output /path/to/output_txdb.sqlite \
#     --organism "Streptococcus pneumoniae" \
#     --genome_version "ATCC700669_v1"  # optional, e.g. reference strain name
#
# Examples:
#   Rscript build_gff_db.R --gff spneu.gff3 --output spneu_txdb.sqlite \
#     --organism "Streptococcus pneumoniae" --genome_version "ATCC700669"
#
#   Rscript build_gff_db.R --gff mtb.gff3 --output mtb_txdb.sqlite \
#     --organism "Mycobacterium tuberculosis" --genome_version "H37Rv"

suppressPackageStartupMessages({
  library(optparse)
  library(GenomicFeatures)
  library(AnnotationDbi)
})

# ---------------------------------------------------------------------------
# CLI arguments
# ---------------------------------------------------------------------------
option_list <- list(
  make_option("--gff",
    type = "character", default = NULL,
    help = "Path to reference GFF3 or GTF annotation file [required]"
  ),
  make_option("--output",
    type = "character", default = NULL,
    help = "Output path for the .sqlite ensembldb file [required]"
  ),
  make_option("--organism",
    type = "character", default = "Unknown organism",
    help = "Organism name string (e.g. 'Streptococcus pneumoniae') [default: 'Unknown organism']"
  ),
  make_option("--genome_version",
    type = "character", default = "v1",
    help = "Genome version/strain label (e.g. 'H37Rv') [default: 'v1']"
  ),
  make_option("--db_version",
    type = "integer", default = 1L,
    help = "Database version integer [default: 1]"
  )
)

opt <- parse_args(OptionParser(
  usage = "%prog [options]",
  option_list = option_list,
  description = paste(
    "Convert a reference GFF3/GTF file to a TxDb SQLite database",
    "for use with locuszoomr bacterial GWAS locus plots."
  )
))

# ---------------------------------------------------------------------------
# Validate inputs
# ---------------------------------------------------------------------------
if (is.null(opt$gff)) stop("--gff is required")
if (is.null(opt$output)) stop("--output is required")
if (!file.exists(opt$gff)) stop("GFF file not found: ", opt$gff)

out_dir <- dirname(opt$output)
if (!dir.exists(out_dir)) {
  message("Creating output directory: ", out_dir)
  dir.create(out_dir, recursive = TRUE)
}

# ---------------------------------------------------------------------------
# Build database
# ---------------------------------------------------------------------------
message("Building TxDb from: ", opt$gff)
message("  Organism    : ", opt$organism)
message("  Genome      : ", opt$genome_version)
message("  Output      : ", opt$output)

tryCatch({
  txdb <- makeTxDbFromGFF(
    file      = opt$gff,
    format    = "gff3",
    organism  = opt$organism,
    dataSource = paste0("NCBI RefSeq ", opt$genome_version)
  )
  saveDb(txdb, file = opt$output)
  message("Done. Database written to: ", opt$output)
  message("\nTo verify, load in R with:")
  message("  library(AnnotationDbi)")
  message("  txdb <- loadDb('", opt$output, "')")
  message("  genes(txdb)")
}, error = function(e) {
  message("\nError building database: ", conditionMessage(e))
  quit(status = 1)
})
