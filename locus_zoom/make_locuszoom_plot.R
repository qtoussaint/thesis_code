#!/usr/bin/env Rscript
# make_locuszoom_plot.R
# LocusZoom-style regional association plots for bacterial GWAS results
# from the gwas_workflow pipeline (RATE values + Bayesian metrics).
#
# Requires a pre-built ensembldb SQLite for gene track annotation.
# Build one with:  Rscript build_gff_db.R --gff ref.gff --output species.sqlite
#
# Usage:
#   Rscript make_locuszoom_plot.R \
#     --rate_file        results/cppRATE_results/RATE_values.txt \
#     --positions_file   data/variant_positions.csv \
#     --genotype_matrix  results/cppRATE_matrices/design_matrix.csv \
#     --ensdb            spneu_ensdb.sqlite \
#     --lead_variant     4521 \
#     --window           500000 \
#     --output           locus_plot.png
#
# OR specify a region directly:
#   --region 1:1200000-1700000
#
# Optional:
#   --variant_effects    results/fitted_model/depruned_variant_effects.csv
#   --annotations        data/snpeff_annotations.tsv
#   --title              "S. pneumoniae GWAS — pbp2x locus"
#   --width 10 --height 7

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
  library(locuszoomr)
  library(AnnotationDbi)
})

# ---------------------------------------------------------------------------
# CLI arguments
# ---------------------------------------------------------------------------
option_list <- list(
  # Required inputs
  make_option("--rate_file",
    type = "character", default = NULL,
    help = "Path to RATE_values.txt or RATE_values_depruned.txt [required]"
  ),
  make_option("--positions_file",
    type = "character", default = NULL,
    help = paste(
      "Path to positions CSV (phandango format: column 2 = base-pair positions,",
      "one row per variant in same order as RATE file) [required]"
    )
  ),
  make_option("--genotype_matrix",
    type = "character", default = NULL,
    help = paste(
      "Path to design_matrix.csv (header row = variant IDs as integers,",
      "subsequent rows = 0/1 genotypes per sample) [required]"
    )
  ),
  make_option("--ensdb",
    type = "character", default = NULL,
    help = "Path to pre-built ensembldb .sqlite file (from build_gff_db.R) [required]"
  ),
  # Region specification (one of these is required)
  make_option("--lead_variant",
    type = "integer", default = NULL,
    help = "Variant index (integer) of the lead/index variant"
  ),
  make_option("--region",
    type = "character", default = NULL,
    help = "Genomic region as 'seqname:start-end' (e.g. '1:1200000-1700000')"
  ),
  make_option("--window",
    type = "integer", default = 500000L,
    help = "Window size (bp) around lead variant [default: 500000]"
  ),
  # Optional inputs
  make_option("--variant_effects",
    type = "character", default = NULL,
    help = paste(
      "Optional: path to depruned_variant_effects.csv — adds significance",
      "markers (filled points) for variants whose 89% CI excludes zero"
    )
  ),
  make_option("--annotations",
    type = "character", default = NULL,
    help = paste(
      "Optional: path to SNPEff annotations TSV with POS and ANN....GENE",
      "columns — labels top RATE variants with gene names"
    )
  ),
  # Output options
  make_option("--output",
    type = "character", default = "locuszoom_plot.png",
    help = "Output file path (PNG or PDF) [default: locuszoom_plot.png]"
  ),
  make_option("--title",
    type = "character", default = NULL,
    help = "Optional plot title (e.g. species name or locus name)"
  ),
  make_option("--width",
    type = "double", default = 10.0,
    help = "Plot width in inches [default: 10]"
  ),
  make_option("--height",
    type = "double", default = 7.0,
    help = "Plot height in inches [default: 7]"
  ),
  make_option("--top_n_labels",
    type = "integer", default = 5L,
    help = "Number of top RATE variants to label with gene names (if --annotations provided) [default: 5]"
  )
)

opt <- parse_args(OptionParser(
  usage = "%prog [options]",
  option_list = option_list,
  description = paste(
    "Create a LocusZoom-style regional association plot for bacterial GWAS",
    "results from the gwas_workflow pipeline. Y-axis = RATE values.",
    "Dots colored by empirical r² with the lead variant."
  )
))

# ---------------------------------------------------------------------------
# Validate inputs
# ---------------------------------------------------------------------------
stop_if_missing <- function(path, flag) {
  if (is.null(path)) stop(flag, " is required")
  if (!file.exists(path)) stop("File not found: ", path, " (", flag, ")")
}

stop_if_missing(opt$rate_file,       "--rate_file")
stop_if_missing(opt$positions_file,  "--positions_file")
stop_if_missing(opt$genotype_matrix, "--genotype_matrix")
stop_if_missing(opt$ensdb,           "--ensdb")

if (is.null(opt$lead_variant) && is.null(opt$region)) {
  stop("Specify either --lead_variant (variant index) or --region (seqname:start-end)")
}
if (!is.null(opt$variant_effects) && !file.exists(opt$variant_effects)) {
  stop("File not found: ", opt$variant_effects, " (--variant_effects)")
}
if (!is.null(opt$annotations) && !file.exists(opt$annotations)) {
  stop("File not found: ", opt$annotations, " (--annotations)")
}

# ---------------------------------------------------------------------------
# Step 1: Load RATE values
# ---------------------------------------------------------------------------
message("Loading RATE values: ", opt$rate_file)
rate_df <- read.table(opt$rate_file, comment.char = "#", header = FALSE)

# ---------------------------------------------------------------------------
# Step 2: Load genomic positions
# ---------------------------------------------------------------------------
message("Loading positions: ", opt$positions_file)
pos_df    <- fread(opt$positions_file, header = TRUE)
positions <- as.integer(pos_df[[2]])
message("  Loaded ", length(positions), " variant positions")

# Original cpprate format: single column (RATE only)
# Depruned format: 3 columns (SNP_ID, RATE, KLD)
rate_col    <- ifelse(ncol(rate_df) >= 2, 2, 1)
rate_values <- as.numeric(rate_df[, rate_col])
if (length(rate_values) != length(positions)) {
  stop("Mismatch: ", length(rate_values), " RATE values but ", length(positions), " positions.")
}
message("  Loaded ", length(rate_values), " RATE values")

# ---------------------------------------------------------------------------
# Step 3: Build full GWAS data frame
# ---------------------------------------------------------------------------
gwas_df <- data.table(
  SNP  = seq_along(rate_values),   # variant index (1-based integer)
  CHR  = 1L,                       # bacteria: single chromosome coded as 1
  BP   = positions,
  RATE = rate_values
)

# ---------------------------------------------------------------------------
# Step 4: Optionally load variant effects for significance annotation
# ---------------------------------------------------------------------------
sig_variants <- NULL
if (!is.null(opt$variant_effects)) {
  message("Loading variant effects: ", opt$variant_effects)
  eff_df <- fread(opt$variant_effects)
  # Expected columns: variant_id, median, signif, cutpoint, ...
  if (!"signif" %in% names(eff_df)) {
    warning("--variant_effects file has no 'signif' column; skipping significance markers")
  } else if (!"variant_id" %in% names(eff_df)) {
    warning("--variant_effects file has no 'variant_id' column; skipping significance markers")
  } else {
    sig_variants <- eff_df[signif == TRUE, variant_id]
    message("  ", length(sig_variants), " variants with 89% CI excluding zero")
  }
}

# ---------------------------------------------------------------------------
# Step 5: Optionally load SNPEff annotations for gene labels
# ---------------------------------------------------------------------------
ann_df <- NULL
if (!is.null(opt$annotations)) {
  message("Loading SNPEff annotations: ", opt$annotations)
  ann_df <- fread(opt$annotations)
  if (!"POS" %in% names(ann_df)) {
    warning("--annotations file has no 'POS' column; skipping gene labels")
    ann_df <- NULL
  }
}

# ---------------------------------------------------------------------------
# Step 6: Determine region
# ---------------------------------------------------------------------------
if (!is.null(opt$region)) {
  # Parse 'seqname:start-end'
  parts   <- strsplit(opt$region, ":")[[1]]
  seqname <- parts[1]
  coords  <- as.integer(strsplit(parts[2], "-")[[1]])
  region_start <- coords[1]
  region_end   <- coords[2]
  # Find lead variant as the one with highest RATE in region
  in_region <- gwas_df[BP >= region_start & BP <= region_end]
  if (nrow(in_region) == 0) stop("No variants found in region: ", opt$region)
  lead_idx     <- in_region[which.max(RATE), SNP]
  lead_pos     <- in_region[which.max(RATE), BP]
} else {
  lead_idx     <- opt$lead_variant
  lead_pos_row <- gwas_df[SNP == lead_idx]
  if (nrow(lead_pos_row) == 0) stop("Lead variant not found: ", lead_idx)
  lead_pos     <- lead_pos_row$BP
  seqname      <- "1"
  region_start <- max(1L, lead_pos - opt$window)
  region_end   <- lead_pos + opt$window
}

message(sprintf(
  "Region: %s:%d-%d  |  Lead variant: %d (pos %d)",
  seqname, region_start, region_end, lead_idx, lead_pos
))

# Filter to region
regional_df <- gwas_df[BP >= region_start & BP <= region_end]
if (nrow(regional_df) == 0) stop("No variants in region after filtering")
message("  ", nrow(regional_df), " variants in region")

# ---------------------------------------------------------------------------
# Step 7: Compute r² with lead variant from genotype matrix
# ---------------------------------------------------------------------------
message("Computing r² from genotype matrix: ", opt$genotype_matrix)
message("  (This may take a moment for large matrices...)")

# Read only the header to get variant IDs, then load selectively
geno_header <- fread(opt$genotype_matrix, nrows = 0)
all_variant_ids <- as.integer(names(geno_header))

# Find columns to load: regional variants + lead variant
regional_ids  <- regional_df$SNP
cols_needed   <- unique(c(lead_idx, regional_ids))
col_positions <- which(all_variant_ids %in% cols_needed)

if (length(col_positions) == 0) {
  stop("No regional variant columns found in genotype matrix. ",
       "Check that variant IDs in design_matrix.csv match those in positions file.")
}

# Load only the needed columns (select= uses 1-based column indices)
geno_sub <- fread(opt$genotype_matrix, select = col_positions)

# Align column names back to variant IDs
sub_ids <- all_variant_ids[col_positions]
names(geno_sub) <- as.character(sub_ids)

lead_col_name <- as.character(lead_idx)
if (!lead_col_name %in% names(geno_sub)) {
  stop("Lead variant column '", lead_col_name, "' not found in genotype matrix")
}

lead_geno <- as.numeric(geno_sub[[lead_col_name]])

# r² = cor(variant, lead)^2 for each regional variant
compute_r2 <- function(x, y) {
  r <- suppressWarnings(cor(x, y, use = "complete.obs"))
  if (is.na(r)) return(0)
  r^2
}

r2_vals <- vapply(
  as.character(regional_ids),
  function(vid) {
    if (!vid %in% names(geno_sub)) return(NA_real_)
    compute_r2(as.numeric(geno_sub[[vid]]), lead_geno)
  },
  FUN.VALUE = numeric(1)
)

regional_df[, r2 := r2_vals]

# Lead variant itself gets r² = 1 by definition
regional_df[SNP == lead_idx, r2 := 1.0]

message(sprintf(
  "  r² range: %.3f – %.3f (median %.3f)",
  min(regional_df$r2, na.rm = TRUE),
  max(regional_df$r2, na.rm = TRUE),
  median(regional_df$r2, na.rm = TRUE)
))

# Warn for M. tuberculosis / highly clonal organisms where all r² ~ 1
high_ld_frac <- mean(regional_df$r2 > 0.8, na.rm = TRUE)
if (high_ld_frac > 0.9) {
  message(sprintf(
    "  NOTE: %.0f%% of variants have r² > 0.8 with the lead variant.",
    high_ld_frac * 100
  ))
  message("  This is expected for M. tuberculosis (non-recombining).")
  message("  LD coloring will have limited discriminatory power in this region.")
}

# ---------------------------------------------------------------------------
# Step 8: Add gene labels from SNPEff annotations (optional)
# ---------------------------------------------------------------------------
if (!is.null(ann_df)) {
  gene_col <- grep("GENE", names(ann_df), value = TRUE)[1]
  if (!is.na(gene_col)) {
    regional_df[, gene := ann_df[[gene_col]][match(BP, ann_df$POS)]]
  }
}

# Mark significance (filled vs open points)
regional_df[, significant := FALSE]
if (!is.null(sig_variants)) {
  regional_df[SNP %in% sig_variants, significant := TRUE]
}

# ---------------------------------------------------------------------------
# Step 9: Build locuszoomr locus object
# ---------------------------------------------------------------------------
message("Loading gene annotation db: ", opt$ensdb)
ensdb <- loadDb(opt$ensdb)

# locuszoomr expects certain column names; create P as a dummy (RATE is yvar)
# We pass a tiny dummy P column — pcutoff line will be suppressed below.
regional_df[, P := 1.0]   # dummy; yvar = "RATE" overrides this on y-axis

message("Creating locus object...")
loc <- tryCatch({
  locus(
    data      = as.data.frame(regional_df),
    seqname   = seqname,
    xrange    = c(region_start, region_end),
    chrom     = "CHR",
    pos       = "BP",
    p         = "P",
    yvar      = "RATE",
    LD        = "r2",
    labs      = if ("gene" %in% names(regional_df)) "gene" else NULL,
    index_snp = lead_idx,
    ens_db    = ensdb
  )
}, error = function(e) {
  message("\nError creating locus object: ", conditionMessage(e))
  message("Trying with seqname auto-detected from ensembldb...")
  # Some GFFs use numeric chromosomes, others use names like "Chromosome"
  available_seqs <- seqlevels(ensdb)
  message("  Available seqnames in ensdb: ", paste(available_seqs, collapse = ", "))
  stop("Adjust --region seqname to match one of the above, then retry.")
})

# ---------------------------------------------------------------------------
# Step 10: Build ggplot2 locus plot
# ---------------------------------------------------------------------------
message("Building LocusZoom plot...")

# Standard LocusZoom LD colour gradient (grey → blue → cyan → green → orange → red → purple)
ld_colours <- c("grey80", "royalblue", "cyan2", "green3", "orange", "red2", "purple")

p <- locus_ggplot(
  loc,
  ylab      = "RATE value",
  xlab      = paste0("Genomic position (", seqname, ")"),
  showLD    = TRUE,
  LD_scheme = ld_colours,
  pcutoff   = NULL,          # no significance threshold line (RATE has no natural cutoff)
  border    = FALSE,
  cex       = 1.2
)

# Add significance shape aesthetic for variants where 89% CI excludes zero
if (!is.null(sig_variants) && any(regional_df$significant)) {
  # Overlay filled triangles for significant variants
  sig_data <- regional_df[significant == TRUE]
  p <- p +
    geom_point(
      data = as.data.frame(sig_data),
      aes(x = BP, y = RATE),
      shape = 17, size = 2.5, colour = "black", alpha = 0.7,
      inherit.aes = FALSE
    )
}

# Title
plot_title <- if (!is.null(opt$title)) {
  opt$title
} else {
  sprintf("Regional association plot — %s:%d-%d", seqname, region_start, region_end)
}

subtitle_parts <- c(
  sprintf("Lead variant: position %d (index %d)", lead_pos, lead_idx),
  sprintf("n = %d variants", nrow(regional_df))
)
if (high_ld_frac > 0.9) {
  subtitle_parts <- c(subtitle_parts,
    sprintf("%.0f%% variants in r² > 0.8 with lead (high clonality)", high_ld_frac * 100)
  )
}

p <- p +
  labs(
    title    = plot_title,
    subtitle = paste(subtitle_parts, collapse = "  |  ")
  ) +
  theme(
    plot.title    = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 9, colour = "grey40")
  )

# ---------------------------------------------------------------------------
# Step 11: Save
# ---------------------------------------------------------------------------
out_dir <- dirname(opt$output)
if (!dir.exists(out_dir) && out_dir != ".") dir.create(out_dir, recursive = TRUE)

message("Saving plot to: ", opt$output)
ggsave(
  filename = opt$output,
  plot     = p,
  width    = opt$width,
  height   = opt$height,
  dpi      = 300
)

message("Done.")
