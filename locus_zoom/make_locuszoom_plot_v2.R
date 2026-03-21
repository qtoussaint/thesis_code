#!/usr/bin/env Rscript
# make_locuszoom_plot_v2.R
# LocusZoom-style regional association plots for bacterial GWAS results.
# Pure ggplot2 + patchwork implementation — no locuszoomr or Bioconductor
# packages required. Gene track is drawn directly from the reference GFF3.
#
# Usage:
#   Rscript make_locuszoom_plot_v2.R \
#     --rate_file        results/cppRATE_results/RATE_values.txt \
#     --positions_file   data/variant_positions.csv \
#     --genotype_matrix  results/cppRATE_matrices/design_matrix.csv \
#     --gff              reference/genomic.gff \
#     --lead_variant     4521 \
#     --window           250000 \
#     --output           locus_plot.png
#
# OR specify a region directly:
#   --region NC_011900.1:1900000-2000000
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
  library(ggrepel)
  library(patchwork)
})

# ---------------------------------------------------------------------------
# CLI arguments
# ---------------------------------------------------------------------------
option_list <- list(
  make_option("--rate_file",
    type = "character", default = NULL,
    help = "Path to RATE_values.txt or RATE_values_depruned.txt [required]"
  ),
  make_option("--positions_file",
    type = "character", default = NULL,
    help = paste(
      "Path to positions CSV (column 2 = base-pair positions,",
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
  make_option("--gff",
    type = "character", default = NULL,
    help = "Path to reference GFF3 annotation file (e.g. genomic.gff) [required]"
  ),
  make_option("--lead_variant",
    type = "integer", default = NULL,
    help = "Variant index (integer) of the lead/index variant"
  ),
  make_option("--region",
    type = "character", default = NULL,
    help = "Genomic region as 'seqname:start-end' (e.g. 'NC_011900.1:1900000-2000000')"
  ),
  make_option("--window",
    type = "integer", default = 500000L,
    help = "Window size (bp) around lead variant [default: 500000]"
  ),
  make_option("--variant_effects",
    type = "character", default = NULL,
    help = paste(
      "Optional: path to depruned_variant_effects.csv — marks variants",
      "whose 89% CI excludes zero with filled triangles"
    )
  ),
  make_option("--annotations",
    type = "character", default = NULL,
    help = paste(
      "Optional: path to SNPEff annotations TSV with POS and ANN....GENE",
      "columns — labels top RATE variants with gene names"
    )
  ),
  make_option("--top_n_labels",
    type = "integer", default = 5L,
    help = "Number of top RATE variants to label with gene names [default: 5]"
  ),
  make_option("--output",
    type = "character", default = "locuszoom_plot.png",
    help = "Output file path (PNG or PDF) [default: locuszoom_plot.png]"
  ),
  make_option("--title",
    type = "character", default = NULL,
    help = "Optional plot title"
  ),
  make_option("--width",
    type = "double", default = 10.0,
    help = "Plot width in inches [default: 10]"
  ),
  make_option("--height",
    type = "double", default = 7.0,
    help = "Plot height in inches [default: 7]"
  )
)

opt <- parse_args(OptionParser(
  usage = "%prog [options]",
  option_list = option_list,
  description = paste(
    "LocusZoom-style regional association plot for bacterial GWAS.",
    "Y-axis = RATE values. Points coloured by r² with lead variant.",
    "Gene track drawn directly from reference GFF3."
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
stop_if_missing(opt$gff,             "--gff")

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
  SNP  = seq_along(rate_values),
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
  if (!"signif" %in% names(eff_df) || !"variant_id" %in% names(eff_df)) {
    warning("--variant_effects file missing 'signif' or 'variant_id' column; skipping")
  } else {
    sig_variants <- eff_df[signif == TRUE, variant_id]
    message("  ", length(sig_variants), " variants with 89% CI excluding zero")
  }
}

# ---------------------------------------------------------------------------
# Step 5: Optionally load SNPEff annotations
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
  parts        <- strsplit(opt$region, ":")[[1]]
  seqname      <- parts[1]
  coords       <- as.integer(strsplit(parts[2], "-")[[1]])
  region_start <- coords[1]
  region_end   <- coords[2]
  in_region    <- gwas_df[BP >= region_start & BP <= region_end]
  if (nrow(in_region) == 0) stop("No variants found in region: ", opt$region)
  lead_idx <- in_region[which.max(RATE), SNP]
  lead_pos <- in_region[which.max(RATE), BP]
} else {
  lead_idx     <- opt$lead_variant
  lead_pos_row <- gwas_df[SNP == lead_idx]
  if (nrow(lead_pos_row) == 0) stop("Lead variant not found: ", lead_idx)
  lead_pos     <- lead_pos_row$BP
  seqname      <- NA_character_   # will be filled after GFF load
  region_start <- max(1L, lead_pos - opt$window)
  region_end   <- lead_pos + opt$window
}

message(sprintf(
  "Region: %s:%d-%d  |  Lead variant: %d (pos %d)",
  if (is.na(seqname)) "?" else seqname,
  region_start, region_end, lead_idx, lead_pos
))

regional_df <- gwas_df[BP >= region_start & BP <= region_end]
if (nrow(regional_df) == 0) stop("No variants in region after filtering")
message("  ", nrow(regional_df), " variants in region")

# ---------------------------------------------------------------------------
# Step 7: Compute r² with lead variant from genotype matrix
# ---------------------------------------------------------------------------
message("Computing r² from genotype matrix: ", opt$genotype_matrix)
message("  (This may take a moment for large matrices...)")

# No header row: column i corresponds to variant index i
first_row       <- fread(opt$genotype_matrix, nrows = 1, header = FALSE)
n_variants      <- ncol(first_row)
all_variant_ids <- seq_len(n_variants)

regional_ids  <- regional_df$SNP
cols_needed   <- unique(c(lead_idx, regional_ids))
col_positions <- which(all_variant_ids %in% cols_needed)

if (length(col_positions) == 0) {
  stop("No regional variant columns found in genotype matrix.")
}

geno_sub      <- fread(opt$genotype_matrix, header = FALSE, select = col_positions)
sub_ids       <- all_variant_ids[col_positions]
names(geno_sub) <- as.character(sub_ids)

lead_col_name <- as.character(lead_idx)
if (!lead_col_name %in% names(geno_sub)) {
  stop("Lead variant column '", lead_col_name, "' not found in genotype matrix")
}

lead_geno <- as.numeric(geno_sub[[lead_col_name]])

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
regional_df[SNP == lead_idx, r2 := 1.0]

message(sprintf(
  "  r² range: %.3f – %.3f (median %.3f)",
  min(regional_df$r2, na.rm = TRUE),
  max(regional_df$r2, na.rm = TRUE),
  median(regional_df$r2, na.rm = TRUE)
))

high_ld_frac <- mean(regional_df$r2 > 0.8, na.rm = TRUE)
if (high_ld_frac > 0.9) {
  message(sprintf(
    "  NOTE: %.0f%% of variants have r² > 0.8 (high clonality / low recombination).",
    high_ld_frac * 100
  ))
}

# ---------------------------------------------------------------------------
# Step 8: Add gene name labels from SNPEff annotations (optional)
# ---------------------------------------------------------------------------
regional_df[, gene_label := NA_character_]
if (!is.null(ann_df)) {
  gene_col <- grep("GENE", names(ann_df), value = TRUE)[1]
  if (!is.na(gene_col)) {
    regional_df[, gene_label := ann_df[[gene_col]][match(BP, ann_df$POS)]]
  }
}

regional_df[, significant := FALSE]
if (!is.null(sig_variants)) {
  regional_df[SNP %in% sig_variants, significant := TRUE]
}

# Keep only top N labels to avoid overplotting
if (!all(is.na(regional_df$gene_label))) {
  top_snps <- regional_df[order(-RATE)][seq_len(min(opt$top_n_labels, .N)), SNP]
  regional_df[!SNP %in% top_snps, gene_label := NA_character_]
}

# ---------------------------------------------------------------------------
# Step 9: Parse GFF3 for gene track
# ---------------------------------------------------------------------------
message("Parsing GFF3: ", opt$gff)

gff_cols <- c("seqname", "source", "feature", "start", "end",
              "score", "strand", "frame", "attributes")

gff_raw <- fread(
  cmd          = paste0("grep -v '^#' ", shQuote(opt$gff)),
  sep          = "\t",
  header       = FALSE,
  col.names    = gff_cols,
  showProgress = FALSE
)

# Detect seqname from GFF if not already known
if (is.na(seqname)) {
  seqname <- gff_raw$seqname[1]
  message("  Auto-detected seqname from GFF: ", seqname)
}

genes_region <- gff_raw[
  feature == "gene" &
  seqname == get("seqname") &
  end   >= region_start &
  start <= region_end
]

if (nrow(genes_region) == 0) {
  message("  WARNING: No gene features found in region — gene track will be empty.")
  message("  Available seqnames in GFF: ",
          paste(unique(gff_raw$seqname), collapse = ", "))
} else {
  message("  ", nrow(genes_region), " genes in region")
}

# Extract gene name: prefer Name= attribute, fall back to gene= then locus_tag=
extract_attr <- function(attrs, key) {
  pattern <- paste0("(?:^|;)", key, "=([^;]+)")
  m <- regmatches(attrs, regexpr(pattern, attrs, perl = TRUE))
  ifelse(nzchar(m),
         sub(paste0(".*", key, "=([^;]+).*"), "\\1", m),
         NA_character_)
}

if (nrow(genes_region) > 0) {
  genes_region[, name := extract_attr(attributes, "Name")]
  genes_region[is.na(name), name := extract_attr(attributes[is.na(name)], "gene")]
  genes_region[is.na(name), name := extract_attr(attributes[is.na(name)], "locus_tag")]
  genes_region[is.na(name), name := paste0("gene_", seq_len(sum(is.na(name))))]

  # Clip gene coordinates to the plot window
  genes_region[, draw_start := pmax(start, region_start)]
  genes_region[, draw_end   := pmin(end,   region_end)]
  genes_region[, mid        := (draw_start + draw_end) / 2]
}

# ---------------------------------------------------------------------------
# Step 10: Build scatter panel
# ---------------------------------------------------------------------------
message("Building plot...")

ld_breaks <- c(0, 0.2, 0.4, 0.6, 0.8, 1.001)
ld_labels <- c("0.0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1.0")
ld_colors <- c("grey75", "royalblue", "cyan3", "green3", "orange", "red2")

regional_df[, ld_bin := cut(r2, breaks = ld_breaks, labels = ld_labels,
                             include.lowest = TRUE, right = TRUE)]
# Lead variant and NAs get their own category
regional_df[is.na(r2), ld_bin := factor(NA, levels = ld_labels)]

plot_title <- if (!is.null(opt$title)) {
  opt$title
} else {
  sprintf("Regional association plot — %s:%d-%d", seqname, region_start, region_end)
}

subtitle_parts <- c(
  sprintf("Lead variant: position %d (index %d)", lead_pos, lead_idx),
  sprintf("n = %d variants in window", nrow(regional_df))
)
if (high_ld_frac > 0.9) {
  subtitle_parts <- c(subtitle_parts,
    sprintf("%.0f%% variants r² > 0.8 with lead (high clonality)", high_ld_frac * 100))
}

# Background (non-lead) points
p_scatter <- ggplot(
    regional_df[SNP != lead_idx],
    aes(x = BP, y = RATE, colour = ld_bin)
  ) +
  geom_point(size = 1.5, alpha = 0.85) +
  scale_colour_manual(
    values  = stats::setNames(ld_colors, ld_labels),
    name    = expression(r^2),
    na.value = "grey60",
    drop    = FALSE
  ) +
  # Lead variant as filled diamond
  geom_point(
    data   = regional_df[SNP == lead_idx],
    aes(x  = BP, y = RATE),
    shape  = 23, size = 3.5, fill = "purple", colour = "black",
    inherit.aes = FALSE
  ) +
  scale_x_continuous(
    limits = c(region_start, region_end),
    expand = expansion(mult = 0.01),
    labels = scales::comma
  ) +
  labs(
    y        = "RATE value",
    title    = plot_title,
    subtitle = paste(subtitle_parts, collapse = "  |  ")
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.title.x     = element_blank(),
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title       = element_text(size = 12, face = "bold"),
    plot.subtitle    = element_text(size = 8.5, colour = "grey40"),
    legend.position  = "right"
  )

# Add gene name labels if available
if (!all(is.na(regional_df$gene_label))) {
  label_data <- regional_df[!is.na(gene_label)]
  p_scatter <- p_scatter +
    geom_text_repel(
      data         = as.data.frame(label_data),
      aes(x = BP, y = RATE, label = gene_label),
      size         = 3,
      colour       = "black",
      inherit.aes  = FALSE,
      max.overlaps = 20,
      box.padding  = 0.4,
      segment.colour = "grey60",
      segment.size   = 0.3
    )
}

# ---------------------------------------------------------------------------
# Step 11: Build gene track panel
# ---------------------------------------------------------------------------
if (nrow(genes_region) > 0) {
  # y position: stagger overlapping genes to reduce label collisions
  # simple approach: alternate between y = 0 and y = 1 for dense loci
  genes_region[, y_level := 0]
  if (nrow(genes_region) > 1) {
    for (i in seq(2, nrow(genes_region))) {
      prev_end <- genes_region[i - 1, draw_end]
      if (genes_region[i, draw_start] < prev_end + 5000) {
        genes_region[i, y_level := 1 - genes_region[i - 1, y_level]]
      }
    }
  }
  gene_height <- 0.35
  genes_region[, ymin := y_level - gene_height]
  genes_region[, ymax := y_level + gene_height]

  # Identify variants that fall within a gene's drawn coordinates
  variant_ticks <- rbindlist(lapply(seq_len(nrow(genes_region)), function(gi) {
    g <- genes_region[gi]
    hits <- regional_df[BP >= g$draw_start & BP <= g$draw_end,
                        .(BP, ld_bin, SNP)]
    if (nrow(hits) == 0) return(NULL)
    hits[, `:=`(tick_ymin = g$ymin, tick_ymax = g$ymax)]
    hits
  }))

  p_genes <- ggplot(genes_region) +
    geom_rect(
      aes(xmin = draw_start, xmax = draw_end,
          ymin = ymin, ymax = ymax, fill = strand),
      colour = "grey30", linewidth = 0.3
    ) +
    # Tick marks showing variant positions within genes, coloured by LD r²
    {if (nrow(variant_ticks) > 0)
      geom_segment(
        data      = as.data.frame(variant_ticks[SNP != lead_idx]),
        aes(x = BP, xend = BP, y = tick_ymin, yend = tick_ymax, colour = ld_bin),
        linewidth = 0.3, alpha = 0.7, inherit.aes = FALSE
      )
    } +
    {if (nrow(variant_ticks) > 0 && lead_idx %in% variant_ticks$SNP) {
      lead_tick <- variant_ticks[SNP == lead_idx]
      geom_segment(
        data      = as.data.frame(lead_tick),
        aes(x = BP, xend = BP, y = tick_ymin, yend = tick_ymax),
        colour    = "purple", linewidth = 0.6, inherit.aes = FALSE
      )
    }} +
    scale_colour_manual(
      values   = stats::setNames(ld_colors, ld_labels),
      name     = expression(r^2),
      na.value = "grey60",
      drop     = FALSE,
      guide    = "none"
    ) +
    geom_text_repel(
      aes(x = mid, y = ymax, label = name),
      size           = 2.8,
      vjust          = 0,
      nudge_y        = 0.2,
      segment.colour = "grey60",
      segment.size   = 0.3,
      max.overlaps   = Inf,
      direction      = "x"
    ) +
    scale_fill_manual(
      values = c("+" = "steelblue3", "-" = "tomato3"),
      labels = c("+" = "Forward (+)", "-" = "Reverse (-)"),
      name   = "Strand"
    ) +
    scale_x_continuous(
      limits = c(region_start, region_end),
      expand = expansion(mult = 0.01),
      labels = scales::comma,
      name   = paste0("Genomic position (", seqname, ")")
    ) +
    scale_y_continuous(expand = expansion(add = 0.5)) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.y      = element_blank(),
      axis.ticks.y     = element_blank(),
      axis.title.y     = element_blank(),
      panel.grid       = element_blank(),
      legend.position  = "right"
    )
} else {
  # Empty gene track placeholder
  p_genes <- ggplot() +
    annotate("text", x = (region_start + region_end) / 2, y = 0,
             label = "No gene features found in region", colour = "grey50") +
    scale_x_continuous(
      limits = c(region_start, region_end),
      labels = scales::comma,
      name   = paste0("Genomic position (", seqname, ")")
    ) +
    theme_bw(base_size = 11)
}

# ---------------------------------------------------------------------------
# Step 12: Assemble and save
# ---------------------------------------------------------------------------
combined <- p_scatter / p_genes + plot_layout(heights = c(3, 1))

out_dir <- dirname(opt$output)
if (!dir.exists(out_dir) && out_dir != ".") dir.create(out_dir, recursive = TRUE)

message("Saving plot to: ", opt$output)
ggsave(
  filename = opt$output,
  plot     = combined,
  width    = opt$width,
  height   = opt$height,
  dpi      = 300
)

message("Done.")
