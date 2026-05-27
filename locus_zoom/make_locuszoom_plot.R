#!/usr/bin/env Rscript
# make_locuszoom_plot.R
# LocusZoom-style regional association plots for bacterial GWAS results.
# Pure ggplot2 + patchwork implementation; gene track is drawn directly
# from the reference GFF3.
#
# The y-axis quantity is selectable via --y_metric:
#   rate            relative centrality (RATE) from cppRATE output
#   abs_median      absolute median posterior effect from variant_effects
#   exp_abs_median  exp() of the absolute median effect
#
# PPOM runs supply per-cutpoint inputs. All cutpoints are overlaid on a
# single plot: color = r² with lead, shape = cutpoint, purple diamond on
# the (variant, cutpoint) pair with the maximum metric value.
#
# Usage (RATE, single cutpoint):
#   Rscript make_locuszoom_plot.R \
#     --y_metric         rate \
#     --rate_files       results/cppRATE_results/RATE_values_depruned.txt \
#     --positions_file   data/variant_positions.csv \
#     --genotype_matrix  results/cppRATE_matrices/design_matrix.csv \
#     --gff              reference/genomic.gff \
#     --lead_variant     4521 \
#     --window           25000 \
#     --output           locus_plot.png
#
# Usage (RATE, multiple cutpoints overlaid):
#   --rate_files RATE_values_cutpoint1_depruned.txt,RATE_values_cutpoint2_depruned.txt,...
#
# Usage (|median| or exp|median| from variant_effects):
#   --y_metric abs_median --variant_effects fitted_model/depruned_variant_effects.csv

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
  make_option("--y_metric",
    type = "character", default = "rate",
    help = "Y-axis metric: rate | abs_median | exp_abs_median [default: rate]"
  ),
  make_option("--rate_files",
    type = "character", default = NULL,
    help = paste(
      "Comma-separated RATE files (one per cutpoint). Required when",
      "--y_metric=rate. Cutpoint is parsed from the filename token",
      "'cutpoint<K>'; if absent, cutpoint = 1."
    )
  ),
  make_option("--rate_file",
    type = "character", default = NULL,
    help = "DEPRECATED alias for --rate_files (single file)."
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
  make_option("--lead_cutpoint",
    type = "integer", default = 1L,
    help = "Cutpoint stratum of the lead variant [default: 1]"
  ),
  make_option("--region",
    type = "character", default = NULL,
    help = "Genomic region as 'seqname:start-end'"
  ),
  make_option("--window",
    type = "integer", default = 25000L,
    help = "Window size (bp) around lead variant [default: 25000]"
  ),
  make_option("--variant_effects",
    type = "character", default = NULL,
    help = paste(
      "Path to depruned_variant_effects.csv. Required when",
      "--y_metric is abs_median or exp_abs_median. Also used (when",
      "available) to mark CI-significant variants with filled triangles."
    )
  ),
  make_option("--annotations",
    type = "character", default = NULL,
    help = "Optional: SNPEff annotations TSV with POS + ANN....GENE columns"
  ),
  make_option("--genes_of_interest",
    type = "character", default = NULL,
    help = paste(
      "Optional 2-col file (annotation_name, display_name) -- e.g.",
      "spn_penicillin_genesofinterest.txt. Gene labels matching col 1 are",
      "rendered using col 2 (e.g. penA -> pbp2b)."
    )
  ),
  make_option("--genes_of_interest_gff",
    type = "character", default = NULL,
    help = paste(
      "Optional 2-col file (gff_name, display_name) used to remap GFF",
      "gene-track labels. Mirrors --genes_of_interest but keyed against the",
      "GFF Name attribute (e.g. pbp2X -> pbp2x, SPN23F_RS01610 -> clpL)."
    )
  ),
  make_option("--top_n_labels",
    type = "integer", default = 40L,
    help = "Number of top variants to label with gene names [default: 40]"
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
  option_list = option_list
))

valid_metrics <- c("rate", "abs_median", "exp_abs_median")
if (!opt$y_metric %in% valid_metrics) {
  stop("--y_metric must be one of: ", paste(valid_metrics, collapse = ", "))
}

# ---------------------------------------------------------------------------
# Validate inputs
# ---------------------------------------------------------------------------
stop_if_missing <- function(path, flag) {
  if (is.null(path)) stop(flag, " is required")
  if (!file.exists(path)) stop("File not found: ", path, " (", flag, ")")
}

stop_if_missing(opt$positions_file,  "--positions_file")
stop_if_missing(opt$genotype_matrix, "--genotype_matrix")
stop_if_missing(opt$gff,             "--gff")

# Fold --rate_file (legacy) into --rate_files
if (is.null(opt$rate_files) && !is.null(opt$rate_file)) {
  opt$rate_files <- opt$rate_file
}

if (opt$y_metric == "rate") {
  if (is.null(opt$rate_files)) stop("--rate_files is required when --y_metric=rate")
} else {
  if (is.null(opt$variant_effects)) {
    stop("--variant_effects is required when --y_metric=", opt$y_metric)
  }
}

if (is.null(opt$lead_variant) && is.null(opt$region)) {
  stop("Specify either --lead_variant (integer) or --region (seqname:start-end)")
}
if (!is.null(opt$variant_effects) && !file.exists(opt$variant_effects)) {
  stop("File not found: ", opt$variant_effects, " (--variant_effects)")
}
if (!is.null(opt$annotations) && !file.exists(opt$annotations)) {
  stop("File not found: ", opt$annotations, " (--annotations)")
}

# ---------------------------------------------------------------------------
# Step 1: Load positions
# ---------------------------------------------------------------------------
message("Loading positions: ", opt$positions_file)
pos_df    <- fread(opt$positions_file, header = TRUE)
positions <- as.integer(pos_df[[2]])
n_variants_pos <- length(positions)
message("  Loaded ", n_variants_pos, " variant positions")

# ---------------------------------------------------------------------------
# Step 2: Build long-format Y values (one row per variant per cutpoint)
# ---------------------------------------------------------------------------
parse_cutpoint <- function(path) {
  m <- regmatches(basename(path),
                  regexpr("cutpoint([0-9]+)", basename(path), perl = TRUE))
  if (!length(m) || !nzchar(m)) return(1L)
  as.integer(sub("cutpoint", "", m))
}

if (opt$y_metric == "rate") {
  rate_paths <- trimws(strsplit(opt$rate_files, ",")[[1]])
  for (p in rate_paths) stop_if_missing(p, "--rate_files entry")

  long_list <- lapply(rate_paths, function(p) {
    cp <- parse_cutpoint(p)
    rate_raw <- read.table(p, comment.char = "#", header = FALSE)
    # depruned files have 3 cols (SNP_ID, RATE, KLD); original has 1 (RATE)
    rate_col <- ifelse(ncol(rate_raw) >= 2, 2, 1)
    rate_values <- as.numeric(rate_raw[, rate_col])
    if (length(rate_values) != n_variants_pos) {
      stop("Mismatch reading ", p, ": ", length(rate_values),
           " RATE values but ", n_variants_pos, " positions")
    }
    data.table(SNP = seq_along(rate_values), BP = positions,
               cutpoint = cp, Y = rate_values)
  })
  gwas_long <- rbindlist(long_list)
  y_label   <- "relative centrality (RATE)"
  message("  Loaded RATE values for ", length(rate_paths),
          " cutpoint(s): ", paste(sort(unique(gwas_long$cutpoint)), collapse = ", "))
} else {
  message("Loading variant effects: ", opt$variant_effects)
  eff_dt <- fread(opt$variant_effects)
  need <- c("variant_id", "median")
  if (!all(need %in% names(eff_dt))) {
    stop("variant_effects file missing required columns: ",
         paste(setdiff(need, names(eff_dt)), collapse = ", "))
  }
  if (!"cutpoint" %in% names(eff_dt)) eff_dt[, cutpoint := 1L]

  abs_med <- abs(eff_dt$median)
  y_vals  <- if (opt$y_metric == "abs_median") abs_med else exp(abs_med)
  gwas_long <- data.table(
    SNP      = as.integer(eff_dt$variant_id),
    BP       = positions[as.integer(eff_dt$variant_id)],
    cutpoint = as.integer(eff_dt$cutpoint),
    Y        = y_vals
  )
  gwas_long <- gwas_long[!is.na(BP)]
  y_label <- if (opt$y_metric == "abs_median") {
    expression(abs(tilde(beta)))
  } else {
    expression("e"^{abs(tilde(beta))})
  }
  message("  Loaded ", nrow(gwas_long), " (variant, cutpoint) rows from variant_effects",
          " — ", length(unique(gwas_long$cutpoint)), " cutpoint(s)")
}

n_cutpoints <- length(unique(gwas_long$cutpoint))

# ---------------------------------------------------------------------------
# Step 3: Optional CI-significance flag (purely for the diamond layer below)
# ---------------------------------------------------------------------------
sig_variants <- NULL
if (!is.null(opt$variant_effects)) {
  eff_sig <- fread(opt$variant_effects)
  if ("signif" %in% names(eff_sig) && "variant_id" %in% names(eff_sig)) {
    sig_variants <- eff_sig[signif == TRUE, unique(variant_id)]
    message("  ", length(sig_variants), " variants with 89% CI excluding zero")
  }
}

# ---------------------------------------------------------------------------
# Step 4: Optionally load SNPEff annotations
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

# Optional genes-of-interest display-name map (col1 = annotation name,
# col2 = display name).
goi_df <- NULL
if (!is.null(opt$genes_of_interest)) {
  stop_if_missing(opt$genes_of_interest, "--genes_of_interest")
  message("Loading genes_of_interest: ", opt$genes_of_interest)
  goi_df <- fread(opt$genes_of_interest, header = FALSE, sep = ",",
                  strip.white = TRUE)
  if (ncol(goi_df) < 2) {
    warning("--genes_of_interest needs >=2 columns; ignoring file")
    goi_df <- NULL
  } else {
    setnames(goi_df, c("annot_name", "display_name")[seq_len(ncol(goi_df))])
    goi_df[, annot_name   := trimws(annot_name)]
    goi_df[, display_name := trimws(display_name)]
  }
}

apply_display_names <- function(x) {
  if (is.null(goi_df) || is.null(x)) return(x)
  m <- match(x, goi_df$annot_name)
  hit <- !is.na(m)
  x[hit] <- goi_df$display_name[m[hit]]
  x
}

# Parallel map keyed against GFF Name (e.g. pbp2X, SPN23F_RS01610). Used to
# remap gene-track labels so they match the SNP-label display names.
goi_gff_df <- NULL
if (!is.null(opt$genes_of_interest_gff)) {
  stop_if_missing(opt$genes_of_interest_gff, "--genes_of_interest_gff")
  message("Loading genes_of_interest_gff: ", opt$genes_of_interest_gff)
  goi_gff_df <- fread(opt$genes_of_interest_gff, header = FALSE, sep = ",",
                      strip.white = TRUE)
  if (ncol(goi_gff_df) < 2) {
    warning("--genes_of_interest_gff needs >=2 columns; ignoring file")
    goi_gff_df <- NULL
  } else {
    setnames(goi_gff_df,
             c("gff_name", "display_name")[seq_len(ncol(goi_gff_df))])
    goi_gff_df[, gff_name     := trimws(gff_name)]
    goi_gff_df[, display_name := trimws(display_name)]
  }
}

apply_gff_display_names <- function(x) {
  if (is.null(x)) return(x)
  if (!is.null(goi_gff_df)) {
    m <- match(x, goi_gff_df$gff_name)
    hit <- !is.na(m)
    x[hit] <- goi_gff_df$display_name[m[hit]]
  }
  # Fall back to the SNPEff-keyed map for any names that happen to overlap
  # (e.g. clpX, ciaH) so a single-file setup still works.
  apply_display_names(x)
}

italic_gene_expr <- function(x) paste0("italic('", gsub("'", "", x), "')")

# ---------------------------------------------------------------------------
# Step 5: Determine region and lead pair (variant, cutpoint)
# ---------------------------------------------------------------------------
if (!is.null(opt$region)) {
  parts        <- strsplit(opt$region, ":")[[1]]
  seqname      <- parts[1]
  coords       <- as.integer(strsplit(parts[2], "-")[[1]])
  region_start <- coords[1]
  region_end   <- coords[2]
  in_region    <- gwas_long[BP >= region_start & BP <= region_end]
  if (nrow(in_region) == 0) stop("No variants found in region: ", opt$region)
  top_row       <- in_region[which.max(Y)]
  lead_idx      <- as.integer(top_row$SNP)
  lead_pos      <- as.integer(top_row$BP)
  lead_cutpoint <- as.integer(top_row$cutpoint)
} else {
  lead_idx      <- opt$lead_variant
  lead_cutpoint <- opt$lead_cutpoint
  if (!(lead_idx %in% gwas_long$SNP)) {
    stop("Lead variant not found in y-axis data: ", lead_idx)
  }
  lead_pos     <- gwas_long[SNP == lead_idx, BP[1]]
  seqname      <- NA_character_   # filled after GFF load
  region_start <- max(1L, lead_pos - opt$window)
  region_end   <- lead_pos + opt$window
}

message(sprintf(
  "Region: %s:%d-%d  |  Lead variant: %d (pos %d, cutpoint %d)",
  if (is.na(seqname)) "?" else seqname,
  region_start, region_end, lead_idx, lead_pos, lead_cutpoint
))

regional_long <- gwas_long[BP >= region_start & BP <= region_end]
if (nrow(regional_long) == 0) stop("No variants in region after filtering")
message("  ", nrow(regional_long), " (variant, cutpoint) rows in region — ",
        length(unique(regional_long$SNP)), " distinct variants")

# ---------------------------------------------------------------------------
# Step 6: Compute r² with lead variant from genotype matrix
# ---------------------------------------------------------------------------
message("Computing r² from genotype matrix: ", opt$genotype_matrix)

# design_matrix.csv is the LD-pruned genotype matrix; column k holds the
# genotype for the k-th *representative* variant, not variant k. Use the
# companion cppRATE files to translate depruned variant IDs to columns:
#   bacprune_rust_results.csv — header row lists rep IDs in column order
#   direction_of_correlation.csv — every variant ID → its representative
geno_dir       <- dirname(opt$genotype_matrix)
bacprune_path  <- file.path(geno_dir, "bacprune_rust_results.csv")
direction_path <- file.path(geno_dir, "direction_of_correlation.csv")
for (p in c(bacprune_path, direction_path)) {
  if (!file.exists(p)) {
    stop("Required cppRATE companion file not found next to genotype matrix: ", p)
  }
}

design_col_ids <- as.integer(strsplit(readLines(bacprune_path, n = 1), ",")[[1]])

dir_dt <- fread(direction_path)
setnames(dir_dt, c("Variant", "Representative Variant"), c("variant_id", "rep_id"))
rep_of_variant <- integer(max(dir_dt$variant_id))
rep_of_variant[dir_dt$variant_id] <- as.integer(dir_dt$rep_id)

regional_ids  <- unique(regional_long$SNP)
lead_rep      <- rep_of_variant[lead_idx]
regional_reps <- rep_of_variant[regional_ids]
needed_reps   <- unique(c(lead_rep, regional_reps))
col_positions <- match(needed_reps, design_col_ids)
if (anyNA(col_positions)) {
  stop("Representative IDs missing from design_matrix header: ",
       paste(needed_reps[is.na(col_positions)], collapse = ", "))
}

# Read selected columns in file order so column-to-rep alignment is unambiguous.
ord           <- order(col_positions)
col_positions <- col_positions[ord]
needed_reps   <- needed_reps[ord]

geno_sub        <- fread(opt$genotype_matrix, header = FALSE, select = col_positions)
names(geno_sub) <- as.character(needed_reps)

lead_geno <- as.numeric(geno_sub[[as.character(lead_rep)]])

compute_r2 <- function(x, y) {
  r <- suppressWarnings(cor(x, y, use = "complete.obs"))
  if (is.na(r)) return(0)
  r^2
}

r2_per_rep <- vapply(
  as.character(needed_reps),
  function(r) compute_r2(as.numeric(geno_sub[[r]]), lead_geno),
  FUN.VALUE = numeric(1)
)
r2_dt <- data.table(
  SNP = regional_ids,
  r2  = unname(r2_per_rep[as.character(regional_reps)])
)
r2_dt[SNP == lead_idx, r2 := 1.0]

regional_long <- merge(regional_long, r2_dt, by = "SNP", all.x = TRUE)

message(sprintf(
  "  r² range (per variant): %.3f – %.3f (median %.3f)",
  min(r2_dt$r2, na.rm = TRUE),
  max(r2_dt$r2, na.rm = TRUE),
  median(r2_dt$r2, na.rm = TRUE)
))

high_ld_frac <- mean(r2_dt$r2 > 0.8, na.rm = TRUE)
if (high_ld_frac > 0.9) {
  message(sprintf(
    "  NOTE: %.0f%% of variants have r² > 0.8 (high clonality / low recombination).",
    high_ld_frac * 100
  ))
}

# ---------------------------------------------------------------------------
# Step 7: Gene-name labels on top-N variants (one label per variant)
# ---------------------------------------------------------------------------
regional_long[, gene_label := NA_character_]
if (!is.null(ann_df)) {
  gene_col <- grep("GENE", names(ann_df), value = TRUE)[1]
  if (!is.na(gene_col)) {
    regional_long[, gene_label := ann_df[[gene_col]][match(BP, ann_df$POS)]]
    regional_long[, gene_label := apply_display_names(gene_label)]
  }
}

regional_long[, significant := FALSE]
if (!is.null(sig_variants)) {
  regional_long[SNP %in% sig_variants, significant := TRUE]
}

# Pick the top-N labels within the lead-cutpoint stratum to avoid the same
# variant being labeled K times.
if (!all(is.na(regional_long$gene_label))) {
  cp_for_ranking <- if (lead_cutpoint %in% regional_long$cutpoint) {
    lead_cutpoint
  } else {
    regional_long$cutpoint[1]
  }
  top_snps <- regional_long[cutpoint == cp_for_ranking
                            ][order(-Y)
                            ][seq_len(min(opt$top_n_labels, .N)), SNP]
  regional_long[!SNP %in% top_snps, gene_label := NA_character_]
  # When a variant is labeled, keep the label only on its lead-cutpoint row
  regional_long[SNP %in% top_snps & cutpoint != cp_for_ranking,
                gene_label := NA_character_]
}

# ---------------------------------------------------------------------------
# Step 8: Parse GFF3 for gene track
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

extract_attr <- function(attrs, key) {
  pattern <- paste0("(?:^|;)", key, "=([^;]+)")
  out <- rep(NA_character_, length(attrs))
  hit <- regexpr(pattern, attrs, perl = TRUE) > 0
  if (any(hit)) {
    out[hit] <- sub(paste0(".*(?:^|;)", key, "=([^;]+).*"), "\\1",
                    attrs[hit], perl = TRUE)
  }
  out
}

if (nrow(genes_region) > 0) {
  genes_region[, name := extract_attr(attributes, "Name")]
  genes_region[is.na(name), name := extract_attr(attributes, "gene")]
  genes_region[is.na(name), name := extract_attr(attributes, "locus_tag")]
  genes_region[is.na(name), name := paste0("gene_", seq_len(sum(is.na(name))))]

  genes_region[, name := apply_gff_display_names(name)]

  genes_region[, draw_start := pmax(start, region_start)]
  genes_region[, draw_end   := pmin(end,   region_end)]
  genes_region[, mid        := (draw_start + draw_end) / 2]
}

# ---------------------------------------------------------------------------
# Step 9: LD palette (locuszoomr-style) + bin assignment
# ---------------------------------------------------------------------------
message("Building plot...")

ld_breaks <- c(0, 0.2, 0.4, 0.6, 0.8, 1.001)
ld_labels <- c("0.0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1.0")
ld_colors <- c("grey80", "lightskyblue", "darkgreen", "orange", "red")

regional_long[, ld_bin := cut(r2, breaks = ld_breaks, labels = ld_labels,
                              include.lowest = TRUE, right = TRUE)]
regional_long[is.na(r2), ld_bin := factor(NA, levels = ld_labels)]

# Stable cutpoint → shape mapping (point shape codes)
cutpoint_levels <- sort(unique(c(regional_long$cutpoint, lead_cutpoint)))
cutpoint_shapes <- c(16, 17, 15, 18, 3, 4, 8, 1, 2, 0)
if (length(cutpoint_levels) > length(cutpoint_shapes)) {
  cutpoint_shapes <- rep_len(cutpoint_shapes, length(cutpoint_levels))
}
shape_values <- setNames(cutpoint_shapes[seq_along(cutpoint_levels)],
                         as.character(cutpoint_levels))
regional_long[, cutpoint_f := factor(cutpoint, levels = cutpoint_levels)]

# ---------------------------------------------------------------------------
# Step 10: Build scatter panel
# ---------------------------------------------------------------------------
plot_title <- if (!is.null(opt$title)) {
  opt$title
} else {
  sprintf("Regional association plot — %s:%d-%d", seqname, region_start, region_end)
}

subtitle_parts <- c(
  sprintf("Lead variant: position %d (index %d, cutpoint %d)",
          lead_pos, lead_idx, lead_cutpoint),
  sprintf("n = %d variants in window", length(unique(regional_long$SNP)))
)
if (high_ld_frac > 0.9) {
  subtitle_parts <- c(subtitle_parts,
    sprintf("%.0f%% variants r² > 0.8 with lead (high clonality)",
            high_ld_frac * 100))
}

lead_row <- regional_long[SNP == lead_idx & cutpoint == lead_cutpoint]

p_scatter_base <- ggplot(
    regional_long[!(SNP == lead_idx & cutpoint == lead_cutpoint)],
    aes(x = BP, y = Y, colour = ld_bin, shape = cutpoint_f)
  ) +
  geom_point(size = 1.7, alpha = 0.85) +
  scale_colour_manual(
    values   = setNames(ld_colors, ld_labels),
    name     = expression(r^2),
    na.value = "grey70",
    drop     = FALSE
  ) +
  scale_shape_manual(
    values = shape_values,
    name   = "cutpoint",
    guide  = if (n_cutpoints > 1) "legend" else "none"
  ) +
  geom_point(
    data   = lead_row,
    aes(x  = BP, y = Y),
    shape  = 23, size = 4, fill = "purple", colour = "black",
    inherit.aes = FALSE
  ) +
  scale_x_continuous(
    limits = c(region_start, region_end),
    expand = expansion(mult = 0.01),
    labels = scales::comma
  ) +
  labs(y = y_label) +
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

p_scatter_labeled <- p_scatter_base +
  labs(
    title    = plot_title,
    subtitle = paste(subtitle_parts, collapse = "  |  ")
  )
if (!all(is.na(regional_long$gene_label))) {
  label_data <- regional_long[!is.na(gene_label)]
  label_data[, gene_expr := italic_gene_expr(gene_label)]
  p_scatter_labeled <- p_scatter_labeled +
    geom_text_repel(
      data         = as.data.frame(label_data),
      aes(x = BP, y = Y, label = gene_expr),
      parse        = TRUE,
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
# Step 11: Gene track panel
# ---------------------------------------------------------------------------
if (nrow(genes_region) > 0) {
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

  # For the gene-track tick layer, use one row per variant (the lead-cutpoint
  # stratum if available, else the first cutpoint) so we don't draw K coloured
  # ticks per variant.
  cp_for_ticks <- if (lead_cutpoint %in% regional_long$cutpoint) lead_cutpoint
                  else regional_long$cutpoint[1]
  tick_rows <- regional_long[cutpoint == cp_for_ticks]

  variant_ticks <- rbindlist(lapply(seq_len(nrow(genes_region)), function(gi) {
    g <- genes_region[gi]
    hits <- tick_rows[BP >= g$draw_start & BP <= g$draw_end,
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
      values   = setNames(ld_colors, ld_labels),
      name     = expression(r^2),
      na.value = "grey70",
      drop     = FALSE,
      guide    = "none"
    ) +
    geom_text_repel(
      aes(x = mid, y = ymax, label = name),
      size           = 2.8,
      vjust          = 0,
      nudge_y        = 0.2,
      fontface       = "italic",
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
      name   = "genomic position"
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
  p_genes <- ggplot() +
    annotate("text", x = (region_start + region_end) / 2, y = 0,
             label = "No gene features found in region", colour = "grey50") +
    scale_x_continuous(
      limits = c(region_start, region_end),
      labels = scales::comma,
      name   = "genomic position"
    ) +
    theme_bw(base_size = 11)
}

# ---------------------------------------------------------------------------
# Step 12: Assemble and save
# ---------------------------------------------------------------------------
combined_labeled   <- p_scatter_labeled / p_genes + plot_layout(heights = c(3, 1))
combined_unlabeled <- p_scatter_base    / p_genes + plot_layout(heights = c(3, 1))

ext  <- tools::file_ext(opt$output)
stem <- tools::file_path_sans_ext(opt$output)
nolabels_path <- paste0(stem, "_nolabels.", if (nzchar(ext)) ext else "png")

out_dir <- dirname(opt$output)
if (!dir.exists(out_dir) && out_dir != ".") dir.create(out_dir, recursive = TRUE)

message("Saving labeled plot to:   ", opt$output)
ggsave(
  filename = opt$output,
  plot     = combined_labeled,
  width    = opt$width,
  height   = opt$height,
  dpi      = 300
)

message("Saving unlabeled plot to: ", nolabels_path)
ggsave(
  filename = nolabels_path,
  plot     = combined_unlabeled,
  width    = opt$width,
  height   = opt$height,
  dpi      = 300
)

message("Done.")
