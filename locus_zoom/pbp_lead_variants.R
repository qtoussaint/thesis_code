#!/usr/bin/env Rscript
# pbp_lead_variants.R
# For each requested gene, find the variant with the largest |median| effect
# in fitted_model/depruned_variant_effects.csv (scanning across all cutpoints)
# whose genomic position falls inside the gene's GFF coordinates. Writes a TSV
# with one row per gene that downstream wrappers feed into make_locuszoom_plot.R.

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

option_list <- list(
  make_option("--variant_effects",
    type = "character", default = NULL,
    help = "Path to fitted_model/depruned_variant_effects.csv [required]"
  ),
  make_option("--positions_file",
    type = "character", default = NULL,
    help = "Path to variant index CSV (column 2 = base-pair position) [required]"
  ),
  make_option("--gff",
    type = "character", default = NULL,
    help = "Path to reference GFF3 [required]"
  ),
  make_option("--genes",
    type = "character", default = "pbp2X,pbp1a,pbp1b,pbp2a,pbp2b",
    help = paste(
      "Comma-separated GFF Name= values to look up.",
      "Default: top 5 pbp genes for spn_penicillin."
    )
  ),
  make_option("--output",
    type = "character", default = "pbp_lead_variants.tsv",
    help = "Output TSV path [default: pbp_lead_variants.tsv]"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

stop_if_missing <- function(path, flag) {
  if (is.null(path)) stop(flag, " is required")
  if (!file.exists(path)) stop("File not found: ", path, " (", flag, ")")
}
stop_if_missing(opt$variant_effects, "--variant_effects")
stop_if_missing(opt$positions_file,  "--positions_file")
stop_if_missing(opt$gff,             "--gff")

requested_genes <- trimws(strsplit(opt$genes, ",")[[1]])
if (!length(requested_genes)) stop("--genes is empty")

# ---------------------------------------------------------------------------
# Parse GFF for gene coordinates (same approach as make_locuszoom_plot.R)
# ---------------------------------------------------------------------------
message("Parsing GFF: ", opt$gff)
gff_cols <- c("seqname", "source", "feature", "start", "end",
              "score", "strand", "frame", "attributes")
gff_raw <- fread(
  cmd          = paste0("grep -v '^#' ", shQuote(opt$gff)),
  sep          = "\t",
  header       = FALSE,
  col.names    = gff_cols,
  showProgress = FALSE
)

extract_attr <- function(attrs, key) {
  pattern <- paste0("(?:^|;)", key, "=([^;]+)")
  m <- regmatches(attrs, regexpr(pattern, attrs, perl = TRUE))
  ifelse(nzchar(m),
         sub(paste0(".*", key, "=([^;]+).*"), "\\1", m),
         NA_character_)
}

genes_dt <- gff_raw[feature == "gene"]
genes_dt[, name := extract_attr(attributes, "Name")]
genes_dt[is.na(name), name := extract_attr(attributes[is.na(name)], "gene")]
genes_dt[is.na(name), name := extract_attr(attributes[is.na(name)], "locus_tag")]

# Case-insensitive match against the request list
genes_dt[, name_lc := tolower(name)]
req_lc <- tolower(requested_genes)
matched <- genes_dt[name_lc %in% req_lc,
                    .(name, seqname, start = as.integer(start), end = as.integer(end))]

missing_genes <- setdiff(req_lc, tolower(matched$name))
if (length(missing_genes)) {
  warning("Gene(s) not found in GFF: ", paste(missing_genes, collapse = ", "))
}

# ---------------------------------------------------------------------------
# Positions and variant effects
# ---------------------------------------------------------------------------
message("Loading positions: ", opt$positions_file)
pos_df <- fread(opt$positions_file, header = TRUE)
positions_dt <- data.table(
  variant_id = seq_len(nrow(pos_df)),
  BP         = as.integer(pos_df[[2]])
)

message("Loading variant effects: ", opt$variant_effects)
eff_dt <- fread(opt$variant_effects)
need <- c("variant_id", "median", "cutpoint")
if (!all(need %in% names(eff_dt))) {
  stop("variant_effects file missing required columns: ",
       paste(setdiff(need, names(eff_dt)), collapse = ", "))
}
eff_dt[, abs_median := abs(median)]
eff_dt <- merge(eff_dt, positions_dt, by = "variant_id", all.x = TRUE)

# ---------------------------------------------------------------------------
# Pick the (variant, cutpoint) with max |median| inside each gene
# ---------------------------------------------------------------------------
results <- rbindlist(lapply(seq_len(nrow(matched)), function(i) {
  g <- matched[i]
  hits <- eff_dt[!is.na(BP) & BP >= g$start & BP <= g$end]
  if (nrow(hits) == 0) {
    warning("No variants inside ", g$name, " [", g$start, "-", g$end, "]")
    return(data.table(
      gene = g$name, lead_variant = NA_integer_, lead_pos = NA_integer_,
      lead_cutpoint = NA_integer_, lead_abs_median = NA_real_
    ))
  }
  top <- hits[which.max(abs_median)]
  data.table(
    gene            = g$name,
    lead_variant    = as.integer(top$variant_id),
    lead_pos        = as.integer(top$BP),
    lead_cutpoint   = as.integer(top$cutpoint),
    lead_abs_median = top$abs_median
  )
}))

# Preserve the order the user requested where possible
results <- results[match(tolower(requested_genes), tolower(gene), nomatch = 0)]

out_dir <- dirname(opt$output)
if (!dir.exists(out_dir) && out_dir != ".") dir.create(out_dir, recursive = TRUE)
fwrite(results, opt$output, sep = "\t")
message("Wrote ", nrow(results), " rows to ", opt$output)
print(results)
