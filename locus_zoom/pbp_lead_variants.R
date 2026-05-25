#!/usr/bin/env Rscript
# pbp_lead_variants.R
# Identify the lead (variant, cutpoint) pair per gene (or genome-wide) by the
# chosen y-axis metric. Writes a TSV consumed by the locus-zoom wrappers.
#
# Per-gene mode (default): supply --gff + --genes; one row per requested gene.
# Genome-wide mode:        supply --whole_genome; one row, the global top hit.

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
    help = "Path to reference GFF3 (required unless --whole_genome)"
  ),
  make_option("--genes",
    type = "character", default = "pbp2X,pbp1a,pbp1b,pbp2a,pbp2b",
    help = paste(
      "Comma-separated GFF Name= values to look up.",
      "Default: top 5 pbp genes for spn_penicillin."
    )
  ),
  make_option("--y_metric",
    type = "character", default = "abs_median",
    help = "Lead-selection metric: rate | abs_median | exp_abs_median [default: abs_median]"
  ),
  make_option("--rate_dir",
    type = "character", default = NULL,
    help = paste(
      "Directory containing RATE_values_cutpoint*_depruned.txt files",
      "(or a single RATE_values_depruned.txt). Required when --y_metric=rate."
    )
  ),
  make_option("--whole_genome",
    action = "store_true", default = FALSE,
    help = "Skip per-gene matching; return the global top (variant, cutpoint) pair."
  ),
  make_option("--window",
    type = "integer", default = 25000L,
    help = paste(
      "Per-gene search half-window in bp. Lead variant is the max-metric SNP",
      "in [gene_start - window, gene_end + window]. [default: 25000]"
    )
  ),
  make_option("--output",
    type = "character", default = "pbp_lead_variants.tsv",
    help = "Output TSV path [default: pbp_lead_variants.tsv]"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

valid_metrics <- c("rate", "abs_median", "exp_abs_median")
if (!opt$y_metric %in% valid_metrics) {
  stop("--y_metric must be one of: ", paste(valid_metrics, collapse = ", "))
}

stop_if_missing <- function(path, flag) {
  if (is.null(path)) stop(flag, " is required")
  if (!file.exists(path)) stop("File not found: ", path, " (", flag, ")")
}
stop_if_missing(opt$positions_file, "--positions_file")
if (opt$y_metric != "rate") stop_if_missing(opt$variant_effects, "--variant_effects")
if (opt$y_metric == "rate") stop_if_missing(opt$rate_dir, "--rate_dir")
if (!opt$whole_genome) stop_if_missing(opt$gff, "--gff")

# ---------------------------------------------------------------------------
# Positions
# ---------------------------------------------------------------------------
message("Loading positions: ", opt$positions_file)
pos_df <- fread(opt$positions_file, header = TRUE)
positions_dt <- data.table(
  variant_id = seq_len(nrow(pos_df)),
  BP         = as.integer(pos_df[[2]])
)

# ---------------------------------------------------------------------------
# Effect table (variant_id, cutpoint, BP, metric_value)
# ---------------------------------------------------------------------------
if (opt$y_metric == "rate") {
  rate_paths <- list.files(opt$rate_dir,
                           pattern = "^RATE_values(_cutpoint[0-9]+)?_depruned\\.txt$",
                           full.names = TRUE)
  if (!length(rate_paths)) {
    rate_paths <- list.files(opt$rate_dir, pattern = "^RATE_values.*\\.txt$",
                             full.names = TRUE)
  }
  if (!length(rate_paths)) stop("No RATE files found in ", opt$rate_dir)
  message("Loading RATE values from ", length(rate_paths), " file(s)")

  parse_cutpoint <- function(path) {
    m <- regmatches(basename(path),
                    regexpr("cutpoint([0-9]+)", basename(path), perl = TRUE))
    if (!length(m) || !nzchar(m)) return(1L)
    as.integer(sub("cutpoint", "", m))
  }

  metric_dt <- rbindlist(lapply(rate_paths, function(p) {
    cp <- parse_cutpoint(p)
    raw <- read.table(p, comment.char = "#", header = FALSE)
    val_col <- ifelse(ncol(raw) >= 2, 2, 1)
    rate_values <- as.numeric(raw[, val_col])
    if (length(rate_values) != nrow(positions_dt)) {
      stop("Mismatch reading ", p, ": ", length(rate_values),
           " RATE values but ", nrow(positions_dt), " positions")
    }
    data.table(
      variant_id   = seq_along(rate_values),
      cutpoint     = cp,
      BP           = positions_dt$BP,
      metric_value = rate_values
    )
  }))
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
  metric_dt <- data.table(
    variant_id   = as.integer(eff_dt$variant_id),
    cutpoint     = as.integer(eff_dt$cutpoint),
    metric_value = if (opt$y_metric == "abs_median") abs_med else exp(abs_med)
  )
  metric_dt <- merge(metric_dt, positions_dt, by = "variant_id", all.x = TRUE)
}

# ---------------------------------------------------------------------------
# Whole-genome mode: single best (variant, cutpoint)
# ---------------------------------------------------------------------------
if (opt$whole_genome) {
  hits <- metric_dt[!is.na(BP)]
  if (!nrow(hits)) stop("No variants with positions in metric table")
  top <- hits[which.max(metric_value)]
  results <- data.table(
    gene             = "genome",
    lead_variant     = as.integer(top$variant_id),
    lead_pos         = as.integer(top$BP),
    lead_cutpoint    = as.integer(top$cutpoint),
    lead_metric      = opt$y_metric,
    lead_metric_value = top$metric_value
  )
} else {
  # ---------------------------------------------------------------------------
  # Parse GFF for gene coordinates
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

  requested_genes <- trimws(strsplit(opt$genes, ",")[[1]])
  if (!length(requested_genes)) stop("--genes is empty")

  genes_dt <- gff_raw[feature == "gene"]
  genes_dt[, name := extract_attr(attributes, "Name")]
  genes_dt[is.na(name), name := extract_attr(attributes[is.na(name)], "gene")]
  genes_dt[is.na(name), name := extract_attr(attributes[is.na(name)], "locus_tag")]

  genes_dt[, name_lc := tolower(name)]
  req_lc <- tolower(requested_genes)
  matched <- genes_dt[name_lc %in% req_lc,
                      .(name, seqname, start = as.integer(start), end = as.integer(end))]

  missing_genes <- setdiff(req_lc, tolower(matched$name))
  if (length(missing_genes)) {
    warning("Gene(s) not found in GFF: ", paste(missing_genes, collapse = ", "))
  }

  # -------------------------------------------------------------------------
  # Pick (variant, cutpoint) with max metric inside each gene
  # -------------------------------------------------------------------------
  if (nrow(matched) == 0) {
    results <- data.table(
      gene = character(), lead_variant = integer(), lead_pos = integer(),
      lead_cutpoint = integer(), lead_metric = character(),
      lead_metric_value = numeric()
    )
    out_dir <- dirname(opt$output)
    if (!dir.exists(out_dir) && out_dir != ".") dir.create(out_dir, recursive = TRUE)
    fwrite(results, opt$output, sep = "\t")
    stop("No requested genes matched in GFF; wrote empty TSV to ", opt$output)
  }

  results <- rbindlist(lapply(seq_len(nrow(matched)), function(i) {
    g <- matched[i]
    win_lo <- g$start - opt$window
    win_hi <- g$end   + opt$window
    hits <- metric_dt[!is.na(BP) & BP >= win_lo & BP <= win_hi]
    if (nrow(hits) == 0) {
      warning("No variants within ", opt$window, " bp of ", g$name,
              " [", g$start, "-", g$end, "]")
      return(data.table(
        gene = g$name, lead_variant = NA_integer_, lead_pos = NA_integer_,
        lead_cutpoint = NA_integer_, lead_metric = opt$y_metric,
        lead_metric_value = NA_real_
      ))
    }
    top <- hits[which.max(metric_value)]
    data.table(
      gene             = g$name,
      lead_variant     = as.integer(top$variant_id),
      lead_pos         = as.integer(top$BP),
      lead_cutpoint    = as.integer(top$cutpoint),
      lead_metric      = opt$y_metric,
      lead_metric_value = top$metric_value
    )
  }))

  results <- results[match(tolower(requested_genes), tolower(gene), nomatch = 0)]
}

out_dir <- dirname(opt$output)
if (!dir.exists(out_dir) && out_dir != ".") dir.create(out_dir, recursive = TRUE)
fwrite(results, opt$output, sep = "\t")
message("Wrote ", nrow(results), " rows to ", opt$output)
print(results)
