#!/usr/bin/env Rscript
# Add eggNOG-derived gene names to the SPN top5 beta/RATE CSVs. For every row whose
# `gene` is a SPN23F locus tag, look up the eggNOG-mapper annotation of that locus tag's
# Spn23F protein and add three columns (right after `gene`):
#   gene_eggnog    -- eggNOG Preferred_name (gene symbol; "" when eggNOG has none)
#   desc_eggnog    -- eggNOG free-text Description
#   eggnog_pident  -- % identity of the diamond seed-ortholog hit
# Rows with a non-locus-tag gene (already a curated name) are left blank in these columns.
# TB CSVs are not touched (no TB proteome was annotated).
#
# Usage:
#   mamba activate gwas_pipeline   # (base R is enough)
#   Rscript gwas_variantsofinterest/add_eggnog_to_csvs.R

EGG  <- "/nfs/research/jlees/jacqueline/thesis_results/gwas_variantsofinterest/eggnog"
CSVD <- "/nfs/research/jlees/jacqueline/thesis_results/gwas_variantsofinterest"
ANNOT <- file.path(EGG, "spn23f.emapper.annotations")
SEEDS <- file.path(EGG, "spn23f.emapper.seed_orthologs")

# Read an emapper TSV: drop the leading/trailing "##" comment lines, use the single
# header line (which begins with "#") for column names.
read_emapper <- function(path) {
  ln <- readLines(path)
  hdr_i <- grep("^#[^#]", ln)[1]                 # the one "#query"/"#qseqid" header line
  cols  <- strsplit(sub("^#", "", ln[hdr_i]), "\t")[[1]]
  body  <- ln[(hdr_i + 1):length(ln)]
  body  <- body[!grepl("^##", body) & nzchar(body)]
  df <- read.delim(text = paste(body, collapse = "\n"), header = FALSE,
                   sep = "\t", quote = "", stringsAsFactors = FALSE)
  names(df) <- cols[seq_len(ncol(df))]
  df
}

ann <- read_emapper(ANNOT)
seed <- read_emapper(SEEDS)

dash_na <- function(x) { x[x == "-" | is.na(x)] <- ""; x }
qcol_a <- names(ann)[1]                          # "query"
qcol_s <- names(seed)[1]                         # "qseqid"
pid_col <- grep("pident", names(seed), ignore.case = TRUE, value = TRUE)[1]

lut <- data.frame(
  locus         = ann[[qcol_a]],
  gene_eggnog   = dash_na(ann[["Preferred_name"]]),
  desc_eggnog   = dash_na(ann[["Description"]]),
  stringsAsFactors = FALSE)
lut$eggnog_pident <- seed[[pid_col]][match(lut$locus, seed[[qcol_s]])]

cat(sprintf("eggNOG annotations: %d loci, %d with a gene symbol\n",
            nrow(lut), sum(nzchar(lut$gene_eggnog))))

# SPN top5 CSVs only (exclude the 3 TB files).
files <- Sys.glob(file.path(CSVD, "*_top5_beta_RATE.csv"))
files <- files[grepl("spn_", basename(files))]

for (f in files) {
  d <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  is_lt <- grepl("^SPN23F", d$gene)
  m <- match(d$gene, lut$locus)
  d$gene_eggnog   <- ifelse(is_lt, lut$gene_eggnog[m],   "")
  d$desc_eggnog   <- ifelse(is_lt, lut$desc_eggnog[m],   "")
  d$eggnog_pident <- ifelse(is_lt, lut$eggnog_pident[m], NA_real_)
  d$gene_eggnog[is.na(d$gene_eggnog)]   <- ""
  d$desc_eggnog[is.na(d$desc_eggnog)]   <- ""

  # place the three new columns right after `gene`
  rest <- setdiff(names(d), c("gene", "gene_eggnog", "desc_eggnog", "eggnog_pident"))
  d <- d[, c("gene", "gene_eggnog", "desc_eggnog", "eggnog_pident", rest)]

  write.csv(d, f, row.names = FALSE)
  cat(sprintf("  %s: %d locus-tag rows annotated (%d got a symbol)\n",
              basename(f), sum(is_lt), sum(is_lt & nzchar(d$gene_eggnog))))
}
cat("Done.\n")
