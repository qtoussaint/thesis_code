#!/usr/bin/env Rscript
# Add eggNOG-derived gene names to a TB variants-of-interest CSV, mirroring
# add_eggnog_to_csvs.R (SPN). The TB CSV `gene` column mixes RefSeq gene symbols
# (gyrA, katG, ...) and Rv locus tags; both are resolved to an Rv locus tag (symbols via
# the H37Rv RefSeq GFF), then looked up in the H37Rv eggNOG-mapper annotation. Three
# columns are inserted right after `gene`:
#   gene_eggnog    -- eggNOG Preferred_name (gene symbol; "" when eggNOG has none)
#   desc_eggnog    -- eggNOG free-text Description
#   eggnog_pident  -- % identity of the diamond seed-ortholog hit
#
# Usage (after run_eggnog_h37rv.sh finishes):
#   mamba activate gwas_pipeline
#   Rscript add_eggnog_to_tb_csv.R [<csv> ...]
# Default CSV: tb_rifampicin_extra_genes_by_position_with_alleles.csv

EGG   <- "/nfs/research/jlees/jacqueline/thesis_results/gwas_variantsofinterest/eggnog_tb"
CSVD  <- "/nfs/research/jlees/jacqueline/thesis_results/gwas_variantsofinterest"
ANNOT <- file.path(EGG, "h37rv.emapper.annotations")
SEEDS <- file.path(EGG, "h37rv.emapper.seed_orthologs")
GFF   <- "/nfs/research/jlees/jacqueline/gwas_data/tuberculosis/genotype/tb_ref_NC_000962.3.gff3"

args  <- commandArgs(trailingOnly = TRUE)
files <- if (length(args)) args else
  file.path(CSVD, "tb_rifampicin_extra_genes_by_position_with_alleles.csv")

stopifnot(file.exists(ANNOT), file.exists(SEEDS))

# emapper TSV reader (drops "##" comments, uses the single "#query"/"#qseqid" header).
read_emapper <- function(path) {
  ln <- readLines(path)
  hdr_i <- grep("^#[^#]", ln)[1]
  cols  <- strsplit(sub("^#", "", ln[hdr_i]), "\t")[[1]]
  body  <- ln[(hdr_i + 1):length(ln)]
  body  <- body[!grepl("^##", body) & nzchar(body)]
  df <- read.delim(text = paste(body, collapse = "\n"), header = FALSE,
                   sep = "\t", quote = "", stringsAsFactors = FALSE)
  names(df) <- cols[seq_len(ncol(df))]
  df
}

ann  <- read_emapper(ANNOT)
seed <- read_emapper(SEEDS)
dash_na <- function(x) { x[x == "-" | is.na(x)] <- ""; x }
pid_col <- grep("pident", names(seed), ignore.case = TRUE, value = TRUE)[1]

lut <- data.frame(
  locus         = ann[[1]],
  gene_eggnog   = dash_na(ann[["Preferred_name"]]),
  desc_eggnog   = dash_na(ann[["Description"]]),
  stringsAsFactors = FALSE)
lut$eggnog_pident <- seed[[pid_col]][match(lut$locus, seed[[1]])]

# gene symbol -> Rv locus tag, from the RefSeq GFF `gene` features.
gff <- readLines(GFF); gff <- gff[!grepl("^#", gff)]
fld <- strsplit(gff, "\t")
is_gene <- vapply(fld, function(x) length(x) >= 3 && x[3] == "gene", logical(1))
attr9 <- vapply(fld[is_gene], function(x) x[9], character(1))
pull <- function(x, k) {                       # element-wise; NA when key absent
  m <- regmatches(x, regexpr(paste0(k, "=[^;]+"), x))
  out <- rep(NA_character_, length(x)); hit <- regexpr(paste0(k, "=[^;]+"), x) > 0
  out[hit] <- sub(paste0(k, "="), "", m); out
}
sym2lt <- stats::setNames(pull(attr9, "locus_tag"), pull(attr9, "gene"))
sym2lt <- sym2lt[!is.na(names(sym2lt)) & !is.na(sym2lt)]

cat(sprintf("eggNOG: %d loci (%d with a symbol); GFF symbol map: %d genes\n",
            nrow(lut), sum(nzchar(lut$gene_eggnog)), length(sym2lt)))

for (f in files) {
  d <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  g <- trimws(d$gene)
  locus <- ifelse(grepl("^Rv", g), g, unname(sym2lt[g]))   # Rv tag direct, else via GFF
  m <- match(locus, lut$locus)
  d$gene_eggnog   <- dash_na(lut$gene_eggnog[m])
  d$desc_eggnog   <- dash_na(lut$desc_eggnog[m])
  d$eggnog_pident <- lut$eggnog_pident[m]

  rest <- setdiff(names(d), c("gene", "gene_eggnog", "desc_eggnog", "eggnog_pident"))
  d <- d[, c("gene", "gene_eggnog", "desc_eggnog", "eggnog_pident", rest)]
  write.csv(d, f, row.names = FALSE)

  na_locus <- unique(g[is.na(locus)])
  cat(sprintf("  %s: %d/%d rows annotated%s\n", basename(f),
              sum(!is.na(m)), nrow(d),
              if (length(na_locus)) paste0(" (no locus for: ", paste(na_locus, collapse = ", "), ")") else ""))
}
cat("Done.\n")
