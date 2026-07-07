#!/usr/bin/env Rscript
# Read a per-variant beta + RATE table (output of extract_gene_beta_RATE.sh) into R.
#
# Usage:
#   Rscript read_gene_beta_RATE.R <CSV>
# Default CSV if none given:
#   ../../thesis_results/gwas_variantsofinterest/spn_penicillin_PPOM_phpA_beta_RATE_by_cutpoint.csv
#
# Returns a data.frame `dat` (wide: one row per variant, beta_cp1..N and RATE_cp1..N
# columns) and also builds `dat_long` (tidy: one row per variant x cutpoint x measure).

args <- commandArgs(trailingOnly = TRUE)
csv <- if (length(args) >= 1) args[1] else
  "/nfs/research/jlees/jacqueline/thesis_results/gwas_variantsofinterest/01_spn_penicillin_binary_logistic_top5_beta_RATE.csv"

# cutpoint -> MIC (ug/mL) for SPN penicillin (7 cutpoints)
cutpoint_mic <- c(0.016, 0.03, 0.06, 0.12, 0.25, 1, 2)

dat <- read.csv(csv, stringsAsFactors = FALSE, check.names = FALSE)

beta_cols <- grep("^beta_cp", names(dat), value = TRUE)
rate_cols <- grep("^RATE_cp", names(dat), value = TRUE)
n_cp <- length(beta_cols)

# Long/tidy form: variant x cutpoint, with beta and RATE side by side
id_cols <- c("variant_id", "POS", "effect", "impact", "HGVS_P")
dat_long <- do.call(rbind, lapply(seq_len(n_cp), function(c) {
  data.frame(
    dat[id_cols],
    cutpoint = c,
    mic      = if (c <= length(cutpoint_mic)) cutpoint_mic[c] else NA_real_,
    beta     = dat[[paste0("beta_cp", c)]],
    RATE     = dat[[paste0("RATE_cp", c)]],
    stringsAsFactors = FALSE
  )
}))

cat(sprintf("Loaded %d variants x %d cutpoints from:\n  %s\n\n",
            nrow(dat), n_cp, csv))

# Quick look: top variants by their best (max) RATE across cutpoints
dat$max_RATE  <- apply(dat[rate_cols], 1, max, na.rm = TRUE)
dat$max_absbeta <- apply(abs(dat[beta_cols]), 1, max, na.rm = TRUE)
top <- head(dat[order(-dat$max_RATE), c("variant_id", "POS", "HGVS_P", "effect", "max_RATE", "max_absbeta")], 10)
cat("Top 10 phpA variants by max RATE across cutpoints:\n")
print(top, row.names = FALSE)

# Top 30 genes by their largest RATE (the gene's single highest-RATE variant)
gene_max_RATE <- aggregate(max_RATE ~ gene, data = dat, FUN = max, na.rm = TRUE)
gene_anno <- unique(dat[c("gene", "gene_eggnog", "desc_eggnog")])
gene_max_RATE <- merge(gene_max_RATE, gene_anno, by = "gene", all.x = TRUE)
# Where an eggnog annotation exists (and differs), append it: SPN23F17610/sstT
has_eggnog <- !is.na(gene_max_RATE$gene_eggnog) & nzchar(gene_max_RATE$gene_eggnog) &
  gene_max_RATE$gene_eggnog != gene_max_RATE$gene
gene_max_RATE$gene <- ifelse(has_eggnog,
                             paste0(gene_max_RATE$gene, "/", gene_max_RATE$gene_eggnog),
                             gene_max_RATE$gene)
top_genes <- head(gene_max_RATE[order(-gene_max_RATE$max_RATE),
                                c("gene", "gene_eggnog", "desc_eggnog", "max_RATE")], 30)
cat("\nTop 30 genes by largest RATE:\n")
print(top_genes, row.names = FALSE)

invisible(list(dat = dat, dat_long = dat_long, top_genes = top_genes))

# ---------------------------------------------------------------------------
# Run the top-30-genes summary over every *beta_RATE*.csv and combine into one.
# ---------------------------------------------------------------------------
results_dir <- "/nfs/research/jlees/jacqueline/thesis_results/gwas_variantsofinterest"
combined_out <- file.path(results_dir, "all_runs_top30_genes_by_RATE.csv")

top30_genes_for_csv <- function(path) {
  d <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  # Some runs (e.g. TB rifampicin) have no eggnog annotation columns
  if (!"gene_eggnog" %in% names(d)) d$gene_eggnog <- NA_character_
  if (!"desc_eggnog" %in% names(d)) d$desc_eggnog <- NA_character_
  d$gene <- trimws(d$gene)
  d$gene_eggnog <- trimws(d$gene_eggnog)
  rc <- grep("^RATE_cp", names(d), value = TRUE)
  d$max_RATE <- apply(d[rc], 1, max, na.rm = TRUE)
  g <- aggregate(max_RATE ~ gene, data = d, FUN = max, na.rm = TRUE)
  anno <- unique(d[c("gene", "gene_eggnog", "desc_eggnog")])
  g <- merge(g, anno, by = "gene", all.x = TRUE)
  he <- !is.na(g$gene_eggnog) & nzchar(g$gene_eggnog) & g$gene_eggnog != g$gene
  g$gene <- ifelse(he, paste0(g$gene, "/", g$gene_eggnog), g$gene)
  g <- head(g[order(-g$max_RATE), c("gene", "gene_eggnog", "desc_eggnog", "max_RATE")], 30)
  # Strip the trailing "_top5_beta_RATE" so source reads like the run name
  run <- sub("_top5_beta_RATE$", "", sub("\\.csv$", "", basename(path)))
  cbind(source = run, g, row.names = NULL)
}

all_csvs <- sort(Sys.glob(file.path(results_dir, "*beta_RATE*.csv")))
combined <- do.call(rbind, lapply(all_csvs, top30_genes_for_csv))
write.csv(combined, combined_out, row.names = FALSE)
cat(sprintf("\nWrote top-30 genes for %d runs (%d rows) to:\n  %s\n",
            length(all_csvs), nrow(combined), combined_out))

gene_groups <- read.csv("/nfs/research/jlees/jacqueline/thesis_results/gwas_variantsofinterest/spn_penicillin_deprune_gene_groups.csv")
