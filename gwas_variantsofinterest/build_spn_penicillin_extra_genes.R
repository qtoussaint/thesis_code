#!/usr/bin/env Rscript
# Strongest-variant-per-gene x run position table for a set of extra SPN penicillin
# cell-wall / cell-division genes (ddlA, mraY, mraW, ftsL, gpsB, mapZ, recU, murF, clpL,
# recR, stk1, fmt), same format as tb_rifampicin_extra_genes_by_position_with_alleles.csv:
# POS + allele context + model + binning, with eggNOG columns after `gene`.
#
# Reuses build_pos_table.R's logic (reads the per-run top5 CSVs; picks the highest-RATE
# variant per gene per run; summarises every allele at that POS from the snpEff annotation),
# then adds the three eggNOG columns like add_eggnog_to_csvs.R.
#
# Usage:
#   mamba activate gwas_pipeline
#   Rscript gwas_variantsofinterest/build_spn_penicillin_extra_genes.R

OUT   <- "/nfs/research/jlees/jacqueline/thesis_results/gwas_variantsofinterest"
ANN_F <- "/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/genotype/fields_filtered_maf05_multiallelic.txt"
EGG   <- file.path(OUT, "eggnog")

runs <- c("01_spn_penicillin_binary_logistic",
          "02_spn_penicillin_MIC_POM","02_spn_penicillin_MIC_PPOM",
          "10_spn_penicillin_MIC_coarse_dilutions_POM","10_spn_penicillin_MIC_coarse_dilutions_PPOM",
          "11_spn_penicillin_MIC_large_minbin_POM","11_spn_penicillin_MIC_large_minbin_PPOM",
          "16_spn_penicillin_MIC_minimabinning_POM","16_spn_penicillin_MIC_minimabinning_PPOM")

# Each extra gene: the label shown in `gene`, the Spn23F locus tag for eggNOG lookup, and
# the name(s) it appears under in the top5 CSV `gene` column (symbol for most, locus tag
# for mapZ, which the reference GFF names only by its tag).
genes <- list(
  list(label = "ddlA", locus = "SPN23F16720", match = "ddlA"),
  list(label = "mraY", locus = "SPN23F03090", match = "mraY"),
  list(label = "mraW", locus = "SPN23F03060", match = "mraW"),
  list(label = "ftsL", locus = "SPN23F03070", match = "ftsL"),
  list(label = "gpsB", locus = "SPN23F03440", match = "gpsB"),
  list(label = "mapZ", locus = "SPN23F03460", match = "SPN23F03460"),
  list(label = "recU", locus = "SPN23F03420", match = "recU"),
  list(label = "murF", locus = "SPN23F16710", match = "murF"),
  list(label = "clpL", locus = "SPN23F03110", match = "clpL"),
  list(label = "recR", locus = "SPN23F16730", match = "recR"),
  list(label = "stk1", locus = "SPN23F17350", match = "stk1"),
  list(label = "fmt",  locus = "SPN23F17380", match = "fmt")
)
GENE_LEVELS <- vapply(genes, `[[`, character(1), "label")

# ---- annotation allele context per POS (same as build_pos_table.R) ----------
ann <- read.delim(ANN_F, stringsAsFactors = FALSE, check.names = FALSE)
names(ann) <- c("POS","REF","ALT","effect","impact","gene","HGVS_C","HGVS_P","LOF_gene")
first_res <- function(x) sub("^(p\\.)?([A-Za-z]+[0-9]+).*", "\\2", x)
pos_summary <- function(p) {
  sub <- ann[ann$POS == p, ]
  rep_hgvs  <- sub("^p\\.", "", sub$HGVS_P[1])
  uniq_hgvs <- sub("^p\\.", "", unique(sub$HGVS_P))
  variants  <- if (length(uniq_hgvs) > 3) "high polymorphism" else paste(uniq_hgvs, collapse = ";")
  data.frame(POS = p, ref_residue = first_res(rep_hgvs), variants = variants,
             n_alleles_nt = nrow(sub), n_alleles_aa = length(uniq_hgvs),
             effects = paste(sort(unique(unlist(strsplit(sub$effect, "&")))), collapse = ";"),
             all_HGVS_P = paste(uniq_hgvs, collapse = ";"), stringsAsFactors = FALSE)
}

# model / binning from the run name
model_of   <- function(run) ifelse(grepl("logistic", run), "logistic",
                            ifelse(grepl("_PPOM$", run), "PPOM", "POM"))
binning_of <- function(run) ifelse(grepl("logistic", run), "--",
                            ifelse(grepl("coarse_dilutions", run), "4-fold (>=5%)",
                             ifelse(grepl("large_minbin", run), "doubling (>=10%)",
                              ifelse(grepl("minimabinning", run), "minima breakpoints (K=5)",
                               "doubling (>=5%)"))))
BIN_LEVELS <- c("--","doubling (>=5%)","doubling (>=10%)","4-fold (>=5%)","minima breakpoints (K=5)")

# strongest variant per gene (highest RATE across cutpoints) for one run
strongest <- function(d, run) {
  d$gene <- trimws(d$gene)
  rc <- grep("^RATE_cp", names(d), value = TRUE); bc <- grep("^beta_cp", names(d), value = TRUE)
  do.call(rbind, lapply(genes, function(g) {
    dg <- d[d$gene %in% g$match, , drop = FALSE]
    if (nrow(dg) == 0) return(NULL)
    rmat <- as.matrix(dg[rc]); rmat[is.na(rmat)] <- -Inf
    ix <- which(rmat == max(rmat), arr.ind = TRUE)[1, ]; r <- ix[1]; cp <- ix[2]
    data.frame(run = run, gene = g$label, POS = dg$POS[r], best_cp = cp,
               RATE = dg[[rc[cp]]][r], beta = dg[[bc[cp]]][r], stringsAsFactors = FALSE)
  }))
}

# ---- eggNOG lookup by locus tag (same source as add_eggnog_to_csvs.R) -------
read_emapper <- function(path) {
  ln <- readLines(path)
  hdr_i <- grep("^#[^#]", ln)[1]
  cols  <- strsplit(sub("^#", "", ln[hdr_i]), "\t")[[1]]
  body  <- ln[(hdr_i + 1):length(ln)]
  body  <- body[!grepl("^##", body) & nzchar(body)]
  df <- read.delim(text = paste(body, collapse = "\n"), header = FALSE,
                   sep = "\t", quote = "", stringsAsFactors = FALSE)
  names(df) <- cols[seq_len(ncol(df))]; df
}
ea <- read_emapper(file.path(EGG, "spn23f.emapper.annotations"))
es <- read_emapper(file.path(EGG, "spn23f.emapper.seed_orthologs"))
dash_na <- function(x) { x[x == "-" | is.na(x)] <- ""; x }
pid_col <- grep("pident", names(es), ignore.case = TRUE, value = TRUE)[1]
egg_lut <- data.frame(locus = ea[[1]],
                      gene_eggnog = dash_na(ea[["Preferred_name"]]),
                      desc_eggnog = dash_na(ea[["Description"]]), stringsAsFactors = FALSE)
egg_lut$eggnog_pident <- es[[pid_col]][match(egg_lut$locus, es[[1]])]
locus_of <- stats::setNames(vapply(genes, `[[`, character(1), "locus"), GENE_LEVELS)

# ---- build --------------------------------------------------------------------
res <- do.call(rbind, lapply(runs, function(run) {
  f <- file.path(OUT, paste0(run, "_top5_beta_RATE.csv"))
  strongest(read.csv(f, check.names = FALSE, stringsAsFactors = FALSE), run)
}))
miss <- setdiff(GENE_LEVELS, unique(res$gene))
if (length(miss)) message("NOT FOUND: ", paste(miss, collapse = ", "))

summ <- do.call(rbind, lapply(unique(res$POS), pos_summary))
out  <- merge(res, summ, by = "POS", all.x = TRUE)
out$model   <- model_of(out$run)
out$binning <- binning_of(out$run)

# eggNOG columns keyed on the gene's locus tag
m <- match(locus_of[out$gene], egg_lut$locus)
out$gene_eggnog   <- egg_lut$gene_eggnog[m]
out$desc_eggnog   <- egg_lut$desc_eggnog[m]
out$eggnog_pident <- egg_lut$eggnog_pident[m]

out$gene    <- factor(out$gene, levels = GENE_LEVELS)     # preserve requested order
out$binning <- factor(out$binning, levels = BIN_LEVELS)
out$model   <- factor(out$model, levels = c("logistic","POM","PPOM"))
out <- out[order(out$gene, out$binning, out$model), ]
out$gene <- as.character(out$gene); out$binning <- as.character(out$binning)
out$model <- as.character(out$model)

out <- out[, c("gene","gene_eggnog","desc_eggnog","eggnog_pident",
               "binning","model","run","POS","best_cp","RATE","beta",
               "ref_residue","variants","n_alleles_nt","n_alleles_aa","effects","all_HGVS_P")]
path <- file.path(OUT, "spn_penicillin_extra_genes_by_position_with_alleles.csv")
write.csv(out, path, row.names = FALSE)
message("wrote ", basename(path), " (", nrow(out), " rows, ",
        length(unique(out$gene)), " genes)")
