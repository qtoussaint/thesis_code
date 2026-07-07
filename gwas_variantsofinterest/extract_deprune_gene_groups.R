#!/usr/bin/env Rscript
# Which genes were depruned together: one CSV per drug/species listing each perfect-LD
# linkage block (the variants the pipeline collapsed into a single representative) and
# the genes they fall in. ld_pruning_summary.csv is byte-identical across a species'
# runs, so one table per drug/species suffices -- taken here from each species' MIC POM
# run.
#
# For the SPN species, SPN23F locus-tag genes are also resolved to eggNOG gene symbols
# (representative_gene_eggnog, member_genes_eggnog), reusing the eggNOG-mapper annotation
# of the Spn23F proteome (eggnog/spn23f.emapper.annotations). Non-locus-tag genes keep
# their curated name; TB has no proteome run, so its file gets no eggNOG columns.
#
# Usage:
#   mamba activate gwas_pipeline
#   Rscript gwas_variantsofinterest/extract_deprune_gene_groups.R
#
# Output: one <key>_deprune_gene_groups.csv per species in
#   /nfs/research/jlees/jacqueline/thesis_results/gwas_variantsofinterest/

# annotate_genes (POS -> gene, with genes-of-interest display remapping) from the figure
# helpers; sourcing does not run its main().
script_dir <- dirname(sub("^--file=", "",
        grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(script_dir, "..", "paper_figures", "replot_rate_faceted_cutpoints.R"))

RES        <- "/nfs/research/jlees/jacqueline/thesis_results"
DATASETS   <- file.path(RES, "gwas_datasets", "inference")
SPN_ANNOT  <- "/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/genotype/fields_filtered_maf05_multiallelic.txt"
TB_ANNOT   <- "/nfs/research/jlees/jacqueline/gwas_data/tuberculosis/cryptic_regeno_snpeff/fields_filtered.txt"
GOI_DIR    <- "/nfs/research/jlees/jacqueline/thesis_code/gwas_genesofinterest"
OUT_DIR    <- file.path(RES, "gwas_variantsofinterest")
EGG_ANNOT  <- file.path(OUT_DIR, "eggnog", "spn23f.emapper.annotations")

# locus_tag -> eggNOG Preferred_name (gene symbol), from emapper's annotations file.
# Drops the leading/trailing "##" comments; the one "#query" line gives column names.
load_eggnog_symbols <- function(path) {
  if (!file.exists(path)) {
    message("  (no eggNOG annotations at ", path, " -- eggNOG columns skipped)")
    return(NULL)
  }
  ln <- readLines(path)
  hdr_i <- grep("^#[^#]", ln)[1]
  cols  <- strsplit(sub("^#", "", ln[hdr_i]), "\t")[[1]]
  body  <- ln[(hdr_i + 1):length(ln)]
  body  <- body[!grepl("^##", body) & nzchar(body)]
  df <- read.delim(text = paste(body, collapse = "\n"), header = FALSE,
                   sep = "\t", quote = "", stringsAsFactors = FALSE)
  names(df) <- cols[seq_len(ncol(df))]
  sym <- df[["Preferred_name"]]
  sym[sym == "-" | is.na(sym)] <- ""
  stats::setNames(sym, df[[1]])
}

# Replace SPN23F locus tags with their eggNOG symbol (falling back to the locus tag when
# eggNOG has none); leave non-locus-tag (curated) names untouched.
resolve_eggnog <- function(g, egg) {
  out <- g
  lt  <- !is.na(g) & grepl("^SPN23F", g)
  if (any(lt)) {
    sym <- unname(egg[g[lt]])
    fb  <- is.na(sym) | !nzchar(sym)
    sym[fb] <- g[lt][fb]
    out[lt] <- sym
  }
  out
}

egg_sym <- load_eggnog_symbols(EGG_ANNOT)

# Per species: the run whose ld_pruning_summary we read, the variant_index dataset for
# index -> position, and the annotation + genes-of-interest list for position -> gene.
species_cfgs <- list(
  list(key = "tb_rifampicin", run = "08_tb_rifampicin_MIC_POM",
       dataset = "08_tb_rifampicin_MIC", annot = TB_ANNOT,
       goi = file.path(GOI_DIR, "tb_rifampicin_genesofinterest.txt")),
  list(key = "spn_penicillin", run = "02_spn_penicillin_MIC_POM",
       dataset = "02_spn_penicillin_MIC", annot = SPN_ANNOT,
       goi = file.path(GOI_DIR, "spn_penicillin_genesofinterest.txt")),
  list(key = "spn_trimethoprim", run = "05_spn_trimethoprim_MIC_POM",
       dataset = "05_spn_trimethoprim_MIC", annot = SPN_ANNOT,
       goi = file.path(GOI_DIR, "spn_trimethoprim_genesofinterest.txt"))
)

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

process_species <- function(cfg) {
  message("[", cfg$key, "]")
  ld <- read.csv(file.path(RES, paste0("gwas_", cfg$key), "inference", cfg$run,
                           "cppRATE_matrices", "ld_pruning_summary.csv"),
                 stringsAsFactors = FALSE, check.names = FALSE)
  vi <- read.csv(file.path(DATASETS, cfg$dataset,
                           paste0(cfg$dataset, "_variant_index.csv")),
                 stringsAsFactors = FALSE)
  pos <- as.numeric(vi$position)                       # indexed by variant_id (row #)

  # gene per variant_id, via annotate_genes (same mapping + display remap as the figures)
  gene_by_id <- annotate_genes(data.frame(pos = pos), cfg$annot, cfg$goi)$gene

  rep_idx    <- as.integer(ld[[1]])
  pruned_str <- ld[[2]]

  rows <- lapply(seq_along(rep_idx), function(i) {
    pruned <- suppressWarnings(as.integer(trimws(strsplit(pruned_str[i], ",")[[1]])))
    pruned <- pruned[!is.na(pruned)]
    if (length(pruned) == 0L) return(NULL)             # no LD partners -> not a block
    members  <- c(rep_idx[i], pruned)
    genes    <- gene_by_id[members]
    real     <- genes[!is.na(genes) & genes != "MODIFIER"]
    data.frame(
      representative_variant_id = rep_idx[i],
      representative_pos        = pos[rep_idx[i]],
      representative_gene       = gene_by_id[rep_idx[i]],
      n_members                 = length(members),
      n_distinct_genes          = length(unique(real)),
      member_variant_ids        = paste(members, collapse = ","),
      member_positions          = paste(pos[members], collapse = ","),
      member_genes              = paste(unique(genes), collapse = ","),
      stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  out <- cbind(block_id = seq_len(nrow(out)), out)

  # SPN species: resolve SPN23F locus tags to eggNOG symbols, placing each new column
  # next to the raw one it derives from. TB (no proteome) keeps the raw columns only.
  egg <- if (grepl("^spn", cfg$key)) egg_sym else NULL
  if (!is.null(egg)) {
    out$representative_gene_eggnog <- resolve_eggnog(out$representative_gene, egg)
    out$member_genes_eggnog <- vapply(strsplit(out$member_genes, ","), function(gs)
      paste(unique(resolve_eggnog(gs, egg)), collapse = ","), character(1))
    out <- out[, c("block_id", "representative_variant_id", "representative_pos",
                   "representative_gene", "representative_gene_eggnog",
                   "n_members", "n_distinct_genes", "member_variant_ids",
                   "member_positions", "member_genes", "member_genes_eggnog")]
  }

  path <- file.path(OUT_DIR, paste0(cfg$key, "_deprune_gene_groups.csv"))
  write.csv(out, path, row.names = FALSE)
  message("    wrote ", basename(path), " (", nrow(out), " blocks, ",
          sum(out$n_distinct_genes >= 2), " span >1 gene)")
}

for (cfg in species_cfgs) process_species(cfg)

message("Done. CSVs in ", OUT_DIR)
