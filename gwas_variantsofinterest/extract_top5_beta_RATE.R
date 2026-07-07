#!/usr/bin/env Rscript
# Per-run beta + RATE tables behind the faceted_rate_beta_summary.R figures, one CSV
# per GWAS run (logistic, POM, PPOM; unitigs excluded). Generalises
# extract_gene_beta_RATE.sh from one gene (phpA) to every mapped gene, keeping the
# top 5 variants in each gene -- the union of the 5 with the largest RATE and the 5
# with the largest |median beta| (across cutpoints for PPOM). Same wide layout as
# read_gene_beta_RATE.R, with a leading `gene` column added.
#
# Usage:
#   mamba activate gwas_pipeline
#   Rscript gwas_variantsofinterest/extract_top5_beta_RATE.R
#
# Output: one <run>_top5_beta_RATE.csv per run in
#   /nfs/research/jlees/jacqueline/thesis_results/gwas_variantsofinterest/

# Reuse read_rate / annotate_genes from the figure helpers (sourcing does not run its
# main()). Path is resolved relative to this script, like faceted_rate_beta_summary.R.
script_dir <- dirname(sub("^--file=", "",
        grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(script_dir, "..", "paper_figures", "replot_rate_faceted_cutpoints.R"))

RES        <- "/nfs/research/jlees/jacqueline/thesis_results"
DATASETS   <- file.path(RES, "gwas_datasets", "inference")
SPN_ANNOT  <- "/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/genotype/fields_filtered_maf05_multiallelic.txt"
TB_ANNOT   <- "/nfs/research/jlees/jacqueline/gwas_data/tuberculosis/cryptic_regeno_snpeff/fields_filtered.txt"
GOI_DIR    <- "/nfs/research/jlees/jacqueline/thesis_code/gwas_genesofinterest"
OUT_DIR    <- file.path(RES, "gwas_variantsofinterest")
TOP_N      <- 5L

# Species config mirrors faceted_rate_beta_summary.R (unitig runs 17/18 excluded). Each
# binning stub already carries its NN prefix; the POM/PPOM run dirs append _POM/_PPOM.
species_cfgs <- list(
  list(key = "tb_rifampicin", annot = TB_ANNOT,
       goi = file.path(GOI_DIR, "tb_rifampicin_genesofinterest.txt"),
       logistic = "07_tb_rifampicin_binary_logistic",
       binnings = c("08_tb_rifampicin_MIC",
                    "14_tb_rifampicin_MIC_coarse_dilutions",
                    "15_tb_rifampicin_MIC_large_minbin")),
  list(key = "spn_penicillin", annot = SPN_ANNOT,
       goi = file.path(GOI_DIR, "spn_penicillin_genesofinterest.txt"),
       logistic = "01_spn_penicillin_binary_logistic",
       binnings = c("02_spn_penicillin_MIC",
                    "10_spn_penicillin_MIC_coarse_dilutions",
                    "11_spn_penicillin_MIC_large_minbin",
                    "16_spn_penicillin_MIC_minimabinning")),
  list(key = "spn_trimethoprim", annot = SPN_ANNOT,
       goi = file.path(GOI_DIR, "spn_trimethoprim_genesofinterest.txt"),
       logistic = "04_spn_trimethoprim_binary_logistic",
       binnings = c("05_spn_trimethoprim_MIC",
                    "12_spn_trimethoprim_MIC_coarse_dilutions",
                    "13_spn_trimethoprim_MIC_large_minbin"))
)

run_dir   <- function(key, run) file.path(RES, paste0("gwas_", key), "inference", run)
# variant_index lives in gwas_datasets under the run name minus its model suffix.
varidx_path <- function(run) {
  ds <- sub("_(POM|PPOM|logistic)$", "", run)
  file.path(DATASETS, ds, paste0(ds, "_variant_index.csv"))
}

# POS for every variant_id (row order = variant_id, 1-based). Used in place of the
# phandango BP column, which is degenerate for the TB POM/PPOM runs.
read_positions_varidx <- function(run) {
  vi <- read.csv(varidx_path(run), stringsAsFactors = FALSE)
  as.numeric(vi$position)
}

# Cutpoint RATE files for a PPOM run, ordered numerically by cutpoint index.
cutpoint_rate_files <- function(rdir) {
  f <- Sys.glob(file.path(rdir, "cppRATE_results", "RATE_values_cutpoint*_depruned.txt"))
  if (length(f) == 0L) return(character(0))
  n <- as.integer(sub(".*cutpoint([0-9]+)_depruned\\.txt$", "\\1", f))
  f[order(n)]
}

# Build the wide per-variant frame (variant_id, pos, beta_cp1..N, RATE_cp1..N) for one
# run. PPOM -> one column per cutpoint; logistic/POM -> a single cp1. Returns NULL when
# the run has no RATE output yet. beta is NA when depruned_variant_effects.csv is absent
# (POM runs frequently skip it), leaving those columns NA.
build_run_wide <- function(rdir, run) {
  pos <- read_positions_varidx(run)
  n_var <- length(pos)
  eff_path <- file.path(rdir, "fitted_model", "depruned_variant_effects.csv")
  eff <- if (file.exists(eff_path)) read.csv(eff_path, stringsAsFactors = FALSE) else NULL

  cp_files <- cutpoint_rate_files(rdir)
  single_rate <- file.path(rdir, "cppRATE_results", "RATE_values_depruned.txt")

  if (length(cp_files) > 0L) {                       # PPOM
    n_cp <- length(cp_files)
    rate_mat <- matrix(NA_real_, n_var, n_cp)
    beta_mat <- matrix(NA_real_, n_var, n_cp)
    mic <- rep(NA_character_, n_cp)
    for (c in seq_len(n_cp)) {
      r <- read_rate(cp_files[c])
      if (length(r) != n_var)
        stop("Cutpoint ", c, " RATE length ", length(r), " != ", n_var, " in ", run)
      rate_mat[, c] <- r
      if (!is.null(eff)) {
        b <- eff$median[eff$cutpoint == c]
        if (length(b) != n_var)
          stop("Cutpoint ", c, " beta length ", length(b), " != ", n_var, " in ", run)
        beta_mat[, c] <- b
        m <- unique(eff$cutpoint_MIC[eff$cutpoint == c])
        if (length(m) == 1L) mic[c] <- as.character(m)
      }
    }
    if (!all(is.na(mic)))
      message("    cutpoint MIC (ug/mL): ",
              paste(sprintf("cp%d=%s", seq_len(n_cp), mic), collapse = ", "))
  } else if (file.exists(single_rate)) {             # logistic / POM
    n_cp <- 1L
    r <- read_rate(single_rate)
    if (length(r) != n_var)
      stop("RATE length ", length(r), " != ", n_var, " in ", run)
    rate_mat <- matrix(r, ncol = 1L)
    beta_mat <- matrix(if (!is.null(eff)) eff$median else NA_real_, ncol = 1L)
    if (!is.null(eff) && nrow(eff) != n_var)
      stop("beta length ", nrow(eff), " != ", n_var, " in ", run)
  } else {
    return(NULL)
  }

  df <- data.frame(variant_id = seq_len(n_var), pos = pos, stringsAsFactors = FALSE)
  for (c in seq_len(n_cp)) df[[paste0("beta_cp", c)]] <- beta_mat[, c]
  for (c in seq_len(n_cp)) df[[paste0("RATE_cp", c)]] <- rate_mat[, c]
  list(df = df, n_cp = n_cp)
}

# Per-POS annotation fields (effect, impact, HGVS_P), first occurrence per POS, matching
# extract_gene_beta_RATE.sh. read.delim mangles the snpEff header (e.g. ANN[*].EFFECT ->
# ANN....EFFECT).
annot_fields <- function(annot) {
  ann <- read.delim(annot, stringsAsFactors = FALSE)
  ann <- ann[!duplicated(ann$POS), ]
  data.frame(POS = ann$POS,
             effect = ann[["ANN....EFFECT"]],
             impact = ann[["ANN....IMPACT"]],
             HGVS_P = ann[["ANN....HGVS_P"]],
             stringsAsFactors = FALSE)
}

# Top-N variant rows per gene: union of the N with the largest max-RATE and the N with
# the largest max-|beta| (NA beta -> RATE only). `df` already has the gene column.
top5_per_gene <- function(df, beta_cols, rate_cols, n = TOP_N) {
  df$max_RATE    <- apply(df[rate_cols], 1, max, na.rm = TRUE)
  absbeta        <- abs(as.matrix(df[beta_cols]))
  df$max_absbeta <- if (all(is.na(absbeta))) NA_real_
                    else apply(absbeta, 1, function(x) if (all(is.na(x))) NA_real_
                                                       else max(x, na.rm = TRUE))
  parts <- lapply(split(df, df$gene), function(g) {
    keep <- head(order(g$max_RATE, decreasing = TRUE), n)
    if (any(is.finite(g$max_absbeta)))
      keep <- union(keep, head(order(g$max_absbeta, decreasing = TRUE), n))
    g[sort(keep), ]
  })
  do.call(rbind, parts)
}

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

process_run <- function(cfg, run) {
  rdir <- run_dir(cfg$key, run)
  message("[", run, "]")
  built <- tryCatch(build_run_wide(rdir, run), error = function(e) {
    message("    ERROR: ", conditionMessage(e)); NULL })
  if (is.null(built)) { message("    skip (no RATE output yet)"); return(invisible()) }

  df <- annotate_genes(built$df, cfg$annot, cfg$goi)   # adds $gene (display-remapped)
  df <- df[!is.na(df$gene) & df$gene != "MODIFIER", ]
  if (nrow(df) == 0L) { message("    skip (no mapped genes)"); return(invisible()) }

  beta_cols <- grep("^beta_cp", names(df), value = TRUE)
  rate_cols <- grep("^RATE_cp", names(df), value = TRUE)
  sel <- top5_per_gene(df, beta_cols, rate_cols)

  fields <- annot_fields(cfg$annot)
  m <- match(sel$pos, fields$POS)
  out <- data.frame(gene = sel$gene, variant_id = sel$variant_id, POS = sel$pos,
                    effect = fields$effect[m], impact = fields$impact[m],
                    HGVS_P = fields$HGVS_P[m], stringsAsFactors = FALSE)
  out <- cbind(out, sel[beta_cols], sel[rate_cols])
  out <- out[order(out$gene, -sel$max_RATE), ]

  path <- file.path(OUT_DIR, paste0(run, "_top5_beta_RATE.csv"))
  write.csv(out, path, row.names = FALSE)
  message("    wrote ", basename(path), " (", nrow(out), " variants, ",
          length(unique(out$gene)), " genes)")
}

for (cfg in species_cfgs) {
  process_run(cfg, cfg$logistic)
  for (stub in cfg$binnings) {
    process_run(cfg, paste0(stub, "_POM"))
    process_run(cfg, paste0(stub, "_PPOM"))
  }
}

message("Done. CSVs in ", OUT_DIR)
