#!/usr/bin/env Rscript
############################################################
## write_klebsiella_jsons.R
##
## Builds inference Stan JSON datasets for the simulated
## K. pneumoniae homoplasic (single causal variant) phenotypes,
## one per effect size:
##
##   kleb_homoplasic_EF1.5_H09
##   kleb_homoplasic_EF2.5_H09
##   kleb_homoplasic_EF10_H09
##   kleb_homoplasic_EF30_H09
##
## All four share one genotype (12,410 simulated SNPs x 4,368
## samples) and one causal variant (rs_26645); only the phenotype
## changes. Output is consumed by continuous_inference.stan via
## gwas_finalruns/generate_klebsiella_run_scripts.R.
##
## Paths and parameters come from config_klebsiella.R.
## Shared builders come from utils.R.
##
## Usage:
##   sbatch gwas_datasets/submit_klebsiella.sh
############################################################

script_args <- commandArgs(trailingOnly = FALSE)
file_arg   <- grep("--file=", script_args, value = TRUE)
script_dir <- if (length(file_arg) > 0) dirname(normalizePath(sub("--file=", "", file_arg))) else "."

# config.R is sourced only for the shared OUT_* conventions and library loads
# that utils.R expects; the Klebsiella-specific paths live in config_klebsiella.R.
source(file.path(script_dir, "config.R"))
source(file.path(script_dir, "utils.R"))
source(file.path(script_dir, "config_klebsiella.R"))

pdf(NULL)


############################################################
## LOADERS
############################################################

#' Load the Klebsiella genotype matrix.
#' The file is variants x samples with an ID column and a sample-ID header,
#' i.e. the same orientation load_spn_genotype() returns, so it can be handed
#' straight to intersect_and_align(geno_in_cols = TRUE).
#' @return integer matrix (variants x samples); rownames = rs_<N>, colnames = sample IDs
load_klebsiella_genotype <- function(path) {
  message("Loading Klebsiella genotype from: ", path)
  dt <- fread(path, sep = "\t", header = TRUE)
  variant_names <- as.character(dt[[1]])
  dt[, 1 := NULL]
  mat <- as.matrix(dt)
  storage.mode(mat) <- "integer"
  rownames(mat) <- variant_names
  message("  Loaded: ", nrow(mat), " variants x ", ncol(mat), " samples")
  mat
}


#' Load Klebsiella lineages (PopPUNK) and sublineages (PopPIPE/fastbaps).
#'
#' Mirrors load_spn_lineages(): every sample gets a sublineage label, and
#' samples in PopPUNK clusters too small for PopPIPE to subcluster fall back to
#' a cluster-level singleton. Two Klebsiella-specific fixes are applied:
#'
#'   1. ID punctuation. The sublineage file writes 573_NNNNN where the genotype
#'      matrix and PopPUNK file write 573.NNNNN. SPARK_*_C1 IDs are left alone.
#'   2. Label scope. fastbaps Cluster_1 numbering restarts inside each PopPUNK
#'      cluster (only 7 distinct values across the whole file), so a raw label
#'      like "3" refers to a different subcluster in every strain. Labels are
#'      made globally unique by prefixing the parent PopPUNK cluster, giving
#'      "<PP>_<Cluster_1>"; the fallback for unsubclustered samples is "<PP>_0".
#'
#' @return list(lineages = data.frame(ID, lineage), sublineages = data.frame(ID, sublineage))
load_klebsiella_lineages <- function(poppunk_path, sublineage_path) {
  message("Loading Klebsiella PopPUNK clusters from: ", poppunk_path)
  pp <- read.delim(poppunk_path, header = FALSE, stringsAsFactors = FALSE,
                   col.names = c("ID", "lineage"))
  message("  ", nrow(pp), " samples, ", length(unique(pp$lineage)), " PopPUNK clusters")

  message("Loading Klebsiella PopPIPE subclusters from: ", sublineage_path)
  sub <- read.csv(sublineage_path, stringsAsFactors = FALSE)
  sub$ID <- sub("^573_", "573.", sub$Isolates)

  n_matched <- sum(sub$ID %in% pp$ID)
  message("  ", nrow(sub), " subclustered isolates, ", n_matched,
          " of which match a PopPUNK cluster")

  # Globally-unique sublineage label, with a cluster-level fallback.
  sub_map <- setNames(as.character(sub$Cluster_1), sub$ID)
  raw_sub <- sub_map[pp$ID]
  has_sub <- !is.na(raw_sub)
  sublabel <- ifelse(has_sub,
                     paste0(pp$lineage, "_", raw_sub),
                     paste0(pp$lineage, "_0"))

  message("  ", sum(has_sub), " samples with a fastbaps subcluster, ",
          sum(!has_sub), " on the <PP>_0 fallback")
  message("  Distinct sublineages: ", length(unique(sublabel)))

  list(
    lineages    = data.frame(ID = pp$ID, lineage    = pp$lineage,  stringsAsFactors = FALSE),
    sublineages = data.frame(ID = pp$ID, sublineage = sublabel,    stringsAsFactors = FALSE)
  )
}


#' Choose the reference (dropped) sublineage, restricted to sublineages with at
#' least `min_size` members.
#'
#' utils.R:select_reference_sublineage() takes the lowest-mean sublineage with no
#' size guard, which is safe for SPN (PopPIPE min_cluster_size = 3) but not for
#' Klebsiella, where 138 of 607 sublineages hold a single isolate. See the
#' KLEB_MIN_REF_SIZE note in config_klebsiella.R for why anchoring alpha on a
#' singleton destroys the fit.
#'
#' @param sublin_ids character vector of sublineage ID per sample
#' @param pheno_vec  numeric phenotype vector
#' @param min_size   minimum members for a sublineage to be eligible
#' @return the sublineage ID to use as reference level
select_reference_sublineage_minsize <- function(sublin_ids, pheno_vec, min_size) {
  counts   <- table(sublin_ids)
  eligible <- names(counts)[counts >= min_size]
  if (length(eligible) == 0L) {
    stop("No sublineage has at least ", min_size, " members; lower KLEB_MIN_REF_SIZE")
  }
  means <- tapply(pheno_vec, sublin_ids, mean, na.rm = TRUE)[eligible]
  ref   <- names(which.min(means))
  message("  Reference (dropped) sublineage: ", ref,
          "  (n = ", counts[[ref]], ", mean phenotype = ",
          round(min(means, na.rm = TRUE), 4), ")")
  message("    eligible sublineages (n >= ", min_size, "): ", length(eligible),
          " of ", length(counts))
  ref
}


#' Encode Klebsiella lineages/sublineages as one-hot matrices.
#'
#' Identical to utils.R:encode_lineages_spn() except that the reference
#' sublineage comes from select_reference_sublineage_minsize(). Kept here rather
#' than patched into utils.R so the spn/tb datasets keep their existing
#' behaviour byte-for-byte.
#'
#' @return list(lineage_matrix, sublineage_matrix, parent_lineage)
encode_lineages_klebsiella <- function(lineages_df, sublin_df, pheno_vec, min_ref_size) {
  ref_sublin <- select_reference_sublineage_minsize(sublin_df[[2]], pheno_vec, min_ref_size)
  sublin_fac <- relevel(factor(sublin_df[[2]]), ref = as.character(ref_sublin))

  sublin_mat <- model.matrix(~ sublin_fac - 1)
  colnames(sublin_mat) <- levels(sublin_fac)

  lin_fac     <- factor(lineages_df[[2]])
  lin_mat_tmp <- model.matrix(~ lin_fac - 1)
  colnames(lin_mat_tmp) <- levels(lin_fac)
  ref_sub_idx <- which(sublin_mat[, 1] == 1)
  ref_lin     <- names(which.max(colSums(lin_mat_tmp[ref_sub_idx, , drop = FALSE])))
  message("  Reference (dropped) lineage: ", ref_lin)

  lin_fac <- relevel(lin_fac, ref = ref_lin)
  lin_mat <- model.matrix(~ lin_fac - 1)
  colnames(lin_mat) <- levels(lin_fac)

  parent_lineage <- integer(ncol(sublin_mat))
  for (k in seq_len(ncol(sublin_mat))) {
    sub_idx <- which(sublin_mat[, k] == 1)
    parent_lineage[k] <- which.max(colSums(lin_mat[sub_idx, , drop = FALSE]))
  }

  list(lineage_matrix    = lin_mat,
       sublineage_matrix = sublin_mat,
       parent_lineage    = parent_lineage)
}


#' Load a simulated .phen file (whitespace-delimited, no header: FID IID value).
#' @return data.frame(ID, phenotype)
load_klebsiella_phenotype <- function(path) {
  ph <- read.table(path, header = FALSE, stringsAsFactors = FALSE)
  data.frame(ID = as.character(ph[[1]]), phenotype = as.numeric(ph[[3]]),
             stringsAsFactors = FALSE)
}


############################################################
## SHARED INPUTS (loaded once, reused across effect sizes)
############################################################

write_inputs_manifest(
  paths = c(
    kleb_genotype    = KLEB_GENOTYPE_PATH,
    kleb_annotations = KLEB_ANNOTATIONS_PATH,
    kleb_poppunk     = KLEB_POPPUNK_PATH,
    kleb_sublineages = KLEB_SUBLINEAGE_PATH,
    kleb_phenotypes  = KLEB_PHENO_DIR
  ),
  out_base = KLEB_OUT_BASE
)

kleb_geno <- load_klebsiella_genotype(KLEB_GENOTYPE_PATH)

# MAF filter (see KLEB_MIN_MAF in config_klebsiella.R -- required for the model
# to fit at all, not a cosmetic threshold).
kleb_af  <- rowMeans(kleb_geno)
kleb_maf <- pmin(kleb_af, 1 - kleb_af)
keep_maf <- kleb_maf > KLEB_MIN_MAF
message("MAF > ", KLEB_MIN_MAF, " keeps ", sum(keep_maf), " of ", length(keep_maf),
        " variants (dropped ", sum(!keep_maf), ")")
stopifnot("causal variant must survive the MAF filter" =
            isTRUE(keep_maf[["rs_26645"]]))
kleb_geno <- kleb_geno[keep_maf, , drop = FALSE]

kleb_lin  <- load_klebsiella_lineages(KLEB_POPPUNK_PATH, KLEB_SUBLINEAGE_PATH)

# Positions are the integer parsed out of each rs_<N> ID. Verified unique across
# all 12,410 variants, which the annotation join in pipeline.R relies on
# (it does match(variant_positions, ann$POS)).
variant_names     <- rownames(kleb_geno)
variant_positions <- as.integer(sub("^rs_", "", variant_names))
stopifnot(!anyNA(variant_positions),
          !anyDuplicated(variant_positions))


############################################################
## ONE DATASET PER EFFECT SIZE
############################################################

# Effect-size label -> directory-safe suffix ("10.0" -> "10")
ef_tag <- function(ef) sub("\\.0$", "", ef)

summary_rows <- list()

for (ef in KLEB_EFFECT_SIZES) {

  dataset_name <- paste0("kleb_homoplasic_EF", ef_tag(ef), "_H",
                         sub("\\.", "", KLEB_HERITABILITY))
  message("\n=== ", dataset_name, " ===")
  out_dir <- file.path(KLEB_OUT_INFER, dataset_name)

  pheno_path <- file.path(
    KLEB_PHENO_DIR,
    sprintf("simul_homoplasic_EF_%s_herit_%s_traitprev_0.1.phen",
            ef, KLEB_HERITABILITY)
  )
  if (!file.exists(pheno_path)) stop("Missing phenotype file: ", pheno_path)
  pheno_df <- load_klebsiella_phenotype(pheno_path)

  # Causal variant for this architecture/effect size
  causal_path <- file.path(KLEB_CAUSAL_DIR, paste0("homoplasic_EffSize_", ef))
  causal <- read.table(causal_path, header = FALSE, stringsAsFactors = FALSE)
  causal_variant <- as.character(causal[[1]][1])
  causal_effect  <- as.numeric(causal[[2]][1])
  message("  Causal variant: ", causal_variant, " (simulated effect ", causal_effect, ")")

  aligned <- intersect_and_align(
    pheno_df     = pheno_df,
    geno         = kleb_geno,
    lineages_df  = kleb_lin$lineages,
    sublin_df    = kleb_lin$sublineages,
    id_col       = "ID",
    geno_in_cols = TRUE
  )

  pheno_vec <- aligned$pheno$phenotype

  # Phenotype transform (see KLEB_PHENO_TRANSFORM in config_klebsiella.R).
  pheno_mean <- mean(pheno_vec)
  pheno_sd   <- sd(pheno_vec)
  if (identical(KLEB_PHENO_TRANSFORM, "zscore")) {
    pheno_fit <- (pheno_vec - pheno_mean) / pheno_sd
    message("  Phenotype z-scored (mean ", round(pheno_mean, 4),
            ", sd ", round(pheno_sd, 4), "); recover raw beta as beta_z * sd")
  } else {
    pheno_fit <- pheno_vec
    message("  Phenotype used raw: sd = ", round(pheno_sd, 3),
            ", implied residual sd at h2=", KLEB_HERITABILITY, " = ",
            round(pheno_sd * sqrt(1 - as.numeric(KLEB_HERITABILITY)), 3))
  }

  enc <- encode_lineages_klebsiella(
    lineages_df  = aligned$lineages,
    sublin_df    = aligned$sublineages,
    pheno_vec    = pheno_fit,
    min_ref_size = KLEB_MIN_REF_SIZE
  )

  stan_list <- build_stan_inference(
    pheno      = pheno_fit,
    geno_mat   = aligned$geno_mat,
    lin_mat    = enc$lineage_matrix,
    sublin_mat = enc$sublineage_matrix,
    parent_lin = enc$parent_lineage
  )

  message("  N = ", stan_list$N, "  V = ", stan_list$V,
          "  L = ", stan_list$L, "  S = ", stan_list$S)

  write_dataset(
    stan_list         = stan_list,
    sample_ids        = aligned$sample_ids,
    variant_names     = colnames(aligned$geno_mat),
    parent_lin        = enc$parent_lineage,
    outdir            = out_dir,
    dataset_name      = dataset_name,
    variant_positions = as.integer(sub("^rs_", "", colnames(aligned$geno_mat)))
  )

  # Record the truth alongside the dataset so the figure step never has to
  # re-derive it. causal_row is the 1-based column index of the causal variant
  # in variant_matrix, i.e. beta_variant[causal_row].
  causal_row <- match(causal_variant, colnames(aligned$geno_mat))
  stopifnot(!is.na(causal_row))
  carriers <- aligned$geno_mat[, causal_row] == 1
  observed_diff <- mean(pheno_vec[carriers]) - mean(pheno_vec[!carriers])

  truth <- data.frame(
    dataset          = dataset_name,
    effect_size      = as.numeric(ef),
    heritability     = as.numeric(KLEB_HERITABILITY),
    causal_variant   = causal_variant,
    causal_position  = as.integer(sub("^rs_", "", causal_variant)),
    causal_row       = causal_row,
    simulated_effect = causal_effect,
    n_carriers       = sum(carriers),
    allele_freq      = mean(carriers),
    observed_diff    = observed_diff,
    pheno_mean       = pheno_mean,
    pheno_sd         = pheno_sd,
    transform        = KLEB_PHENO_TRANSFORM,
    stringsAsFactors = FALSE
  )
  write.csv(truth, file.path(out_dir, paste0(dataset_name, "_truth.csv")),
            row.names = FALSE)
  message("  Causal variant is column ", causal_row, " (beta_variant[", causal_row, "]); ",
          "observed carrier-vs-non difference = ", round(observed_diff, 3))

  summary_rows[[dataset_name]] <- truth
}

summary_df <- do.call(rbind, summary_rows)
write.csv(summary_df,
          file.path(KLEB_OUT_INFER, "klebsiella_homoplasic_truth_summary.csv"),
          row.names = FALSE)

message("\nWrote ", nrow(summary_df), " Klebsiella inference datasets to ", KLEB_OUT_INFER)
print(summary_df[, c("dataset", "effect_size", "causal_row",
                     "observed_diff", "pheno_sd")])
