#!/usr/bin/env Rscript
# Repair all TB *_variant_index.csv files that were truncated to a single row
# by the previous fread(...)[[1]] reader in load_tb_genotype.
#
# Reads the source variant index TSV (one row, tab-separated), re-derives
# positions the same way write_dataset() does, and overwrites each TB
# dataset's variant_index.csv. The JSONs are untouched.

TB_VARIANT_IDX_PATH <- "/nfs/research/jlees/jacqueline/gwas_data/tuberculosis/cryptic_regeno_snpeff/subsampled_presabs/variant_index_presence_absence_RIF_K8.tsv"

INFER_DIR <- "/nfs/research/jlees/jacqueline/thesis_results/gwas_datasets/inference"
PRED_DIR  <- "/nfs/research/jlees/jacqueline/thesis_results/gwas_datasets/prediction"

TB_DATASETS_INFER <- c(
  "07_tb_rifampicin_binary",
  "08_tb_rifampicin_MIC",
  "09_tb_rifampicin_continuous",
  "14_tb_rifampicin_MIC_coarse_dilutions",
  "15_tb_rifampicin_MIC_large_minbin"
)
TB_DATASETS_PRED <- c(
  "07_tb_rifampicin_binary",       "07_tb_rifampicin_binary_loso",
  "08_tb_rifampicin_MIC",          "08_tb_rifampicin_MIC_loso",
  "09_tb_rifampicin_continuous",   "09_tb_rifampicin_continuous_loso",
  "14_tb_rifampicin_MIC_coarse_dilutions",
  "14_tb_rifampicin_MIC_coarse_dilutions_loso",
  "15_tb_rifampicin_MIC_large_minbin",
  "15_tb_rifampicin_MIC_large_minbin_loso"
)

message("Reading source variant index: ", TB_VARIANT_IDX_PATH)
vnames <- scan(TB_VARIANT_IDX_PATH, what = character(), sep = "\t",
               quiet = TRUE, na.strings = "")
message("  ", length(vnames), " variants")

# Same position-parsing logic as write_dataset() in utils.R
vnames_clean <- sub("^Chromosome_", "", vnames)
pos_str      <- sub("_[^_]+$", "", vnames_clean)
pos_str      <- sub("_[^_]+$", "", pos_str)
variant_pos  <- suppressWarnings(as.numeric(pos_str))
n_na <- sum(is.na(variant_pos))
if (n_na > 0) {
  message("  WARNING: ", n_na, " variant names did not parse to a position")
}

vi <- data.frame(variant_name = vnames, position = variant_pos,
                 stringsAsFactors = FALSE)

repair_one <- function(base_dir, dataset_name) {
  out_path <- file.path(base_dir, dataset_name,
                        paste0(dataset_name, "_variant_index.csv"))
  if (!file.exists(out_path)) {
    message("  skip (no existing file): ", out_path)
    return(invisible(NULL))
  }
  old_n <- length(readLines(out_path)) - 1L
  write.csv(vi, out_path, row.names = FALSE)
  new_n <- length(readLines(out_path)) - 1L
  message(sprintf("  wrote %s (%d -> %d rows)", out_path, old_n, new_n))
}

message("\n== Inference datasets ==")
for (d in TB_DATASETS_INFER) repair_one(INFER_DIR, d)

message("\n== Prediction datasets ==")
for (d in TB_DATASETS_PRED) repair_one(PRED_DIR, d)

message("\nDone.")
