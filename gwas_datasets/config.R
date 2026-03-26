############################################################
## config.R
## Single source of truth for all input paths and parameters
## used in write_inference_jsons.R and write_prediction_jsons.R
############################################################

# === INPUT DATA PATHS ===

# SPN genotype (MAF > 0.05, no modifier annotations, multiallelic; last modified Jul 14 2025)
SPN_GENOTYPE_PATH <- "/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/genotype/genotype_maf05_NOMODIFIERS_multiallelic.tsv"

# SPN lineages — PopPIPE/PopPunk sublineages, min_cluster_size = 3
SPN_LINEAGES_PATH <- "/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/lineages/sublineages_from_poppipe/min_cluster_3/all_clusters.txt"

# SPN phenotype sources
SPN_METADATA_PATH    <- "/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/phenotype/sparc_metadata.csv"

# TB genotype — presence/absence matrix filtered to samples with RIF MIC data (Jan 20 2025)
TB_PRESABS_PATH     <- "/nfs/research/jlees/jacqueline/gwas_data/tuberculosis/cryptic_regeno_snpeff/subsampled_presabs/presence_absence_final_RIF_K8.tsv"
TB_VARIANT_IDX_PATH  <- "/nfs/research/jlees/jacqueline/gwas_data/tuberculosis/cryptic_regeno_snpeff/subsampled_presabs/variant_index_presence_absence_RIF_K8.tsv"
TB_SAMPLE_INDEX_PATH <- "/nfs/research/jlees/jacqueline/gwas_data/tuberculosis/cryptic_regeno_snpeff/sample_index_RIF_K8.txt"

# TB phenotype — CRyPTIC reuse table (Sep 17 2024; most recent version)
# Relevant columns: col 1 = ENA_RUN, col 27 = RIF_MIC, col 44 = REGENOTYPED_VCF
TB_PHENOTYPE_PATH <- "/nfs/research/jlees/jacqueline/gwas_data/tuberculosis/CRyPTIC_reuse_table_20240917.csv"

# TB lineages — fastlin barcode assignments
TB_LINEAGES_PATH <- "/nfs/research/jlees/jacqueline/gwas_data/code/gwas_data/tuberculosis/01_lineages/output_fastlin.txt"


# === PHENOTYPE PARAMETERS ===

# Standard microbiological doubling-dilution breakpoints (needed by load_breakpoints below)
MIC_STANDARD_DILUTIONS <- c(0, 0.008, 0.016, 0.03, 0.06, 0.12, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64)

# Load clinical breakpoints from CSV and derive S_MAX / R_MIN from standard dilutions.
#   type = "r_min"     : breakpoint is R threshold; S_MAX = next lower dilution
#   type = "s_max"     : breakpoint is S threshold; R_MIN = next higher dilution
#   type = "threshold" : simple cutoff (S if <=); no intermediate zone
load_breakpoints <- function(csv_path, dilutions = MIC_STANDARD_DILUTIONS) {
  bp <- read.csv(csv_path, stringsAsFactors = FALSE)
  stopifnot(all(c("species_drug", "breakpoint", "type", "source") %in% names(bp)))

  result <- list()
  for (i in seq_len(nrow(bp))) {
    row <- bp[i, ]
    val <- row$breakpoint
    entry <- list(breakpoint = val, source = row$source)

    if (row$type == "r_min") {
      entry$r_min <- val
      lower <- dilutions[dilutions < val]
      if (length(lower) == 0) stop("No standard dilution below ", val, " for ", row$species_drug)
      entry$s_max <- max(lower)
    } else if (row$type == "s_max") {
      entry$s_max <- val
      upper <- dilutions[dilutions > val]
      if (length(upper) == 0) stop("No standard dilution above ", val, " for ", row$species_drug)
      entry$r_min <- min(upper)
    } else if (row$type == "threshold") {
      entry$s_max <- val
      entry$r_min <- NULL
    } else {
      stop("Unknown breakpoint type: ", row$type)
    }
    result[[row$species_drug]] <- entry
  }
  result
}

BREAKPOINTS <- load_breakpoints(file.path(script_dir, "EUCAST_breakpoints.csv"))

SPN_PEN_BINARY_S_MAX    <- BREAKPOINTS$spn_benzylpenicillin$s_max
SPN_PEN_BINARY_R_MIN    <- BREAKPOINTS$spn_benzylpenicillin$r_min
SPN_TMP_BINARY_S_MAX    <- BREAKPOINTS$spn_trimethoprim$s_max
SPN_TMP_BINARY_R_MIN    <- BREAKPOINTS$spn_trimethoprim$r_min
TB_RIF_BINARY_THRESHOLD <- BREAKPOINTS$tb_rifampicin$s_max

# MIC cleaning replacement maps (raw string -> numeric value)
# Note: for TB RIF, "<=0.03" -> 0.03 (v1 scripts had a bug mapping this to 8)
SPN_PEN_MIC_REPLACEMENTS <- c(
  ">4.00" = 4, ">4" = 4, "<=0.03" = 0.03, "<0.03" = 0.03,
  "<0.016" = 0.016, "<.016" = 0.016, "< 0.01" = 0.01
)

SPN_TMP_MIC_REPLACEMENTS <- c(
  ">2.00" = 2, "<=0.25" = 0.25, ">32" = 32,
  "4/76" = 4, ">4/76" = 4, "0.25/4" = 0.25,
  "0.5/9." = 0.5, "0.12/2" = 0.12, "1/19" = 1,
  "<0.06/" = 0.06, "2/38" = 2
)

TB_RIF_MIC_REPLACEMENTS <- c(
  "<=0.06" = 0.06, ">2" = 2, ">4" = 4, ">8" = 8, "<=0.03" = 0.03
)

# Minimum fraction of samples required per bin after auto-binning
MIC_MIN_BIN_FRAC        <- 0.05   # 5%  — standard MIC datasets
MIC_MIN_BIN_FRAC_COARSE <- 0.10   # 10% — coarser-binned MIC datasets

# 4-fold (every-other-step) dilution grid for coarser MIC datasets
MIC_COARSE_DILUTIONS <- MIC_STANDARD_DILUTIONS[c(TRUE, FALSE)]


# === OUTPUT PATHS ===

OUT_BASE  <- "/nfs/research/jlees/jacqueline/thesis_results/gwas_datasets"
OUT_INFER <- file.path(OUT_BASE, "inference")
OUT_PRED  <- file.path(OUT_BASE, "prediction")
OUT_HIST  <- file.path(OUT_BASE, "MIC_bin_histograms")

# Prediction train/test split parameters
PRED_SEED       <- 12345L
PRED_TRAIN_PROP <- 0.8
