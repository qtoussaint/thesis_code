############################################################
## write_inference_jsons.R
##
## Creates 15 inference Stan JSON datasets:
##   01 SPN penicillin    binary
##   02 SPN penicillin    MIC (ordinal, dilution-merge >= 30/bin)
##   03 SPN penicillin    continuous (log2)
##   04 SPN trimethoprim  binary
##   05 SPN trimethoprim  MIC (ordinal, dilution-merge >= 30/bin)
##   06 SPN trimethoprim  continuous (log2)
##   07 TB  rifampicin    binary
##   08 TB  rifampicin    MIC (ordinal, dilution-merge >= 30/bin)
##   09 TB  rifampicin    continuous (log2)
##   10 SPN penicillin    MIC (ordinal, equal-frequency bins)
##   11 SPN penicillin    MIC (ordinal, peak-valley bins)
##   12 SPN trimethoprim  MIC (ordinal, equal-frequency bins)
##   13 SPN trimethoprim  MIC (ordinal, peak-valley bins)
##   14 TB  rifampicin    MIC (ordinal, equal-frequency bins)
##   15 TB  rifampicin    MIC (ordinal, peak-valley bins)
##
## All input paths and parameters come from config.R.
## Shared functions come from utils.R.
############################################################

script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("--file=", script_args, value = TRUE)
script_dir <- if (length(file_arg) > 0) dirname(normalizePath(sub("--file=", "", file_arg))) else "."
source(file.path(script_dir, "config.R"))
source(file.path(script_dir, "utils.R"))

# Write manifest of all input files used
write_inputs_manifest(
  paths = c(
    spn_genotype     = SPN_GENOTYPE_PATH,
    spn_lineages     = SPN_LINEAGES_PATH,
    spn_metadata     = SPN_METADATA_PATH,
    spn_pen_binary   = SPN_PEN_BINARY_PATH,
    tb_presabs       = TB_PRESABS_PATH,
    tb_variant_index = TB_VARIANT_IDX_PATH,
    tb_phenotype     = TB_PHENOTYPE_PATH,
    tb_lineages      = TB_LINEAGES_PATH
  ),
  out_base = OUT_BASE
)


############################################################
## SPN DATA LOADING (shared across datasets 01-06)
############################################################

message("\n=== Loading SPN data ===")
spn_geno <- load_spn_genotype(SPN_GENOTYPE_PATH)
spn_lin  <- load_spn_lineages(SPN_LINEAGES_PATH)

# Load SPN MIC metadata (used for penicillin and trimethoprim MIC/continuous)
spn_meta <- read.csv(SPN_METADATA_PATH, na.strings = "")

# Load raw MIC columns and derive phenotype data frames
pen_raw <- data.frame(
  ID  = spn_meta$Lane.name,
  MIC = as.character(spn_meta$Benzylpenicillin.MIC..ug.mL.)
)
pen_raw <- pen_raw[complete.cases(pen_raw), ]
pen_raw$MIC_num <- clean_mic_values(pen_raw$MIC, SPN_PEN_MIC_REPLACEMENTS)
pen_raw <- pen_raw[!is.na(pen_raw$MIC_num), ]

tmp_raw <- data.frame(
  ID  = spn_meta$Lane.name,
  MIC = as.character(spn_meta$Trimethoprim.MIC..ug.mL.)
)
tmp_raw <- tmp_raw[complete.cases(tmp_raw), ]
tmp_raw$MIC_num <- clean_mic_values(tmp_raw$MIC, SPN_TMP_MIC_REPLACEMENTS)
tmp_raw <- tmp_raw[!is.na(tmp_raw$MIC_num), ]

# Load pre-computed EUCAST binary for penicillin
spn_pen_bin_raw <- read.csv(SPN_PEN_BINARY_PATH, na.strings = "")


############################################################
## 01: SPN PENICILLIN BINARY
############################################################

message("\n=== 01: SPN penicillin binary ===")
dataset_name <- "01_spn_penicillin_binary"
out_dir <- file.path(OUT_INFER, dataset_name)

pheno_df <- spn_pen_bin_raw[complete.cases(spn_pen_bin_raw), ]

aligned <- intersect_and_align(
  pheno_df    = pheno_df,
  geno        = spn_geno,
  lineages_df = spn_lin$lineages,
  sublin_df   = spn_lin$sublineages,
  id_col      = "ID",
  geno_in_cols = TRUE
)

enc <- encode_lineages_spn(
  lineages_df = aligned$lineages,
  sublin_df   = aligned$sublineages,
  pheno_vec   = aligned$pheno$resistance
)

stan_list <- build_stan_inference(
  pheno      = aligned$pheno$resistance,
  geno_mat   = aligned$geno_mat,
  lin_mat    = enc$lineage_matrix,
  sublin_mat = enc$sublineage_matrix,
  parent_lin = enc$parent_lineage
)

write_dataset(
  stan_list    = stan_list,
  sample_ids   = aligned$sample_ids,
  variant_names = rownames(spn_geno),
  parent_lin   = enc$parent_lineage,
  outdir       = out_dir,
  dataset_name = dataset_name
)


############################################################
## 02: SPN PENICILLIN MIC (ordinal, auto-binned)
############################################################

message("\n=== 02: SPN penicillin MIC (ordinal) ===")
dataset_name <- "02_spn_penicillin_MIC"
out_dir <- file.path(OUT_INFER, dataset_name)

aligned <- intersect_and_align(
  pheno_df    = pen_raw,
  geno        = spn_geno,
  lineages_df = spn_lin$lineages,
  sublin_df   = spn_lin$sublineages,
  id_col      = "ID",
  geno_in_cols = TRUE
)

binning <- bin_mic_auto(
  mic_numeric   = aligned$pheno$MIC_num,
  min_bin_size  = MIC_MIN_BIN_SIZE,
  dilutions     = MIC_STANDARD_DILUTIONS,
  hist_path     = file.path(OUT_HIST, paste0(dataset_name, "_bins.png")),
  dataset_label = "SPN Penicillin"
)

enc <- encode_lineages_spn(
  lineages_df = aligned$lineages,
  sublin_df   = aligned$sublineages,
  pheno_vec   = binning$bins
)

stan_list <- build_stan_inference(
  pheno      = binning$bins,
  geno_mat   = aligned$geno_mat,
  lin_mat    = enc$lineage_matrix,
  sublin_mat = enc$sublineage_matrix,
  parent_lin = enc$parent_lineage,
  K          = binning$K,
  mic_bkpts  = binning$mic_breakpoints
)

write_dataset(
  stan_list    = stan_list,
  sample_ids   = aligned$sample_ids,
  variant_names = rownames(spn_geno),
  parent_lin   = enc$parent_lineage,
  outdir       = out_dir,
  dataset_name = dataset_name
)


############################################################
## 03: SPN PENICILLIN CONTINUOUS (log2 MIC)
############################################################

message("\n=== 03: SPN penicillin continuous ===")
dataset_name <- "03_spn_penicillin_continuous"
out_dir <- file.path(OUT_INFER, dataset_name)

aligned <- intersect_and_align(
  pheno_df    = pen_raw,
  geno        = spn_geno,
  lineages_df = spn_lin$lineages,
  sublin_df   = spn_lin$sublineages,
  id_col      = "ID",
  geno_in_cols = TRUE
)

log2_mic <- log2(aligned$pheno$MIC_num)

enc <- encode_lineages_spn(
  lineages_df = aligned$lineages,
  sublin_df   = aligned$sublineages,
  pheno_vec   = log2_mic
)

stan_list <- build_stan_inference(
  pheno      = log2_mic,
  geno_mat   = aligned$geno_mat,
  lin_mat    = enc$lineage_matrix,
  sublin_mat = enc$sublineage_matrix,
  parent_lin = enc$parent_lineage
)

write_dataset(
  stan_list    = stan_list,
  sample_ids   = aligned$sample_ids,
  variant_names = rownames(spn_geno),
  parent_lin   = enc$parent_lineage,
  outdir       = out_dir,
  dataset_name = dataset_name
)


############################################################
## 04: SPN TRIMETHOPRIM BINARY
############################################################

message("\n=== 04: SPN trimethoprim binary ===")
dataset_name <- "04_spn_trimethoprim_binary"
out_dir <- file.path(OUT_INFER, dataset_name)

tmp_bin <- tmp_raw
tmp_bin$binary <- mic_to_binary(
  mic_numeric = tmp_bin$MIC_num,
  s_max       = SPN_TMP_BINARY_S_MAX,
  r_min       = SPN_TMP_BINARY_R_MIN
)
tmp_bin <- tmp_bin[!is.na(tmp_bin$binary), ]
tmp_bin$pheno <- tmp_bin$binary

aligned <- intersect_and_align(
  pheno_df    = tmp_bin,
  geno        = spn_geno,
  lineages_df = spn_lin$lineages,
  sublin_df   = spn_lin$sublineages,
  id_col      = "ID",
  geno_in_cols = TRUE
)

enc <- encode_lineages_spn(
  lineages_df = aligned$lineages,
  sublin_df   = aligned$sublineages,
  pheno_vec   = aligned$pheno$pheno
)

stan_list <- build_stan_inference(
  pheno      = aligned$pheno$pheno,
  geno_mat   = aligned$geno_mat,
  lin_mat    = enc$lineage_matrix,
  sublin_mat = enc$sublineage_matrix,
  parent_lin = enc$parent_lineage
)

write_dataset(
  stan_list    = stan_list,
  sample_ids   = aligned$sample_ids,
  variant_names = rownames(spn_geno),
  parent_lin   = enc$parent_lineage,
  outdir       = out_dir,
  dataset_name = dataset_name
)


############################################################
## 05: SPN TRIMETHOPRIM MIC (ordinal, auto-binned)
############################################################

message("\n=== 05: SPN trimethoprim MIC (ordinal) ===")
dataset_name <- "05_spn_trimethoprim_MIC"
out_dir <- file.path(OUT_INFER, dataset_name)

aligned <- intersect_and_align(
  pheno_df    = tmp_raw,
  geno        = spn_geno,
  lineages_df = spn_lin$lineages,
  sublin_df   = spn_lin$sublineages,
  id_col      = "ID",
  geno_in_cols = TRUE
)

binning <- bin_mic_auto(
  mic_numeric   = aligned$pheno$MIC_num,
  min_bin_size  = MIC_MIN_BIN_SIZE,
  dilutions     = MIC_STANDARD_DILUTIONS,
  hist_path     = file.path(OUT_HIST, paste0(dataset_name, "_bins.png")),
  dataset_label = "SPN Trimethoprim"
)

enc <- encode_lineages_spn(
  lineages_df = aligned$lineages,
  sublin_df   = aligned$sublineages,
  pheno_vec   = binning$bins
)

stan_list <- build_stan_inference(
  pheno      = binning$bins,
  geno_mat   = aligned$geno_mat,
  lin_mat    = enc$lineage_matrix,
  sublin_mat = enc$sublineage_matrix,
  parent_lin = enc$parent_lineage,
  K          = binning$K,
  mic_bkpts  = binning$mic_breakpoints
)

write_dataset(
  stan_list    = stan_list,
  sample_ids   = aligned$sample_ids,
  variant_names = rownames(spn_geno),
  parent_lin   = enc$parent_lineage,
  outdir       = out_dir,
  dataset_name = dataset_name
)


############################################################
## 06: SPN TRIMETHOPRIM CONTINUOUS (log2 MIC)
############################################################

message("\n=== 06: SPN trimethoprim continuous ===")
dataset_name <- "06_spn_trimethoprim_continuous"
out_dir <- file.path(OUT_INFER, dataset_name)

aligned <- intersect_and_align(
  pheno_df    = tmp_raw,
  geno        = spn_geno,
  lineages_df = spn_lin$lineages,
  sublin_df   = spn_lin$sublineages,
  id_col      = "ID",
  geno_in_cols = TRUE
)

log2_mic <- log2(aligned$pheno$MIC_num)

enc <- encode_lineages_spn(
  lineages_df = aligned$lineages,
  sublin_df   = aligned$sublineages,
  pheno_vec   = log2_mic
)

stan_list <- build_stan_inference(
  pheno      = log2_mic,
  geno_mat   = aligned$geno_mat,
  lin_mat    = enc$lineage_matrix,
  sublin_mat = enc$sublineage_matrix,
  parent_lin = enc$parent_lineage
)

write_dataset(
  stan_list    = stan_list,
  sample_ids   = aligned$sample_ids,
  variant_names = rownames(spn_geno),
  parent_lin   = enc$parent_lineage,
  outdir       = out_dir,
  dataset_name = dataset_name
)


############################################################
## TB DATA LOADING (shared across datasets 07-09)
############################################################

message("\n=== Loading TB data ===")
tb_geno_raw <- load_tb_genotype(TB_PRESABS_PATH, TB_VARIANT_IDX_PATH)
tb_geno     <- tb_geno_raw$genotype
tb_var_names <- tb_geno_raw$variant_names
tb_lin_raw  <- load_tb_lineages(TB_LINEAGES_PATH)

# Load and clean TB phenotype table (cols 1=ENA_RUN, 27=RIF_MIC, 44=REGENOTYPED_VCF)
tb_pheno_raw <- read.csv(TB_PHENOTYPE_PATH)
tb_pheno_raw <- tb_pheno_raw[, c(1, 27, 44)]
colnames(tb_pheno_raw) <- c("ENA_RUN", "RIF_MIC", "REGENOTYPED_VCF")

# Remove NAs
tb_pheno_raw$RIF_MIC[tb_pheno_raw$RIF_MIC == ""] <- NA
tb_pheno_raw <- na.omit(tb_pheno_raw)

# Clean MIC strings
tb_pheno_raw$RIF_MIC_num <- clean_mic_values(
  tb_pheno_raw$RIF_MIC, TB_RIF_MIC_REPLACEMENTS
)
tb_pheno_raw <- tb_pheno_raw[!is.na(tb_pheno_raw$RIF_MIC_num), ]

message("TB phenotype samples with RIF MIC: ", nrow(tb_pheno_raw))

# Align phenotype and lineages (lineages is the "sublin", sample IDs for genotype
# come from the order in the presence/absence matrix which matches phenotype ENA_RUN)
# The presence_absence matrix rows correspond to the sample order in the phenotype
# table (see creation script). We use ENA_RUN to link.
# Genotype row indices are in the same order as phenotype rows after the RIF_K8 filter.

# Align lineages to phenotype
sample_in_lin <- match(tb_pheno_raw$ENA_RUN, tb_lin_raw$X.sample)
tb_pheno_with_lin <- tb_pheno_raw[!is.na(sample_in_lin), ]
tb_lin_aligned    <- tb_lin_raw[na.omit(sample_in_lin), ]

# Now align genotype rows to this order
# Genotype sample order = tb_pheno_raw$ENA_RUN (rows of presence_absence_final)
geno_row_idx <- match(tb_pheno_with_lin$ENA_RUN, tb_pheno_raw$ENA_RUN)
geno_row_idx <- geno_row_idx[!is.na(geno_row_idx)]
tb_geno_mat  <- as.matrix(tb_geno[geno_row_idx, ])
storage.mode(tb_geno_mat) <- "integer"

# Final aligned sample IDs and MIC values
tb_sample_ids <- tb_pheno_with_lin$ENA_RUN
tb_mic_num    <- tb_pheno_with_lin$RIF_MIC_num

message("TB aligned samples: ", length(tb_sample_ids))
message("TB unique MIC values (cleaned): ", paste(sort(unique(tb_mic_num)), collapse = ", "))


############################################################
## 07: TB RIFAMPICIN BINARY
############################################################

message("\n=== 07: TB rifampicin binary ===")
dataset_name <- "07_tb_rifampicin_binary"
out_dir <- file.path(OUT_INFER, dataset_name)

tb_binary <- mic_to_binary(tb_mic_num, s_max = TB_RIF_BINARY_THRESHOLD)
keep_idx  <- !is.na(tb_binary)
message("  Binary class counts: 0=", sum(tb_binary[keep_idx]==0),
        "  1=", sum(tb_binary[keep_idx]==1))

enc <- encode_lineages_tb(
  lineages_df = tb_lin_aligned[keep_idx, ],
  pheno_vec   = tb_binary[keep_idx]
)

stan_list <- build_stan_inference(
  pheno      = tb_binary[keep_idx],
  geno_mat   = tb_geno_mat[keep_idx, , drop = FALSE],
  lin_mat    = enc$lineage_matrix,
  sublin_mat = enc$sublineage_matrix,
  parent_lin = enc$parent_lineage
)

write_dataset(
  stan_list    = stan_list,
  sample_ids   = tb_sample_ids[keep_idx],
  variant_names = tb_var_names,
  parent_lin   = enc$parent_lineage,
  outdir       = out_dir,
  dataset_name = dataset_name
)


############################################################
## 08: TB RIFAMPICIN MIC (ordinal, auto-binned)
############################################################

message("\n=== 08: TB rifampicin MIC (ordinal) ===")
dataset_name <- "08_tb_rifampicin_MIC"
out_dir <- file.path(OUT_INFER, dataset_name)

binning <- bin_mic_auto(
  mic_numeric   = tb_mic_num,
  min_bin_size  = MIC_MIN_BIN_SIZE,
  dilutions     = MIC_STANDARD_DILUTIONS,
  hist_path     = file.path(OUT_HIST, paste0(dataset_name, "_bins.png")),
  dataset_label = "TB Rifampicin"
)

enc <- encode_lineages_tb(
  lineages_df = tb_lin_aligned,
  pheno_vec   = binning$bins
)

stan_list <- build_stan_inference(
  pheno      = binning$bins,
  geno_mat   = tb_geno_mat,
  lin_mat    = enc$lineage_matrix,
  sublin_mat = enc$sublineage_matrix,
  parent_lin = enc$parent_lineage,
  K          = binning$K,
  mic_bkpts  = binning$mic_breakpoints
)

write_dataset(
  stan_list    = stan_list,
  sample_ids   = tb_sample_ids,
  variant_names = tb_var_names,
  parent_lin   = enc$parent_lineage,
  outdir       = out_dir,
  dataset_name = dataset_name
)


############################################################
## 09: TB RIFAMPICIN CONTINUOUS (log2 MIC)
############################################################

message("\n=== 09: TB rifampicin continuous ===")
dataset_name <- "09_tb_rifampicin_continuous"
out_dir <- file.path(OUT_INFER, dataset_name)

log2_mic <- log2(tb_mic_num)

enc <- encode_lineages_tb(
  lineages_df = tb_lin_aligned,
  pheno_vec   = log2_mic
)

stan_list <- build_stan_inference(
  pheno      = log2_mic,
  geno_mat   = tb_geno_mat,
  lin_mat    = enc$lineage_matrix,
  sublin_mat = enc$sublineage_matrix,
  parent_lin = enc$parent_lineage
)

write_dataset(
  stan_list    = stan_list,
  sample_ids   = tb_sample_ids,
  variant_names = tb_var_names,
  parent_lin   = enc$parent_lineage,
  outdir       = out_dir,
  dataset_name = dataset_name
)

############################################################
## 10: SPN PENICILLIN MIC (ordinal, equal-frequency bins)
############################################################

message("\n=== 10: SPN penicillin MIC (equal-freq) ===")
dataset_name <- "10_spn_penicillin_MIC_equalfreq"
out_dir <- file.path(OUT_INFER, dataset_name)

aligned <- intersect_and_align(
  pheno_df    = pen_raw,
  geno        = spn_geno,
  lineages_df = spn_lin$lineages,
  sublin_df   = spn_lin$sublineages,
  id_col      = "ID",
  geno_in_cols = TRUE
)

binning <- bin_mic_equalfreq(
  mic_numeric   = aligned$pheno$MIC_num,
  min_bin_size  = MIC_MIN_BIN_SIZE,
  dilutions     = MIC_STANDARD_DILUTIONS,
  hist_path     = file.path(OUT_HIST, paste0(dataset_name, "_bins.png")),
  dataset_label = "SPN Penicillin"
)

enc <- encode_lineages_spn(
  lineages_df = aligned$lineages,
  sublin_df   = aligned$sublineages,
  pheno_vec   = binning$bins
)

stan_list <- build_stan_inference(
  pheno      = binning$bins,
  geno_mat   = aligned$geno_mat,
  lin_mat    = enc$lineage_matrix,
  sublin_mat = enc$sublineage_matrix,
  parent_lin = enc$parent_lineage,
  K          = binning$K,
  mic_bkpts  = binning$mic_breakpoints
)

write_dataset(
  stan_list    = stan_list,
  sample_ids   = aligned$sample_ids,
  variant_names = rownames(spn_geno),
  parent_lin   = enc$parent_lineage,
  outdir       = out_dir,
  dataset_name = dataset_name
)


############################################################
## 11: SPN PENICILLIN MIC (ordinal, peak-valley bins)
############################################################

message("\n=== 11: SPN penicillin MIC (peaks) ===")
dataset_name <- "11_spn_penicillin_MIC_peaks"
out_dir <- file.path(OUT_INFER, dataset_name)

aligned <- intersect_and_align(
  pheno_df    = pen_raw,
  geno        = spn_geno,
  lineages_df = spn_lin$lineages,
  sublin_df   = spn_lin$sublineages,
  id_col      = "ID",
  geno_in_cols = TRUE
)

binning <- bin_mic_peaks(
  mic_numeric   = aligned$pheno$MIC_num,
  min_bin_size  = MIC_MIN_BIN_SIZE,
  dilutions     = MIC_STANDARD_DILUTIONS,
  hist_path     = file.path(OUT_HIST, paste0(dataset_name, "_bins.png")),
  dataset_label = "SPN Penicillin"
)

enc <- encode_lineages_spn(
  lineages_df = aligned$lineages,
  sublin_df   = aligned$sublineages,
  pheno_vec   = binning$bins
)

stan_list <- build_stan_inference(
  pheno      = binning$bins,
  geno_mat   = aligned$geno_mat,
  lin_mat    = enc$lineage_matrix,
  sublin_mat = enc$sublineage_matrix,
  parent_lin = enc$parent_lineage,
  K          = binning$K,
  mic_bkpts  = binning$mic_breakpoints
)

write_dataset(
  stan_list    = stan_list,
  sample_ids   = aligned$sample_ids,
  variant_names = rownames(spn_geno),
  parent_lin   = enc$parent_lineage,
  outdir       = out_dir,
  dataset_name = dataset_name
)


############################################################
## 12: SPN TRIMETHOPRIM MIC (ordinal, equal-frequency bins)
############################################################

message("\n=== 12: SPN trimethoprim MIC (equal-freq) ===")
dataset_name <- "12_spn_trimethoprim_MIC_equalfreq"
out_dir <- file.path(OUT_INFER, dataset_name)

aligned <- intersect_and_align(
  pheno_df    = tmp_raw,
  geno        = spn_geno,
  lineages_df = spn_lin$lineages,
  sublin_df   = spn_lin$sublineages,
  id_col      = "ID",
  geno_in_cols = TRUE
)

binning <- bin_mic_equalfreq(
  mic_numeric   = aligned$pheno$MIC_num,
  min_bin_size  = MIC_MIN_BIN_SIZE,
  dilutions     = MIC_STANDARD_DILUTIONS,
  hist_path     = file.path(OUT_HIST, paste0(dataset_name, "_bins.png")),
  dataset_label = "SPN Trimethoprim"
)

enc <- encode_lineages_spn(
  lineages_df = aligned$lineages,
  sublin_df   = aligned$sublineages,
  pheno_vec   = binning$bins
)

stan_list <- build_stan_inference(
  pheno      = binning$bins,
  geno_mat   = aligned$geno_mat,
  lin_mat    = enc$lineage_matrix,
  sublin_mat = enc$sublineage_matrix,
  parent_lin = enc$parent_lineage,
  K          = binning$K,
  mic_bkpts  = binning$mic_breakpoints
)

write_dataset(
  stan_list    = stan_list,
  sample_ids   = aligned$sample_ids,
  variant_names = rownames(spn_geno),
  parent_lin   = enc$parent_lineage,
  outdir       = out_dir,
  dataset_name = dataset_name
)


############################################################
## 13: SPN TRIMETHOPRIM MIC (ordinal, peak-valley bins)
############################################################

message("\n=== 13: SPN trimethoprim MIC (peaks) ===")
dataset_name <- "13_spn_trimethoprim_MIC_peaks"
out_dir <- file.path(OUT_INFER, dataset_name)

aligned <- intersect_and_align(
  pheno_df    = tmp_raw,
  geno        = spn_geno,
  lineages_df = spn_lin$lineages,
  sublin_df   = spn_lin$sublineages,
  id_col      = "ID",
  geno_in_cols = TRUE
)

binning <- bin_mic_peaks(
  mic_numeric   = aligned$pheno$MIC_num,
  min_bin_size  = MIC_MIN_BIN_SIZE,
  dilutions     = MIC_STANDARD_DILUTIONS,
  hist_path     = file.path(OUT_HIST, paste0(dataset_name, "_bins.png")),
  dataset_label = "SPN Trimethoprim"
)

enc <- encode_lineages_spn(
  lineages_df = aligned$lineages,
  sublin_df   = aligned$sublineages,
  pheno_vec   = binning$bins
)

stan_list <- build_stan_inference(
  pheno      = binning$bins,
  geno_mat   = aligned$geno_mat,
  lin_mat    = enc$lineage_matrix,
  sublin_mat = enc$sublineage_matrix,
  parent_lin = enc$parent_lineage,
  K          = binning$K,
  mic_bkpts  = binning$mic_breakpoints
)

write_dataset(
  stan_list    = stan_list,
  sample_ids   = aligned$sample_ids,
  variant_names = rownames(spn_geno),
  parent_lin   = enc$parent_lineage,
  outdir       = out_dir,
  dataset_name = dataset_name
)


############################################################
## 14: TB RIFAMPICIN MIC (ordinal, equal-frequency bins)
############################################################

message("\n=== 14: TB rifampicin MIC (equal-freq) ===")
dataset_name <- "14_tb_rifampicin_MIC_equalfreq"
out_dir <- file.path(OUT_INFER, dataset_name)

binning <- bin_mic_equalfreq(
  mic_numeric   = tb_mic_num,
  min_bin_size  = MIC_MIN_BIN_SIZE,
  dilutions     = MIC_STANDARD_DILUTIONS,
  hist_path     = file.path(OUT_HIST, paste0(dataset_name, "_bins.png")),
  dataset_label = "TB Rifampicin"
)

enc <- encode_lineages_tb(
  lineages_df = tb_lin_aligned,
  pheno_vec   = binning$bins
)

stan_list <- build_stan_inference(
  pheno      = binning$bins,
  geno_mat   = tb_geno_mat,
  lin_mat    = enc$lineage_matrix,
  sublin_mat = enc$sublineage_matrix,
  parent_lin = enc$parent_lineage,
  K          = binning$K,
  mic_bkpts  = binning$mic_breakpoints
)

write_dataset(
  stan_list    = stan_list,
  sample_ids   = tb_sample_ids,
  variant_names = tb_var_names,
  parent_lin   = enc$parent_lineage,
  outdir       = out_dir,
  dataset_name = dataset_name
)


############################################################
## 15: TB RIFAMPICIN MIC (ordinal, peak-valley bins)
############################################################

message("\n=== 15: TB rifampicin MIC (peaks) ===")
dataset_name <- "15_tb_rifampicin_MIC_peaks"
out_dir <- file.path(OUT_INFER, dataset_name)

binning <- bin_mic_peaks(
  mic_numeric   = tb_mic_num,
  min_bin_size  = MIC_MIN_BIN_SIZE,
  dilutions     = MIC_STANDARD_DILUTIONS,
  hist_path     = file.path(OUT_HIST, paste0(dataset_name, "_bins.png")),
  dataset_label = "TB Rifampicin"
)

enc <- encode_lineages_tb(
  lineages_df = tb_lin_aligned,
  pheno_vec   = binning$bins
)

stan_list <- build_stan_inference(
  pheno      = binning$bins,
  geno_mat   = tb_geno_mat,
  lin_mat    = enc$lineage_matrix,
  sublin_mat = enc$sublineage_matrix,
  parent_lin = enc$parent_lineage,
  K          = binning$K,
  mic_bkpts  = binning$mic_breakpoints
)

write_dataset(
  stan_list    = stan_list,
  sample_ids   = tb_sample_ids,
  variant_names = tb_var_names,
  parent_lin   = enc$parent_lineage,
  outdir       = out_dir,
  dataset_name = dataset_name
)


message("\n=== All 15 inference datasets written to: ", OUT_INFER, " ===")
