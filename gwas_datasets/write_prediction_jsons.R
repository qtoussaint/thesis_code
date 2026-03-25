############################################################
## write_prediction_jsons.R
##
## Creates 15 prediction Stan JSON datasets (80/20 train/test):
##
##   SPN PENICILLIN
##   01 SPN penicillin    binary
##   02 SPN penicillin    MIC (ordinal, >= 5% per bin)
##   10 SPN penicillin    MIC (ordinal, coarse dilutions, >= 5% per bin)
##   11 SPN penicillin    MIC (ordinal, >= 10% per bin)
##   03 SPN penicillin    continuous (log2)
##
##   SPN TRIMETHOPRIM
##   04 SPN trimethoprim  binary
##   05 SPN trimethoprim  MIC (ordinal, >= 5% per bin)
##   12 SPN trimethoprim  MIC (ordinal, coarse dilutions, >= 5% per bin)
##   13 SPN trimethoprim  MIC (ordinal, >= 10% per bin)
##   06 SPN trimethoprim  continuous (log2)
##
##   TB RIFAMPICIN
##   07 TB  rifampicin    binary
##   08 TB  rifampicin    MIC (ordinal, >= 5% per bin)
##   14 TB  rifampicin    MIC (ordinal, coarse dilutions, >= 5% per bin)
##   15 TB  rifampicin    MIC (ordinal, >= 10% per bin)
##   09 TB  rifampicin    continuous (log2)
##
## MIC binning is deterministic given the data and MIC_MIN_BIN_FRAC.
## Histograms are saved by the inference script (hist_path = NULL here).
############################################################

script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("--file=", script_args, value = TRUE)
script_dir <- if (length(file_arg) > 0) dirname(normalizePath(sub("--file=", "", file_arg))) else "."
source(file.path(script_dir, "config.R"))
source(file.path(script_dir, "utils.R"))


############################################################
## SPN DATA LOADING
############################################################

message("\n=== Loading SPN data ===")
spn_geno <- load_spn_genotype(SPN_GENOTYPE_PATH)
spn_lin  <- load_spn_lineages(SPN_LINEAGES_PATH)

spn_meta <- read.csv(SPN_METADATA_PATH, na.strings = "")

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

pen_bin <- pen_raw
pen_bin$binary <- mic_to_binary(
  mic_numeric = pen_bin$MIC_num,
  s_max       = SPN_PEN_BINARY_S_MAX,
  r_min       = SPN_PEN_BINARY_R_MIN
)
pen_bin <- pen_bin[!is.na(pen_bin$binary), ]
pen_bin$pheno <- pen_bin$binary

tmp_bin <- tmp_raw
tmp_bin$binary <- mic_to_binary(
  mic_numeric = tmp_bin$MIC_num,
  s_max       = SPN_TMP_BINARY_S_MAX,
  r_min       = SPN_TMP_BINARY_R_MIN
)
tmp_bin <- tmp_bin[!is.na(tmp_bin$binary), ]
tmp_bin$pheno <- tmp_bin$binary


############################################################
## 01: SPN PENICILLIN BINARY
############################################################

message("\n=== 01-pred: SPN penicillin binary ===")
dataset_name <- "01_spn_penicillin_binary"
out_dir <- file.path(OUT_PRED, dataset_name)

aligned <- intersect_and_align(
  pheno_df    = pen_bin,
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

pred <- build_stan_prediction(
  pheno      = aligned$pheno$pheno,
  geno_mat   = aligned$geno_mat,
  lin_mat    = enc$lineage_matrix,
  sublin_mat = enc$sublineage_matrix,
  parent_lin = enc$parent_lineage,
  sample_ids = aligned$sample_ids
)

write_dataset(
  stan_list    = pred$stan_list,
  sample_ids   = pred$train_ids,
  variant_names = rownames(spn_geno),
  parent_lin   = enc$parent_lineage,
  outdir       = out_dir,
  dataset_name = dataset_name,
  test_ids     = pred$test_ids,
  test_pheno   = pred$stan_list$test_phenotype
)


############################################################
## 02: SPN PENICILLIN MIC (ordinal, auto-binned)
############################################################

message("\n=== 02-pred: SPN penicillin MIC (ordinal) ===")
dataset_name <- "02_spn_penicillin_MIC"
out_dir <- file.path(OUT_PRED, dataset_name)

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
  min_bin_frac  = MIC_MIN_BIN_FRAC,
  dilutions     = MIC_STANDARD_DILUTIONS,
  hist_path     = NULL,   # histograms already saved by inference script
  dataset_label = "SPN Penicillin"
)

enc <- encode_lineages_spn(
  lineages_df = aligned$lineages,
  sublin_df   = aligned$sublineages,
  pheno_vec   = binning$bins
)

pred <- build_stan_prediction(
  pheno      = binning$bins,
  geno_mat   = aligned$geno_mat,
  lin_mat    = enc$lineage_matrix,
  sublin_mat = enc$sublineage_matrix,
  parent_lin = enc$parent_lineage,
  sample_ids = aligned$sample_ids,
  K          = binning$K,
  mic_bkpts  = binning$mic_breakpoints
)

write_dataset(
  stan_list    = pred$stan_list,
  sample_ids   = pred$train_ids,
  variant_names = rownames(spn_geno),
  parent_lin   = enc$parent_lineage,
  outdir       = out_dir,
  dataset_name = dataset_name,
  test_ids     = pred$test_ids,
  test_pheno   = pred$stan_list$test_phenotype
)


############################################################
## 10: SPN PENICILLIN MIC (ordinal, coarse dilutions)
############################################################

message("\n=== 10-pred: SPN penicillin MIC coarse dilutions (ordinal) ===")
dataset_name <- "10_spn_penicillin_MIC_coarse_dilutions"
out_dir <- file.path(OUT_PRED, dataset_name)

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
  min_bin_frac  = MIC_MIN_BIN_FRAC,
  dilutions     = MIC_COARSE_DILUTIONS,
  hist_path     = NULL,
  dataset_label = "SPN Penicillin (coarse dilutions)"
)

enc <- encode_lineages_spn(
  lineages_df = aligned$lineages,
  sublin_df   = aligned$sublineages,
  pheno_vec   = binning$bins
)

pred <- build_stan_prediction(
  pheno      = binning$bins,
  geno_mat   = aligned$geno_mat,
  lin_mat    = enc$lineage_matrix,
  sublin_mat = enc$sublineage_matrix,
  parent_lin = enc$parent_lineage,
  sample_ids = aligned$sample_ids,
  K          = binning$K,
  mic_bkpts  = binning$mic_breakpoints
)

write_dataset(
  stan_list    = pred$stan_list,
  sample_ids   = pred$train_ids,
  variant_names = rownames(spn_geno),
  parent_lin   = enc$parent_lineage,
  outdir       = out_dir,
  dataset_name = dataset_name,
  test_ids     = pred$test_ids,
  test_pheno   = pred$stan_list$test_phenotype
)


############################################################
## 11: SPN PENICILLIN MIC (ordinal, large min bin)
############################################################

message("\n=== 11-pred: SPN penicillin MIC large minbin (ordinal) ===")
dataset_name <- "11_spn_penicillin_MIC_large_minbin"
out_dir <- file.path(OUT_PRED, dataset_name)

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
  min_bin_frac  = MIC_MIN_BIN_FRAC_COARSE,
  dilutions     = MIC_STANDARD_DILUTIONS,
  hist_path     = NULL,
  dataset_label = "SPN Penicillin (large minbin)"
)

enc <- encode_lineages_spn(
  lineages_df = aligned$lineages,
  sublin_df   = aligned$sublineages,
  pheno_vec   = binning$bins
)

pred <- build_stan_prediction(
  pheno      = binning$bins,
  geno_mat   = aligned$geno_mat,
  lin_mat    = enc$lineage_matrix,
  sublin_mat = enc$sublineage_matrix,
  parent_lin = enc$parent_lineage,
  sample_ids = aligned$sample_ids,
  K          = binning$K,
  mic_bkpts  = binning$mic_breakpoints
)

write_dataset(
  stan_list    = pred$stan_list,
  sample_ids   = pred$train_ids,
  variant_names = rownames(spn_geno),
  parent_lin   = enc$parent_lineage,
  outdir       = out_dir,
  dataset_name = dataset_name,
  test_ids     = pred$test_ids,
  test_pheno   = pred$stan_list$test_phenotype
)


############################################################
## 03: SPN PENICILLIN CONTINUOUS (log2 MIC)
############################################################

message("\n=== 03-pred: SPN penicillin continuous ===")
dataset_name <- "03_spn_penicillin_continuous"
out_dir <- file.path(OUT_PRED, dataset_name)

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

pred <- build_stan_prediction(
  pheno      = log2_mic,
  geno_mat   = aligned$geno_mat,
  lin_mat    = enc$lineage_matrix,
  sublin_mat = enc$sublineage_matrix,
  parent_lin = enc$parent_lineage,
  sample_ids = aligned$sample_ids
)

write_dataset(
  stan_list    = pred$stan_list,
  sample_ids   = pred$train_ids,
  variant_names = rownames(spn_geno),
  parent_lin   = enc$parent_lineage,
  outdir       = out_dir,
  dataset_name = dataset_name,
  test_ids     = pred$test_ids,
  test_pheno   = pred$stan_list$test_phenotype
)


############################################################
## 04: SPN TRIMETHOPRIM BINARY
############################################################

message("\n=== 04-pred: SPN trimethoprim binary ===")
dataset_name <- "04_spn_trimethoprim_binary"
out_dir <- file.path(OUT_PRED, dataset_name)

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

pred <- build_stan_prediction(
  pheno      = aligned$pheno$pheno,
  geno_mat   = aligned$geno_mat,
  lin_mat    = enc$lineage_matrix,
  sublin_mat = enc$sublineage_matrix,
  parent_lin = enc$parent_lineage,
  sample_ids = aligned$sample_ids
)

write_dataset(
  stan_list    = pred$stan_list,
  sample_ids   = pred$train_ids,
  variant_names = rownames(spn_geno),
  parent_lin   = enc$parent_lineage,
  outdir       = out_dir,
  dataset_name = dataset_name,
  test_ids     = pred$test_ids,
  test_pheno   = pred$stan_list$test_phenotype
)

save_binary_histogram(
  mic_numeric   = aligned$pheno$MIC_num,
  binary_vec    = aligned$pheno$pheno,
  dataset_label = "SPN Trimethoprim",
  hist_path     = file.path(OUT_HIST, paste0(dataset_name, "_dist.png")),
  s_max         = SPN_TMP_BINARY_S_MAX,
  r_min         = SPN_TMP_BINARY_R_MIN
)


############################################################
## 05: SPN TRIMETHOPRIM MIC (ordinal, auto-binned)
############################################################

message("\n=== 05-pred: SPN trimethoprim MIC (ordinal) ===")
dataset_name <- "05_spn_trimethoprim_MIC"
out_dir <- file.path(OUT_PRED, dataset_name)

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
  min_bin_frac  = MIC_MIN_BIN_FRAC,
  dilutions     = MIC_STANDARD_DILUTIONS,
  hist_path     = NULL,
  dataset_label = "SPN Trimethoprim"
)

enc <- encode_lineages_spn(
  lineages_df = aligned$lineages,
  sublin_df   = aligned$sublineages,
  pheno_vec   = binning$bins
)

pred <- build_stan_prediction(
  pheno      = binning$bins,
  geno_mat   = aligned$geno_mat,
  lin_mat    = enc$lineage_matrix,
  sublin_mat = enc$sublineage_matrix,
  parent_lin = enc$parent_lineage,
  sample_ids = aligned$sample_ids,
  K          = binning$K,
  mic_bkpts  = binning$mic_breakpoints
)

write_dataset(
  stan_list    = pred$stan_list,
  sample_ids   = pred$train_ids,
  variant_names = rownames(spn_geno),
  parent_lin   = enc$parent_lineage,
  outdir       = out_dir,
  dataset_name = dataset_name,
  test_ids     = pred$test_ids,
  test_pheno   = pred$stan_list$test_phenotype
)


############################################################
## 12: SPN TRIMETHOPRIM MIC (ordinal, coarse dilutions)
############################################################

message("\n=== 12-pred: SPN trimethoprim MIC coarse dilutions (ordinal) ===")
dataset_name <- "12_spn_trimethoprim_MIC_coarse_dilutions"
out_dir <- file.path(OUT_PRED, dataset_name)

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
  min_bin_frac  = MIC_MIN_BIN_FRAC,
  dilutions     = MIC_COARSE_DILUTIONS,
  hist_path     = NULL,
  dataset_label = "SPN Trimethoprim (coarse dilutions)"
)

enc <- encode_lineages_spn(
  lineages_df = aligned$lineages,
  sublin_df   = aligned$sublineages,
  pheno_vec   = binning$bins
)

pred <- build_stan_prediction(
  pheno      = binning$bins,
  geno_mat   = aligned$geno_mat,
  lin_mat    = enc$lineage_matrix,
  sublin_mat = enc$sublineage_matrix,
  parent_lin = enc$parent_lineage,
  sample_ids = aligned$sample_ids,
  K          = binning$K,
  mic_bkpts  = binning$mic_breakpoints
)

write_dataset(
  stan_list    = pred$stan_list,
  sample_ids   = pred$train_ids,
  variant_names = rownames(spn_geno),
  parent_lin   = enc$parent_lineage,
  outdir       = out_dir,
  dataset_name = dataset_name,
  test_ids     = pred$test_ids,
  test_pheno   = pred$stan_list$test_phenotype
)


############################################################
## 13: SPN TRIMETHOPRIM MIC (ordinal, large min bin)
############################################################

message("\n=== 13-pred: SPN trimethoprim MIC large minbin (ordinal) ===")
dataset_name <- "13_spn_trimethoprim_MIC_large_minbin"
out_dir <- file.path(OUT_PRED, dataset_name)

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
  min_bin_frac  = MIC_MIN_BIN_FRAC_COARSE,
  dilutions     = MIC_STANDARD_DILUTIONS,
  hist_path     = NULL,
  dataset_label = "SPN Trimethoprim (large minbin)"
)

enc <- encode_lineages_spn(
  lineages_df = aligned$lineages,
  sublin_df   = aligned$sublineages,
  pheno_vec   = binning$bins
)

pred <- build_stan_prediction(
  pheno      = binning$bins,
  geno_mat   = aligned$geno_mat,
  lin_mat    = enc$lineage_matrix,
  sublin_mat = enc$sublineage_matrix,
  parent_lin = enc$parent_lineage,
  sample_ids = aligned$sample_ids,
  K          = binning$K,
  mic_bkpts  = binning$mic_breakpoints
)

write_dataset(
  stan_list    = pred$stan_list,
  sample_ids   = pred$train_ids,
  variant_names = rownames(spn_geno),
  parent_lin   = enc$parent_lineage,
  outdir       = out_dir,
  dataset_name = dataset_name,
  test_ids     = pred$test_ids,
  test_pheno   = pred$stan_list$test_phenotype
)


############################################################
## 06: SPN TRIMETHOPRIM CONTINUOUS (log2 MIC)
############################################################

message("\n=== 06-pred: SPN trimethoprim continuous ===")
dataset_name <- "06_spn_trimethoprim_continuous"
out_dir <- file.path(OUT_PRED, dataset_name)

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

pred <- build_stan_prediction(
  pheno      = log2_mic,
  geno_mat   = aligned$geno_mat,
  lin_mat    = enc$lineage_matrix,
  sublin_mat = enc$sublineage_matrix,
  parent_lin = enc$parent_lineage,
  sample_ids = aligned$sample_ids
)

write_dataset(
  stan_list    = pred$stan_list,
  sample_ids   = pred$train_ids,
  variant_names = rownames(spn_geno),
  parent_lin   = enc$parent_lineage,
  outdir       = out_dir,
  dataset_name = dataset_name,
  test_ids     = pred$test_ids,
  test_pheno   = pred$stan_list$test_phenotype
)


############################################################
## TB DATA LOADING
############################################################

message("\n=== Loading TB data ===")
tb_geno_raw  <- load_tb_genotype(TB_PRESABS_PATH, TB_VARIANT_IDX_PATH)
tb_geno      <- tb_geno_raw$genotype
tb_var_names <- tb_geno_raw$variant_names
tb_lin_raw   <- load_tb_lineages(TB_LINEAGES_PATH)

tb_pheno_raw <- read.csv(TB_PHENOTYPE_PATH)
tb_pheno_raw <- tb_pheno_raw[, c(1, 27, 44)]
colnames(tb_pheno_raw) <- c("ENA_RUN", "RIF_MIC", "REGENOTYPED_VCF")
tb_pheno_raw$RIF_MIC[tb_pheno_raw$RIF_MIC == ""] <- NA
tb_pheno_raw <- na.omit(tb_pheno_raw)
tb_pheno_raw$RIF_MIC_num <- clean_mic_values(
  tb_pheno_raw$RIF_MIC, TB_RIF_MIC_REPLACEMENTS
)
tb_pheno_raw <- tb_pheno_raw[!is.na(tb_pheno_raw$RIF_MIC_num), ]

sample_in_lin     <- match(tb_pheno_raw$ENA_RUN, tb_lin_raw$X.sample)
tb_pheno_with_lin <- tb_pheno_raw[!is.na(sample_in_lin), ]
tb_lin_aligned    <- tb_lin_raw[na.omit(sample_in_lin), ]
geno_row_idx      <- match(tb_pheno_with_lin$ENA_RUN, tb_pheno_raw$ENA_RUN)
geno_row_idx      <- geno_row_idx[!is.na(geno_row_idx)]
tb_geno_mat       <- as.matrix(tb_geno[geno_row_idx, ])
storage.mode(tb_geno_mat) <- "integer"
tb_sample_ids <- tb_pheno_with_lin$ENA_RUN
tb_mic_num    <- tb_pheno_with_lin$RIF_MIC_num
message("TB aligned samples: ", length(tb_sample_ids))


############################################################
## 07: TB RIFAMPICIN BINARY
############################################################

message("\n=== 07-pred: TB rifampicin binary ===")
dataset_name <- "07_tb_rifampicin_binary"
out_dir <- file.path(OUT_PRED, dataset_name)

tb_binary <- mic_to_binary(tb_mic_num, s_max = TB_RIF_BINARY_THRESHOLD)
keep_idx  <- !is.na(tb_binary)

enc <- encode_lineages_tb(
  lineages_df = tb_lin_aligned[keep_idx, ],
  pheno_vec   = tb_binary[keep_idx]
)

pred <- build_stan_prediction(
  pheno      = tb_binary[keep_idx],
  geno_mat   = tb_geno_mat[keep_idx, , drop = FALSE],
  lin_mat    = enc$lineage_matrix,
  sublin_mat = enc$sublineage_matrix,
  parent_lin = enc$parent_lineage,
  sample_ids = tb_sample_ids[keep_idx]
)

write_dataset(
  stan_list    = pred$stan_list,
  sample_ids   = pred$train_ids,
  variant_names = tb_var_names,
  parent_lin   = enc$parent_lineage,
  outdir       = out_dir,
  dataset_name = dataset_name,
  test_ids     = pred$test_ids,
  test_pheno   = pred$stan_list$test_phenotype
)

save_binary_histogram(
  mic_numeric   = tb_mic_num[keep_idx],
  binary_vec    = tb_binary[keep_idx],
  dataset_label = "TB Rifampicin",
  hist_path     = file.path(OUT_HIST, paste0(dataset_name, "_dist.png")),
  s_max         = TB_RIF_BINARY_THRESHOLD
)


############################################################
## 08: TB RIFAMPICIN MIC (ordinal, auto-binned)
############################################################

message("\n=== 08-pred: TB rifampicin MIC (ordinal) ===")
dataset_name <- "08_tb_rifampicin_MIC"
out_dir <- file.path(OUT_PRED, dataset_name)

binning <- bin_mic_auto(
  mic_numeric   = tb_mic_num,
  min_bin_frac  = MIC_MIN_BIN_FRAC,
  dilutions     = MIC_STANDARD_DILUTIONS,
  hist_path     = NULL,
  dataset_label = "TB Rifampicin"
)

enc <- encode_lineages_tb(
  lineages_df = tb_lin_aligned,
  pheno_vec   = binning$bins
)

pred <- build_stan_prediction(
  pheno      = binning$bins,
  geno_mat   = tb_geno_mat,
  lin_mat    = enc$lineage_matrix,
  sublin_mat = enc$sublineage_matrix,
  parent_lin = enc$parent_lineage,
  sample_ids = tb_sample_ids,
  K          = binning$K,
  mic_bkpts  = binning$mic_breakpoints
)

write_dataset(
  stan_list    = pred$stan_list,
  sample_ids   = pred$train_ids,
  variant_names = tb_var_names,
  parent_lin   = enc$parent_lineage,
  outdir       = out_dir,
  dataset_name = dataset_name,
  test_ids     = pred$test_ids,
  test_pheno   = pred$stan_list$test_phenotype
)


############################################################
## 14: TB RIFAMPICIN MIC (ordinal, coarse dilutions)
############################################################

message("\n=== 14-pred: TB rifampicin MIC coarse dilutions (ordinal) ===")
dataset_name <- "14_tb_rifampicin_MIC_coarse_dilutions"
out_dir <- file.path(OUT_PRED, dataset_name)

binning <- bin_mic_auto(
  mic_numeric   = tb_mic_num,
  min_bin_frac  = MIC_MIN_BIN_FRAC,
  dilutions     = MIC_COARSE_DILUTIONS,
  hist_path     = NULL,
  dataset_label = "TB Rifampicin (coarse dilutions)"
)

enc <- encode_lineages_tb(
  lineages_df = tb_lin_aligned,
  pheno_vec   = binning$bins
)

pred <- build_stan_prediction(
  pheno      = binning$bins,
  geno_mat   = tb_geno_mat,
  lin_mat    = enc$lineage_matrix,
  sublin_mat = enc$sublineage_matrix,
  parent_lin = enc$parent_lineage,
  sample_ids = tb_sample_ids,
  K          = binning$K,
  mic_bkpts  = binning$mic_breakpoints
)

write_dataset(
  stan_list    = pred$stan_list,
  sample_ids   = pred$train_ids,
  variant_names = tb_var_names,
  parent_lin   = enc$parent_lineage,
  outdir       = out_dir,
  dataset_name = dataset_name,
  test_ids     = pred$test_ids,
  test_pheno   = pred$stan_list$test_phenotype
)


############################################################
## 15: TB RIFAMPICIN MIC (ordinal, large min bin)
############################################################

message("\n=== 15-pred: TB rifampicin MIC large minbin (ordinal) ===")
dataset_name <- "15_tb_rifampicin_MIC_large_minbin"
out_dir <- file.path(OUT_PRED, dataset_name)

binning <- bin_mic_auto(
  mic_numeric   = tb_mic_num,
  min_bin_frac  = MIC_MIN_BIN_FRAC_COARSE,
  dilutions     = MIC_STANDARD_DILUTIONS,
  hist_path     = NULL,
  dataset_label = "TB Rifampicin (large minbin)"
)

enc <- encode_lineages_tb(
  lineages_df = tb_lin_aligned,
  pheno_vec   = binning$bins
)

pred <- build_stan_prediction(
  pheno      = binning$bins,
  geno_mat   = tb_geno_mat,
  lin_mat    = enc$lineage_matrix,
  sublin_mat = enc$sublineage_matrix,
  parent_lin = enc$parent_lineage,
  sample_ids = tb_sample_ids,
  K          = binning$K,
  mic_bkpts  = binning$mic_breakpoints
)

write_dataset(
  stan_list    = pred$stan_list,
  sample_ids   = pred$train_ids,
  variant_names = tb_var_names,
  parent_lin   = enc$parent_lineage,
  outdir       = out_dir,
  dataset_name = dataset_name,
  test_ids     = pred$test_ids,
  test_pheno   = pred$stan_list$test_phenotype
)


############################################################
## 09: TB RIFAMPICIN CONTINUOUS (log2 MIC)
############################################################

message("\n=== 09-pred: TB rifampicin continuous ===")
dataset_name <- "09_tb_rifampicin_continuous"
out_dir <- file.path(OUT_PRED, dataset_name)

log2_mic <- log2(tb_mic_num)

enc <- encode_lineages_tb(
  lineages_df = tb_lin_aligned,
  pheno_vec   = log2_mic
)

pred <- build_stan_prediction(
  pheno      = log2_mic,
  geno_mat   = tb_geno_mat,
  lin_mat    = enc$lineage_matrix,
  sublin_mat = enc$sublineage_matrix,
  parent_lin = enc$parent_lineage,
  sample_ids = tb_sample_ids
)

write_dataset(
  stan_list    = pred$stan_list,
  sample_ids   = pred$train_ids,
  variant_names = tb_var_names,
  parent_lin   = enc$parent_lineage,
  outdir       = out_dir,
  dataset_name = dataset_name,
  test_ids     = pred$test_ids,
  test_pheno   = pred$stan_list$test_phenotype
)


message("\n=== All 15 prediction datasets written to: ", OUT_PRED, " ===")
