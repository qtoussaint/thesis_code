############################################################
## write_unitig_datasets.R
##
## Builds ONLY the unitig-genotype SPN penicillin datasets, so they can be
## (re)generated without re-running the full inference/prediction pipelines:
##
##   17 SPN penicillin MIC (ordinal, >= 5% per bin, UNITIG genotype)
##   18 SPN penicillin MIC (ordinal, breakpoint minima binning, K=5, UNITIG genotype)
##
## For each it writes the inference dataset plus the two prediction splits
## (80/20 random and LOSO). The unitig matrix is loaded once and reused.
##
## These same blocks also live in write_inference_jsons.R / write_prediction_jsons.R;
## this script is a convenience runner and writes to the identical OUT_INFER / OUT_PRED
## locations. It deliberately does NOT write inputs_manifest.csv (that belongs to the
## full pipeline run and would be clobbered here).
##
## All input paths and parameters come from config.R; shared functions from utils.R.
############################################################

script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("--file=", script_args, value = TRUE)
script_dir <- if (length(file_arg) > 0) dirname(normalizePath(sub("--file=", "", file_arg))) else "."
source(file.path(script_dir, "config.R"))
source(file.path(script_dir, "utils.R"))

# Open a null graphics device so gridExtra::arrangeGrob() does not create Rplots.pdf.
pdf(NULL)


############################################################
## SPN DATA LOADING (penicillin MIC + lineages + unitigs)
############################################################

message("\n=== Loading SPN data ===")
spn_lin  <- load_spn_lineages(SPN_LINEAGES_PATH)
spn_meta <- read.csv(SPN_METADATA_PATH, na.strings = "")

pen_raw <- data.frame(
  ID  = spn_meta$Lane.name,
  MIC = as.character(spn_meta$Benzylpenicillin.MIC..ug.mL.)
)
pen_raw <- pen_raw[complete.cases(pen_raw), ]
pen_raw$MIC_num <- clean_mic_values(pen_raw$MIC, SPN_PEN_MIC_REPLACEMENTS)
pen_raw <- pen_raw[!is.na(pen_raw$MIC_num), ]

message("\n=== Loading SPN unitigs ===")
spn_unitigs <- load_spn_unitigs(SPN_UNITIGS_PATH, SPN_UNITIG_MIN_AF, SPN_UNITIG_MAX_AF)


############################################################
## 17 (inference): SPN PENICILLIN MIC (ordinal, auto-binned, unitig genotype)
############################################################

message("\n=== 17-infer: SPN penicillin MIC unitigs (ordinal) ===")
dataset_name <- "17_spn_penicillin_MIC_unitigs"
out_dir <- file.path(OUT_INFER, dataset_name)

aligned <- intersect_and_align(
  pheno_df    = pen_raw,
  geno        = spn_unitigs,
  lineages_df = spn_lin$lineages,
  sublin_df   = spn_lin$sublineages,
  id_col      = "ID",
  geno_in_cols = TRUE
)

binning <- bin_mic_auto(
  mic_numeric    = aligned$pheno$MIC_num,
  min_bin_frac   = MIC_MIN_BIN_FRAC,
  dilutions      = MIC_STANDARD_DILUTIONS,
  hist_path      = file.path(OUT_HIST, paste0(dataset_name, "_bins.png")),
  dataset_label  = "SPN Penicillin (unitigs)",
  drug_name      = "benzylpenicillin",
  species_name   = "S. pneumoniae",
  strategy_label = "doubling dilutions with 5% minimum frequency per bin, unitig genotype"
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
  variant_names = rownames(spn_unitigs),
  variant_positions = unitig_variant_positions(rownames(spn_unitigs), SPN_UNITIG_VARINDEX_PATH),
  parent_lin   = enc$parent_lineage,
  outdir       = out_dir,
  dataset_name = dataset_name
)


############################################################
## 18 (inference): SPN PENICILLIN MIC (ordinal, breakpoint minima binning, unitig genotype)
############################################################

message("\n=== 18-infer: SPN penicillin MIC minima binning unitigs (ordinal) ===")
dataset_name <- "18_spn_penicillin_MIC_minimabinning_unitigs"
out_dir <- file.path(OUT_INFER, dataset_name)

aligned <- intersect_and_align(
  pheno_df    = pen_raw,
  geno        = spn_unitigs,
  lineages_df = spn_lin$lineages,
  sublin_df   = spn_lin$sublineages,
  id_col      = "ID",
  geno_in_cols = TRUE
)

# Fixed breakpoints placed at the natural minima of the SPN penicillin
# MIC distribution, giving K = 5 ordered categories.
mic_bkpts_full <- c(0, 0.032, 0.065, 0.2, 2, 6)
bins <- as.integer(cut(aligned$pheno$MIC_num,
                       breaks = mic_bkpts_full,
                       include.lowest = TRUE))
K <- length(mic_bkpts_full) - 1L
counts <- tabulate(bins, nbins = K)
bin_labels <- paste0("(", mic_bkpts_full[-length(mic_bkpts_full)],
                     ", ", mic_bkpts_full[-1], "]")
cat("\n--- 18_spn_penicillin_MIC_minimabinning_unitigs ---\n")
cat("K =", K, "bins (fixed breakpoint minima)\n")
print(data.frame(bin = 1:K, interval = bin_labels, count = counts))

.save_bin_histogram(
  mic_numeric    = aligned$pheno$MIC_num,
  bins_after     = bins,
  breaks_after   = mic_bkpts_full,
  dataset_label  = "SPN Penicillin (breakpoint minima binning, unitigs)",
  hist_path      = file.path(OUT_HIST, paste0(dataset_name, "_bins.png")),
  drug_name      = "benzylpenicillin",
  species_name   = "S. pneumoniae",
  strategy_label = "fixed breakpoints at natural minima"
)

enc <- encode_lineages_spn(
  lineages_df = aligned$lineages,
  sublin_df   = aligned$sublineages,
  pheno_vec   = bins
)

stan_list <- build_stan_inference(
  pheno      = bins,
  geno_mat   = aligned$geno_mat,
  lin_mat    = enc$lineage_matrix,
  sublin_mat = enc$sublineage_matrix,
  parent_lin = enc$parent_lineage,
  K          = K,
  mic_bkpts  = mic_bkpts_full[-c(1, length(mic_bkpts_full))]
)

write_dataset(
  stan_list    = stan_list,
  sample_ids   = aligned$sample_ids,
  variant_names = rownames(spn_unitigs),
  variant_positions = unitig_variant_positions(rownames(spn_unitigs), SPN_UNITIG_VARINDEX_PATH),
  parent_lin   = enc$parent_lineage,
  outdir       = out_dir,
  dataset_name = dataset_name
)


############################################################
## 17 (prediction): SPN PENICILLIN MIC (ordinal, auto-binned, unitig genotype)
############################################################

message("\n=== 17-pred: SPN penicillin MIC unitigs (ordinal) ===")
dataset_name <- "17_spn_penicillin_MIC_unitigs"
out_dir <- file.path(OUT_PRED, dataset_name)

aligned <- intersect_and_align(
  pheno_df    = pen_raw,
  geno        = spn_unitigs,
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
  dataset_label = "SPN Penicillin (unitigs)"
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
  variant_names = rownames(spn_unitigs),
  variant_positions = unitig_variant_positions(rownames(spn_unitigs), SPN_UNITIG_VARINDEX_PATH),
  parent_lin   = enc$parent_lineage,
  outdir       = out_dir,
  dataset_name = dataset_name,
  test_ids     = pred$test_ids,
  test_pheno   = pred$stan_list$test_phenotype
)

save_prediction_ordinal_histogram(
  train_pheno  = pred$stan_list$training_phenotype,
  test_pheno   = pred$stan_list$test_phenotype,
  K            = binning$K,
  breakpoints  = binning$breakpoints,
  dataset_label = "SPN Penicillin (unitigs)",
  hist_path    = file.path(OUT_HIST, paste0(dataset_name, "_pred_dist.png")),
  drug_name    = "benzylpenicillin",
  species_name = "S. pneumoniae",
  split_label  = "80/20 random"
)


# ---- LOSO ----
message("\n=== 17-pred LOSO: SPN penicillin MIC unitigs (ordinal) ===")
dataset_name <- "17_spn_penicillin_MIC_unitigs_loso"
out_dir <- file.path(OUT_PRED, dataset_name)

pred <- build_stan_prediction_loso(
  pheno         = binning$bins,
  geno_mat      = aligned$geno_mat,
  lin_mat       = enc$lineage_matrix,
  sublin_mat    = enc$sublineage_matrix,
  parent_lin    = enc$parent_lineage,
  sample_ids    = aligned$sample_ids,
  sublineage_vec = aligned$sublineages[[2]],
  K             = binning$K,
  mic_bkpts     = binning$mic_breakpoints
)

write_dataset(
  stan_list    = pred$stan_list,
  sample_ids   = pred$train_ids,
  variant_names = rownames(spn_unitigs),
  variant_positions = unitig_variant_positions(rownames(spn_unitigs), SPN_UNITIG_VARINDEX_PATH),
  parent_lin   = enc$parent_lineage,
  outdir       = out_dir,
  dataset_name = dataset_name,
  test_ids     = pred$test_ids,
  test_pheno   = pred$stan_list$test_phenotype
)

save_prediction_ordinal_histogram(
  train_pheno    = pred$stan_list$training_phenotype,
  test_pheno     = pred$stan_list$test_phenotype,
  K              = binning$K,
  breakpoints    = binning$breakpoints,
  dataset_label  = "SPN Penicillin (unitigs)",
  hist_path      = file.path(OUT_HIST, paste0(dataset_name, "_pred_dist.png")),
  drug_name      = "benzylpenicillin",
  species_name   = "S. pneumoniae",
  split_label    = "LOSO",
  held_out_name  = pred$held_out_sublineage
)


############################################################
## 18 (prediction): SPN PENICILLIN MIC (ordinal, breakpoint minima binning, unitig genotype)
############################################################

message("\n=== 18-pred: SPN penicillin MIC minima binning unitigs (ordinal) ===")
dataset_name <- "18_spn_penicillin_MIC_minimabinning_unitigs"
out_dir <- file.path(OUT_PRED, dataset_name)

aligned <- intersect_and_align(
  pheno_df    = pen_raw,
  geno        = spn_unitigs,
  lineages_df = spn_lin$lineages,
  sublin_df   = spn_lin$sublineages,
  id_col      = "ID",
  geno_in_cols = TRUE
)

# Fixed breakpoints placed at the natural minima of the SPN penicillin
# MIC distribution, giving K = 5 ordered categories.
mic_bkpts_full <- c(0, 0.032, 0.065, 0.2, 2, 6)
bins <- as.integer(cut(aligned$pheno$MIC_num,
                       breaks = mic_bkpts_full,
                       include.lowest = TRUE))
K <- length(mic_bkpts_full) - 1L

enc <- encode_lineages_spn(
  lineages_df = aligned$lineages,
  sublin_df   = aligned$sublineages,
  pheno_vec   = bins
)

pred <- build_stan_prediction(
  pheno      = bins,
  geno_mat   = aligned$geno_mat,
  lin_mat    = enc$lineage_matrix,
  sublin_mat = enc$sublineage_matrix,
  parent_lin = enc$parent_lineage,
  sample_ids = aligned$sample_ids,
  K          = K,
  mic_bkpts  = mic_bkpts_full[-c(1, length(mic_bkpts_full))]
)

write_dataset(
  stan_list    = pred$stan_list,
  sample_ids   = pred$train_ids,
  variant_names = rownames(spn_unitigs),
  variant_positions = unitig_variant_positions(rownames(spn_unitigs), SPN_UNITIG_VARINDEX_PATH),
  parent_lin   = enc$parent_lineage,
  outdir       = out_dir,
  dataset_name = dataset_name,
  test_ids     = pred$test_ids,
  test_pheno   = pred$stan_list$test_phenotype
)

save_prediction_ordinal_histogram(
  train_pheno  = pred$stan_list$training_phenotype,
  test_pheno   = pred$stan_list$test_phenotype,
  K            = K,
  breakpoints  = mic_bkpts_full,
  dataset_label = "SPN Penicillin (breakpoint minima binning, unitigs)",
  hist_path    = file.path(OUT_HIST, paste0(dataset_name, "_pred_dist.png")),
  drug_name    = "benzylpenicillin",
  species_name = "S. pneumoniae",
  split_label  = "80/20 random"
)


# ---- LOSO ----
message("\n=== 18-pred LOSO: SPN penicillin MIC minima binning unitigs (ordinal) ===")
dataset_name <- "18_spn_penicillin_MIC_minimabinning_unitigs_loso"
out_dir <- file.path(OUT_PRED, dataset_name)

pred <- build_stan_prediction_loso(
  pheno         = bins,
  geno_mat      = aligned$geno_mat,
  lin_mat       = enc$lineage_matrix,
  sublin_mat    = enc$sublineage_matrix,
  parent_lin    = enc$parent_lineage,
  sample_ids    = aligned$sample_ids,
  sublineage_vec = aligned$sublineages[[2]],
  K             = K,
  mic_bkpts     = mic_bkpts_full[-c(1, length(mic_bkpts_full))]
)

write_dataset(
  stan_list    = pred$stan_list,
  sample_ids   = pred$train_ids,
  variant_names = rownames(spn_unitigs),
  variant_positions = unitig_variant_positions(rownames(spn_unitigs), SPN_UNITIG_VARINDEX_PATH),
  parent_lin   = enc$parent_lineage,
  outdir       = out_dir,
  dataset_name = dataset_name,
  test_ids     = pred$test_ids,
  test_pheno   = pred$stan_list$test_phenotype
)

save_prediction_ordinal_histogram(
  train_pheno    = pred$stan_list$training_phenotype,
  test_pheno     = pred$stan_list$test_phenotype,
  K              = K,
  breakpoints    = mic_bkpts_full,
  dataset_label  = "SPN Penicillin (breakpoint minima binning, unitigs)",
  hist_path      = file.path(OUT_HIST, paste0(dataset_name, "_pred_dist.png")),
  drug_name      = "benzylpenicillin",
  species_name   = "S. pneumoniae",
  split_label    = "LOSO",
  held_out_name  = pred$held_out_sublineage
)


message("\n=== Unitig datasets written: inference 17,18 to ", OUT_INFER,
        " and prediction 17,18 (+loso) to ", OUT_PRED, " ===")
