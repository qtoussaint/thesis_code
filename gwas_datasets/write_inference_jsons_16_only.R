############################################################
## write_inference_jsons_16_only.R
## Generates ONLY dataset 16 (16_spn_penicillin_MIC_minimabinning).
## Mirrors the corresponding block in write_inference_jsons.R but
## skips loading TB data and all other SPN datasets.
############################################################

script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("--file=", script_args, value = TRUE)
script_dir <- if (length(file_arg) > 0) dirname(normalizePath(sub("--file=", "", file_arg))) else "."
source(file.path(script_dir, "config.R"))
source(file.path(script_dir, "utils.R"))

pdf(NULL)

############################################################
## SPN DATA LOADING (subset needed by dataset 16)
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

############################################################
## 16: SPN PENICILLIN MIC (ordinal, breakpoint minima binning)
############################################################

message("\n=== 16: SPN penicillin MIC breakpoint minima binning (ordinal) ===")
dataset_name <- "16_spn_penicillin_MIC_minimabinning"
out_dir <- file.path(OUT_INFER, dataset_name)

aligned <- intersect_and_align(
  pheno_df    = pen_raw,
  geno        = spn_geno,
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
cat("\n--- 16_spn_penicillin_MIC_minimabinning ---\n")
cat("K =", K, "bins (fixed breakpoint minima)\n")
print(data.frame(bin = 1:K, interval = bin_labels, count = counts))

.save_bin_histogram(
  mic_numeric    = aligned$pheno$MIC_num,
  bins_after     = bins,
  breaks_after   = mic_bkpts_full,
  dataset_label  = "SPN Penicillin (breakpoint minima binning)",
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
  variant_names = rownames(spn_geno),
  parent_lin   = enc$parent_lineage,
  outdir       = out_dir,
  dataset_name = dataset_name
)

message("\n=== Done: dataset 16 only ===")
