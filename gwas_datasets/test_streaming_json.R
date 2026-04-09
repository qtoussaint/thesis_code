## Quick test: compare write_stan_json_streaming() vs cmdstanr::write_stan_json()
## on the SPN penicillin binary dataset (dataset 01).

script_dir <- "/nfs/research/jlees/jacqueline/thesis_code/gwas_datasets"
source(file.path(script_dir, "config.R"))
source(file.path(script_dir, "utils.R"))

pdf(NULL)

# ---- Load SPN data (same as write_inference_jsons.R lines 62-75) ----
message("Loading SPN data...")
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

# ---- Build dataset 01 stan_list ----
message("Building dataset 01 stan_list...")
pen_bin <- pen_raw
pen_bin$binary <- mic_to_binary(
  mic_numeric = pen_bin$MIC_num,
  s_max       = SPN_PEN_BINARY_S_MAX,
  r_min       = SPN_PEN_BINARY_R_MIN
)
pen_bin <- pen_bin[!is.na(pen_bin$binary), ]
pen_bin$pheno <- pen_bin$binary

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

stan_list <- build_stan_inference(
  pheno      = aligned$pheno$pheno,
  geno_mat   = aligned$geno_mat,
  lin_mat    = enc$lineage_matrix,
  sublin_mat = enc$sublineage_matrix,
  parent_lin = enc$parent_lineage
)

# ---- Write with both methods ----
tmp_dir <- tempdir()
ref_path    <- file.path(tmp_dir, "ref_cmdstanr.json")
stream_path <- file.path(tmp_dir, "streaming.json")

message("Writing with cmdstanr::write_stan_json()...")
write_stan_json(data = stan_list, file = ref_path)

message("Writing with write_stan_json_streaming()...")
write_stan_json_streaming(data = stan_list, file = stream_path)

# ---- Compare files ----
message("Comparing output files...")

# First check: MD5 checksums
md5_ref    <- tools::md5sum(ref_path)
md5_stream <- tools::md5sum(stream_path)
message("  ref MD5:    ", md5_ref)
message("  stream MD5: ", md5_stream)

if (md5_ref == md5_stream) {
  message("\nPASSED: files are byte-identical.")
} else {
  message("\nMD5 differs — checking line-by-line (with and without whitespace)...")
  ref_con    <- file(ref_path, open = "rt")
  stream_con <- file(stream_path, open = "rt")
  on.exit({ close(ref_con); close(stream_con) }, add = TRUE)

  line_num      <- 0L
  n_diff_raw    <- 0L
  n_diff_normed <- 0L
  max_show      <- 5L
  repeat {
    r <- readLines(ref_con, n = 1L)
    s <- readLines(stream_con, n = 1L)
    line_num <- line_num + 1L
    if (length(r) == 0L && length(s) == 0L) break
    if (length(r) == 0L || length(s) == 0L || r != s) {
      n_diff_raw <- n_diff_raw + 1L
      # Check if difference is whitespace-only
      r_normed <- gsub("\\s", "", if (length(r)) r else "")
      s_normed <- gsub("\\s", "", if (length(s)) s else "")
      if (r_normed != s_normed) {
        n_diff_normed <- n_diff_normed + 1L
        if (n_diff_normed <= max_show) {
          message(sprintf("  NON-WHITESPACE diff at line %d:", line_num))
          message("    ref:    ", substr(r_normed, 1, 120))
          message("    stream: ", substr(s_normed, 1, 120))
        }
      }
    }
  }
  message(sprintf("\nTotal differing lines (raw): %d", n_diff_raw))
  message(sprintf("Total differing lines (ignoring whitespace): %d", n_diff_normed))
  if (n_diff_normed == 0L) {
    message("PASSED: content identical (whitespace-only formatting differences).")
  } else {
    message("FAILED: non-whitespace differences detected.")
  }
}

# File sizes
message(sprintf("  ref size:    %.1f MB", file.info(ref_path)$size / 1e6))
message(sprintf("  stream size: %.1f MB", file.info(stream_path)$size / 1e6))

# Clean up
unlink(c(ref_path, stream_path))
