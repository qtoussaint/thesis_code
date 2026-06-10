#!/usr/bin/env Rscript
# Combined train/test phenotype-distribution figure across all prediction datasets.
#
# Mirrors the per-dataset *_pred_dist.png histograms in
# thesis_results/gwas_datasets/MIC_bin_histograms/ (see save_prediction_*_histogram
# in gwas_datasets/utils.R), but combines every binary + ordinal dataset into
# faceted figures: dodged Train (blue) vs Test (orange) bars over the phenotype
# categories, one panel per dataset, with counts labelled above the bars.
#
# Two figures are produced — one per split strategy (80/20 random, LOSO).
#
# Phenotype vectors are pulled straight from each dataset JSON. The JSONs are huge
# (TB ~1.7 GB) because of the genotype arrays, but training_phenotype/test_phenotype
# sit in the first ~1 KB and K/mic_breakpoints in the last ~100 bytes, so we read only
# those edges.
#
# Run with:
#   mamba activate gwas_pipeline
#   Rscript paper_figures/dataset_split_summary.R
#
# Output: <output_dir>/dataset_split_summary_{random,loso}.png + dataset_split_summary.csv

suppressPackageStartupMessages({
  library(ggplot2)
})

RESULTS_ROOT <- "/nfs/research/jlees/jacqueline/thesis_results"
PRED_ROOT    <- file.path(RESULTS_ROOT, "gwas_datasets", "prediction")
OUTPUT_DIR   <- file.path(RESULTS_ROOT, "paper_figures", "dataset_binning")

# Consistent colours: training = blue, test = orange (matches the ordinal reference
# histograms; binary reference uses firebrick but we keep blue:orange throughout).
SET_COLOURS <- c(Train = "steelblue", Test = "darkorange")

# -----------------------------------------------------------------------------
# Dataset specs: binary + ordinal datasets (continuous excluded — no categories).
# Numbering from write_prediction_jsons.R; ordered drug -> dataset number.
# species/drug are the properly cased names used in the panel strip labels
# (species italicised via label_parsed).
# -----------------------------------------------------------------------------

specs <- list(
  list(dir_stub = "01_spn_penicillin_binary",                species = "S. pneumoniae",   drug = "benzylpenicillin", pheno = "binary",          type = "binary"),
  list(dir_stub = "02_spn_penicillin_MIC",                   species = "S. pneumoniae",   drug = "benzylpenicillin", pheno = "MIC doubling (≥5%)",        type = "ordinal"),
  list(dir_stub = "10_spn_penicillin_MIC_coarse_dilutions",  species = "S. pneumoniae",   drug = "benzylpenicillin", pheno = "MIC 4-fold (≥5%)", type = "ordinal"),
  list(dir_stub = "11_spn_penicillin_MIC_large_minbin",      species = "S. pneumoniae",   drug = "benzylpenicillin", pheno = "MIC doubling (≥10%)",       type = "ordinal"),
  list(dir_stub = "16_spn_penicillin_MIC_minimabinning",     species = "S. pneumoniae",   drug = "benzylpenicillin", pheno = "MIC minima",                type = "ordinal"),

  list(dir_stub = "04_spn_trimethoprim_binary",              species = "S. pneumoniae",   drug = "trimethoprim",     pheno = "binary",          type = "binary"),
  list(dir_stub = "05_spn_trimethoprim_MIC",                 species = "S. pneumoniae",   drug = "trimethoprim",     pheno = "MIC doubling (≥5%)",        type = "ordinal"),
  list(dir_stub = "12_spn_trimethoprim_MIC_coarse_dilutions",species = "S. pneumoniae",   drug = "trimethoprim",     pheno = "MIC 4-fold (≥5%)", type = "ordinal"),
  list(dir_stub = "13_spn_trimethoprim_MIC_large_minbin",    species = "S. pneumoniae",   drug = "trimethoprim",     pheno = "MIC doubling (≥10%)",       type = "ordinal"),

  list(dir_stub = "07_tb_rifampicin_binary",                 species = "M. tuberculosis", drug = "rifampicin",       pheno = "binary",          type = "binary"),
  list(dir_stub = "08_tb_rifampicin_MIC",                    species = "M. tuberculosis", drug = "rifampicin",       pheno = "MIC doubling (≥5%)",        type = "ordinal"),
  list(dir_stub = "14_tb_rifampicin_MIC_coarse_dilutions",   species = "M. tuberculosis", drug = "rifampicin",       pheno = "MIC 4-fold (≥5%)", type = "ordinal"),
  list(dir_stub = "15_tb_rifampicin_MIC_large_minbin",       species = "M. tuberculosis", drug = "rifampicin",       pheno = "MIC doubling (≥10%)",       type = "ordinal")
)

# Strip label as a plotmath expression string (parsed by label_parsed): species
# in italics, then drug and binning, e.g. italic("S. pneumoniae")*", trimethoprim, MIC doubling (≥5%)".
panel_expr <- function(species, drug, pheno) {
  paste0('italic("', species, '")*", ', drug, ', ', pheno, '"')
}

SPLITS <- c(Random = "", LOSO = "_loso")

# -----------------------------------------------------------------------------
# Cheap JSON edge readers (avoid loading the multi-GB genotype arrays)
# -----------------------------------------------------------------------------

problems <- character(0)

# Parse a "key":[ints] array out of an already-read text chunk. Non-greedy to the
# first closing bracket, since phenotype arrays hold only integers and commas.
extract_int_array <- function(txt, key) {
  m <- regmatches(txt, regexpr(paste0('"', key, '":\\s*\\[[^]]*\\]'), txt, perl = TRUE))
  if (length(m) == 0) return(NULL)
  inner <- sub(paste0('.*?\\['), "", m)   # drop key + "[" prefix
  inner <- sub("\\]$", "", inner)
  as.integer(strsplit(inner, ",")[[1]])
}

# training_phenotype and test_phenotype live in the first ~1 KB; 512 KB is a safe
# prefix even for the largest TB set (~11.6k integers).
read_phenotypes <- function(json_path) {
  con <- file(json_path, "rb")
  on.exit(close(con))
  txt <- readChar(con, 524288L, useBytes = TRUE)
  list(train = extract_int_array(txt, "training_phenotype"),
       test  = extract_int_array(txt, "test_phenotype"))
}

# mic_breakpoints (inner breakpoints) sits in the last ~100 bytes.
read_breakpoints <- function(json_path) {
  sz  <- file.size(json_path)
  con <- file(json_path, "rb")
  on.exit(close(con))
  seek(con, where = max(0, sz - 1024L))
  tail_txt <- readChar(con, 1024L, useBytes = TRUE)
  m <- regmatches(tail_txt, regexpr('"mic_breakpoints":\\s*\\[[^]]*\\]', tail_txt, perl = TRUE))
  if (length(m) == 0) return(NULL)
  inner <- sub('.*?\\[', "", m)
  inner <- sub("\\]$", "", inner)
  as.numeric(strsplit(inner, ",")[[1]])
}

# -----------------------------------------------------------------------------
# Category labels
# -----------------------------------------------------------------------------

binary_labels <- c("Susceptible (0)", "Resistant (1)")

# Compact ordinal bin labels from inner breakpoints b (length K-1), no outer edges:
#   1 -> "<=b1", i -> "(b[i-1], b[i]]", K -> ">b[K-1]"
ordinal_labels <- function(b, K) {
  b <- signif(b, 3)
  if (K == 1) return("all")
  labs <- character(K)
  labs[1] <- paste0("≤", b[1])
  if (K > 2) for (i in 2:(K - 1)) labs[i] <- paste0("(", b[i - 1], ", ", b[i], "]")
  labs[K] <- paste0(">", b[K - 1])
  labs
}

# -----------------------------------------------------------------------------
# Build long data frame: one row per (dataset, split, category, set)
# -----------------------------------------------------------------------------

rows <- list()
panel_order <- character(0)

for (s in specs) {
  for (split in names(SPLITS)) {
    dir_name  <- paste0(s$dir_stub, SPLITS[[split]])
    json_path <- file.path(PRED_ROOT, dir_name, paste0(dir_name, ".json"))
    panel     <- panel_expr(s$species, s$drug, s$pheno)
    panel_order <- c(panel_order, panel)

    if (!file.exists(json_path)) {
      problems <- c(problems, paste0("missing JSON: ", json_path)); next
    }
    ph <- read_phenotypes(json_path)
    if (is.null(ph$train) || is.null(ph$test) ||
        anyNA(ph$train) || anyNA(ph$test)) {
      problems <- c(problems, paste0("unparseable phenotypes: ", json_path)); next
    }

    if (s$type == "binary") {
      levels_int <- c(0, 1)
      labs       <- binary_labels
      train_counts <- tabulate(factor(ph$train, levels = levels_int), nbins = 2)
      test_counts  <- tabulate(factor(ph$test,  levels = levels_int), nbins = 2)
    } else {
      K <- max(c(ph$train, ph$test))
      b <- read_breakpoints(json_path)
      if (is.null(b) || length(b) != K - 1) {
        problems <- c(problems, paste0("bad breakpoints (K=", K, "): ", json_path))
        labs <- as.character(seq_len(K))           # fall back to category index
      } else {
        labs <- ordinal_labels(b, K)
      }
      train_counts <- tabulate(ph$train, nbins = K)
      test_counts  <- tabulate(ph$test,  nbins = K)
    }

    n_cat <- length(labs)
    rows[[length(rows) + 1]] <- data.frame(
      species  = s$species,
      drug     = s$drug,
      pheno    = s$pheno,
      split    = split,
      panel    = panel,
      category = factor(rep(labs, 2), levels = labs),
      set      = factor(rep(c("Train", "Test"), each = n_cat), levels = c("Train", "Test")),
      count    = c(train_counts, test_counts),
      stringsAsFactors = FALSE)
  }
}

df <- do.call(rbind, rows)

if (length(problems) > 0) {
  warning(length(problems), " dataset issue(s):\n", paste(" ", problems, collapse = "\n"))
}

# Panel order follows the spec order (drug -> dataset number).
df$panel <- factor(df$panel, levels = unique(panel_order))

# -----------------------------------------------------------------------------
# Plot: one figure per split strategy
# -----------------------------------------------------------------------------

make_figure <- function(split_name) {
  d <- df[df$split == split_name, ]
  d$panel <- factor(d$panel, levels = intersect(levels(df$panel), unique(d$panel)))

  ggplot(d, aes(x = category, y = count, fill = set)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6,
             colour = "white", linewidth = 0.2) +
    geom_text(aes(label = ifelse(count > 0, count, "")),
              position = position_dodge(width = 0.7), vjust = -0.3, size = 2.4) +
    facet_wrap(~ panel, scales = "free", ncol = 3, labeller = label_parsed) +
    scale_fill_manual(values = SET_COLOURS) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(x = NULL, y = "Number of samples", fill = NULL) +
    theme_bw(base_size = 11) +
    theme(legend.position = "top",
          panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "grey92", colour = NA),
          axis.text.x = element_text(angle = 35, hjust = 1, size = 7))
}

# -----------------------------------------------------------------------------
# Save outputs
# -----------------------------------------------------------------------------

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

figures <- list(
  Random = list(file = "dataset_split_summary_random.png"),
  LOSO   = list(file = "dataset_split_summary_loso.png")
)

for (split_name in names(figures)) {
  spec     <- figures[[split_name]]
  png_path <- file.path(OUTPUT_DIR, spec$file)
  ggsave(png_path, make_figure(split_name),
         width = 14, height = 15, dpi = 300, bg = "white", limitsize = FALSE)
  message("Wrote ", png_path)
}

csv_path <- file.path(OUTPUT_DIR, "dataset_split_summary.csv")
write.csv(df, csv_path, row.names = FALSE)
message("Wrote ", csv_path)
