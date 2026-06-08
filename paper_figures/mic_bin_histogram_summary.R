#!/usr/bin/env Rscript
# Compile the per-dataset MIC bin histogram inference figures into one image per dataset.
#
# Each dataset's individual PNGs are stacked vertically in filename order. The panels are
# reused exactly as-is: every panel is padded to the common width on a white background
# (centred) so widths match for stacking, but no panel is rescaled or relabelled.
#
# Usage:
#   Rscript mic_bin_histogram_summary.R

suppressPackageStartupMessages({
  library(magick)
})

INPUT_DIR  <- "/nfs/research/jlees/jacqueline/thesis_results/gwas_datasets/MIC_bin_histograms"
OUTPUT_DIR <- "/nfs/research/jlees/jacqueline/thesis_results/paper_figures/dataset_binning"

# dataset key -> regex selecting that dataset's PNGs. Penicillin excludes the unitigs
# variants (prefixes 17 and 18).
datasets <- list(
  spn_penicillin   = function(f) grepl("spn_penicillin", f) & !grepl("^(17|18)_", f),
  spn_trimethoprim = function(f) grepl("spn_trimethoprim", f),
  tb_rifampicin    = function(f) grepl("tb_rifampicin", f)
)

# Keep only the inference distribution / bin histograms; drop prediction panels.
all_pngs <- sort(list.files(INPUT_DIR, pattern = "\\.png$"))
all_pngs <- all_pngs[!grepl("pred", all_pngs)]

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

for (key in names(datasets)) {
  files <- all_pngs[datasets[[key]](all_pngs)]
  if (length(files) == 0) stop("No PNGs matched dataset: ", key)

  message(key, ": ", length(files), " panels")

  imgs  <- image_read(file.path(INPUT_DIR, files))
  info  <- image_info(imgs)
  max_w <- max(info$width)

  # Pad each panel to the common width on white, keeping native pixels, then stack.
  padded  <- image_extent(imgs,
                          geometry = geometry_size_pixels(width = max_w, height = max(info$height)),
                          gravity = "center", color = "white")
  stacked <- image_append(padded, stack = TRUE)

  out_path <- file.path(OUTPUT_DIR, paste0("mic_bin_histogram_summary_", key, ".png"))
  image_write(stacked, out_path)
  message("Wrote ", out_path)

  # Free magick bitmaps before the next dataset to keep peak memory down.
  image_destroy(imgs); image_destroy(padded); image_destroy(stacked)
  gc()
}
