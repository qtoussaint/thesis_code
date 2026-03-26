############################################################
## utils.R
## Shared helper functions for GWAS dataset preparation
############################################################

library(dplyr)
library(data.table)
library(cmdstanr)
library(ggplot2)


# ============================================================
#  DATA LOADING
# ============================================================

#' Load SPN genotype matrix
#' @return matrix (variants x samples) with integer 0/1 values;
#'         rownames = variant IDs, colnames = sample IDs
load_spn_genotype <- function(path) {
  message("Loading SPN genotype from: ", path)
  geno <- read.delim(path, sep = "", header = FALSE, row.names = 1)
  colnames(geno) <- geno[1, ]
  geno <- geno[-1, ]
  colnames(geno) <- sub("^results_", "", colnames(geno))
  geno[] <- lapply(geno, function(x) as.integer(as.character(x)))
  message("  Loaded: ", nrow(geno), " variants x ", ncol(geno), " samples")
  geno
}


#' Load TB genotype presence/absence matrix and variant index
#' @return list(genotype = data.table N_samples x N_variants,
#'              variant_names = character vector)
load_tb_genotype <- function(presabs_path, variant_idx_path) {
  message("Loading TB genotype from: ", presabs_path)
  geno <- fread(presabs_path, sep = "\t", header = FALSE,
                na.strings = "", colClasses = "numeric")
  variant_names <- fread(variant_idx_path, sep = "\t", header = FALSE,
                         na.strings = "")[[1]]
  message("  Loaded: ", nrow(geno), " samples x ", ncol(geno), " variants")
  list(genotype = geno, variant_names = as.character(variant_names))
}


#' Load SPN poppipe lineages/sublineages
#' @return list(lineages = data.frame(Taxon, Strain),
#'              sublineages = data.frame(Taxon, Subcluster_1))
#'         NA subclusters filled with paste(coarse_cluster, "0000")
load_spn_lineages <- function(path) {
  message("Loading SPN lineages from: ", path)
  raw <- read.csv(path)
  lineages    <- raw[, c(1, 2)]  # Taxon, Strain
  sublineages <- raw[, c(1, 3)]  # Taxon, Subcluster_1
  na_idx <- which(is.na(sublineages[, 2]))
  if (length(na_idx) > 0) {
    sublineages[na_idx, 2] <- as.numeric(
      paste0(lineages[na_idx, 2], "0000")
    )
  }
  list(lineages = lineages, sublineages = sublineages)
}


#' Load TB fastlin lineage assignments
#' Filters out samples with multiple lineage designations (contamination)
#' and cleans parenthetical k-mer counts from lineage strings.
#' Derives coarse lineage (first element of dot-delimited string).
#' @return data.frame with columns: X.sample, lineages, lineages_coarse
load_tb_lineages <- function(path) {
  message("Loading TB lineages from: ", path)
  lin <- read.delim(path)
  n_before <- nrow(lin)
  lin <- lin %>%
    filter(!grepl(",", lineages)) %>%
    mutate(lineages = sub("\\s+\\(.*\\)$", "", lineages)) %>%
    filter(lineages != "")
  message("  Kept ", nrow(lin), " / ", n_before,
          " samples after filtering multi-lineage and empty designations")
  lin$lineages_coarse <- sapply(
    strsplit(lin$lineages, split = "[.]"),
    function(x) x[[1]]
  )
  lin
}


# ============================================================
#  PHENOTYPE PREPARATION
# ============================================================

#' Clean MIC string values by applying a named replacement map,
#' then coercing to numeric.
#' @param mic_vec character vector of raw MIC values
#' @param replacements named numeric vector: names are raw strings, values are numerics
#' @return numeric vector (NAs where conversion failed)
clean_mic_values <- function(mic_vec, replacements) {
  mic_vec <- as.character(mic_vec)
  for (raw in names(replacements)) {
    mic_vec[mic_vec == raw] <- as.character(replacements[[raw]])
  }
  as.numeric(mic_vec)
}


#' Convert numeric MIC to binary (0 = susceptible, 1 = resistant).
#' Samples with MIC strictly between s_max and r_min are excluded (returned as NA).
#' @param mic_numeric numeric vector
#' @param s_max  susceptible if MIC <= s_max
#' @param r_min  resistant  if MIC >= r_min (default: s_max + epsilon, i.e. no gap)
#' @return integer vector (0, 1, or NA)
mic_to_binary <- function(mic_numeric, s_max, r_min = NULL) {
  if (is.null(r_min)) r_min <- s_max + .Machine$double.eps
  result <- rep(NA_integer_, length(mic_numeric))
  result[!is.na(mic_numeric) & mic_numeric <= s_max] <- 0L
  result[!is.na(mic_numeric) & mic_numeric >= r_min] <- 1L
  n_excl <- sum(is.na(result) & !is.na(mic_numeric))
  if (n_excl > 0) {
    message("  mic_to_binary: excluded ", n_excl,
            " samples with MIC in intermediate zone (", s_max, ", ", r_min, ")")
  }
  result
}


#' Auto-bin MIC values using the standard microbiological doubling-dilution series.
#' Merges adjacent bins iteratively until all bins have >= min_bin_size samples.
#' Saves a two-panel before/after histogram PNG.
#'
#' @param mic_numeric  numeric vector of cleaned MIC values
#' @param min_bin_frac  minimum fraction of samples per bin (default 0.05 = 5%)
#' @param dilutions    candidate breakpoints (default: MIC_STANDARD_DILUTIONS from config)
#' @param hist_path    path to save histogram PNG (NULL = no save)
#' @param dataset_label label used in histogram title
#' @return list(
#'   bins            = integer vector of ordinal category (1..K),
#'   breakpoints     = final numeric breakpoints used,
#'   K               = number of bins,
#'   bin_counts      = table of counts per bin,
#'   mic_breakpoints = inner breakpoints (for PPOM cutpoints)
#' )
bin_mic_auto <- function(mic_numeric,
                         min_bin_frac  = 0.05,
                         dilutions     = MIC_STANDARD_DILUTIONS,
                         hist_path     = NULL,
                         dataset_label = "",
                         drug_name     = "",
                         species_name  = "",
                         strategy_label = "") {

  min_bin_size <- max(1L, round(min_bin_frac * length(mic_numeric)))

  # 1. Restrict breakpoints to span observed range
  lo <- min(mic_numeric, na.rm = TRUE)
  hi <- max(mic_numeric, na.rm = TRUE)
  breaks <- dilutions[dilutions < hi]
  breaks <- c(breaks[breaks <= lo][max(1, sum(breaks <= lo))],
              breaks[breaks > lo],
              max(dilutions[dilutions >= hi], hi * 1.01))
  breaks <- unique(sort(breaks))

  # Helper: assign bins from breakpoints
  assign_bins <- function(br) {
    as.integer(cut(mic_numeric, breaks = br, include.lowest = TRUE))
  }

  # 2. Initial binning
  bins_before <- assign_bins(breaks)
  counts_before <- tabulate(bins_before, nbins = length(breaks) - 1)
  breaks_before <- breaks

  # 3. Iteratively merge smallest bin with its smaller neighbour
  while (any(counts_before < min_bin_size, na.rm = TRUE) &&
         length(breaks) > 2) {
    counts_tmp <- tabulate(assign_bins(breaks), nbins = length(breaks) - 1)
    small_idx <- which.min(counts_tmp)
    # decide merge direction: merge with left or right neighbour (whichever is smaller)
    merge_left  <- if (small_idx > 1)                 counts_tmp[small_idx - 1] else Inf
    merge_right <- if (small_idx < length(counts_tmp)) counts_tmp[small_idx + 1] else Inf
    if (merge_left <= merge_right) {
      breaks <- breaks[-(small_idx)]      # remove left boundary of small bin
    } else {
      breaks <- breaks[-(small_idx + 1)]  # remove right boundary of small bin
    }
    counts_before <- tabulate(assign_bins(breaks), nbins = length(breaks) - 1)
  }

  bins_after  <- assign_bins(breaks)
  counts_after <- tabulate(bins_after, nbins = length(breaks) - 1)
  K <- length(breaks) - 1

  # 4. Print summary
  cat("\n--- bin_mic_auto:", dataset_label, "---\n")
  cat("Final K =", K, "bins\n")
  bin_labels <- paste0("(", breaks[-length(breaks)], ", ", breaks[-1], "]")
  print(data.frame(bin = 1:K, interval = bin_labels, count = counts_after))

  # 5. Save histogram PNG
  if (!is.null(hist_path)) {
    .save_bin_histogram(
      mic_numeric    = mic_numeric,
      bins_after     = bins_after,
      breaks_after   = breaks,
      dataset_label  = dataset_label,
      hist_path      = hist_path,
      drug_name      = drug_name,
      species_name   = species_name,
      strategy_label = strategy_label
    )
  }

  list(
    bins            = bins_after,
    breakpoints     = breaks,
    K               = K,
    bin_counts      = counts_after,
    mic_breakpoints = breaks[-c(1, length(breaks))]  # inner cuts only (for PPOM)
  )
}


#' Bin MIC values using equal-frequency binning constrained to the doubling-dilution grid.
#' Finds the subset of dilution breakpoints that minimises the coefficient of variation
#' (CV) of bin sizes, subject to every bin having >= min_bin_size samples.
#' Brute-force over all subsets of interior dilution breakpoints; feasible because
#' the dilution series spans at most ~10 interior points.
#'
#' @param mic_numeric  numeric vector of cleaned MIC values
#' @param min_bin_size minimum samples per bin (default 30)
#' @param dilutions    candidate breakpoints (default: MIC_STANDARD_DILUTIONS from config)
#' @param hist_path    path to save histogram PNG (NULL = no save)
#' @param dataset_label label used in histogram title
#' @return list(bins, breakpoints, K, bin_counts, mic_breakpoints) — same structure as bin_mic_auto
bin_mic_equalfreq <- function(mic_numeric,
                               min_bin_size  = 30,
                               dilutions     = MIC_STANDARD_DILUTIONS,
                               hist_path     = NULL,
                               dataset_label = "",
                               drug_name     = "",
                               species_name  = "",
                               strategy_label = "") {

  # 1. Same initial breakpoints as bin_mic_auto (restrict to observed range)
  lo <- min(mic_numeric, na.rm = TRUE)
  hi <- max(mic_numeric, na.rm = TRUE)
  breaks_full <- dilutions[dilutions < hi]
  breaks_full <- c(breaks_full[breaks_full <= lo][max(1, sum(breaks_full <= lo))],
                   breaks_full[breaks_full > lo],
                   max(dilutions[dilutions >= hi], hi * 1.01))
  breaks_full <- unique(sort(breaks_full))

  assign_bins <- function(br) as.integer(cut(mic_numeric, breaks = br, include.lowest = TRUE))

  bins_before  <- assign_bins(breaks_full)
  breaks_before <- breaks_full

  interior <- breaks_full[-c(1, length(breaks_full))]
  n_int    <- length(interior)

  best_cv     <- Inf
  best_breaks <- NULL

  # 2. Enumerate all non-empty subsets of interior breakpoints (K >= 2 bins)
  for (k in seq_len(n_int)) {
    combos <- combn(n_int, k, simplify = FALSE)
    for (idx in combos) {
      br_try  <- c(breaks_full[1], interior[idx], breaks_full[length(breaks_full)])
      br_try  <- unique(sort(br_try))
      counts  <- tabulate(assign_bins(br_try), nbins = length(br_try) - 1)
      if (any(counts < min_bin_size)) next
      cv_val  <- sd(counts) / mean(counts)
      if (cv_val < best_cv) {
        best_cv     <- cv_val
        best_breaks <- br_try
      }
    }
  }

  # 3. Fallback: if no valid subset meets min_bin_size, split at the middle interior point
  if (is.null(best_breaks)) {
    warning("bin_mic_equalfreq: no subset meets min_bin_size=", min_bin_size,
            " for ", dataset_label, "; falling back to midpoint split.")
    mid_idx     <- ceiling(n_int / 2)
    best_breaks <- c(breaks_full[1], interior[mid_idx], breaks_full[length(breaks_full)])
    best_cv     <- NA_real_
  }

  bins_after   <- assign_bins(best_breaks)
  counts_after <- tabulate(bins_after, nbins = length(best_breaks) - 1)
  K <- length(best_breaks) - 1

  cat("\n--- bin_mic_equalfreq:", dataset_label, "---\n")
  cat("Final K =", K, "bins  (CV =", round(best_cv, 3), ")\n")
  bin_labels <- paste0("(", best_breaks[-length(best_breaks)], ", ", best_breaks[-1], "]")
  print(data.frame(bin = 1:K, interval = bin_labels, count = counts_after))

  if (!is.null(hist_path)) {
    .save_bin_histogram(
      mic_numeric    = mic_numeric,
      bins_after     = bins_after,
      breaks_after   = best_breaks,
      dataset_label  = paste(dataset_label, "(equal-freq)"),
      hist_path      = hist_path,
      drug_name      = drug_name,
      species_name   = species_name,
      strategy_label = strategy_label
    )
  }

  list(
    bins            = bins_after,
    breakpoints     = best_breaks,
    K               = K,
    bin_counts      = counts_after,
    mic_breakpoints = best_breaks[-c(1, length(best_breaks))]
  )
}


#' Bin MIC values by placing breakpoints at valleys between natural peaks in the MIC
#' distribution, then applying a min-size merge pass identical to bin_mic_auto.
#' Valley bins are identified from a moving-average-smoothed count histogram.
#'
#' @param mic_numeric   numeric vector of cleaned MIC values
#' @param min_bin_size  minimum samples per bin (default 30)
#' @param dilutions     candidate breakpoints (default: MIC_STANDARD_DILUTIONS from config)
#' @param smooth_span   moving-average window for count smoothing (default 3)
#' @param hist_path     path to save histogram PNG (NULL = no save)
#' @param dataset_label label used in histogram title
#' @return list(bins, breakpoints, K, bin_counts, mic_breakpoints) — same structure as bin_mic_auto
bin_mic_peaks <- function(mic_numeric,
                           min_bin_size  = 30,
                           dilutions     = MIC_STANDARD_DILUTIONS,
                           smooth_span   = 3,
                           hist_path     = NULL,
                           dataset_label = "",
                           drug_name     = "",
                           species_name  = "",
                           strategy_label = "") {

  # 1. Same initial breakpoints as bin_mic_auto
  lo <- min(mic_numeric, na.rm = TRUE)
  hi <- max(mic_numeric, na.rm = TRUE)
  breaks_full <- dilutions[dilutions < hi]
  breaks_full <- c(breaks_full[breaks_full <= lo][max(1, sum(breaks_full <= lo))],
                   breaks_full[breaks_full > lo],
                   max(dilutions[dilutions >= hi], hi * 1.01))
  breaks_full <- unique(sort(breaks_full))

  assign_bins <- function(br) as.integer(cut(mic_numeric, breaks = br, include.lowest = TRUE))

  bins_before   <- assign_bins(breaks_full)
  breaks_before <- breaks_full
  counts_full   <- tabulate(bins_before, nbins = length(breaks_full) - 1)
  n_bins_full   <- length(counts_full)

  # 2. Smooth counts with a moving average; replace edge NAs with original counts
  smooth_raw <- as.numeric(stats::filter(counts_full, rep(1/smooth_span, smooth_span), sides = 2))
  smooth_counts <- ifelse(is.na(smooth_raw), counts_full, smooth_raw)

  # 3. Find interior bins that are local minima (valleys between peaks)
  valley_idx <- integer(0)
  if (n_bins_full >= 3) {
    interior_idx <- seq(2, n_bins_full - 1)
    valley_idx   <- interior_idx[
      smooth_counts[interior_idx] < smooth_counts[interior_idx - 1] &
      smooth_counts[interior_idx] < smooth_counts[interior_idx + 1]
    ]
  }

  # 4. Keep only breakpoints adjacent to valley bins (left and right edges of each valley)
  #    Bin i spans (breaks_full[i], breaks_full[i+1]], so its edges are
  #    breaks_full[valley_idx] (left) and breaks_full[valley_idx + 1] (right).
  if (length(valley_idx) > 0) {
    near_valley <- sort(unique(c(breaks_full[valley_idx], breaks_full[valley_idx + 1])))
    breaks      <- unique(sort(c(breaks_full[1], near_valley, breaks_full[length(breaks_full)])))
  } else {
    # No valleys detected (unimodal / monotone): fall back to bin_mic_auto behaviour
    message("  bin_mic_peaks: no valleys detected for ", dataset_label,
            "; falling back to min-size merge.")
    breaks <- breaks_full
  }

  # 5. Min-size merge pass (identical to bin_mic_auto)
  counts_cur <- tabulate(assign_bins(breaks), nbins = length(breaks) - 1)
  while (any(counts_cur < min_bin_size, na.rm = TRUE) && length(breaks) > 2) {
    counts_tmp  <- tabulate(assign_bins(breaks), nbins = length(breaks) - 1)
    small_idx   <- which.min(counts_tmp)
    merge_left  <- if (small_idx > 1)                   counts_tmp[small_idx - 1] else Inf
    merge_right <- if (small_idx < length(counts_tmp))  counts_tmp[small_idx + 1] else Inf
    if (merge_left <= merge_right) {
      breaks <- breaks[-(small_idx)]
    } else {
      breaks <- breaks[-(small_idx + 1)]
    }
    counts_cur <- tabulate(assign_bins(breaks), nbins = length(breaks) - 1)
  }

  bins_after   <- assign_bins(breaks)
  counts_after <- tabulate(bins_after, nbins = length(breaks) - 1)
  K <- length(breaks) - 1

  cat("\n--- bin_mic_peaks:", dataset_label, "---\n")
  cat("Final K =", K, "bins  (", length(valley_idx), "valleys detected)\n")
  bin_labels <- paste0("(", breaks[-length(breaks)], ", ", breaks[-1], "]")
  print(data.frame(bin = 1:K, interval = bin_labels, count = counts_after))

  if (!is.null(hist_path)) {
    .save_bin_histogram(
      mic_numeric    = mic_numeric,
      bins_after     = bins_after,
      breaks_after   = breaks,
      dataset_label  = paste(dataset_label, "(peaks)"),
      hist_path      = hist_path,
      drug_name      = drug_name,
      species_name   = species_name,
      strategy_label = strategy_label
    )
  }

  list(
    bins            = bins_after,
    breakpoints     = breaks,
    K               = K,
    bin_counts      = counts_after,
    mic_breakpoints = breaks[-c(1, length(breaks))]
  )
}


# Internal helper: save three-panel histogram
# Panel 0 (left):   raw MIC histogram with vertical lines at bin boundaries
# Panel 1 (middle): log2(MIC) histogram with vertical lines at bin boundaries
# Panel 2 (right):  categorical bar chart of sample counts per final bin
.save_bin_histogram <- function(mic_numeric, bins_after, breaks_after,
                                dataset_label, hist_path,
                                drug_name = "", species_name = "",
                                strategy_label = "") {
  K            <- length(breaks_after) - 1
  counts_after <- tabulate(bins_after, nbins = K)
  log2_mic     <- log2(mic_numeric)
  inner_cuts   <- log2(breaks_after[-c(1, length(breaks_after))])
  inner_cuts_raw <- breaks_after[-c(1, length(breaks_after))]
  bin_labels   <- paste0("(", round(breaks_after[-length(breaks_after)], 4),
                         ", ", round(breaks_after[-1], 4), "]")

  # Panel 0: raw (non-log) MIC distribution with breakpoint lines
  df_raw_mic <- data.frame(mic = mic_numeric)
  sorted_mics <- sort(unique(mic_numeric))
  mic_bw <- if (length(sorted_mics) > 1) diff(range(sorted_mics)) / 30 else 0.01
  p0 <- ggplot(df_raw_mic, aes(x = mic)) +
    geom_histogram(binwidth = mic_bw, fill = "steelblue", colour = "white", linewidth = 0.2) +
    geom_vline(xintercept = inner_cuts_raw,
               colour = "red", linetype = "dashed", linewidth = 0.7) +
    labs(
      x = expression(MIC~(mu*g%.%mL^{-1})), y = "Number of samples"
    ) +
    theme_bw(base_size = 12)

  # Panel 1: log2(MIC) distribution with breakpoint lines
  # Bars use binwidth=0.25, boundary=0; breakpoint lines at inner_cuts+0.5 fall
  # cleanly between bars.
  df_raw <- data.frame(log2_mic = log2_mic)
  p1 <- ggplot(df_raw, aes(x = log2_mic)) +
    geom_histogram(binwidth = 0.25, boundary = 0, fill = "steelblue", colour = "white", linewidth = 0.2) +
    geom_vline(xintercept = inner_cuts + 0.5, colour = "red", linetype = "dashed", linewidth = 0.7) +
    labs(
      x = expression(log[2](MIC)), y = "Number of samples"
    ) +
    theme_bw(base_size = 12)

  # Panel 2: categorical binned bar chart (equal spacing, no log2 x axis)
  df_bins <- data.frame(
    bin   = factor(seq_len(K), levels = seq_len(K), labels = bin_labels),
    count = counts_after
  )
  p2 <- ggplot(df_bins, aes(x = bin, y = count, fill = bin)) +
    geom_col(colour = "white", linewidth = 0.2) +
    labs(
      x = "Binned MICs", y = "Number of samples"
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 35, hjust = 1, size = 8))

  # Build overall title with italic species name
  overall_title <- if (nzchar(drug_name) && nzchar(species_name)) {
    grid::textGrob(
      bquote("Distribution of ordered categories using " * .(strategy_label) * ", "
             * .(drug_name) * " resistance in " * italic(.(species_name))),
      gp = grid::gpar(fontsize = 14)
    )
  } else {
    grid::textGrob(dataset_label, gp = grid::gpar(fontsize = 14))
  }

  dir.create(dirname(hist_path), showWarnings = FALSE, recursive = TRUE)
  ggsave(hist_path, gridExtra::arrangeGrob(p0, p1, p2, ncol = 3,
                                           top = overall_title),
         width = 15, height = 5.5, dpi = 150)
  message("  Saved MIC bin histogram to: ", hist_path)
}


#' Save a two-panel phenotype distribution plot for a binary (S/R) dataset.
#'
#' Panel 1 (if mic_numeric provided): plain log2(MIC) histogram of all samples.
#'   If s_max is given, draws a red dashed breakpoint line shifted right by 0.5 in
#'   log2 space to sit cleanly between histogram bars.
#' Panel 2: bar chart of susceptible (0) vs resistant (1) counts.
#'
#' @param mic_numeric  numeric MIC vector (same length as binary_vec), or NULL to
#'                     omit Panel 1.
#' @param binary_vec   integer/numeric vector of 0/1 phenotype values.
#' @param dataset_label  character label used in plot titles.
#' @param hist_path    file path for the output PNG.
#' @param s_max        susceptible threshold (MIC <= s_max => 0); draws a threshold
#'                     line when provided.
#' @param r_min        resistant threshold (optional, for annotation only).
save_binary_histogram <- function(mic_numeric = NULL, binary_vec,
                                  dataset_label, hist_path,
                                  s_max = NULL, r_min = NULL,
                                  drug_name = "", species_name = "",
                                  strategy_label = "") {
  binary_fac <- factor(binary_vec, levels = c(0, 1),
                       labels = c("Susceptible (0)", "Resistant (1)"))
  counts     <- table(binary_fac)

  panels <- list()

  if (!is.null(mic_numeric)) {
    # Panel 0 (leftmost): raw (non-log) MIC histogram with breakpoint line
    df_raw_mic <- data.frame(mic = mic_numeric)
    sorted_mics <- sort(unique(mic_numeric))
    mic_bw <- if (length(sorted_mics) > 1) diff(range(sorted_mics)) / 30 else 0.01
    p0 <- ggplot(df_raw_mic, aes(x = mic)) +
      geom_histogram(binwidth = mic_bw, fill = "steelblue", colour = "white", linewidth = 0.2) +
      labs(
        x = expression(MIC~(mu*g%.%mL^{-1})), y = "Number of samples"
      ) +
      theme_bw(base_size = 12)
    if (!is.null(s_max)) {
      p0 <- p0 + geom_vline(xintercept = s_max, colour = "red",
                             linetype = "dashed", linewidth = 0.7)
    }
    panels[["p0"]] <- p0

    # Panel 1 (middle): log2(MIC) histogram with breakpoint line
    df_raw <- data.frame(log2_mic = log2(mic_numeric))
    p1 <- ggplot(df_raw, aes(x = log2_mic)) +
      geom_histogram(binwidth = 0.25, boundary = 0, fill = "steelblue",
                     colour = "white", linewidth = 0.2) +
      labs(
        x = expression(log[2](MIC)), y = "Number of samples"
      ) +
      theme_bw(base_size = 12)
    if (!is.null(s_max)) {
      threshold_cut <- log2(s_max) + 0.5
      p1 <- p1 +
        geom_vline(xintercept = threshold_cut, colour = "red",
                   linetype = "dashed", linewidth = 0.7)
    }
    panels[["p1"]] <- p1
  }

  # Panel 2 (rightmost): bar chart of S/R counts
  df_counts <- data.frame(
    class = names(counts),
    count = as.integer(counts)
  )
  p2 <- ggplot(df_counts, aes(x = class, y = count,
                               fill = class)) +
    geom_col(colour = "white", linewidth = 0.2) +
    geom_text(aes(label = count), vjust = -0.3, size = 3.5) +
    scale_fill_manual(values = c("Susceptible (0)" = "steelblue",
                                 "Resistant (1)"   = "firebrick")) +
    labs(
      x = NULL, y = "Number of samples"
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "none")
  panels[["p2"]] <- p2

  # Build overall title with italic species name
  overall_title <- if (nzchar(drug_name) && nzchar(species_name)) {
    grid::textGrob(
      bquote("Distribution of binary categories using " * .(strategy_label) * ", "
             * .(drug_name) * " resistance in " * italic(.(species_name))),
      gp = grid::gpar(fontsize = 14)
    )
  } else {
    grid::textGrob(dataset_label, gp = grid::gpar(fontsize = 14))
  }

  grob <- gridExtra::arrangeGrob(grobs = panels, ncol = length(panels),
                                 top = overall_title)

  dir.create(dirname(hist_path), showWarnings = FALSE, recursive = TRUE)
  ggsave(hist_path, grob,
         width = 5 * length(panels), height = 5.5, dpi = 150)
  message("  Saved binary histogram to: ", hist_path)
}


# ============================================================
#  SAMPLE ALIGNMENT
# ============================================================

#' Intersect samples present in all four data sources and reorder
#' all objects to the same sample order (defined by phenotype order).
#'
#' @param pheno_df     data.frame with a column named by id_col
#' @param geno         matrix or data.table (variants x samples, colnames = IDs)
#'                     OR (samples x variants, rownames = IDs) — see transpose_geno
#' @param lineages_df  data.frame with column 1 = sample ID
#' @param sublin_df    data.frame with column 1 = sample ID (NULL for TB)
#' @param id_col       name of ID column in pheno_df
#' @param geno_in_cols TRUE if sample IDs are colnames of geno (SPN orientation)
#' @return list(pheno, geno_matrix [samples x variants], lineages, sublineages)
intersect_and_align <- function(pheno_df, geno, lineages_df, sublin_df = NULL,
                                id_col = "ID", geno_in_cols = TRUE) {

  pheno_ids  <- pheno_df[[id_col]]
  lin_ids    <- lineages_df[[1]]
  sublin_ids <- if (!is.null(sublin_df)) sublin_df[[1]] else lin_ids

  if (geno_in_cols) {
    geno_ids <- colnames(geno)
  } else {
    geno_ids <- rownames(geno)
    if (is.null(geno_ids)) geno_ids <- as.character(seq_len(nrow(geno)))
  }

  id_lists <- list(pheno_ids, geno_ids, lin_ids, sublin_ids)
  common   <- Reduce(intersect, id_lists)
  message("  Common samples across all sources: ", length(common))

  ref_order <- pheno_ids[pheno_ids %in% common]

  # Align phenotype
  pheno_out <- pheno_df[pheno_df[[id_col]] %in% common, ]
  pheno_out <- pheno_out[match(ref_order, pheno_out[[id_col]]), ]

  # Align genotype -> samples x variants
  if (geno_in_cols) {
    geno_mat <- t(as.matrix(geno))
    geno_mat <- geno_mat[match(ref_order, rownames(geno_mat)), , drop = FALSE]
  } else {
    geno_mat <- as.matrix(geno)[match(ref_order, geno_ids), , drop = FALSE]
  }
  storage.mode(geno_mat) <- "integer"

  # Align lineages
  lin_out <- lineages_df[lineages_df[[1]] %in% common, ]
  lin_out <- lin_out[match(ref_order, lin_out[[1]]), ]

  # Align sublineages
  sublin_out <- if (!is.null(sublin_df)) {
    s <- sublin_df[sublin_df[[1]] %in% common, ]
    s[match(ref_order, s[[1]]), ]
  } else NULL

  list(pheno     = pheno_out,
       geno_mat  = geno_mat,
       lineages  = lin_out,
       sublineages = sublin_out,
       sample_ids  = ref_order)
}


# ============================================================
#  LINEAGE ENCODING
# ============================================================

#' Choose the reference (dropped) sublineage: the one with the
#' smallest mean phenotype value (lowest resistance / lowest MIC).
#' @param sublin_ids  character/numeric vector of sublineage ID per sample
#' @param pheno_vec   numeric phenotype vector (binary 0/1 or continuous)
#' @return the sublineage ID to use as reference level
select_reference_sublineage <- function(sublin_ids, pheno_vec) {
  means <- tapply(pheno_vec, sublin_ids, mean, na.rm = TRUE)
  ref   <- names(which.min(means))
  message("  Reference (dropped) sublineage: ", ref,
          "  (mean phenotype = ", round(min(means, na.rm = TRUE), 4), ")")
  ref
}


#' Encode SPN lineages and sublineages as one-hot matrices,
#' with the reference sublineage set to the least-resistant cluster.
#' Also computes the parent_lineage mapping vector.
#'
#' @param lineages_df   data.frame(Taxon, Strain)
#' @param sublin_df     data.frame(Taxon, Subcluster_1)
#' @param pheno_vec     numeric phenotype (for reference selection)
#' @return list(lineage_matrix, sublineage_matrix, parent_lineage)
encode_lineages_spn <- function(lineages_df, sublin_df, pheno_vec) {
  ref_sublin <- select_reference_sublineage(sublin_df[[2]], pheno_vec)
  sublin_fac <- relevel(factor(sublin_df[[2]]), ref = as.character(ref_sublin))

  lin_mat   <- model.matrix(~ as.factor(lineages_df[[2]]) - 1)
  sublin_mat <- model.matrix(~ sublin_fac - 1)

  # Shorten column names for readability
  colnames(lin_mat)    <- levels(factor(lineages_df[[2]]))
  colnames(sublin_mat) <- levels(sublin_fac)

  # parent_lineage[k] = index of the majority lineage for sublineage k
  parent_lineage <- integer(ncol(sublin_mat))
  for (k in seq_len(ncol(sublin_mat))) {
    sub_idx <- which(sublin_mat[, k] == 1)
    lin_counts <- colSums(lin_mat[sub_idx, , drop = FALSE])
    parent_lineage[k] <- which.max(lin_counts)
  }

  list(lineage_matrix    = lin_mat,
       sublineage_matrix = sublin_mat,
       parent_lineage    = parent_lineage)
}


#' Encode TB lineages and sublineages as one-hot matrices,
#' with the reference sublineage set to the least-resistant cluster.
#'
#' @param lineages_df  data.frame with columns: X.sample, lineages, lineages_coarse
#' @param pheno_vec    numeric phenotype
#' @return list(lineage_matrix, sublineage_matrix, parent_lineage)
encode_lineages_tb <- function(lineages_df, pheno_vec) {
  ref_sublin <- select_reference_sublineage(lineages_df$lineages, pheno_vec)
  sublin_fac <- relevel(factor(lineages_df$lineages), ref = ref_sublin)
  coarse_fac <- factor(lineages_df$lineages_coarse)

  lin_mat    <- model.matrix(~ coarse_fac - 1)
  sublin_mat <- model.matrix(~ sublin_fac - 1)

  colnames(lin_mat)    <- levels(coarse_fac)
  colnames(sublin_mat) <- levels(sublin_fac)

  parent_lineage <- integer(ncol(sublin_mat))
  for (k in seq_len(ncol(sublin_mat))) {
    sub_idx <- which(sublin_mat[, k] == 1)
    lin_counts <- colSums(lin_mat[sub_idx, , drop = FALSE])
    parent_lineage[k] <- which.max(lin_counts)
  }

  list(lineage_matrix    = lin_mat,
       sublineage_matrix = sublin_mat,
       parent_lineage    = parent_lineage)
}


# ============================================================
#  STAN LIST BUILDERS
# ============================================================

#' Build Stan data list for inference
#' @param pheno        numeric/integer phenotype vector (length N)
#' @param geno_mat     integer matrix (N x V)
#' @param lin_mat      one-hot lineage matrix (N x L)
#' @param sublin_mat   one-hot sublineage matrix (N x S)
#' @param parent_lin   integer vector (length S)
#' @param K            number of ordered categories (NULL for binary/continuous)
#' @param mic_bkpts    inner MIC breakpoints for PPOM (NULL if not needed)
build_stan_inference <- function(pheno, geno_mat, lin_mat, sublin_mat,
                                 parent_lin, K = NULL, mic_bkpts = NULL) {
  N <- nrow(geno_mat)
  V <- ncol(geno_mat)
  L <- ncol(lin_mat)
  S <- ncol(sublin_mat)

  d <- list(
    phenotype         = pheno,
    variant_matrix    = geno_mat,
    lineage_matrix    = lin_mat,
    sublineage_matrix = sublin_mat,
    parent_lineage    = parent_lin,
    N = N, V = V, L = L, S = S
  )
  if (!is.null(K))        d$K               <- as.integer(K)
  if (!is.null(mic_bkpts)) d$mic_breakpoints <- mic_bkpts
  d
}


#' Build Stan data list for prediction (80/20 train/test split)
#' @param pheno      numeric/integer phenotype vector (full dataset)
#' @param geno_mat   integer matrix N x V (full dataset)
#' @param lin_mat    one-hot lineage matrix (full)
#' @param sublin_mat one-hot sublineage matrix (full)
#' @param parent_lin integer vector (length S)
#' @param sample_ids character vector of sample IDs (full dataset)
#' @param K          number of ordered categories (NULL for binary/continuous)
#' @param mic_bkpts  inner MIC breakpoints for PPOM (NULL if not needed)
#' @param train_prop proportion for training set (default 0.8)
#' @param seed       random seed
#' @return list(stan_list, train_ids, test_ids)
build_stan_prediction <- function(pheno, geno_mat, lin_mat, sublin_mat,
                                  parent_lin, sample_ids,
                                  K = NULL, mic_bkpts = NULL,
                                  train_prop = PRED_TRAIN_PROP,
                                  seed = PRED_SEED) {
  N <- nrow(geno_mat)
  set.seed(seed)
  train_idx <- sort(sample(N, floor(N * train_prop)))
  test_idx  <- setdiff(seq_len(N), train_idx)

  N_train <- length(train_idx)
  N_test  <- length(test_idx)
  V <- ncol(geno_mat)
  S <- ncol(sublin_mat)

  d <- list(
    training_phenotype  = pheno[train_idx],
    test_phenotype      = pheno[test_idx],
    training_variants   = geno_mat[train_idx, , drop = FALSE],
    test_variants       = geno_mat[test_idx,  , drop = FALSE],
    lineage_matrix      = lin_mat,          # full — pipeline subsets internally
    sublineage_matrix   = sublin_mat,
    training_lineages   = lin_mat[train_idx, , drop = FALSE],
    test_lineages       = lin_mat[test_idx,  , drop = FALSE],
    training_sublineages = sublin_mat[train_idx, , drop = FALSE],
    test_sublineages     = sublin_mat[test_idx,  , drop = FALSE],
    parent_lineage      = parent_lin,
    N_train = N_train, N_test = N_test, V = V, S = S,
    L = ncol(lin_mat)
  )
  if (!is.null(K))        d$K               <- as.integer(K)
  if (!is.null(mic_bkpts)) d$mic_breakpoints <- mic_bkpts

  list(
    stan_list  = d,
    train_ids  = sample_ids[train_idx],
    test_ids   = sample_ids[test_idx]
  )
}


# ============================================================
#  OUTPUT WRITING
# ============================================================

#' Write all output files for one dataset.
#' Creates the output directory if needed.
#'
#' @param stan_list     named list ready for cmdstanr::write_stan_json
#' @param sample_ids    character vector of sample IDs (in same order as phenotype)
#' @param variant_names character vector of variant IDs
#' @param parent_lin    integer vector
#' @param outdir        path to output directory (created if needed)
#' @param dataset_name  prefix used for file names
#' @param test_ids      character vector of test sample IDs (prediction only, NULL otherwise)
#' @param test_pheno    numeric vector of test phenotypes (prediction only, NULL otherwise)
write_dataset <- function(stan_list, sample_ids, variant_names, parent_lin,
                          outdir, dataset_name,
                          test_ids = NULL, test_pheno = NULL) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  # JSON
  json_path <- file.path(outdir, paste0(dataset_name, ".json"))
  write_stan_json(data = stan_list, file = json_path)
  message("  Wrote JSON: ", json_path)

  # Variant index
  vnames <- as.character(variant_names)
  vnames_clean <- sub("^Chromosome_", "", vnames)
  pos_str <- sub("_[^_]+$", "", vnames_clean)
  pos_str <- sub("_[^_]+$", "", pos_str)
  variant_pos <- suppressWarnings(as.numeric(pos_str))
  vi <- data.frame(variant_name = vnames, position = variant_pos,
                   stringsAsFactors = FALSE)
  write.csv(vi, file.path(outdir, paste0(dataset_name, "_variant_index.csv")),
            row.names = FALSE)

  # Sample index
  write.csv(data.frame(row_index = seq_along(sample_ids), sample_id = sample_ids),
            file.path(outdir, paste0(dataset_name, "_sample_index.csv")),
            row.names = FALSE)

  # Parent lineages
  write.csv(data.frame(sublineage_index = seq_along(parent_lin),
                       parent_lineage_index = parent_lin),
            file.path(outdir, paste0(dataset_name, "_parent_lineages.csv")),
            row.names = FALSE)

  # Test phenotypes (prediction only)
  if (!is.null(test_ids) && !is.null(test_pheno)) {
    write.csv(data.frame(sample_id = test_ids, true_phenotype = test_pheno),
              file.path(outdir, paste0(dataset_name, "_test_phenotypes.csv")),
              row.names = FALSE)
  }
}


#' Write inputs manifest CSV recording all input files and their metadata.
#' @param paths    named character vector of input file paths
#' @param out_base base output directory
write_inputs_manifest <- function(paths, out_base) {
  rows <- lapply(names(paths), function(nm) {
    p  <- paths[[nm]]
    fi <- tryCatch(file.info(p), error = function(e) NULL)
    data.frame(
      input_name = nm,
      path       = p,
      size_bytes = if (!is.null(fi)) fi$size  else NA,
      mtime      = if (!is.null(fi)) as.character(fi$mtime) else NA,
      stringsAsFactors = FALSE
    )
  })
  manifest <- do.call(rbind, rows)
  manifest$created_at <- as.character(Sys.time())
  out_path <- file.path(out_base, "inputs_manifest.csv")
  write.csv(manifest, out_path, row.names = FALSE)
  message("Wrote inputs manifest to: ", out_path)
}
