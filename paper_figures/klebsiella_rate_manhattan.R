#!/usr/bin/env Rscript
# Faceted RATE Manhattan across the simulated Klebsiella effect sizes.
#
# Companion to klebsiella_effectsize_faceted_manhattan.R, same design (one panel
# per effect size, causal variant rs_26645 in red and labelled, no h2 box
# and no title) but plotting cppRATE relative centrality
# instead of |median beta|.
#
# NO observed-effect reference line is drawn here, unlike the effect-size
# Manhattan. RATE is a unitless relative centrality that sums to 1 across
# variants; the carrier vs non-carrier phenotype difference is on the phenotype
# scale, so there is no meaningful horizontal line to draw.
#
# RATE SOURCE. Read from the pipeline's own cppRATE output rather than
# recomputed:
#   no-pruning arm : cppRATE_results/RATE_values.txt          (snp_id 0-based)
#   pruned arm     : cppRATE_results/RATE_values_depruned.txt (snp_id 1-based)
# Both carry one row per variant in variant-index order, verified by checking
# that the causal variant's row is the maximum in every run. The pruned arm's
# de-pruned file is written by depruning_and_write_rates(), so representatives'
# values are already mapped back onto their pruned partners.
#
# Usage:
#   mamba activate gwas_pipeline
#   Rscript paper_figures/klebsiella_rate_manhattan.R
#
# Outputs (to <RES>/gwas_klebsiella_homoplasic/figures/):
#   klebsiella_effectsize_faceted_rate_pruned.png
#   klebsiella_effectsize_faceted_rate_nopruning.png
#   klebsiella_effectsize_rate_summary.csv

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
})
grDevices::pdf(NULL)

RES       <- "/nfs/research/jlees/jacqueline/thesis_results"
DATASETS  <- file.path(RES, "gwas_datasets", "inference")
RUNS      <- file.path(RES, "gwas_klebsiella_homoplasic", "inference")
EXTRACTED <- file.path(RES, "gwas_klebsiella_homoplasic", "extracted")
OUT_DIR   <- file.path(RES, "gwas_klebsiella_homoplasic", "figures")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

COL_NULL     <- "#6b7280"
COL_CAUSAL   <- "#c1121f"
CAUSAL_LABEL <- "simulated causal variant"
EF_LEVELS    <- c("1.5", "2.5", "10", "30")

MODE_TITLES <- c(pruned = "with LD pruning", nopruning = "without LD pruning")

ef_facet <- function(x) paste0("effect size = ", x)

LOG_BREAKS <- scales::trans_breaks("log10", function(x) 10^x, n = 8)
log10_labels <- function(x) {
  ifelse(is.na(x), "", format(log10(x), trim = TRUE, drop0trailing = TRUE))
}

truth <- read.csv(file.path(DATASETS, "klebsiella_homoplasic_truth_summary.csv"),
                  stringsAsFactors = FALSE)
truth$ef_chr <- sub("\\.0$", "", format(truth$effect_size, trim = TRUE))
causal_name  <- unique(truth$causal_variant)
stopifnot(length(causal_name) == 1L)

extract_summary <- {
  f <- file.path(EXTRACTED, "klebsiella_extracted_summary.csv")
  if (file.exists(f)) data.table::fread(f) else NULL
}

# Base-R subsetting: data.table's i-expression would let the arguments shadow
# the identically named columns.
lookup <- function(ds_name, mode_name, column) {
  if (is.null(extract_summary)) return(NA)
  hit <- extract_summary$dataset == ds_name & extract_summary$mode == mode_name
  if (!any(hit)) return(NA)
  extract_summary[[column]][which(hit)[1]]
}

# Same h2 validity gate as the other figures. Worth knowing: cppRATE is
# scale-free, so a run that landed in a bad ADVI optimum can still rank the
# causal variant top on RATE even though its effect estimates collapsed to
# ~0 (EF 2.5 pruned does exactly that). The gate is kept anyway -- the
# underlying posterior is degenerate, so the values are not trustworthy just
# because the ranking happens to survive.
is_valid_run <- function(ds_name, mode_name) {
  if (is.null(extract_summary)) return(TRUE)
  v <- lookup(ds_name, mode_name, "valid")
  if (is.na(v)) return(FALSE)
  isTRUE(as.logical(v))
}

#' Read one cppRATE output. The two arms differ only in filename and in whether
#' snp_id is 0- or 1-based; both list every variant in variant-index order.
read_rate <- function(run_dir, mode, n_variants) {
  f <- file.path(run_dir, "cppRATE_results",
                 if (identical(mode, "pruned")) "RATE_values_depruned.txt"
                 else "RATE_values.txt")
  if (!file.exists(f)) return(NULL)
  d <- data.table::fread(f, skip = 3, header = FALSE, blank.lines.skip = TRUE,
                         col.names = c("snp_id", "rate", "kld"))
  if (nrow(d) != n_variants) {
    warning("RATE rows (", nrow(d), ") != variants (", n_variants, ") in ", f,
            call. = FALSE)
    return(NULL)
  }
  d
}

build_arm <- function(mode) {
  rows <- list(); excluded <- character(0); h2 <- list()
  for (i in seq_len(nrow(truth))) {
    ds       <- truth$dataset[i]
    nickname <- paste0(ds, "_continuous_", mode)
    vi       <- data.table::fread(file.path(DATASETS, ds,
                                            paste0(ds, "_variant_index.csv")))
    r <- if (is_valid_run(ds, mode)) {
      read_rate(file.path(RUNS, nickname), mode, nrow(vi))
    } else NULL

    if (is.null(r)) {
      excluded <- c(excluded, truth$ef_chr[i])
      message("[", mode, "] EF ", truth$ef_chr[i],
              " excluded: no RATE output, or fit did not pass the h2 gate")
      next
    }
    rows[[length(rows) + 1L]] <- data.table(
      ef_chr       = truth$ef_chr[i],
      variant_name = vi$variant_name,
      position     = vi$position,
      rate         = r$rate,
      is_causal    = vi$variant_name == causal_name)
    h2[[length(h2) + 1L]] <- data.table(
      ef_chr = truth$ef_chr[i], h2 = as.numeric(lookup(ds, mode, "h2_narrow")))
  }
  if (length(rows) == 0L) {
    warning("No usable RATE output for arm '", mode, "'", call. = FALSE)
    return(NULL)
  }
  df <- data.table::rbindlist(rows)
  df[, ef_fac := factor(ef_chr, levels = EF_LEVELS, labels = ef_facet(EF_LEVELS))]

  # A log axis cannot show RATE == 0; count rather than silently dropping.
  n_zero <- df[rate <= 0, .N]
  if (n_zero > 0) message("[", mode, "] dropping ", n_zero, " variants with RATE == 0")

  h2_df <- data.table::rbindlist(h2)
  h2_df[, ef_fac := factor(ef_chr, levels = EF_LEVELS, labels = ef_facet(EF_LEVELS))]
  h2_df <- h2_df[!is.na(h2)]

  list(df = df, plot_df = df[rate > 0], causal_df = df[rate > 0 & is_causal],
       excluded = excluded, h2_df = h2_df)
}

make_rate_manhattan <- function(plot_df, causal_df, excluded = character(0)) {
  excl_df <- NULL
  if (length(excluded) > 0) {
    y_mid <- mean(range(plot_df$rate))
    excl_df <- data.table(
      ef_fac   = factor(ef_facet(excluded), levels = ef_facet(EF_LEVELS)),
      position = mean(range(plot_df$position)),
      rate     = y_mid)
  }

  # RATE is a relative centrality that sums to 1 over the V variants, so 1/V is
  # the value every variant would take if importance were spread evenly -- the
  # natural null reference. V is taken from the data rather than hardcoded.
  V <- plot_df[, data.table::uniqueN(variant_name)]
  message("  reference line at 1/V = 1/", V, " = ", signif(1 / V, 4))

  p <- ggplot(plot_df[is_causal == FALSE], aes(x = position, y = rate)) +
    geom_point(colour = COL_NULL, alpha = 0.35, size = 0.7, stroke = 0) +
    geom_hline(yintercept = 1 / V, linetype = "dashed",
               colour = "grey45", linewidth = 0.4) +
    geom_point(data = causal_df, colour = "white", size = 3.6, stroke = 0) +
    geom_point(data = causal_df, colour = COL_CAUSAL, size = 2.6, stroke = 0) +
    geom_text(data = causal_df, label = CAUSAL_LABEL, colour = COL_CAUSAL,
              vjust = -1.1, size = 3.6, fontface = "bold", show.legend = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.22))) +
    scale_x_continuous(labels = scales::comma) +
    facet_wrap(~ ef_fac, nrow = 1) +
    labs(x = "variants, positioned by rs number", y = "RATE") +
    theme_minimal(base_size = 13) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(linewidth = 0.3, colour = "grey88"),
          strip.text = element_text(face = "bold"),
          axis.text.x = element_text(size = 8))

  if (!is.null(excl_df)) {
    p <- p +
      geom_blank(data = excl_df, aes(x = position, y = rate)) +
      geom_text(data = excl_df, aes(x = position, y = rate),
                label = "fit failed\n(no usable result)", colour = "grey35",
                size = 3.6, fontface = "italic", lineheight = 0.95)
  }
  p
}

arms <- list()
for (mode in names(MODE_TITLES)) {
  a <- build_arm(mode)
  if (is.null(a)) next
  arms[[mode]] <- a
  p <- make_rate_manhattan(a$plot_df, a$causal_df, a$excluded)
  png_path <- file.path(
    OUT_DIR, paste0("klebsiella_effectsize_faceted_rate_", mode, ".png"))
  ggsave(png_path, p, width = 15, height = 4.6, dpi = 600, bg = "white")
  message("Wrote ", png_path)
}
if (length(arms) == 0L) stop("No RATE output found under ", RUNS)

rate_summary <- data.table::rbindlist(lapply(names(arms), function(arm) {
  d <- arms[[arm]]$df
  d[, .(mode = arm,
        n_variants   = .N,
        causal_rate  = rate[is_causal],
        causal_rank  = rank(-rate, ties.method = "min")[is_causal],
        max_null_rate = max(rate[!is_causal]),
        rate_ratio   = rate[is_causal] / max(rate[!is_causal])), by = ef_chr]
}))
rate_summary[, ef_num := as.numeric(ef_chr)]
data.table::setorder(rate_summary, mode, ef_num)
rate_summary[, ef_num := NULL]

summary_csv <- file.path(OUT_DIR, "klebsiella_effectsize_rate_summary.csv")
data.table::fwrite(rate_summary, summary_csv)
message("Wrote ", summary_csv)
print(rate_summary)
