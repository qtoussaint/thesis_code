#!/usr/bin/env Rscript
# Reconstruct per-variant median + 89% CI significance from the posterior draws
# a run has ALREADY written. Two reasons to use this rather than the pipeline's
# own fitted_model/depruned_variant_effects.csv:
#
#  1. It is available as soon as the fit finishes, whereas the pipeline writes
#     its CSV only after cppRATE.
#  2. It reports the UNSTANDARDIZED effect. The pipeline's CSV does not:
#     cpprate.R:304 selects draws with grepl("beta_variant", colnames), which
#     also matches beta_variant_std and beta_variant_std_prior. beta_variant_std
#     is a transformed parameter and so is ordered ahead of the beta_variant
#     generated quantity, so the first V columns -- the ones that survive R's
#     recycling in total_med[-pruned] <- medians_kept, with a "number of items to
#     replace" warning -- are the STANDARDIZED effects. Confirmed empirically:
#     canonical/from_draws is positive for every variant and spans 0.0093 to
#     0.5000, exactly the range of sd(x_v) for a binary column.
#
# cppRATE_matrices/coefficients.csv is written by write_coefficient_matrix()
# using fitted_model$draws(variables = "beta_variant"), an exact match, so it
# holds beta_variant and nothing else. Reading it also avoids the 21 GB RDS.
# Mirrors paper_figures/regen_pom_variant_effects.R, but reduces in a single read
# rather than repeated column-chunk passes (one 5.4 GB parse instead of ~17).
#
# Runs fitted with --ld_pruning false are handled too: they write no pruning
# summary, their draws already cover every variant, and the de-pruning step is
# skipped. That is how the no-pruning arm of compare_ld_pruning gets its betas,
# since that branch of the pipeline writes no variant effects CSV of its own.
#
# Usage:
#   Rscript regen_effects_from_draws.R --run-dir <inference dir> [--out <path>]

suppressPackageStartupMessages({
  library(data.table)
})

GWAS_WORKFLOW_R <- "/nfs/research/jlees/jacqueline/gwas_workflow/code/gwas_workflow/R"
# Both files define only functions (no top-level code), so sourcing is safe and
# reuses the pipeline's de-pruning logic verbatim.
source(file.path(GWAS_WORKFLOW_R, "ld_pruning.R"))
source(file.path(GWAS_WORKFLOW_R, "cpprate.R"))

parse_args <- function(argv) {
  out <- list(run_dir = NULL, out = NULL)
  i <- 1
  while (i <= length(argv)) {
    a <- argv[i]
    if (a == "--run-dir") { out$run_dir <- argv[i + 1]; i <- i + 2 }
    else if (a == "--out") { out$out <- argv[i + 1]; i <- i + 2 }
    else stop("Unknown argument: ", a)
  }
  if (is.null(out$run_dir)) stop("Required arg: --run-dir <dir>")
  out
}

main <- function() {
  args    <- parse_args(commandArgs(trailingOnly = TRUE))
  run_dir <- normalizePath(args$run_dir, mustWork = TRUE)
  mat_dir <- file.path(run_dir, "cppRATE_matrices")

  coef_file <- file.path(mat_dir, "coefficients.csv")
  ld_file   <- file.path(mat_dir, "ld_pruning_summary.csv")
  if (!file.exists(coef_file)) stop("Missing: ", coef_file)

  # A run fitted with --ld_pruning false writes no pruning summary. Its columns
  # already span every variant, so there is nothing to expand and the reduced
  # medians are the answer.
  pruned_run <- file.exists(ld_file)
  if (!pruned_run)
    message("No ", basename(ld_file), " found; treating this as an unpruned run.")

  out_csv <- if (!is.null(args$out)) args$out else
    file.path(run_dir, "fitted_model", "depruned_variant_effects_from_draws.csv")

  message("Run dir: ", run_dir)
  message("Reading posterior draws: ", coef_file,
          " (", round(file.size(coef_file) / 1e9, 2), " GB)")
  dt <- data.table::fread(coef_file, header = FALSE, showProgress = FALSE)
  n_kept <- ncol(dt)
  message("  draws = ", nrow(dt), ", kept variants = ", n_kept)

  # Same reduction the pipeline applies in deprun_effects_from_draws():
  # posterior median, and significant iff the 89% CI excludes zero.
  medians_kept <- vapply(dt, stats::median, numeric(1))
  ci           <- vapply(dt, stats::quantile, numeric(2), probs = c(0.055, 0.945))
  signif_kept  <- !(ci[1, ] <= 0 & ci[2, ] >= 0)
  rm(dt); invisible(gc())
  message("  reduced to medians + 89% CI; ", sum(signif_kept),
          " of ", n_kept, " kept variants significant")

  if (pruned_run) {
    final_rep_variants <- parse_ld_summary(ld_file)
    pruned <- unique(na.omit(
      as.vector(as.matrix(final_rep_variants[, 2:ncol(final_rep_variants)]))
    ))
    V <- n_kept + length(pruned)
    message("Total variants V = ", V, "  (kept ", n_kept,
            " + pruned ", length(pruned), ")")

    # Pruned negative-correlation complements get their copied beta negated.
    flip_sign <- read_correlation_directions(run_dir, V)
    message("Sign-flips applied to ",
            if (is.null(flip_sign)) 0 else sum(flip_sign), " pruned complements")

    eff <- deprun_effects_from_draws_reduced(medians_kept, signif_kept,
                                             final_rep_variants, V, flip_sign)

    n_na <- sum(is.na(eff$medians))
    message("Depruned medians: length ", length(eff$medians), ", NA = ", n_na)
    if (n_na > 0L) warning(n_na, " variants have NA median (unmapped after deprune)")
    message("Significant after de-pruning: ", sum(eff$signif, na.rm = TRUE))
  } else {
    V   <- n_kept
    eff <- list(medians = unname(medians_kept), signif = unname(signif_kept))
    message("Total variants V = ", V, "  (nothing pruned)")
  }

  dir.create(dirname(out_csv), showWarnings = FALSE, recursive = TRUE)
  data.table::fwrite(
    data.frame(variant_id = seq_len(V),
               median     = eff$medians,
               signif     = eff$signif),
    out_csv)
  message("Wrote ", out_csv)
}

# The expansion half of cpprate.R::deprun_effects_from_draws(), taking the
# already-reduced per-column medians/signif instead of the full draws matrix.
# Identical copy semantics: each representative's median and significance are
# copied onto every variant pruned under it, negating the median where the
# pruned variant is the representative's complement.
deprun_effects_from_draws_reduced <- function(medians_kept, signif_kept,
                                              final_rep_variants, V,
                                              flip_sign = NULL) {
  pruned <- unique(na.omit(
    as.vector(as.matrix(final_rep_variants[, 2:ncol(final_rep_variants)]))
  ))
  total_med <- rep(NA_real_, V)
  total_sig <- rep(NA, V)
  total_med[-pruned] <- medians_kept
  total_sig[-pruned] <- signif_kept
  for (i in seq_len(nrow(final_rep_variants))) {
    n_p    <- sum(!is.na(final_rep_variants[i, ])) - 1L
    rep_id <- final_rep_variants[i, 1L]
    for (j in seq_len(n_p)) {
      idx        <- final_rep_variants[i, j + 1L]
      rep_median <- total_med[rep_id]
      total_med[idx] <- if (!is.null(flip_sign) && isTRUE(flip_sign[idx]))
        -rep_median else rep_median
      total_sig[idx] <- total_sig[rep_id]
    }
  }
  list(medians = total_med, signif = as.logical(total_sig))
}

main()
