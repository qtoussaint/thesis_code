#!/usr/bin/env Rscript
# Regenerate fitted_model/depruned_variant_effects.csv for a POM run that never
# wrote one (the pipeline only writes it for POM when --phandango was passed).
# We reconstruct it from saved artifacts WITHOUT touching the (14 GB) fitted RDS:
# the posterior beta_variant draws are already on disk as
# cppRATE_matrices/coefficients.csv, so we feed that matrix straight into
# gwas_workflow's own deprun_effects_from_draws() — the exact call the pipeline
# makes for POM effects (pipeline.R POM block ~287-319). Output columns match the
# POM CSV the pipeline writes: variant_id, median, signif.
#
# Usage:
#   Rscript regen_pom_variant_effects.R --run-dir <POM inference dir> [--n-total V]
#
# n-total defaults to (kept variants in coefficients.csv) + (pruned variants in
# ld_pruning_summary.csv), which is exactly the V that deprun_effects_from_draws
# expects.

suppressPackageStartupMessages({
  library(data.table)
})

GWAS_WORKFLOW_R <- "/nfs/research/jlees/jacqueline/gwas_workflow/code/gwas_workflow/R"
# Both files define only functions (no top-level code), so sourcing is safe and
# lets us reuse the pipeline's de-pruning logic verbatim.
source(file.path(GWAS_WORKFLOW_R, "ld_pruning.R"))
source(file.path(GWAS_WORKFLOW_R, "cpprate.R"))

parse_args <- function(argv) {
  out <- list(run_dir = NULL, n_total = NULL)
  i <- 1
  while (i <= length(argv)) {
    a <- argv[i]
    if (a == "--run-dir") {
      out$run_dir <- argv[i + 1]; i <- i + 2
    } else if (a == "--n-total") {
      out$n_total <- as.integer(argv[i + 1]); i <- i + 2
    } else {
      stop("Unknown argument: ", a)
    }
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
  if (!file.exists(ld_file))   stop("Missing: ", ld_file)

  message("Run dir: ", run_dir)
  # coefficients.csv is rows = draws, cols = kept variants (beta_variant[i]); see
  # cpprate.R::write_coefficient_matrix(). The full matrix is ~2 GB, so read it in
  # column chunks and reduce each column to its median + 89% CI on the fly rather
  # than holding all draws in memory at once.
  n_kept <- ncol(data.table::fread(coef_file, header = FALSE, nrows = 1L))
  message("Reducing posterior draws in column chunks: ", coef_file)
  message("  kept variants = ", n_kept)

  medians_kept <- numeric(n_kept)
  signif_kept  <- logical(n_kept)
  chunk <- 4000L
  starts <- seq.int(1L, n_kept, by = chunk)
  for (s in starts) {
    cols <- s:min(s + chunk - 1L, n_kept)
    dt <- data.table::fread(coef_file, header = FALSE, select = cols,
                            showProgress = FALSE)
    medians_kept[cols] <- vapply(dt, stats::median, numeric(1))
    ci <- vapply(dt, stats::quantile, numeric(2), probs = c(0.055, 0.945))
    signif_kept[cols]  <- !(ci[1, ] <= 0 & ci[2, ] >= 0)
    rm(dt)
    message("  reduced columns ", cols[1], "-", cols[length(cols)])
  }

  final_rep_variants <- parse_ld_summary(ld_file)
  n_pruned <- length(unique(na.omit(
    as.vector(as.matrix(final_rep_variants[, 2:ncol(final_rep_variants)]))
  )))

  V <- if (!is.null(args$n_total)) args$n_total else n_kept + n_pruned
  message("Total variants V = ", V, "  (kept ", n_kept, " + pruned ", n_pruned, ")")
  if (V != n_kept + n_pruned) {
    warning("V (", V, ") != kept+pruned (", n_kept + n_pruned,
            "); deprune expects kept == V - pruned.")
  }

  # Pruned negative-correlation complements get their copied beta negated.
  flip_sign <- read_correlation_directions(run_dir, V)
  message("Sign-flips applied to ", if (is.null(flip_sign)) 0 else sum(flip_sign),
          " pruned complement variants")

  # De-prune expansion, replicated from cpprate.R deprun_effects_from_draws()
  # lines 180-204 (the only part that needs the per-column medians/signif, which
  # we computed chunk-wise above instead of from the full draws matrix).
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
      if (!is.null(flip_sign) && isTRUE(flip_sign[idx])) {
        total_med[idx] <- -rep_median
      } else {
        total_med[idx] <- rep_median
      }
      total_sig[idx] <- total_sig[rep_id]
    }
  }
  eff <- list(medians = total_med, signif = as.logical(total_sig))

  n_na <- sum(is.na(eff$medians))
  message("Depruned medians: length ", length(eff$medians), ", NA = ", n_na)
  if (n_na > 0L) warning(n_na, " variants have NA median (unmapped after deprune)")

  out_csv <- file.path(run_dir, "fitted_model", "depruned_variant_effects.csv")
  dir.create(dirname(out_csv), showWarnings = FALSE, recursive = TRUE)
  data.table::fwrite(
    data.frame(variant_id = seq_len(V),
               median     = eff$medians,
               signif     = eff$signif),
    out_csv
  )
  message("Wrote ", out_csv)
}

main()
