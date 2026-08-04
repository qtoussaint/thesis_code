#!/usr/bin/env Rscript
# Per-variant median beta, 89% CI significance and cppRATE on the ORIGINAL
# variant index, de-pruned by following BacPrune's representative chains.
#
# WHY THIS EXISTS RATHER THAN THE PIPELINE'S OWN DE-PRUNING
#
# depruning_rates() and deprun_effects_from_draws() both expand a pruned run by
# reading ld_pruning_summary.csv and copying each representative's value onto
# the variants listed against it. That assumes every representative is itself a
# kept variant. Under r2 >= 1 pruning it is: perfect correlation is transitive,
# so the representative relation collapses into disjoint equivalence classes.
#
# Under |D'| >= 1 it is not. |D'| reaches 1 whenever one of the four haplotypes
# is absent, which is not a transitive relation, and BacPrune's greedy sweep
# then emits chains: variant 1 is represented by 2, which is itself represented
# by 9, and only the end of the chain survives. On the TB rifampicin |D'| arm,
# 3,397 of 3,511 representatives are themselves pruned, chains run up to 8 deep,
# and 69,870 of 75,272 variants sit deeper than one hop. The pipeline's
# de-pruning gives all of those nothing: its RATE_values_depruned.txt has 6,087
# rows where it should have 75,272.
#
# ld_pruning_summary.csv is also the wrong shape to read here. parse_ld_summary()
# widens it to one column per pruned variant, and on this arm a single
# representative prunes 33,025 variants, so the frame is 3,511 x 33,025 and the
# parse alone exhausts memory.
#
# direction_of_correlation.csv is used instead. It has exactly one row per
# variant -- Variant, Status, Representative Variant -- so it needs no widening,
# and it carries the sign relation as well as the link. Representatives are
# self-pointing with status "representative" and are exactly the kept set, both
# checked below.
#
# SIGNS COMPOSE ALONG THE CHAIN. A "negative_correlation" variant is the
# genotype complement of its representative, so its copied beta is negated. Over
# a chain the flips XOR: a complement of a complement is back in the original
# encoding. RATE is unsigned and so is never flipped, and significance survives
# negation because a CI excluding zero still excludes zero after a sign change.
#
# Runs fitted with --ld_pruning false have no direction file, nothing to expand,
# and pass straight through.
#
# Usage:
#   Rscript deprune_from_directions.R --run-dir <inference dir> [--out <path>]
#
# Output CSV columns:
#   variant_id      1..V on the original variant index
#   median, signif  unstandardized posterior median and 89% CI significance
#   rate            cppRATE relative centrality, NA if the run has no RATE yet
#   representative  the kept variant this row's values came from
#   chain_depth     hops to that representative; 0 for a kept variant
#   sign_flipped    whether the median was negated on the way

suppressPackageStartupMessages({
  library(data.table)
})

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

# Posterior median and 89% CI significance per FITTED variant, straight from the
# draws the run already wrote. cppRATE_matrices/coefficients.csv is written with
# fitted_model$draws(variables = "beta_variant"), so it holds beta_variant and
# nothing else -- unlike the pipeline's depruned_variant_effects.csv, whose
# grepl("beta_variant") also matches beta_variant_std and therefore reports the
# STANDARDIZED effect. Reading the CSV also avoids the multi-GB fitted RDS.
reduce_draws <- function(coef_file) {
  message("Reading posterior draws: ", coef_file,
          " (", round(file.size(coef_file) / 1e9, 2), " GB)")
  dt <- data.table::fread(coef_file, header = FALSE, showProgress = FALSE)
  message("  draws = ", nrow(dt), ", fitted variants = ", ncol(dt))
  medians <- vapply(dt, stats::median, numeric(1))
  ci      <- vapply(dt, stats::quantile, numeric(2), probs = c(0.055, 0.945))
  signif  <- !(ci[1, ] <= 0 & ci[2, ] >= 0)
  rm(dt); invisible(gc())
  message("  ", sum(signif), " of ", length(medians),
          " fitted variants significant at the 89% CI")
  list(median = unname(medians), signif = unname(signif))
}

# cppRATE output for the variants the run actually fitted, in fitted order.
# RATE_values.txt is cppRATE's own output, three header lines then one row per
# fitted variant. The pipeline's RATE_values_depruned.txt is deliberately NOT
# read: on a chained arm it is the truncated file this script exists to replace.
read_rate_kept <- function(run_dir, n_kept) {
  f <- file.path(run_dir, "cppRATE_results", "RATE_values.txt")
  if (!file.exists(f)) {
    message("No RATE_values.txt yet; rate column will be NA.")
    return(NULL)
  }
  d <- data.table::fread(f, skip = 3, header = FALSE, blank.lines.skip = TRUE,
                         col.names = c("snp_id", "rate", "kld"))
  if (nrow(d) != n_kept) {
    warning("RATE rows (", nrow(d), ") != fitted variants (", n_kept,
            "); rate column will be NA.", call. = FALSE)
    return(NULL)
  }
  d$rate
}

# Resolve every variant to the kept representative at the end of its chain, and
# to the net sign relation with it, by pointer doubling: repeatedly replace each
# variant's representative with its representative's representative, XOR-ing the
# flips as they compose. Representatives are self-pointing fixed points, so the
# iteration converges once every chain has bottomed out -- in log2(depth) passes
# rather than one pass per hop.
resolve_chains <- function(rep_of, flip_of) {
  root <- rep_of
  flip <- flip_of

  for (iter in seq_len(64L)) {
    parent <- root[root]
    if (identical(parent, root)) break
    flip <- xor(flip, flip[root])
    root <- parent
  }
  if (any(root[root] != root))
    stop("Representative chains did not resolve; the direction file may contain a cycle.")

  list(root = root, flip = flip)
}

# Exact chain depth, walked one hop at a time. resolve_chains() is the fast path
# for the values; this is a small integer pass kept separate so the reported
# depth is the true hop count rather than a doubling artefact.
chain_depths <- function(rep_of) {
  V     <- length(rep_of)
  cur   <- rep_of
  depth <- as.integer(rep_of != seq_len(V))
  for (iter in seq_len(1000L)) {
    moving <- cur != rep_of[cur]
    if (!any(moving)) break
    cur[moving]   <- rep_of[cur[moving]]
    depth[moving] <- depth[moving] + 1L
  }
  depth
}

main <- function() {
  args    <- parse_args(commandArgs(trailingOnly = TRUE))
  run_dir <- normalizePath(args$run_dir, mustWork = TRUE)
  mat_dir <- file.path(run_dir, "cppRATE_matrices")
  message("Run dir: ", run_dir)

  coef_file <- file.path(mat_dir, "coefficients.csv")
  if (!file.exists(coef_file)) stop("Missing: ", coef_file)
  dir_file  <- file.path(mat_dir, "direction_of_correlation.csv")

  out_csv <- if (!is.null(args$out)) args$out else
    file.path(run_dir, "fitted_model", "depruned_effects_and_rate.csv")

  red    <- reduce_draws(coef_file)
  n_kept <- length(red$median)
  rate_kept <- read_rate_kept(run_dir, n_kept)

  if (!file.exists(dir_file)) {
    # Unpruned run: the fitted variants are already every variant.
    message("No direction_of_correlation.csv; treating this as an unpruned run.")
    V   <- n_kept
    out <- data.table(
      variant_id     = seq_len(V),
      median         = red$median,
      signif         = red$signif,
      rate           = if (is.null(rate_kept)) NA_real_ else rate_kept,
      representative = seq_len(V),
      chain_depth    = 0L,
      sign_flipped   = FALSE)
  } else {
    dirs <- data.table::fread(dir_file)
    setnames(dirs, c("Variant", "Status", "Representative Variant"),
             c("variant", "status", "rep"), skip_absent = TRUE)
    V <- nrow(dirs)
    message("Direction file: ", V, " variants")

    data.table::setorder(dirs, variant)
    stopifnot(identical(dirs$variant, seq_len(V)))

    kept_ids <- dirs$variant[dirs$status == "representative"]
    if (length(kept_ids) != n_kept)
      stop("Direction file lists ", length(kept_ids), " representatives but the ",
           "run fitted ", n_kept, " variants.")
    # Representatives must be self-pointing, or the chain walk has no fixed point.
    stopifnot(all(dirs$rep[kept_ids] == kept_ids))

    res   <- resolve_chains(dirs$rep, dirs$status == "negative_correlation")
    depth <- chain_depths(dirs$rep)
    if (!all(res$root %in% kept_ids))
      stop("Some variants resolved to a representative that was not fitted.")

    # kept_ids is ascending and matches the column order of the draws matrix and
    # the row order of RATE_values.txt, so position in kept_ids is the index into
    # both. Verified against bacprune_rust_results.csv's header order.
    pos <- match(res$root, kept_ids)

    out <- data.table(
      variant_id     = seq_len(V),
      median         = ifelse(res$flip, -red$median[pos], red$median[pos]),
      signif         = red$signif[pos],
      rate           = if (is.null(rate_kept)) NA_real_ else rate_kept[pos],
      representative = res$root,
      chain_depth    = depth,
      sign_flipped   = res$flip)

    message("Chain depth: max ", max(depth), ", ",
            sum(depth > 1L), " of ", V, " variants deeper than one hop")
    message("Sign flips applied to ", sum(res$flip), " variants")
  }

  message("Variants covered: ", nrow(out), "; significant: ", sum(out$signif))
  dir.create(dirname(out_csv), showWarnings = FALSE, recursive = TRUE)
  data.table::fwrite(out, out_csv)
  message("Wrote ", out_csv)
}

main()
