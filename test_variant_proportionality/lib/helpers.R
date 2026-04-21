# Helpers for comparing POM vs PPOM variant effect posteriors.
# Mirrors conventions used in gwaspipeline (R/pipeline.R, R/cpprate.R, R/ld_pruning.R).

suppressPackageStartupMessages({
  library(posterior)
})

parse_ld_summary <- function(summary_path) {
  index <- read.csv(file = summary_path)

  split_pruned <- strsplit(index[, 2], ",\\s*")
  max_len <- max(sapply(split_pruned, length))

  expanded <- do.call(rbind, lapply(split_pruned, function(x) {
    length(x) <- max_len
    x
  }))
  expanded <- as.data.frame(expanded, stringsAsFactors = FALSE)
  colnames(expanded) <- paste0("PrunedSNP_", seq_len(max_len))

  result <- cbind(index[, 1], expanded)
  colnames(result)[1] <- "RepresentativeSNP"
  result[, -1] <- lapply(result[, -1], as.numeric)

  result[, -1] <- lapply(result[, -1], function(x) ifelse(is.na(x), NA, x + 1L))
  result[, 1] <- result[, 1] + 1L

  result <- result[order(result[, 1]), ]
  rownames(result) <- seq_len(nrow(result))
  result
}

read_kept_variants <- function(output_dir) {
  pruned_matrix <- read.csv(
    file = file.path(output_dir, "cppRATE_matrices", "bacprune_rust_results.csv"),
    header = TRUE
  )
  as.numeric(sub("X", "", colnames(pruned_matrix)))
}

find_fit_rds <- function(output_dir) {
  candidates <- Sys.glob(file.path(output_dir, "fitted_model", "*.RDS"))
  candidates <- candidates[!grepl("\\.RDS\\.", candidates)]
  if (length(candidates) == 0L) {
    stop("No fitted model RDS found under ", file.path(output_dir, "fitted_model"))
  }
  if (length(candidates) > 1L) {
    warning("Multiple RDS files found in ", output_dir,
            "; using the most recently modified.")
    candidates <- candidates[order(file.info(candidates)$mtime, decreasing = TRUE)]
  }
  candidates[1]
}

load_beta_variant_draws <- function(rds_path) {
  fit <- readRDS(rds_path)
  draws <- fit$draws(variables = "beta_variant")
  posterior::as_draws_matrix(draws)
}

# Expand a [draws × V_kept] matrix to [draws × V_total] by copying each
# representative's draws to its pruned neighbours. Uses the same convention as
# gwaspipeline::deprun_effects: variants not in `final_rep_variants[, -1]` get
# the kept draws directly; pruned variants inherit from their representative.
deprun_draws_single <- function(draws_kept, kept_variants, final_rep_variants,
                                n_total_variants) {
  stopifnot(ncol(draws_kept) == length(kept_variants))

  pruned_all <- unique(na.omit(
    as.vector(as.matrix(final_rep_variants[, 2:ncol(final_rep_variants)]))
  ))

  n_draws <- nrow(draws_kept)
  out <- matrix(NA_real_, nrow = n_draws, ncol = n_total_variants)

  nonpruned_slots <- setdiff(seq_len(n_total_variants), pruned_all)
  if (length(nonpruned_slots) != ncol(draws_kept)) {
    stop("Mismatch: kept variants (", ncol(draws_kept),
         ") vs non-pruned slots (", length(nonpruned_slots), ").")
  }
  out[, nonpruned_slots] <- draws_kept

  # For each representative, copy its draws to its pruned neighbours.
  kept_lookup <- match(final_rep_variants$RepresentativeSNP, kept_variants)
  for (i in seq_len(nrow(final_rep_variants))) {
    rep_kept_col <- kept_lookup[i]
    if (is.na(rep_kept_col)) next
    n_p <- sum(!is.na(final_rep_variants[i, ])) - 1L
    if (n_p == 0L) next
    for (j in seq_len(n_p)) {
      pruned_id <- final_rep_variants[i, j + 1L]
      out[, pruned_id] <- draws_kept[, rep_kept_col]
    }
  }

  out
}

# PPOM draws come back as [draws × (V_kept * (K-1))] organised cutpoint-major.
# Returns a list of length (K-1), each entry a [draws × V_total] matrix.
deprun_draws_ppom <- function(ppom_draws, kept_variants, final_rep_variants,
                              n_total_variants) {
  V_kept <- length(kept_variants)
  total_cols <- ncol(ppom_draws)
  if (total_cols %% V_kept != 0L) {
    stop("PPOM draws column count (", total_cols,
         ") is not a multiple of the kept-variant count (", V_kept, ").")
  }
  n_cutpoints <- total_cols %/% V_kept

  lapply(seq_len(n_cutpoints), function(c) {
    col_start <- V_kept * (c - 1L) + 1L
    col_end   <- V_kept * c
    deprun_draws_single(
      draws_kept         = ppom_draws[, col_start:col_end, drop = FALSE],
      kept_variants      = kept_variants,
      final_rep_variants = final_rep_variants,
      n_total_variants   = n_total_variants
    )
  })
}

# Per-variant median and CI across draws.
column_quantiles <- function(draws_matrix, probs) {
  q <- matrixStats::colQuantiles(draws_matrix, probs = probs, na.rm = TRUE)
  colnames(q) <- paste0("q", probs)
  q
}

# Fallback if matrixStats is not available.
column_quantiles_fallback <- function(draws_matrix, probs) {
  out <- t(apply(draws_matrix, 2, stats::quantile, probs = probs, na.rm = TRUE))
  colnames(out) <- paste0("q", probs)
  out
}

safe_column_quantiles <- function(draws_matrix, probs) {
  if (requireNamespace("matrixStats", quietly = TRUE)) {
    column_quantiles(draws_matrix, probs)
  } else {
    column_quantiles_fallback(draws_matrix, probs)
  }
}

ci_tail_probs <- function(ci_level) {
  alpha <- 1 - ci_level
  c(alpha / 2, 0.5, 1 - alpha / 2)
}

# Map each variant_id to its representative (for metadata).
build_representative_map <- function(final_rep_variants, n_total_variants) {
  rep_of <- seq_len(n_total_variants)  # default: each variant represents itself
  is_rep <- rep(TRUE, n_total_variants)
  for (i in seq_len(nrow(final_rep_variants))) {
    rep_id <- final_rep_variants[i, 1L]
    n_p <- sum(!is.na(final_rep_variants[i, ])) - 1L
    if (n_p == 0L) next
    for (j in seq_len(n_p)) {
      pruned_id <- final_rep_variants[i, j + 1L]
      rep_of[pruned_id] <- rep_id
      is_rep[pruned_id] <- FALSE
    }
  }
  data.frame(
    variant_id        = seq_len(n_total_variants),
    is_representative = is_rep,
    representative_id = rep_of
  )
}
