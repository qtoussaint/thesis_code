#!/usr/bin/env Rscript
# Compute ROPE-region Bayes factor (BF_10), Savage-Dickey BF at point null
# beta = 0, and probability of direction (pd) for each variant in a
# significance list produced by rope_significance.R.
#
# Output: <input dir>/<basename>_with_bf_pd.csv with bf_rope, bf_point, pd
# columns appended (or per-cutpoint variants for PPOM union lists).

suppressPackageStartupMessages({
  library(optparse)
  library(posterior)
})

has_bayestestR <- requireNamespace("bayestestR", quietly = TRUE)
has_logspline  <- requireNamespace("logspline",  quietly = TRUE)

# ---- CLI --------------------------------------------------------------------

option_list <- list(
  make_option("--variant-list",     type = "character", default = NULL,
              help = "Significance list CSV from rope_significance.R (required)"),
  make_option("--run-dir",          type = "character", default = NULL,
              help = "Parent of fitted_model/ (required unless --rds is given)"),
  make_option("--variant-index",    type = "character", default = NULL,
              help = "<dataset>_variant_index.csv (required)"),
  make_option("--rds",              type = "character", default = NULL,
              help = "Explicit path to fitted RDS (overrides run-dir auto-pick)"),
  make_option("--rope-bounds",      type = "character", default = NULL,
              help = "rope_bounds.txt path; default = <variant-list dir>/rope_bounds.txt"),
  make_option("--n-prior-samples",  type = "integer",   default = 200000L,
              help = "Total prior null samples [%default]"),
  make_option("--output",           type = "character", default = NULL,
              help = "Output CSV; default = <input dir>/<basename>_with_bf_pd.csv"),
  make_option("--seed",             type = "integer",   default = 42L,
              help = "RNG seed for prior null sampling [%default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt[["variant-list"]]) ||
    is.null(opt[["variant-index"]]) ||
    (is.null(opt[["run-dir"]]) && is.null(opt[["rds"]]))) {
  stop("Required: --variant-list, --variant-index, and (--run-dir or --rds)")
}

variant_list_path  <- opt[["variant-list"]]
variant_index_path <- opt[["variant-index"]]
run_dir            <- opt[["run-dir"]]
rds_path           <- opt[["rds"]]
rope_bounds_path   <- opt[["rope-bounds"]]
n_prior            <- opt[["n-prior-samples"]]
out_path           <- opt[["output"]]

set.seed(opt[["seed"]])

# ---- Resolve paths ----------------------------------------------------------

list_dir <- dirname(normalizePath(variant_list_path))
if (is.null(rope_bounds_path)) {
  rope_bounds_path <- file.path(list_dir, "rope_bounds.txt")
}
if (is.null(out_path)) {
  base <- tools::file_path_sans_ext(basename(variant_list_path))
  out_path <- file.path(list_dir, paste0(base, "_with_bf_pd.csv"))
}
if (is.null(rds_path)) {
  fitted_dir <- file.path(run_dir, "fitted_model")
  rds_files <- list.files(fitted_dir, pattern = "\\.RDS$", full.names = TRUE)
  if (length(rds_files) == 0L) stop("No .RDS in ", fitted_dir)
  rds_path <- rds_files[which.max(file.mtime(rds_files))]
}

# ---- Load list & detect mode ------------------------------------------------

vlist <- read.csv(variant_list_path, stringsAsFactors = FALSE)
if (!"variant_name" %in% names(vlist)) {
  stop("Input list must contain a 'variant_name' column")
}

has_cutpoint_col <- "cutpoint" %in% names(vlist)
union_cols <- grep("^signif_(median|ci)_k\\d+$", names(vlist), value = TRUE)
is_union   <- length(union_cols) > 0L && !has_cutpoint_col

mode <- if (is_union) "union" else if (has_cutpoint_col) "per_cutpoint" else "vector"
message("[bf_pd] mode: ", mode, "  (", nrow(vlist), " variants)")

# ---- Load fit + posterior draws --------------------------------------------

message("[bf_pd] loading ", rds_path)
fit <- readRDS(rds_path)

draws <- posterior::as_draws_df(fit$draws(variables = c("tau", "c2", "beta_variant_std")))
draws <- as.data.frame(draws)

is_matrix_param <- any(grepl("^beta_variant_std\\[\\d+,\\d+\\]$", names(draws)))
if (mode == "union" && !is_matrix_param) {
  stop("Input list looks like a PPOM union but the fit's beta_variant_std is scalar-indexed")
}
if (mode == "per_cutpoint" && !is_matrix_param) {
  stop("Input list has a 'cutpoint' column but the fit's beta_variant_std is scalar-indexed")
}

# Determine K-1 from draws (PPOM) or list union columns
if (is_matrix_param) {
  bvs_cols <- grep("^beta_variant_std\\[\\d+,\\d+\\]$", names(draws), value = TRUE)
  vk <- do.call(rbind,
                lapply(regmatches(bvs_cols, regexec("^beta_variant_std\\[(\\d+),(\\d+)\\]$", bvs_cols)),
                       function(m) as.integer(m[-1])))
  K_m1 <- max(vk[, 2])
} else {
  K_m1 <- 1L
}

# ---- Variant index lookup ---------------------------------------------------

vidx <- read.csv(variant_index_path, stringsAsFactors = FALSE)
if (!"variant_name" %in% names(vidx)) stop("variant-index missing 'variant_name'")

v_idx <- match(vlist$variant_name, vidx$variant_name)
if (any(is.na(v_idx))) {
  missing <- vlist$variant_name[is.na(v_idx)]
  stop("These variants are not in variant-index: ", paste(head(missing, 5), collapse=", "))
}

# ---- ROPE bounds ------------------------------------------------------------

rope <- NULL
if (file.exists(rope_bounds_path)) {
  lines <- readLines(rope_bounds_path)
  pull_num <- function(prefix) {
    line <- grep(paste0("^", prefix), lines, value = TRUE)
    if (length(line) == 0L) return(NA_real_)
    as.numeric(sub(".*:\\s*", "", line[1]))
  }
  rope <- c(pull_num("rope_lo \\(standardized\\)"),
            pull_num("rope_hi \\(standardized\\)"))
  message(sprintf("[bf_pd] ROPE (standardized): [%.5g, %.5g]", rope[1], rope[2]))
} else {
  warning("rope_bounds.txt not found at ", rope_bounds_path,
          " - bf_rope will be NA; bf_point and pd still computed.")
}

# ---- Build prior null sample ------------------------------------------------

tau_post <- draws$tau
c2_post  <- draws$c2
S <- length(tau_post)
M <- max(1L, ceiling(n_prior / S))
message(sprintf("[bf_pd] prior null: %d posterior x %d prior = %d samples", S, M, S * M))

lambda_null <- abs(rcauchy(S * M, 0, 2))
z_null      <- rnorm(S * M)
tau_rep <- rep(tau_post, each = M)
c2_rep  <- rep(c2_post,  each = M)
lambda_tilde_null <- sqrt((c2_rep * lambda_null^2) /
                          (c2_rep + tau_rep^2 * lambda_null^2))
beta_null_std <- z_null * tau_rep * lambda_tilde_null
rm(lambda_null, z_null, tau_rep, c2_rep, lambda_tilde_null); gc(verbose = FALSE)

# ---- BF / pd helpers --------------------------------------------------------

# BF_10 against ROPE region null
bf_rope_10 <- function(post, prior, rope) {
  if (is.null(rope) || any(is.na(rope))) return(NA_real_)
  if (has_bayestestR) {
    res <- tryCatch(
      bayestestR::bayesfactor_parameters(post, prior = prior,
                                         null = c(rope[1], rope[2]),
                                         verbose = FALSE),
      error = function(e) NULL
    )
    if (!is.null(res)) {
      # bayestestR returns log_BF (natural log) in favor of the alternative
      logbf <- res$log_BF[1]
      return(exp(logbf))
    }
  }
  # Frequency-ratio fallback (asymptotically equivalent to Savage-Dickey for large N)
  p_post_out  <- mean(post  < rope[1] | post  > rope[2])
  p_prior_out <- mean(prior < rope[1] | prior > rope[2])
  if (p_post_out %in% c(0, 1) || p_prior_out %in% c(0, 1)) return(Inf)
  (p_post_out / (1 - p_post_out)) / (p_prior_out / (1 - p_prior_out))
}

# Savage-Dickey BF_10 at point null beta = 0
bf_point_10 <- function(post, prior) {
  if (has_bayestestR) {
    res <- tryCatch(
      bayestestR::bayesfactor_pointnull(post, prior = prior, verbose = FALSE),
      error = function(e) NULL
    )
    if (!is.null(res)) return(exp(res$log_BF[1]))
  }
  if (has_logspline) {
    fp <- tryCatch(logspline::logspline(post),  error = function(e) NULL)
    fq <- tryCatch(logspline::logspline(prior), error = function(e) NULL)
    if (!is.null(fp) && !is.null(fq)) {
      d_post  <- logspline::dlogspline(0, fp)
      d_prior <- logspline::dlogspline(0, fq)
      if (d_post > 0) return(d_prior / d_post)
    }
  }
  # density() fallback
  rng <- range(c(post, prior, 0))
  d_post  <- approx(density(post,  from = rng[1], to = rng[2])$x,
                    density(post,  from = rng[1], to = rng[2])$y, xout = 0)$y
  d_prior <- approx(density(prior, from = rng[1], to = rng[2])$x,
                    density(prior, from = rng[1], to = rng[2])$y, xout = 0)$y
  if (is.na(d_post) || d_post <= 0) return(NA_real_)
  d_prior / d_post
}

pd_compute <- function(post) {
  if (has_bayestestR) {
    res <- tryCatch(bayestestR::p_direction(post), error = function(e) NULL)
    if (!is.null(res)) return(as.numeric(res$pd[1]))
  }
  max(mean(post > 0), mean(post < 0))
}

# ---- Extract beta draws for one (variant, cutpoint) ------------------------

get_beta_draws <- function(v, k = NULL) {
  col <- if (is.null(k)) sprintf("beta_variant_std[%d]", v)
         else            sprintf("beta_variant_std[%d,%d]", v, k)
  if (!col %in% names(draws)) stop("Missing draws column: ", col)
  draws[[col]]
}

# ---- Compute per variant ----------------------------------------------------

out <- vlist
if (mode == "vector") {
  out$bf_rope  <- NA_real_
  out$bf_point <- NA_real_
  out$pd       <- NA_real_
  for (i in seq_len(nrow(vlist))) {
    post <- get_beta_draws(v_idx[i])
    out$bf_rope[i]  <- bf_rope_10(post, beta_null_std, rope)
    out$bf_point[i] <- bf_point_10(post, beta_null_std)
    out$pd[i]       <- pd_compute(post)
  }

} else if (mode == "per_cutpoint") {
  out$bf_rope  <- NA_real_
  out$bf_point <- NA_real_
  out$pd       <- NA_real_
  for (i in seq_len(nrow(vlist))) {
    k_i  <- as.integer(vlist$cutpoint[i])
    post <- get_beta_draws(v_idx[i], k_i)
    out$bf_rope[i]  <- bf_rope_10(post, beta_null_std, rope)
    out$bf_point[i] <- bf_point_10(post, beta_null_std)
    out$pd[i]       <- pd_compute(post)
  }

} else {  # union
  for (k in seq_len(K_m1)) {
    out[[sprintf("bf_rope_k%d",  k)]] <- NA_real_
    out[[sprintf("bf_point_k%d", k)]] <- NA_real_
    out[[sprintf("pd_k%d",       k)]] <- NA_real_
  }
  for (i in seq_len(nrow(vlist))) {
    for (k in seq_len(K_m1)) {
      post <- get_beta_draws(v_idx[i], k)
      out[[sprintf("bf_rope_k%d",  k)]][i] <- bf_rope_10(post, beta_null_std, rope)
      out[[sprintf("bf_point_k%d", k)]][i] <- bf_point_10(post, beta_null_std)
      out[[sprintf("pd_k%d",       k)]][i] <- pd_compute(post)
    }
  }
}

# ---- Console summary --------------------------------------------------------

summarise_col <- function(col_name) {
  if (!col_name %in% names(out)) return(invisible())
  v <- out[[col_name]]
  message(sprintf("  %-14s  median=%.3g  max=%.3g  n_NA=%d",
                  col_name, median(v, na.rm = TRUE), max(v, na.rm = TRUE),
                  sum(is.na(v))))
}

message("[bf_pd] summary:")
for (cn in grep("^(bf_rope|bf_point|pd)", names(out), value = TRUE)) summarise_col(cn)

# ---- Write output -----------------------------------------------------------

write.csv(out, out_path, row.names = FALSE)
message("[bf_pd] wrote ", out_path)
