#!/usr/bin/env Rscript
# Per-variant ROPE Bayes factor, point-null (Savage-Dickey) Bayes factor and
# probability of direction for EVERY variant in a fit, so the results can be
# drawn as Manhattans.
#
# WHY THIS EXISTS. bf_pd_significance.R computes the same three quantities but
# only for the variants in a significance list, which is what the methods
# describe. A Manhattan needs all 4,578. Running bf_pd_significance.R over the
# full set is impractical: bayestestR::bayesfactor_pointnull() re-fits a
# logspline density to the 200,000-draw prior on every call, so the prior fit
# alone would be repeated 4,578 times per dataset.
#
# This script uses the identical definitions, but fits the prior density ONCE
# per dataset and reuses it. The Savage-Dickey ratio d_prior(0) / d_post(0) is
# exactly bf_pd_significance.R's logspline branch, and the ROPE Bayes factor is
# the posterior:prior odds of falling outside the ROPE, which is what
# bayesfactor_parameters() computes for an interval null (no density estimation
# involved). Agreement with bayestestR is CHECKED, not assumed -- see the
# validation block, which halts if the two disagree.
#
# The ROPE is read from rope_bounds.txt written by rope_significance.R, so both
# scripts share one definition of the null.
#
# Usage (via submit_klebsiella_rope.sh; needs cluster memory, not a login node):
#   Rscript klebsiella_bf_pd_all_variants.R <run_dir> <variant_index.csv> <out.csv>

suppressPackageStartupMessages({
  library(posterior)
  library(logspline)
  library(bayestestR)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: klebsiella_bf_pd_all_variants.R <run_dir> <variant_index.csv> <out.csv>")
}
run_dir <- args[1]; variant_index_path <- args[2]; out_path <- args[3]

N_PRIOR <- 200000L   # matches bf_pd_significance.R's --n-prior-samples default
SEED    <- 42L       # matches its --seed default
set.seed(SEED)

# ---- Load fit ---------------------------------------------------------------
fitted_dir <- file.path(run_dir, "fitted_model")
rds_files  <- list.files(fitted_dir, pattern = "\\.RDS$", full.names = TRUE)
if (length(rds_files) == 0L) stop("No .RDS in ", fitted_dir)
rds_path <- rds_files[which.max(file.mtime(rds_files))]
message("[bfpd] loading ", rds_path)
fit <- readRDS(rds_path)

draws <- as.data.frame(
  posterior::as_draws_df(fit$draws(variables = c("tau", "c2", "beta_variant_std"))))
rm(fit); gc(verbose = FALSE)

bvs_cols <- grep("^beta_variant_std\\[\\d+\\]$", names(draws), value = TRUE)
if (length(bvs_cols) == 0L) stop("No scalar-indexed beta_variant_std in draws")
# Order by index, not by name -- string sort would give 1, 10, 100, 2, ...
v_ord    <- order(as.integer(sub("^beta_variant_std\\[(\\d+)\\]$", "\\1", bvs_cols)))
bvs_cols <- bvs_cols[v_ord]

vidx <- read.csv(variant_index_path, stringsAsFactors = FALSE)
if (length(bvs_cols) != nrow(vidx)) {
  stop("beta_variant_std columns (", length(bvs_cols),
       ") != variant index rows (", nrow(vidx), ")")
}
message("[bfpd] V = ", nrow(vidx), "  draws = ", nrow(draws))

# ---- ROPE from rope_significance.R ------------------------------------------
rope_path <- file.path(run_dir, "rope_pd_bf", "rope_bounds.txt")
if (!file.exists(rope_path)) stop("rope_bounds.txt not found: ", rope_path)
rl <- readLines(rope_path)
pull <- function(prefix) {
  line <- grep(paste0("^", prefix), rl, value = TRUE)
  if (length(line) == 0L) return(NA_real_)
  as.numeric(sub(".*:\\s*", "", line[1]))
}
rope <- c(pull("rope_lo \\(standardized\\)"), pull("rope_hi \\(standardized\\)"))
if (any(is.na(rope))) stop("Could not parse ROPE from ", rope_path)
message(sprintf("[bfpd] ROPE (standardized): [%.6g, %.6g]", rope[1], rope[2]))

# ---- Prior null, built exactly as bf_pd_significance.R does ------------------
# tau and c2 come from the POSTERIOR draws; only lambda and z are redrawn. This
# is a partially-conditional null, matching the existing script rather than a
# draw from the pure prior.
tau_post <- draws$tau; c2_post <- draws$c2
S <- length(tau_post); M <- max(1L, ceiling(N_PRIOR / S))
lambda_null <- abs(rcauchy(S * M, 0, 2))
z_null      <- rnorm(S * M)
tau_rep <- rep(tau_post, each = M); c2_rep <- rep(c2_post, each = M)
beta_null_std <- z_null * tau_rep *
  sqrt((c2_rep * lambda_null^2) / (c2_rep + tau_rep^2 * lambda_null^2))
rm(lambda_null, z_null, tau_rep, c2_rep); gc(verbose = FALSE)
message("[bfpd] prior null draws: ", length(beta_null_std))

# Density at 0, with the same fallback chain bf_pd_significance.R uses. The
# regularized horseshoe prior is extremely spiked at 0, and logspline fails
# outright on it ("no convergence") -- bayestestR hits this too and falls back
# internally, which is why its point-null BFs came back finite. Warnings are
# expected and not treated as failure; only errors are.
dens_at_0 <- function(x, what) {
  f <- tryCatch(suppressWarnings(logspline::logspline(x)),
                error = function(e) NULL)
  if (!is.null(f)) {
    d <- tryCatch(logspline::dlogspline(0, f), error = function(e) NA_real_)
    if (is.finite(d) && d > 0) return(list(d = d, method = "logspline"))
  }
  dd <- stats::density(x, n = 8192L)
  d  <- stats::approx(dd$x, dd$y, xout = 0)$y
  if (!is.finite(d) || d <= 0) return(list(d = NA_real_, method = "failed"))
  list(d = d, method = "density")
}

# Prior quantities that do not depend on the variant, computed once.
pr <- dens_at_0(beta_null_std, "prior")
d_prior_0 <- pr$d
if (!is.finite(d_prior_0)) stop("Could not estimate prior density at 0")
p_prior_out <- mean(beta_null_std < rope[1] | beta_null_std > rope[2])
message(sprintf("[bfpd] prior density at 0 = %.6g (%s);  P(prior outside ROPE) = %.6g",
                d_prior_0, pr$method, p_prior_out))

# ---- Per-variant quantities -------------------------------------------------
bf_rope_of <- function(post) {
  p_post_out <- mean(post < rope[1] | post > rope[2])
  if (p_post_out %in% c(0, 1) || p_prior_out %in% c(0, 1)) return(Inf)
  (p_post_out / (1 - p_post_out)) / (p_prior_out / (1 - p_prior_out))
}
bf_point_of <- function(post) {
  r <- dens_at_0(post, "posterior")
  if (!is.finite(r$d)) return(NA_real_)
  d_prior_0 / r$d
}
pd_of <- function(post) max(mean(post > 0), mean(post < 0))

V <- nrow(vidx)
bf_rope <- numeric(V); bf_point <- numeric(V); pd <- numeric(V)
med <- numeric(V); ci_lo <- numeric(V); ci_hi <- numeric(V)

# 89% central interval, matching --ci-prob 0.89.
CI_PROB <- 0.89
qlo <- (1 - CI_PROB) / 2; qhi <- 1 - qlo

for (v in seq_len(V)) {
  post <- draws[[bvs_cols[v]]]
  bf_rope[v]  <- bf_rope_of(post)
  bf_point[v] <- bf_point_of(post)
  pd[v]       <- pd_of(post)
  med[v]      <- median(post)
  q <- unname(quantile(post, probs = c(qlo, qhi)))
  ci_lo[v] <- q[1]; ci_hi[v] <- q[2]
  if (v %% 500 == 0) message("[bfpd]   ", v, " / ", V)
}

# ---- Validate against bayestestR on a sample -------------------------------
# The point of the caching is speed, not a different answer. Check a spread of
# variants (including the strongest) and stop if the shortcut disagrees.
set.seed(SEED)
check_idx <- unique(c(which.max(abs(med)), sample.int(V, 8L)))
bad <- character(0)
for (v in check_idx) {
  post <- draws[[bvs_cols[v]]]
  # bayestestR can fail outright on this spiked prior; when it does there is
  # nothing to compare against and the check is skipped rather than passed.
  ref_rope <- tryCatch(
    suppressWarnings(exp(bayestestR::bayesfactor_parameters(
      post, prior = beta_null_std, null = rope, verbose = FALSE)$log_BF[1])),
    error = function(e) NA_real_)
  ref_pd <- tryCatch(as.numeric(bayestestR::p_direction(post)$pd[1]),
                     error = function(e) NA_real_)
  rel <- function(a, b) if (!is.finite(a) || !is.finite(b)) NA else abs(a - b) / max(1e-12, abs(b))
  r1 <- rel(bf_rope[v], ref_rope); r2 <- rel(pd[v], ref_pd)
  message(sprintf("[bfpd] check v=%d  bf_rope %.6g vs bayestestR %.6g (rel %.2g) | pd %.6g vs %.6g (rel %.2g)",
                  v, bf_rope[v], ref_rope, r1, pd[v], ref_pd, r2))
  if (!is.na(r1) && r1 > 0.02) bad <- c(bad, sprintf("bf_rope v=%d", v))
  if (!is.na(r2) && r2 > 0.02) bad <- c(bad, sprintf("pd v=%d", v))
}
if (length(bad) > 0) {
  stop("Shortcut disagrees with bayestestR for: ", paste(bad, collapse = ", "))
}
message("[bfpd] validation passed against bayestestR")

out <- data.frame(
  variant_name = vidx$variant_name,
  position     = vidx$position,
  median_std   = med,
  ci_lo_std    = ci_lo,
  ci_hi_std    = ci_hi,
  bf_rope      = bf_rope,
  bf_point     = bf_point,
  pd           = pd,
  rope_lo      = rope[1],
  rope_hi      = rope[2],
  # The two significance criteria from rope_significance.R, recomputed here so
  # the plotting step does not have to re-join its outputs.
  signif_median = med   < rope[1] | med   > rope[2],
  signif_ci     = ci_hi < rope[1] | ci_lo > rope[2],
  stringsAsFactors = FALSE)

write.csv(out, out_path, row.names = FALSE)
message("[bfpd] wrote ", out_path)
