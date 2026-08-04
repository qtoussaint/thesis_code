#!/usr/bin/env Rscript
# Extract per-variant effects and Bayes factors from the simulated Klebsiella
# fits, once per run, into small CSVs the figure scripts can read cheaply.
#
# The runs use --ld_pruning false (see generate_klebsiella_run_scripts.R for
# why), so there is no depruned_variant_effects.csv and no de-pruning mapping:
# beta_variant columns correspond 1:1 to the rows of the dataset variant index.
#
# Each fitted RDS holds 10,000 draws x 4,578 variants across several parameter
# blocks and is several GB, so it is loaded exactly once here and never by the
# plotting scripts.
#
# BAYES FACTORS. The prior is NOT the N(0,5) used by the older
# bayesian_bridge.R -- that matched BBGWAS-cont-association.stan, whereas
# continuous_inference.stan puts a regularized horseshoe on the standardized
# effects. continuous_inference.stan does declare beta_variant_std_prior in
# generated quantities "to expose priors for post-hoc verification", but it
# evaluates z_variant[v] * tau * lambda_tilde_variant[v], the exact expression
# used for beta_variant_std in transformed parameters -- so it is a copy of the
# posterior, not a prior draw. That is checked against the real draws below.
# The prior is therefore simulated from the model's own generative structure.
#
# The null is a ROPE of +/- 0.1 * sd(phenotype) per dataset (the Kruschke /
# bayestestR convention) rather than the fixed +/- 0.1 of the older script,
# because these four phenotypes differ in scale by a factor of ~20.
#
# Usage:
#   mamba activate gwas_pipeline
#   Rscript paper_figures/klebsiella_extract_effects.R

suppressPackageStartupMessages({
  library(data.table)
  library(cmdstanr)
  library(jsonlite)
  library(dplyr)
})

# De-pruning helpers for the LD-pruned arm: read_pruned_variants(),
# parse_ld_summary() and depruning_rates(). Sourced from the pipeline package so
# the mapping from representative to pruned variant is identical to the one the
# pipeline applies itself, rather than reimplemented here.
PKG <- "/nfs/research/jlees/jacqueline/gwas_workflow/code/gwas_workflow/R"
invisible(lapply(file.path(PKG, c("ld_pruning.R", "cpprate.R")), source))

RES      <- "/nfs/research/jlees/jacqueline/thesis_results"
DATASETS <- file.path(RES, "gwas_datasets", "inference")
RUNS     <- file.path(RES, "gwas_klebsiella_homoplasic", "inference")
OUT_DIR  <- file.path(RES, "gwas_klebsiella_homoplasic", "extracted")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Hyperparameters as hardcoded in continuous_inference.stan.
SLAB_SCALE <- 200; NU <- 4; TAU_0 <- 0.1
N_PRIOR    <- 40000L
BF_CORRECTION <- 0.5   # continuity correction, keeps BF finite

truth <- read.csv(file.path(DATASETS, "klebsiella_homoplasic_truth_summary.csv"),
                  stringsAsFactors = FALSE)

simulate_beta_std_prior <- function(n, seed = 20260804L) {
  set.seed(seed)
  tau    <- abs(rcauchy(n, 0, TAU_0))
  lambda <- abs(rcauchy(n, 0, 2))
  c2     <- 1 / rgamma(n, shape = 0.5 * NU, rate = 0.5 * NU * SLAB_SCALE^2)
  z      <- rnorm(n)
  z * tau * sqrt((c2 * lambda^2) / (c2 + tau^2 * lambda^2))
}

# Interval-null Bayes factor: ratio of posterior to prior odds of lying outside
# the ROPE. This is exactly what bayestestR::bayesfactor_parameters() computes
# for an interval null -- no density estimation involved, so logspline and
# bayestestR (absent from this environment) are not needed. Only a POINT null
# would require the Savage-Dickey density ratio.
bf_interval <- function(post, prior, rope) {
  n1 <- length(post); n0 <- length(prior)
  i1 <- sum(post  >= rope[1] & post  <= rope[2])
  i0 <- sum(prior >= rope[1] & prior <= rope[2])
  ((n1 - i1 + BF_CORRECTION) / (i1 + BF_CORRECTION)) /
  ((n0 - i0 + BF_CORRECTION) / (i0 + BF_CORRECTION))
}

summaries <- list()

# Both LD-pruning arms, matching the nicknames emitted by
# gwas_finalruns/generate_klebsiella_run_scripts.R.
MODES <- c("pruned", "nopruning")

grid <- expand.grid(row = seq_len(nrow(truth)), mode = MODES,
                    stringsAsFactors = FALSE)

for (g in seq_len(nrow(grid))) {
  i        <- grid$row[g]
  mode     <- grid$mode[g]
  ds       <- truth$dataset[i]
  ef       <- sub("\\.0$", "", format(truth$effect_size[i], trim = TRUE))
  nickname <- paste0(ds, "_continuous_", mode)
  run_dir  <- file.path(RUNS, nickname)
  message("\n=== ", nickname, " ===")

  rds <- Sys.glob(file.path(run_dir, "fitted_model", paste0(nickname, "_*.RDS")))
  if (length(rds) == 0L) { warning("No fit for ", nickname, "; skipping", call. = FALSE); next }
  fit <- readRDS(rds[1])

  # Evidence for the beta_variant_std_prior duplication, checked not assumed.
  a <- tryCatch(as.numeric(fit$draws("beta_variant_std[1]", format = "draws_matrix")), error = function(e) NULL)
  b <- tryCatch(as.numeric(fit$draws("beta_variant_std_prior[1]", format = "draws_matrix")), error = function(e) NULL)
  if (!is.null(a) && !is.null(b)) {
    message("  beta_variant_std_prior == posterior beta_variant_std: ",
            isTRUE(all.equal(a, b)), " (max abs diff ", signif(max(abs(a - b)), 3), ")")
  }

  vi <- data.table::fread(file.path(DATASETS, ds, paste0(ds, "_variant_index.csv")))
  post <- fit$draws("beta_variant", format = "draws_matrix")

  # In the pruned arm the fit covers only the kept variants, so every per-variant
  # quantity has to be mapped back onto the full variant set before the two arms
  # can be compared or plotted on the same axis. Each pruned variant inherits its
  # representative's value, which is what the pipeline's own de-pruning does.
  # In the no-pruning arm the columns already line up 1:1.
  ld_csv <- file.path(run_dir, "cppRATE_matrices", "ld_pruning_summary.csv")
  is_pruned_run <- file.exists(ld_csv)

  if (is_pruned_run) {
    kept <- read_pruned_variants(run_dir)
    if (ncol(post) != length(kept)) {
      stop("beta_variant columns (", ncol(post), ") != kept variants (", length(kept), ")")
    }
    ld_summary <- parse_ld_summary(ld_csv)
    expand <- function(vals) {
      full <- depruning_rates(data.frame(id = kept, v = vals), ld_summary, kept)
      as.numeric(full[[2]])
    }
    message("  LD-pruned run: ", length(kept), " kept variants expanded to ", nrow(vi))
  } else {
    if (ncol(post) != nrow(vi)) {
      stop("beta_variant columns (", ncol(post), ") != variant index rows (", nrow(vi), ")")
    }
    kept   <- seq_len(nrow(vi))
    expand <- function(vals) vals
    message("  no-pruning run: ", ncol(post), " variants, 1:1 with the index")
  }

  med <- expand(apply(post, 2, median))
  q05 <- expand(apply(post, 2, quantile, probs = 0.05))
  q95 <- expand(apply(post, 2, quantile, probs = 0.95))

  # Genotype sd per variant, matching Stan's sample sd, to put the standardized
  # prior on the same raw scale as beta_variant. Subset to the fitted columns so
  # it aligns with the posterior in both arms.
  vm <- jsonlite::fromJSON(file.path(DATASETS, ds, paste0(ds, ".json")))$variant_matrix
  sd_all     <- apply(vm, 2, sd); rm(vm); gc()
  sd_variant <- sd_all[kept]

  rope_half <- 0.1 * truth$pheno_sd[i]
  prior_std <- simulate_beta_std_prior(N_PRIOR)
  message("  ROPE = +/-", signif(rope_half, 4), " (0.1 x phenotype sd)")

  bf_fitted <- numeric(ncol(post))
  for (v in seq_len(ncol(post))) {
    pv <- if (sd_variant[v] > 0) prior_std / sd_variant[v] else prior_std
    bf_fitted[v] <- bf_interval(as.numeric(post[, v]), pv, c(-rope_half, rope_half))
  }
  bf <- expand(bf_fitted)
  rm(post); gc()

  dt <- data.table(dataset = ds, mode = mode, ef_chr = ef,
                   variant_name = vi$variant_name, position = vi$position,
                   median = med, q05 = q05, q95 = q95, bf = bf,
                   allele_sd = sd_all)
  dt[, is_causal := variant_name == truth$causal_variant[i]]
  data.table::fwrite(dt, file.path(OUT_DIR, paste0(ds, "_", mode, "_variant_effects.csv")))

  h2 <- tryCatch(fit$summary("h2_narrow", "median")$median, error = function(e) NA_real_)
  sg <- tryCatch(fit$summary("sigma", "median")$median,     error = function(e) NA_real_)
  cz <- dt[is_causal == TRUE]

  summaries[[length(summaries) + 1L]] <- data.table(
    dataset = ds, mode = mode, ef_chr = ef, effect_size = truth$effect_size[i],
    n_fitted = length(kept), n_total = nrow(vi),
    causal_median = cz$median[1], simulated_effect = truth$observed_diff[i],
    ratio = cz$median[1] / truth$observed_diff[i],
    causal_rank = rank(-abs(dt$median), ties.method = "min")[which(dt$is_causal)][1],
    causal_bf = cz$bf[1],
    bf_rank = rank(-dt$bf, ties.method = "min")[which(dt$is_causal)][1],
    n_bf_gt10 = sum(dt$bf > 10, na.rm = TRUE),
    n_variants = nrow(dt),
    h2_narrow = h2, sigma = sg,
    true_residual_sd = truth$pheno_sd[i] * sqrt(1 - truth$heritability[i]),
    rope_half = rope_half,
    # Validity gate. Every phenotype was simulated at h2 = 0.9, so a fit that
    # recovers almost no variant-explained variance has landed in a bad ADVI
    # optimum regardless of what the ELBO reported -- observed failures come
    # back at h2 ~1e-06 to 1e-08, not marginally low, so the threshold is not
    # delicate. Panels failing this must not be plotted as if they were results.
    valid = !is.na(h2) && h2 >= 0.5 * truth$heritability[i])
  message("  causal median = ", signif(cz$median[1], 5),
          " (true ", signif(truth$observed_diff[i], 5), ")",
          "  rank ", summaries[[length(summaries)]]$causal_rank,
          "  BF = ", signif(cz$bf[1], 4))
}

if (length(summaries) == 0L) stop("No fits found.")
s <- data.table::rbindlist(summaries)
data.table::fwrite(s, file.path(OUT_DIR, "klebsiella_extracted_summary.csv"))
message("\nWrote extracts to ", OUT_DIR)
print(s)

# Say plainly which panels are unusable, so a bad optimum cannot slip into a
# figure unnoticed.
bad <- s[valid == FALSE]
if (nrow(bad) > 0) {
  message("\n!! ", nrow(bad), " run(s) failed the h2 validity gate and must not be plotted as results:")
  for (k in seq_len(nrow(bad))) {
    message("   ", bad$dataset[k], " [", bad$mode[k], "]  h2_narrow = ",
            signif(bad$h2_narrow[k], 3), " (simulated 0.9)")
  }
}
expected <- nrow(truth) * length(MODES)
if (nrow(s) < expected) {
  message("\n!! only ", nrow(s), " of ", expected,
          " runs produced a fit at all; the rest failed outright")
}
