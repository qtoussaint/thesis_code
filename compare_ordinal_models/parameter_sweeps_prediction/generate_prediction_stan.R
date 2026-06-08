#!/usr/bin/env Rscript
# Generate out-of-sample (training/test split) PREDICTION Stan models for the
# PPOM parameter-sweep variants (12 regularization-prior + 20 tau-sweep models),
# by mechanically transforming each inference .stan into its prediction twin.
#
# The inference -> prediction transform is the one verified by diffing the
# existing hand-written pair
#   final_ordered_categorical_PPOM_free_cutpoints_wide_drift{,_prediction}.stan:
#   * data:             N -> N_train + N_test; phenotype -> training_phenotype;
#                       variant_matrix/sublineage_matrix -> training_*/test_*;
#                       drop N_ppc / ppc_idx.
#   * transformed data: standardize on TRAINING mean/sd, apply to test
#                       (X_std_train / X_std_test); baseline over training only.
#   * parameters / transformed parameters: copied verbatim (this is where the
#                       ridge/lasso/est-sigma/ncp/fixed-tau prior differences live).
#   * model:            likelihood loop on training_* / X_std_train.
#   * generated quants: drop the PPC block; emit matrix[N_test,K] predicted_phenotype
#                       (Rao-Blackwellized); heritability over the training partition.
#
# Two archetypes: "centered" (standardized X_std, beta_variant_std) and
# "no_centering" (raw genotype, beta_variant_raw acting on variant_matrix). The
# large data / transformed-data blocks for the centered case are lifted byte-exact
# from the base inference/prediction pair at runtime; the rest are exact-string
# block replacements asserted to occur exactly once. A residual-token assertion
# (no bare N / N_ppc / variant_matrix / categorical_rng left) is the semantic net.
#
# Run with:
#   mamba activate gwas_pipeline   # (or conda activate gwas_pipeline)
#   Rscript generate_prediction_stan.R
#
# Validation: regenerates the existing wide_drift / tau5 prediction twins and
# checks the functional blocks match, then stanc-syntax-checks every output.

PPOM_DIR <- "/nfs/research/jlees/jacqueline/thesis_code/compare_ordinal_models/PPOM_models"
OUT_DIR  <- "/nfs/research/jlees/jacqueline/thesis_code/compare_ordinal_models/parameter_sweeps_prediction/stan_models"
STANC    <- "/homes/lilyjacqueline/.cmdstan/cmdstan-2.37.0/bin/stanc"
PFX      <- "final_ordered_categorical_PPOM_free_cutpoints_wide_drift"

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

read_file  <- function(p) paste(readLines(p, warn = FALSE), collapse = "\n")
write_file <- function(p, x) writeLines(x, p)

# Extract a top-level Stan block (`<name> {` .. lone `}`) byte-exact.
extract_block <- function(text, name) {
  lines <- strsplit(text, "\n", fixed = TRUE)[[1]]
  start <- which(grepl(paste0("^", name, " \\{"), lines))[1]
  if (is.na(start)) stop("block not found: ", name)
  end <- which(lines == "}" & seq_along(lines) > start)[1]
  paste(lines[start:end], collapse = "\n")
}

# Replace `find` with `repl` in `text`, asserting `find` occurs exactly `n` times.
sub1 <- function(text, find, repl, n = 1L, tag = "") {
  hits <- length(gregexpr(find, text, fixed = TRUE)[[1]])
  if (identical(gregexpr(find, text, fixed = TRUE)[[1]][1], -1L)) hits <- 0L
  if (hits != n)
    stop(sprintf("expected %d occurrence(s) of [%s] but found %d", n, tag, hits))
  gsub(find, repl, text, fixed = TRUE)
}

# ---------------------------------------------------------------------------
# Lift the large blocks byte-exact from the base inference/prediction pair.
# ---------------------------------------------------------------------------
inf_base  <- read_file(file.path(PPOM_DIR, paste0(PFX, ".stan")))
pred_base <- read_file(file.path(PPOM_DIR, paste0(PFX, "_prediction.stan")))

DATA_INF  <- extract_block(inf_base,  "data")   # comment-free, identical for all 32
DATA_PRED <- extract_block(pred_base, "data")

# transformed-data transform via targeted, comment-free code replacements (whole-
# block replacement is brittle because per-variant comments differ). These three
# code chunks are identical across variants of each archetype.

# Sublineage treatment-contrast split (both archetypes).
SUB_INF  <- "  matrix[N, S-1] X_sublineage = block(sublineage_matrix, 1, 2, N, S-1);"
SUB_PRED <- r"---(  matrix[N_train, S-1] X_sublineage_train =
    block(training_sublineages, 1, 2, N_train, S-1);
  matrix[N_test,  S-1] X_sublineage_test  =
    block(test_sublineages,     1, 2, N_test,  S-1);)---"

# Genotype standardization (centered archetype only).
XSTD_INF <- r"---(  matrix[N, V] X_std;
  vector[V] mean_variant;
  vector[V] sd_variant;
  for (v in 1:V) {
    mean_variant[v] = mean(variant_matrix[, v]);
    vector[N] centred = variant_matrix[, v] - mean_variant[v];
    sd_variant[v] = sd(centred);
    for (n in 1:N)
      X_std[n, v] = (sd_variant[v] > 0) ? centred[n] / sd_variant[v] : 0;
  })---"
XSTD_PRED <- r"---(  matrix[N_train, V] X_std_train;
  matrix[N_test,  V] X_std_test;
  vector[V] mean_variant;
  vector[V] sd_variant;
  for (v in 1:V) {
    mean_variant[v] = mean(training_variants[, v]);
    vector[N_train] centred_train = training_variants[, v] - mean_variant[v];
    sd_variant[v] = sd(centred_train);
    for (n in 1:N_train)
      X_std_train[n, v] = (sd_variant[v] > 0) ? centred_train[n] / sd_variant[v] : 0;
    for (n in 1:N_test)
      X_std_test[n, v] = (sd_variant[v] > 0)
        ? (test_variants[n, v] - mean_variant[v]) / sd_variant[v] : 0;
  })---"

# Reference-sublineage baseline tally (both archetypes).
BASE_INF <- r"---(  for (n in 1:N) {
    if (sublineage_matrix[n, 1] > 0.5) {
      n_ref += 1;
      if (phenotype[n] == 1) n_ref_cat1 += 1;
    }
  })---"
BASE_PRED <- r"---(  for (n in 1:N_train) {
    if (training_sublineages[n, 1] > 0.5) {
      n_ref += 1;
      if (training_phenotype[n] == 1) n_ref_cat1 += 1;
    }
  })---"

# ---------------------------------------------------------------------------
# Small sub-block find/replace pairs (asserted single-occurrence).
# ---------------------------------------------------------------------------
# Centered likelihood loop (model block).
LIK_C_INF <- r"---(  {
    vector[N] sub_eta = X_sublineage * beta_sublineage;
    for (n in 1:N) {
      for (k in 1:(K-1)) {
        real mu_nk = dot_product(X_std[n], beta_variant_std[, k]) + sub_eta[n];
        if (phenotype[n] <= k) {
          target += bernoulli_logit_lpmf(1 | cutpoints[k] - mu_nk);
        } else {
          target += bernoulli_logit_lpmf(0 | cutpoints[k] - mu_nk);
        }
      }
    }
  })---"
LIK_C_PRED <- r"---(  {
    vector[N_train] sub_eta = X_sublineage_train * beta_sublineage;
    for (n in 1:N_train) {
      for (k in 1:(K-1)) {
        real mu_nk = dot_product(X_std_train[n], beta_variant_std[, k]) + sub_eta[n];
        if (training_phenotype[n] <= k) {
          target += bernoulli_logit_lpmf(1 | cutpoints[k] - mu_nk);
        } else {
          target += bernoulli_logit_lpmf(0 | cutpoints[k] - mu_nk);
        }
      }
    }
  })---"

# No_centering likelihood loop (raw genotype).
LIK_NC_INF <- r"---(  {
    vector[N] sub_eta = X_sublineage * beta_sublineage;
    for (n in 1:N) {
      for (k in 1:(K-1)) {
        real mu_nk = dot_product(variant_matrix[n], beta_variant_raw[, k]) + sub_eta[n];
        if (phenotype[n] <= k) {
          target += bernoulli_logit_lpmf(1 | cutpoints[k] - mu_nk);
        } else {
          target += bernoulli_logit_lpmf(0 | cutpoints[k] - mu_nk);
        }
      }
    }
  })---"
LIK_NC_PRED <- r"---(  {
    vector[N_train] sub_eta = X_sublineage_train * beta_sublineage;
    for (n in 1:N_train) {
      for (k in 1:(K-1)) {
        real mu_nk = dot_product(training_variants[n], beta_variant_raw[, k]) + sub_eta[n];
        if (training_phenotype[n] <= k) {
          target += bernoulli_logit_lpmf(1 | cutpoints[k] - mu_nk);
        } else {
          target += bernoulli_logit_lpmf(0 | cutpoints[k] - mu_nk);
        }
      }
    }
  })---"

# Centered PPC block (generated quantities) -> predicted_phenotype.
PPC_C_INF <- r"---(  array[N_ppc] int y_rep_ppc;
  array[N_ppc] int y_true_ppc;
  {
    vector[N] sub_eta_gen = X_sublineage * beta_sublineage;
    for (i in 1:N_ppc) {
      int n = ppc_idx[i];
      vector[K-1] cdf_n;
      for (k in 1:(K-1)) {
        real mu_nk = dot_product(X_std[n], beta_variant_std[, k]) + sub_eta_gen[n];
        cdf_n[k] = inv_logit(cutpoints[k] - mu_nk);
      }
      vector[K] probs;
      probs[1] = cdf_n[1];
      for (k in 2:(K-1)) probs[k] = cdf_n[k] - cdf_n[k-1];
      probs[K] = 1 - cdf_n[K-1];
      for (k in 1:K) if (probs[k] < 0) probs[k] = 0;
      probs /= sum(probs);
      y_rep_ppc[i]  = categorical_rng(probs);
      y_true_ppc[i] = phenotype[n];
    }
  })---"
PRED_PHENO_C <- r"---(  matrix[N_test, K] predicted_phenotype;
  {
    vector[N_test] sub_eta_test = X_sublineage_test * beta_sublineage;
    for (n in 1:N_test) {
      vector[K-1] cdf_n;
      for (k in 1:(K-1)) {
        real mu_nk = dot_product(X_std_test[n], beta_variant_std[, k]) + sub_eta_test[n];
        cdf_n[k] = inv_logit(cutpoints[k] - mu_nk);
      }
      vector[K] probs;
      probs[1] = cdf_n[1];
      for (k in 2:(K-1)) probs[k] = cdf_n[k] - cdf_n[k-1];
      probs[K] = 1 - cdf_n[K-1];
      for (k in 1:K) if (probs[k] < 1e-12) probs[k] = 1e-12;
      probs /= sum(probs);
      for (k in 1:K) predicted_phenotype[n, k] = probs[k];
    }
  })---"

# No_centering PPC block -> predicted_phenotype (raw genotype on test).
PPC_NC_INF <- r"---(  array[N_ppc] int y_rep_ppc;
  array[N_ppc] int y_true_ppc;
  {
    vector[N] sub_eta_gen = X_sublineage * beta_sublineage;
    for (i in 1:N_ppc) {
      int n = ppc_idx[i];
      vector[K-1] cdf_n;
      for (k in 1:(K-1)) {
        real mu_nk = dot_product(variant_matrix[n], beta_variant_raw[, k]) + sub_eta_gen[n];
        cdf_n[k] = inv_logit(cutpoints[k] - mu_nk);
      }
      vector[K] probs;
      probs[1] = cdf_n[1];
      for (k in 2:(K-1)) probs[k] = cdf_n[k] - cdf_n[k-1];
      probs[K] = 1 - cdf_n[K-1];
      for (k in 1:K) if (probs[k] < 0) probs[k] = 0;
      probs /= sum(probs);
      y_rep_ppc[i]  = categorical_rng(probs);
      y_true_ppc[i] = phenotype[n];
    }
  })---"
PRED_PHENO_NC <- r"---(  matrix[N_test, K] predicted_phenotype;
  {
    vector[N_test] sub_eta_test = X_sublineage_test * beta_sublineage;
    for (n in 1:N_test) {
      vector[K-1] cdf_n;
      for (k in 1:(K-1)) {
        real mu_nk = dot_product(test_variants[n], beta_variant_raw[, k]) + sub_eta_test[n];
        cdf_n[k] = inv_logit(cutpoints[k] - mu_nk);
      }
      vector[K] probs;
      probs[1] = cdf_n[1];
      for (k in 2:(K-1)) probs[k] = cdf_n[k] - cdf_n[k-1];
      probs[K] = 1 - cdf_n[K-1];
      for (k in 1:K) if (probs[k] < 1e-12) probs[k] = 1e-12;
      probs /= sum(probs);
      for (k in 1:K) predicted_phenotype[n, k] = probs[k];
    }
  })---"

# Heritability lines (generated quantities).
HER_C_INF   <- "    vector[N] g_pop = X_sublineage * beta_sublineage;"
HER_C_PRED  <- "    vector[N_train] g_pop = X_sublineage_train * beta_sublineage;"
HER_CV_INF  <- "      vector[N] g_variant_k = X_std * beta_variant_std[, k];"
HER_CV_PRED <- "      vector[N_train] g_variant_k = X_std_train * beta_variant_std[, k];"
HER_NCV_INF  <- "      vector[N] g_variant_k = variant_matrix * beta_variant_raw[, k];"
HER_NCV_PRED <- "      vector[N_train] g_variant_k = training_variants * beta_variant_raw[, k];"

# ---------------------------------------------------------------------------
# Transform one inference model into its prediction twin.
# ---------------------------------------------------------------------------
header <- function(name) sprintf(paste0(
  "// AUTO-GENERATED prediction (training/test split) twin of\n",
  "// %s.stan\n",
  "// Produced by parameter_sweeps_prediction/generate_prediction_stan.R; do not\n",
  "// edit by hand. Fits on N_train and emits an N_test x K predicted_phenotype\n",
  "// matrix (Rao-Blackwellized) for out-of-sample prediction scoring.\n"), name)

transform_stan <- function(text, name) {
  centered <- grepl("matrix[N, V] X_std;", text, fixed = TRUE)

  t <- sub1(text, DATA_INF, DATA_PRED, tag = "data")
  t <- sub1(t, SUB_INF,  SUB_PRED,  tag = "sublineage split")
  t <- sub1(t, BASE_INF, BASE_PRED, tag = "baseline tally")
  if (centered) {
    t <- sub1(t, XSTD_INF, XSTD_PRED, tag = "standardization (centered)")
    t <- sub1(t, LIK_C_INF, LIK_C_PRED, tag = "likelihood (centered)")
    t <- sub1(t, PPC_C_INF, PRED_PHENO_C, tag = "ppc->pred (centered)")
    t <- sub1(t, HER_C_INF,  HER_C_PRED,  tag = "her g_pop")
    t <- sub1(t, HER_CV_INF, HER_CV_PRED, tag = "her g_variant (centered)")
  } else {
    t <- sub1(t, LIK_NC_INF, LIK_NC_PRED, tag = "likelihood (no_centering)")
    t <- sub1(t, PPC_NC_INF, PRED_PHENO_NC, tag = "ppc->pred (no_centering)")
    t <- sub1(t, HER_C_INF,   HER_C_PRED,   tag = "her g_pop")
    t <- sub1(t, HER_NCV_INF, HER_NCV_PRED, tag = "her g_variant (no_centering)")
  }

  # Replace the preamble comment with a generated header.
  anchor <- regexpr("data {\n  int<lower=1> N_train;", t, fixed = TRUE)
  if (anchor < 0) stop("could not locate transformed data block start in ", name)
  t <- paste0(header(name), "\n", substring(t, anchor))

  # Residual-token safety net: nothing in-sample should survive the transform.
  bad <- c("\\bN\\b", "\\bN_ppc\\b", "\\bppc_idx\\b", "\\bvariant_matrix\\b",
           "\\bsublineage_matrix\\b", "\\bcategorical_rng\\b", "y_rep_ppc",
           "y_true_ppc", "\\bphenotype\\b")
  for (b in bad)
    if (grepl(b, t, perl = TRUE))
      stop(sprintf("residual in-sample token [%s] survived transform of %s", b, name))
  t
}

# ---------------------------------------------------------------------------
# Input lists.
# ---------------------------------------------------------------------------
REG_VARIANTS <- c(
  "lasso", "lasso_no_centering", "lasso_estscale", "lasso_estscale_no_centering",
  "lasso_estscale_mixture", "lasso_estscale_mixture_no_centering",
  "ridge", "ridge_no_centering", "ridge_estscale", "ridge_estscale_no_centering",
  "ridge_estscale_ncp", "ridge_estscale_ncp_no_centering")

# Tau sweep: the 20 (series x tau) the tau figure uses. Excludes the stray
# fixedtau05_slab3 (no figure slot).
TAU_TOK <- c("0p001", "0p01", "0p05", "1", "1p5", "2", "3", "5")
TAU_VARIANTS <- c(
  paste0("fixedtau", c("0p001", "0p01", "0p05", "1"), "_slab3"),
  paste0("fixedtau", c("0p001", "0p01", "0p05", "1"), "_slab5"),
  paste0("fixedtau", TAU_TOK, "_slab3_lambda2"),
  paste0("fixedtau", c("1p5", "2", "3", "5"), "_slab5_lambda2"))

ALL_VARIANTS <- c(REG_VARIANTS, TAU_VARIANTS)

# ---------------------------------------------------------------------------
# Validation: regenerate the existing hand-written prediction twins and confirm
# the functional blocks match (ignoring header comments + extra diagnostics).
# ---------------------------------------------------------------------------
# Normalise for structural comparison: drop comments, blank lines, and the
# beta_variant_std_prior diagnostic block that the hand-written twins add but the
# generator (rightly) does not, then trim each line.
normalise <- function(x) {
  # Drop the beta_variant_std_prior diagnostic block the hand-written twins add
  # (loop headers in it carry no token, so strip the whole block first).
  x <- gsub("(?s)\\n\\s*matrix\\[V, K-1\\] beta_variant_std_prior;.*?lambda_tilde_variant\\[v, k\\];",
            "", x, perl = TRUE)
  lines <- strsplit(x, "\n", fixed = TRUE)[[1]]
  lines <- sub("//.*$", "", lines)                    # strip inline + full-line comments
  lines <- trimws(lines)
  lines <- lines[nzchar(lines)]
  paste(lines, collapse = "\n")
}

validate_against <- function(variant) {
  gen <- transform_stan(read_file(file.path(PPOM_DIR, paste0(PFX, "_", variant, ".stan"))),
                        paste0(PFX, "_", variant))
  ref <- read_file(file.path(PPOM_DIR, paste0(PFX, "_", variant, "_prediction.stan")))
  ng <- normalise(gen)
  nr <- normalise(ref)
  if (!identical(ng, nr)) {
    # Report first differing line to aid debugging.
    lg <- strsplit(ng, "\n")[[1]]; lr <- strsplit(nr, "\n")[[1]]
    d <- which(lg[seq_len(min(length(lg), length(lr)))] !=
               lr[seq_len(min(length(lg), length(lr)))])[1]
    stop(sprintf("VALIDATION FAILED for %s at line %s:\n  gen: %s\n  ref: %s",
                 variant, d, if (is.na(d)) "(length)" else lg[d],
                 if (is.na(d)) "(length)" else lr[d]))
  }
  message("  validation OK: regenerated ", variant, "_prediction matches hand-written twin")
}

message("== Validation against existing hand-written prediction twins ==")
validate_against("tau5")   # centered, default slab; existing twin

# ---------------------------------------------------------------------------
# Generate all variants.
# ---------------------------------------------------------------------------
message("== Generating ", length(ALL_VARIANTS), " PPOM prediction models ==")
written <- character(0)
for (v in ALL_VARIANTS) {
  src <- file.path(PPOM_DIR, paste0(PFX, "_", v, ".stan"))
  if (!file.exists(src)) stop("missing input: ", src)
  out_name <- paste0(PFX, "_", v, "_prediction.stan")
  out <- file.path(OUT_DIR, out_name)
  write_file(out, transform_stan(read_file(src), paste0(PFX, "_", v)))
  written <- c(written, out)
  message("  wrote ", out_name)
}

# ---------------------------------------------------------------------------
# stanc syntax check on every generated model.
# ---------------------------------------------------------------------------
message("== stanc syntax check ==")
fail <- character(0)
for (f in written) {
  res <- suppressWarnings(system2(STANC, c("--warn-pedantic", shQuote(f)),
                                  stdout = TRUE, stderr = TRUE))
  status <- attr(res, "status")
  if (!is.null(status) && status != 0) {
    fail <- c(fail, f)
    message("  FAIL ", basename(f))
    message(paste0("    ", res, collapse = "\n"))
  }
}
if (length(fail) > 0) stop(length(fail), " model(s) failed stanc syntax check")
message("All ", length(written), " PPOM prediction models passed stanc.")
