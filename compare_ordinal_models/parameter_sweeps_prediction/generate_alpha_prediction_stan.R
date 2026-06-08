#!/usr/bin/env Rscript
# Generate out-of-sample PREDICTION Stan models for the logistic + continuous
# alpha-prior-SD sweeps.
#
# The sweep variants differ from their baseline by ONLY the `alpha_prior_sd`
# literal (verified: diff vs gwas_finalmodels/{logistic,continuous}_inference.stan
# is just that line + a comment). Since hand-written, correct train/test
# prediction twins already exist in gwas_finalmodels, each prediction variant is
# that twin with the one literal set to the swept value.
#
# Run with:
#   mamba activate gwas_pipeline
#   Rscript generate_alpha_prediction_stan.R

FM     <- "/nfs/research/jlees/jacqueline/thesis_code/gwas_finalmodels"
OUT_DIR <- "/nfs/research/jlees/jacqueline/thesis_code/compare_ordinal_models/parameter_sweeps_prediction/stan_models"
STANC  <- "/homes/lilyjacqueline/.cmdstan/cmdstan-2.37.0/bin/stanc"

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
read_file <- function(p) paste(readLines(p, warn = FALSE), collapse = "\n")

# token -> numeric literal as it should appear in the .stan
SD <- c("0p075" = "0.075", "1" = "1", "1p25" = "1.25", "1p5" = "1.5",
        "2" = "2", "2p5" = "2.5", "3" = "3", "5" = "5")

FAMILIES <- list(
  list(name = "logistic",   base = file.path(FM, "logistic_prediction.stan"),
       find = "real alpha_prior_sd   = 1.5;",  # note: aligned spacing in this file
       tokens = c("0p075", "1", "1p25", "1p5", "2", "2p5", "3", "5")),
  list(name = "continuous", base = file.path(FM, "continuous_prediction.stan"),
       find = "real alpha_prior_sd = 0.5;",
       tokens = c("1p5", "3", "5"))
)

header <- function(fam, val) sprintf(paste0(
  "// AUTO-GENERATED alpha-SD sweep prediction model: %s, alpha_prior_sd = %s.\n",
  "// Produced by parameter_sweeps_prediction/generate_alpha_prediction_stan.R\n",
  "// from gwas_finalmodels/%s_prediction.stan; do not edit by hand.\n"),
  fam, val, fam)

written <- character(0)
for (fam in FAMILIES) {
  base_txt <- read_file(fam$base)
  hits <- length(gregexpr(fam$find, base_txt, fixed = TRUE)[[1]])
  if (gregexpr(fam$find, base_txt, fixed = TRUE)[[1]][1] == -1L) hits <- 0L
  if (hits != 1L)
    stop(sprintf("[%s] expected exactly 1 occurrence of '%s', found %d",
                 fam$name, fam$find, hits))
  for (tok in fam$tokens) {
    val  <- SD[[tok]]
    repl <- sub("= [0-9.]+;", paste0("= ", val, ";"), fam$find)
    body <- gsub(fam$find, repl, base_txt, fixed = TRUE)
    body <- sub("(?s)^.*?\\ndata \\{", "data {", body, perl = TRUE)  # drop preamble comment
    out_txt <- paste0(header(fam$name, val), "\n", body)
    out <- file.path(OUT_DIR, sprintf("%s_prediction_alphasd%s.stan", fam$name, tok))
    writeLines(out_txt, out)
    written <- c(written, out)
    message("  wrote ", basename(out), "  (alpha_prior_sd = ", val, ")")
  }
}

message("== stanc syntax check ==")
fail <- character(0)
for (f in written) {
  res <- suppressWarnings(system2(STANC, c("--warn-pedantic", shQuote(f)),
                                  stdout = TRUE, stderr = TRUE))
  st <- attr(res, "status")
  if (!is.null(st) && st != 0) { fail <- c(fail, f); message("  FAIL ", basename(f), "\n",
                                                             paste0("    ", res, collapse = "\n")) }
}
if (length(fail) > 0) stop(length(fail), " model(s) failed stanc syntax check")
message("All ", length(written), " alpha-sweep prediction models passed stanc.")
