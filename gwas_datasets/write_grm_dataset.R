############################################################
## write_grm_dataset.R
##
## Adds a genetic relatedness matrix (GRM) to an existing binary inference
## dataset:
##
##   19 SPN penicillin binary  (from 01)
##   20 TB  rifampicin binary  (from 07)
##
## This is the input for gwas_finalmodels/logistic_grm_inference.stan, the GRM
## twin of logistic_inference.stan. The two models are compared to check that
## the lineage/subcluster encoding of population structure gives the same
## variant effects and RATE values as a continuous GRM.
##
## The dataset is derived from the existing dataset rather than rebuilt from
## the raw genotype, so phenotype, variant_matrix, sample order and the PPC
## subset are identical by construction and the two GWAS runs differ only by
## the GRM.
##
## Usage:
##   Rscript write_grm_dataset.R [source_dataset] [output_dataset]
## Defaults to the SPN penicillin pair.
##
## Kept separate from write_inference_jsons.R because that script rebuilds all
## 18 datasets in one long job.
############################################################

script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("--file=", script_args, value = TRUE)
script_dir <- if (length(file_arg) > 0) dirname(normalizePath(sub("--file=", "", file_arg))) else "."
source(file.path(script_dir, "config.R"))
source(file.path(script_dir, "utils.R"))

cli <- commandArgs(trailingOnly = TRUE)
SOURCE_DATASET <- if (length(cli) >= 1) cli[1] else "01_spn_penicillin_binary"
DATASET_NAME   <- if (length(cli) >= 2) cli[2] else "19_spn_penicillin_binary_grm"
# Optional: round the GRM to this many significant figures to shrink the JSON.
# NA writes full precision.
SIGFIGS <- if (length(cli) >= 3) as.integer(cli[3]) else NA_integer_

# Rounding perturbs entries by up to 10^-SIGFIGS relative, which for a matrix
# scaled to unit mean diagonal swamps a 1e-6 jitter, so the jitter needed to
# keep the matrix positive definite is raised alongside.
JITTER <- if (is.na(SIGFIGS)) 1e-6 else 1e-3

src_dir <- file.path(OUT_INFER, SOURCE_DATASET)
out_dir <- file.path(OUT_INFER, DATASET_NAME)

src_json <- file.path(src_dir, paste0(SOURCE_DATASET, ".json"))
stopifnot(file.exists(src_json))

message("\n=== ", DATASET_NAME, " ===")
message("Reading source dataset: ", src_json)
d <- jsonlite::fromJSON(src_json)

message("  N = ", d$N, "  V = ", d$V, "  L = ", d$L, "  S = ", d$S)
stopifnot(nrow(d$variant_matrix) == d$N, ncol(d$variant_matrix) == d$V)


############################################################
## Build the GRM
############################################################

# VanRaden-style GRM on haploid 0/1 genotypes: standardise each variant by its
# allele frequency, then K = ZZ'/V. Monomorphic columns are dropped because
# their standardisation divides by zero. K is normalised by its mean diagonal
# so that sigma_g in the Stan model is on the same scale as the sublineage
# effect SD in logistic_inference.stan, and the diagonal is jittered so
# cholesky_decompose() succeeds.

# Each full copy of the genotype matrix as doubles is ~7 GB for the TB
# dataset, so intermediates are dropped as soon as they are consumed.
message("Building GRM...")
X <- d$variant_matrix
storage.mode(X) <- "double"

p    <- colMeans(X)
keep <- which(p > 0 & p < 1)
message("  Variants used for the GRM: ", length(keep), " of ", ncol(X),
        " (", ncol(X) - length(keep), " monomorphic dropped)")

Z <- if (length(keep) < ncol(X)) X[, keep, drop = FALSE] else X
rm(X); invisible(gc())

Z <- sweep(Z, 2, p[keep], "-")
Z <- sweep(Z, 2, sqrt(p[keep] * (1 - p[keep])), "/")

message("  Cross-product over ", nrow(Z), " samples x ", ncol(Z), " variants...")
K <- tcrossprod(Z) / ncol(Z)
rm(Z); invisible(gc())

K <- (K + t(K)) / 2
K <- K / mean(diag(K))

if (!is.na(SIGFIGS)) {
  message("  Rounding GRM to ", SIGFIGS, " significant figures")
  K <- signif(K, SIGFIGS)
  K <- (K + t(K)) / 2   # signif is elementwise, but re-symmetrise defensively
}

message("  Diagonal jitter: ", format(JITTER, scientific = TRUE))
diag(K) <- diag(K) + JITTER

message("  diag range:     ", paste(signif(range(diag(K)), 4), collapse = " to "))
message("  off-diag range: ", paste(signif(range(K[upper.tri(K)]), 4), collapse = " to "))

stopifnot(nrow(K) == d$N, ncol(K) == d$N)
stopifnot(isSymmetric(K))
chol_ok <- tryCatch({ chol(K); TRUE }, error = function(e) FALSE)
if (!chol_ok) stop("GRM is not positive definite; cholesky_decompose() would fail in Stan")
message("  Cholesky check: OK")


############################################################
## Write the dataset
############################################################

# Insert GRM after sublineage_matrix (field order is cosmetic for CmdStan).
stan_list <- append(d, list(GRM = K), after = which(names(d) == "sublineage_matrix"))
stopifnot(identical(names(d), setdiff(names(stan_list), "GRM")))

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

json_path <- file.path(out_dir, paste0(DATASET_NAME, ".json"))
# digits = NA keeps the GRM at full precision; the jsonlite default of 4 dp
# would round most off-diagonal entries away. When SIGFIGS is set, I(SIGFIGS)
# tells jsonlite to emit that many significant figures, which is what actually
# shrinks the file (signif() alone leaves values like 0.7498 printing as
# 0.7498000000000000265).
write_stan_json_streaming(data = stan_list, file = json_path,
                          digits = if (is.na(SIGFIGS)) NA else I(SIGFIGS))
message("  Wrote JSON: ", json_path)

# Copy the sidecars from dataset 01 so the variant/sample indices are identical
# to the reference run. The variant index is what the pipeline reads for
# --phandango.
for (suffix in c("_variant_index.csv", "_sample_index.csv", "_parent_lineages.csv")) {
  from <- file.path(src_dir, paste0(SOURCE_DATASET, suffix))
  to   <- file.path(out_dir, paste0(DATASET_NAME, suffix))
  stopifnot(file.exists(from))
  ok <- file.copy(from, to, overwrite = TRUE)
  stopifnot(ok)
  message("  Copied: ", basename(to))
}

message("\nDone. Dataset written to: ", out_dir)
