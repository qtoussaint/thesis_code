#!/usr/bin/env Rscript
# Emits one SLURM script per simulated Klebsiella effect size, running
# continuous_inference.stan through the standard gwas_workflow pipeline.
#
# Kept separate from generate_run_scripts.R so re-running that script (which
# regenerates every spn/tb script from its own dataset table) never clobbers
# these, and so the simulated runs are not mixed into the real-data results
# tree. Run once after editing; re-run to regenerate.
#
# Usage:
#   mamba activate gwas_pipeline
#   Rscript gwas_finalruns/generate_klebsiella_run_scripts.R
#   bash    gwas_finalruns/submit_klebsiella.sh

suppressPackageStartupMessages({
  library(glue)
})

# ----------------------------- static paths -------------------------------- #
ROOT             <- "/nfs/research/jlees/jacqueline/thesis_code/gwas_finalruns"
MODELS_DIR       <- "/nfs/research/jlees/jacqueline/thesis_code/gwas_finalmodels"
DATASETS_DIR     <- "/nfs/research/jlees/jacqueline/thesis_results/gwas_datasets"
RESULTS_BASE     <- "/nfs/research/jlees/jacqueline/thesis_results"
GENES_DIR        <- "/nfs/research/jlees/jacqueline/thesis_code/gwas_genesofinterest"
PIPELINE_RSCRIPT <- "/nfs/research/jlees/jacqueline/gwas_workflow/code/gwas_workflow/inst/scripts/run_pipeline.R"
PRUNING_BIN      <- "/hps/software/users/jlees/jacqueline/manual_installs/bin/BacPrune-Rust/"
CPPRATE_BIN      <- "/hps/software/users/jlees/jacqueline/manual_installs/bin/cpprate-0.2.0/build/bin/cpprate"

# Synthetic annotation file: POS is the integer from each rs_<N> ID and
# ANN[*].GENE is the variant's own name, so Manhattan gene labels read "rs_N".
KLEB_ANNOT <- "/nfs/research/jlees/jacqueline/gwas_data/klebsiella/genotype/klebsiella_rs_annotations.txt"
KLEB_GENES <- file.path(GENES_DIR, "klebsiella_homoplasic_genesofinterest.txt")

RESULTS_DIR <- file.path(RESULTS_BASE, "gwas_klebsiella_homoplasic")
OUT_DIR     <- file.path(ROOT, "inference", "klebsiella_homoplasic")

# ----------------------------- dataset table ------------------------------- #
# One row per effect size. Dataset names must match those written by
# gwas_datasets/write_klebsiella_jsons.R.
DATASETS <- data.frame(
  ef      = c("1.5", "2.5", "10", "30"),
  dataset = c("kleb_homoplasic_EF1.5_H09",
              "kleb_homoplasic_EF2.5_H09",
              "kleb_homoplasic_EF10_H09",
              "kleb_homoplasic_EF30_H09"),
  stringsAsFactors = FALSE
)

# Resource request. The genotype is 4,368 x 12,410; the design matrix is
# standardized inside transformed data, so Stan holds a dense 4,368 x 12,410
# real matrix (~434 MB) on top of the JSON parse. cppRATE is the memory peak,
# which is why the older Klebsiella runs asked for 100G with 30 CPUs.
CPUS <- 32L
MEM  <- "250G"
TIME <- "24:00:00"

# Both LD-pruning arms are emitted, so every effect size is run twice and the
# two can be compared side by side. Threshold 1 prunes only perfectly-correlated
# variants.
#
# The arms are not merely cosmetic: pruning changes which variants are fitted
# (4,369 of 4,578) and therefore the optimization problem the sampler faces. In
# an earlier round at grad_samples = 10 the pruned datasets failed for EF 1.5
# and EF 2.5 (h2_narrow 2e-06 and 2e-08) while the unpruned ones recovered the
# causal variant at rank 1 for all four effect sizes. Expect the pruned arm to
# be the more fragile of the two; check h2_narrow per run before trusting a
# panel.
#
# The arms also differ in what the pipeline writes. run_inference_pipeline()
# only builds Manhattans, the annotation join and depruned_variant_effects.csv
# inside the ld_pruning == "true" branch; the no-pruning branch writes just the
# fit, PPC, heritability and cppRATE. paper_figures/klebsiella_extract_effects.R
# therefore reads effects from the fitted model in both arms, de-pruning where
# needed, so the two are treated identically downstream.
MODES <- list(
  list(key = "pruned",    ld = "true",  tag = "prn"),
  list(key = "nopruning", ld = "false", tag = "npr")
)
LD_THRESHOLD <- 1

# Resume behaviour. With --resume the pipeline reuses any fitted model, pruning
# summary, effects CSV and plots already on disk. That must stay off whenever the
# underlying dataset has been rebuilt, otherwise a stale fit is silently kept.
# With it off, run_inference_pipeline() deletes its own stale RDS files and
# regenerates every downstream artefact, so nothing has to be removed by hand.
USE_RESUME <- FALSE

# ADVI gradient samples per iteration. This is not a tuning nicety -- at the
# default of 1 the fit is unstable on these datasets and the outcome depends on
# the seed: across a 4-seed sweep the causal variant landed at ranks 1, 2, 57,
# 359, 467, 1133, 1286, 1864, 4433 or the fit aborted outright, all while
# reporting "MEDIAN ELBO CONVERGED". At 10, on the unpruned datasets these runs
# use, it recovered rs_26645 at rank 1 for every effect size with effects within
# 2% of the simulated truth and sigma matching the true residual SD (0.701 vs
# 0.705 at EF 1.5 through 13.72 vs 14.03 at EF 30).
#
# Caveat worth keeping in mind: that check covered one seed across all four
# effect sizes, not several seeds. If a future run of one of these datasets
# comes back with h2_narrow near zero, suspect a bad ADVI optimum rather than
# the data, and re-fit before believing it.
GRAD_SAMPLES <- 10L

# Per-run overrides, keyed by nickname, for fits that will not run at the
# defaults. These are retries of HARD FAILURES (ADVI aborting with "dropped
# evaluations ... severely ill-conditioned" and returning no draws at all), not
# a search across successful fits for a preferred answer -- picking among
# converged fits would be cherry-picking, whereas a run that produced nothing
# has to be re-attempted somehow.
#
# kleb_homoplasic_EF1.5_H09_continuous_pruned aborted after 60s at
# grad_samples = 10 / seed 987654321. The pruned arm is the fragile one, and
# EF 1.5 is its weakest case.
#
# Anything listed here differs from its sibling runs, so note it when comparing
# arms.
VI_SEED_DEFAULT <- 987654321L   # cmdstanr/pipeline default

OVERRIDES <- list(
  "kleb_homoplasic_EF1.5_H09_continuous_pruned" = list(grad_samples = 25L, vi_seed = 42L)
)

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

TEMPLATE <- "#!/usr/bin/env bash
#SBATCH --job-name={job_name}
#SBATCH --nodes=1
#SBATCH --cpus-per-task={cpus}
#SBATCH --mem={mem}
#SBATCH --time={time}
#SBATCH --error={run_dir}/logs/{nickname}_%j.err
#SBATCH --output={run_dir}/logs/{nickname}_%j.out

#################################################################################

source ~/.bashrc
mamba activate gwas_pipeline

mkdir -p {run_dir}/logs

RSCRIPT_PATH=\"{pipeline}\"

DATA=\"--data {dataset_dir}/{dataset}.json\"
STAN_MODEL=\"--stan_model {models_dir}/continuous_inference.stan\"
ANALYSIS_TYPE=\"--analysis_type inference\"
ANALYSIS_NICKNAME=\"--analysis_nickname {nickname}\"
OUTPUT_DIR=\"--output_directory {run_dir}\"
THREADS=\"--threads {cpus}\"

LD_PRUNING=\"--ld_pruning {ld_pruning}\"
PRUNING_SOFTWARE=\"--pruning_software {pruning_bin}\"
MAF_CUTOFF=\"--maf_cutoff 0\"
LD_THRESHOLD=\"--ld_threshold {ld_threshold}\"

PHANDANGO=\"--phandango {dataset_dir}/{dataset}_variant_index.csv\"
ANNOTATIONS=\"--annotations {annot}\"
MODEL_TYPE=\"--model_type continuous\"
GRAD_SAMPLES=\"--grad_samples {grad_samples}\"
VI_SEED=\"--vi_seed {vi_seed}\"
GENES_OF_INTEREST=\"--genes_of_interest {genes}\"
RESUME=\"{resume_flag}\"
CPPRATE=\"--cpprate_bin {cpprate_bin}\"

Rscript $RSCRIPT_PATH \\
$DATA \\
$STAN_MODEL \\
$ANALYSIS_TYPE \\
$ANALYSIS_NICKNAME \\
$OUTPUT_DIR \\
$THREADS \\
$LD_PRUNING \\
$PRUNING_SOFTWARE \\
$MAF_CUTOFF \\
$LD_THRESHOLD \\
$PHANDANGO \\
$ANNOTATIONS \\
$MODEL_TYPE \\
$GRAD_SAMPLES \\
$VI_SEED \\
$GENES_OF_INTEREST \\
$RESUME \\
$CPPRATE
"

written <- character(0)

for (m in MODES) {
  for (i in seq_len(nrow(DATASETS))) {
    dataset  <- DATASETS$dataset[i]
    # The arm is part of the nickname, so the two never share an output
    # directory and neither can silently overwrite the other's fit.
    nickname <- paste0(dataset, "_continuous_", m$key)
    run_dir  <- file.path(RESULTS_DIR, "inference", nickname)

    ov <- OVERRIDES[[nickname]]
    gs <- if (!is.null(ov$grad_samples)) ov$grad_samples else GRAD_SAMPLES
    sd_ <- if (!is.null(ov$vi_seed))     ov$vi_seed      else VI_SEED_DEFAULT
    if (!is.null(ov)) {
      message("  override for ", nickname, ": grad_samples=", gs, " vi_seed=", sd_)
    }

    body <- glue(
      TEMPLATE,
      job_name     = paste0("kEF", DATASETS$ef[i], "_", m$tag),
      cpus         = CPUS,
      mem          = MEM,
      time         = TIME,
      run_dir      = run_dir,
      nickname     = nickname,
      pipeline     = PIPELINE_RSCRIPT,
      dataset_dir  = file.path(DATASETS_DIR, "inference", dataset),
      dataset      = dataset,
      models_dir   = MODELS_DIR,
      pruning_bin  = PRUNING_BIN,
      ld_threshold = LD_THRESHOLD,
      annot        = KLEB_ANNOT,
      genes        = KLEB_GENES,
      resume_flag  = if (USE_RESUME) "--resume" else "",
      grad_samples = gs,
      vi_seed      = sd_,
      ld_pruning   = m$ld,
      cpprate_bin  = CPPRATE_BIN,
      .open = "{", .close = "}"
    )

    path <- file.path(OUT_DIR, paste0(nickname, ".sh"))
    writeLines(body, path)
    Sys.chmod(path, "0755")
    written <- c(written, path)
    message("Wrote ", path)
  }
}

# Convenience submitter
submit <- c(
  "#!/usr/bin/env bash",
  "# Submit every simulated Klebsiella inference run (both LD-pruning arms).",
  "set -euo pipefail",
  paste0("for f in ", OUT_DIR, "/*.sh; do"),
  "  echo \"Submitting $f\"",
  "  sbatch \"$f\"",
  "done"
)
submit_path <- file.path(ROOT, "submit_klebsiella.sh")
writeLines(submit, submit_path)
Sys.chmod(submit_path, "0755")
message("Wrote ", submit_path)

message("\nGenerated ", length(written), " run scripts.")
