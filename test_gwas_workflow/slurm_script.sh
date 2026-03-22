#!/usr/bin/env bash

#SBATCH --job-name=test_workflow               # job name
#SBATCH --nodes=1                                  # number of nodes
#SBATCH --cpus-per-task=48                         # CPUs/task (cppRATE uses --threads minus 1)
#SBATCH --mem=90G                                 # memory
#SBATCH --time=5:00:00                            # time limit
#SBATCH --error=/nfs/research/jlees/jacqueline/thesis_code/test_gwas_workflow/test_workflow.err      # error file
#SBATCH --output=/nfs/research/jlees/jacqueline/thesis_code/test_gwas_workflow/test_workflow.out     # output file

#################################################################################

#conda activate gwas_pipeline

RSCRIPT_PATH="/nfs/research/jlees/jacqueline/gwas_workflow/code/gwas_workflow/inst/scripts/run_pipeline.R"

# ---------------------------------------------------------------------------
# Required arguments
# ---------------------------------------------------------------------------

DATA="--data /nfs/research/jlees/jacqueline/thesis_results/test_gwas_workflow/01_spn_penicillin_subclusters_K5_3f7c4f1ce71de0.json"
STAN_MODEL="--stan_model /nfs/research/jlees/jacqueline/thesis_code/test_gwas_workflow/ordinal-subcluster-standard-PPOM-association_SPNPENcutpoints.stan"
ANALYSIS_TYPE="--analysis_type inference"           # inference | prediction
ANALYSIS_NICKNAME="--analysis_nickname test_worklow"
OUTPUT_DIR="--output_directory /nfs/research/jlees/jacqueline/thesis_results/test_gwas_workflow"
THREADS="--threads 48"
CPPRATE_BIN="--cpprate_bin /hps/software/users/jlees/jacqueline/manual_installs/bin/cpprate-0.2.0/build/bin/cpprate"

# ---------------------------------------------------------------------------
# LD pruning arguments (required when --ld_pruning true)
# ---------------------------------------------------------------------------

LD_PRUNING="--ld_pruning true"                      # true | false
PRUNING_SOFTWARE="--pruning_software /nfs/research/jlees/jacqueline/gwas_code/ld_pruning/BacPrune-Rust/"
MAF_CUTOFF="--maf_cutoff 0"                      # 0 disables MAF filtering
#CARGO_BIN="--cargo_bin /path/to/.cargo/bin/cargo"   # defaults to ~/.cargo/bin/cargo if omitted

# ---------------------------------------------------------------------------
# Optional arguments
# ---------------------------------------------------------------------------

PHANDANGO="--phandango /nfs/research/jlees/jacqueline/bayesian_gwas_paper/00_data/inference/01_spn_penicillin_subclusters_K5_variant_index.csv"
#CMDSTAN_PATH="--cmdstan_path /hps/software/users/jlees/jacqueline/manual_installs/bin/cmdstan-2.34.1"  # omit to use auto-detected CmdStan
ANNOTATIONS="--annotations /nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/genotype/fields_filtered_maf05_multiallelic.txt"  # snpEff tab-delimited file; enables gene labels on plots
GENES_OF_INTEREST="--genes_of_interest /nfs/research/jlees/jacqueline/thesis_code/test_gwas_workflow/genes_of_interest.csv"  # two-col CSV (gene name, display label); PPOM only; requires --phandango and --annotations
MODEL_TYPE="--model_type ppom"               # binary | pom | ppom | continuous; for prediction accuracy
#TRUE_PHENOTYPES="--true_phenotypes /path/to/true_phenotypes.csv"  # columns: sample_id, true_phenotype; for prediction accuracy
RESUME="--resume"                              # resume from last completed step
NORATE="--norate"                                       # skip cppRATE and RATE plots; omit --cpprate_bin when using this

# ---------------------------------------------------------------------------

Rscript $RSCRIPT_PATH \
$DATA \
$STAN_MODEL \
$ANALYSIS_TYPE \
$ANALYSIS_NICKNAME \
$OUTPUT_DIR \
$THREADS \
$CPPRATE_BIN \
$LD_PRUNING \
$PRUNING_SOFTWARE \
$MAF_CUTOFF \
$CARGO_BIN \
$PHANDANGO \
$CMDSTAN_PATH \
$ANNOTATIONS \
$GENES_OF_INTEREST \
$MODEL_TYPE \
$RESUME \
#$NORATE \
#$TRUE_PHENOTYPES \
