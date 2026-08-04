############################################################
## config_klebsiella.R
## Input paths and parameters for the simulated K. pneumoniae
## homoplasic-phenotype inference datasets.
##
## Kept separate from config.R so the real-data (spn/tb) config
## stays untouched. Sourced by write_klebsiella_jsons.R.
############################################################

# === INPUT DATA PATHS ===

# Klebsiella genotype — SNP presence/absence, variants x samples.
# Row names are simulated variant IDs of the form rs_<N>; the header row holds
# sample IDs (Ref_NJST258_1, 573.*, SPARK_*_C1). 12,410 variants x 4,368 samples.
# Byte-identical to Kpneu_MIC_prediction/data/SNPs/matrix_and_lineages/pres_abs_matrix.txt.
KLEB_GENOTYPE_PATH <- "/nfs/research/jlees/jacqueline/gwas_data/klebsiella/genotype/genotype_pres_abs.txt"

# Klebsiella annotations — synthetic, one row per variant. These are simulated
# variants with no genome coordinates and no gene assignment, so POS is the
# integer parsed out of the rs_<N> ID and ANN[*].GENE is the variant's own name.
# Generated alongside this config; regenerate with make_klebsiella_annotations.py
# if the genotype matrix ever changes.
KLEB_ANNOTATIONS_PATH <- "/nfs/research/jlees/jacqueline/gwas_data/klebsiella/genotype/klebsiella_rs_annotations.txt"

# Klebsiella lineages — PopPUNK clusters, all 4,368 samples, 347 clusters.
# Tab-delimited, no header: <sample_id> <PP cluster>.
KLEB_POPPUNK_PATH <- "/nfs/research/jlees/jacqueline/gwas_data/klebsiella/lineages/klebsiella_poppunk_clusters.txt"

# Klebsiella sublineages — PopPIPE/fastbaps subclusters, 3,417 of 4,368 samples.
# CSV with header: Isolates,Cluster_1,Cluster_2.
#
# Two quirks handled in load_klebsiella_lineages():
#   1. Isolate IDs use 573_NNNNN where the genotype matrix uses 573.NNNNN.
#   2. Cluster_1 numbering RESTARTS within each PopPUNK cluster (only 7 distinct
#      values across the whole file), so the raw labels are not globally unique
#      and must be prefixed with the parent PopPUNK cluster.
KLEB_SUBLINEAGE_PATH <- "/nfs/research/jlees/jacqueline/gwas_data/klebsiella/lineages/cluster_summary_poppipe_sublineages.txt"

# Simulated phenotype directory. Files are named
#   simul_homoplasic_EF_<EF>_herit_<H>_traitprev_0.1.phen
# and are whitespace-delimited with no header: FID IID value.
KLEB_PHENO_DIR <- "/nfs/research/jlees/jacqueline/gwas_data/klebsiella/phenotypes"

# Causal variant record for the homoplasic architecture: <variant> <effect size>.
# The same variant (rs_26645) is causal at every effect size.
KLEB_CAUSAL_DIR <- "/nfs/research/jlees/jacqueline/gwas_data/klebsiella/causal_phenotypes"


# === DATASET PARAMETERS ===

# Effect sizes to build. Note the "EffSize" label is a multiplier, not the raw
# beta: the observed carrier vs non-carrier phenotype difference is ~4.1x the
# label at every effect size (6.19 / 10.25 / 41.43 / 123.23 for EF 1.5 / 2.5 /
# 10 / 30). Expected recovered beta is therefore ~4.1 * EF, not EF.
KLEB_EFFECT_SIZES <- c("1.5", "2.5", "10.0", "30.0")

# Heritability. Only 0.6 and 0.9 were simulated; 0.9 matches all previous
# Klebsiella runs.
KLEB_HERITABILITY <- "0.9"

# Phenotype transform applied before fitting.
#   "none"   -- raw simulated phenotype, model priors untouched (current choice)
#   "zscore" -- centre and scale to unit variance; recover raw-scale effects
#               afterwards as beta_raw = beta_z * sd(y), which is exact for
#               this linear-Gaussian model
#
# Why this switch exists: continuous_inference.stan is calibrated for log2-MIC
# (sigma ~ half-N(0,1), alpha_prior_sd = 0.5, beta_lineage ~ N(0, 0.1)) but the
# simulated phenotype SD is 2.23 / 3.72 / 14.90 / 44.36 across these four
# effect sizes, implying a residual SD of roughly 0.71 / 1.18 / 4.71 / 14.03 at
# h2 = 0.9. EF 10 and especially EF 30 therefore sit well outside the sigma
# prior. Flip to "zscore" and re-run if the fitted sigma comes back pinned.
KLEB_PHENO_TRANSFORM <- "none"


# Minimum sample count for a sublineage to be eligible as the dropped reference
# level.
#
# Why this exists: continuous_inference.stan anchors its intercept prior on the
# reference sublineage, alpha ~ normal(mean phenotype of reference, 0.5). SPN
# can pick the plain lowest-mean subcluster because its PopPIPE input was built
# with min_cluster_size = 3, so every candidate is a real group. Klebsiella has
# 138 singleton sublineages out of 607 (many PopPUNK clusters hold one isolate),
# so the unguarded lowest-mean rule selects the single most extreme sample and
# pins alpha ~1.3 SD away from the data mean with a tight sd of 0.5. Nothing in
# the model can absorb that offset -- beta_lineage ~ N(0, 0.1) is tight and the
# within-lineage centering forces beta_sublineage to exactly 0 for the 265
# lineages that have only one sublineage -- so sigma inflates and the horseshoe
# shrinks every variant effect to zero. That is what made the first run miss a
# causal variant correlated 0.95 with the phenotype.
KLEB_MIN_REF_SIZE <- 20L

# Minimum minor allele frequency for a variant to enter the dataset.
#
# Why this exists: with all 12,410 variants (the genotype is pre-filtered only to
# MAF >= 0.005) ADVI aborts with "dropped evaluations has reached its maximum ...
# severely ill-conditioned or misspecified" and returns no draws at all. Raising
# --grad_samples to 10, the documented remedy in args.R, does not rescue it.
# The cause is the standardization in transformed data: the rarest variants have
# sd ~0.07, so X_std reaches 14.1 against 4.4 once MAF > 0.05, and the horseshoe
# geometry on top of that is more than mean-field VI can handle.
#
# 0.05 matches the SPN genotype the model was calibrated on
# (genotype_maf05_NOMODIFIERS_multiallelic.tsv), which is the configuration known
# to work -- SPN penicillin recovers pbpX at RATE rank 1. It drops 7,832 of
# 12,410 variants (63%). The causal variant rs_26645 has MAF 0.135 and is
# retained. Intermediate thresholds are untested: MAF > 0.02 keeps 6,826 and
# MAF > 0.01 keeps 9,077.
#
# Note this differs from the older Klebsiella runs in bayesian_gwas_paper, which
# passed --maf_cutoff 0 to a much simpler model with no horseshoe and no
# genotype standardization. Filtering happens here, at dataset build time, so
# the pipeline is still invoked with --maf_cutoff 0 and does not filter twice.
KLEB_MIN_MAF <- 0.05


# === OUTPUT PATHS ===

KLEB_OUT_BASE  <- "/nfs/research/jlees/jacqueline/thesis_results/gwas_datasets"
KLEB_OUT_INFER <- file.path(KLEB_OUT_BASE, "inference")
