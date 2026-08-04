#!/usr/bin/env bash
# Run rope_significance.R + bf_pd_significance.R on the four unpruned
# simulated Klebsiella fits, at 89% for both the ROPE and the posterior CI.
set -eo pipefail   # not -u: /etc/bashrc references unbound BASHRCSOURCED

source ~/.bashrc
mamba activate gwas_pipeline
# /hps is 100% full, so bayestestR/logspline live on NFS.
export R_LIBS_USER=/nfs/research/jlees/jacqueline/Rlibs

CODE=/nfs/research/jlees/jacqueline/thesis_code/rope_pd_bf
RUNS=/nfs/research/jlees/jacqueline/thesis_results/gwas_klebsiella_homoplasic/inference
DATA=/nfs/research/jlees/jacqueline/thesis_results/gwas_datasets/inference
ANNOT=/nfs/research/jlees/jacqueline/gwas_data/klebsiella/genotype/klebsiella_rs_annotations.txt
GOI=/nfs/research/jlees/jacqueline/thesis_code/gwas_genesofinterest/klebsiella_homoplasic_genesofinterest.txt

for DS in kleb_homoplasic_EF1.5_H09 kleb_homoplasic_EF2.5_H09 \
          kleb_homoplasic_EF10_H09 kleb_homoplasic_EF30_H09; do
  NICK="${DS}_continuous_nopruning"
  echo "=================== $NICK ==================="
  Rscript "$CODE/rope_significance.R" \
    --run-dir       "$RUNS/$NICK" \
    --snpeff        "$ANNOT" \
    --variant-index "$DATA/$DS/${DS}_variant_index.csv" \
    --rope-prob 0.89 --ci-prob 0.89 \
    --genes-of-interest "$GOI"
done
