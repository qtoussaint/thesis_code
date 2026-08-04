#!/usr/bin/env bash
# Submit one job per unpruned simulated Klebsiella fit: prior-ROPE significance,
# then BF/pd for the significance list, then BF/pd for every variant (for the
# Manhattans). Needs cluster memory -- the fits are 2.2 GB RDS and the login
# node OOMs silently.
CODE=/nfs/research/jlees/jacqueline/thesis_code/rope_pd_bf
RUNS=/nfs/research/jlees/jacqueline/thesis_results/gwas_klebsiella_homoplasic/inference
DATA=/nfs/research/jlees/jacqueline/thesis_results/gwas_datasets/inference
OUT=/nfs/research/jlees/jacqueline/thesis_results/gwas_klebsiella_homoplasic/rope_pd_bf
ANNOT=/nfs/research/jlees/jacqueline/gwas_data/klebsiella/genotype/klebsiella_rs_annotations.txt
GOI=/nfs/research/jlees/jacqueline/thesis_code/gwas_genesofinterest/klebsiella_homoplasic_genesofinterest.txt

for DS in kleb_homoplasic_EF1.5_H09 kleb_homoplasic_EF2.5_H09 \
          kleb_homoplasic_EF10_H09 kleb_homoplasic_EF30_H09; do
  NICK="${DS}_continuous_nopruning"
  sbatch --job-name="rope_${DS}" --nodes=1 --cpus-per-task=4 --mem=100G --time=12:00:00 \
    --output="$OUT/logs/${NICK}_%j.out" --error="$OUT/logs/${NICK}_%j.err" \
    --wrap="source ~/.bashrc; mamba activate gwas_pipeline; \
export R_LIBS_USER=/nfs/research/jlees/jacqueline/Rlibs; \
Rscript $CODE/rope_significance.R --run-dir $RUNS/$NICK --snpeff $ANNOT \
  --variant-index $DATA/$DS/${DS}_variant_index.csv \
  --rope-prob 0.89 --ci-prob 0.89 --genes-of-interest $GOI && \
for L in \$(ls $RUNS/$NICK/rope_pd_bf/*signif*.csv 2>/dev/null); do \
  Rscript $CODE/bf_pd_significance.R --variant-list \$L --run-dir $RUNS/$NICK \
    --variant-index $DATA/$DS/${DS}_variant_index.csv || true; \
done; \
Rscript $CODE/klebsiella_bf_pd_all_variants.R $RUNS/$NICK \
  $DATA/$DS/${DS}_variant_index.csv $OUT/${DS}_all_variants_bf_pd.csv"
done
