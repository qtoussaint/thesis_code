#!/usr/bin/env bash
#SBATCH --job-name=eggnog_spn23f
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=200G
#SBATCH --time=12:00:00
#SBATCH --error=/nfs/research/jlees/jacqueline/thesis_results/gwas_variantsofinterest/eggnog/logs/eggnog_spn23f_%j.err
#SBATCH --output=/nfs/research/jlees/jacqueline/thesis_results/gwas_variantsofinterest/eggnog/logs/eggnog_spn23f_%j.out

#################################################################################
# Annotate the Spn23F (ATCC 700669) proteome with eggNOG-mapper so the SPN23F
# locus tags in the top5 beta/RATE CSVs can be resolved to gene names. Runs diamond
# against the local bacteria.dmnd and annotates from the local eggNOG 5.0.2 DB.
# The earlier run failed on a small interactive node (OOM building diamond's seed
# array); a compute node with 200G handles it, and --dbmem loads eggnog.db (41G) to
# RAM for a fast annotation step.
#################################################################################

source ~/.bashrc
mamba activate eggnog

EGGNOG_DB=/nfs/research/jlees/jacqueline/eggnog_dbs
CODE=/nfs/research/jlees/jacqueline/thesis_code/gwas_variantsofinterest
GBFF=/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/genotype/ATCC_700669.gbff
OUT=/nfs/research/jlees/jacqueline/thesis_results/gwas_variantsofinterest/eggnog
QUERY=$OUT/spn23f_proteins.faa

mkdir -p "$OUT/logs"

# Query = every CDS protein translation in the Spn23F GenBank file, keyed by locus_tag.
python3 "$CODE/gbff_to_protein_fasta.py" "$GBFF" "$QUERY"

# Diamond search (bacteria scope) + eggNOG annotation. Full sensitivity / iterate, as
# emapper defaults, now that memory is not the limit.
emapper.py \
  -i "$QUERY" \
  -o spn23f \
  --output_dir "$OUT" \
  --data_dir "$EGGNOG_DB" \
  -m diamond \
  --dmnd_db "$EGGNOG_DB/bacteria.dmnd" \
  --itype proteins \
  --sensmode sensitive \
  --dmnd_iterate yes \
  --dbmem \
  --cpu "${SLURM_CPUS_PER_TASK:-32}" \
  --override

echo "emapper finished with exit code $? -- annotations: $OUT/spn23f.emapper.annotations"
