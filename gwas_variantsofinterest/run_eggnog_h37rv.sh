#!/usr/bin/env bash
#SBATCH --job-name=eggnog_h37rv
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=200G
#SBATCH --time=12:00:00
#SBATCH --error=/nfs/research/jlees/jacqueline/thesis_results/gwas_variantsofinterest/eggnog_tb/logs/eggnog_h37rv_%j.err
#SBATCH --output=/nfs/research/jlees/jacqueline/thesis_results/gwas_variantsofinterest/eggnog_tb/logs/eggnog_h37rv_%j.out

#################################################################################
# Annotate the M. tuberculosis H37Rv (GCF_000195955.2) proteome with eggNOG-mapper
# so the Rv locus tags in the TB variants-of-interest CSVs can be resolved to gene
# names / functional descriptions, mirroring run_eggnog_spn23f.sh. Runs diamond
# against the local bacteria.dmnd and annotates from the local eggNOG 5.0.2 DB;
# --dbmem loads eggnog.db (41G) to RAM for a fast annotation step.
#################################################################################

source ~/.bashrc
mamba activate eggnog

EGGNOG_DB=/nfs/research/jlees/jacqueline/eggnog_dbs
CODE=/nfs/research/jlees/jacqueline/thesis_code/gwas_variantsofinterest
GBFF=/nfs/research/jlees/jacqueline/gwas_data/tuberculosis/pyseer_data/trees/reference/GCF_000195955.2/genomic.gbff
OUT=/nfs/research/jlees/jacqueline/thesis_results/gwas_variantsofinterest/eggnog_tb
QUERY=$OUT/h37rv_proteins.faa

mkdir -p "$OUT/logs"

# Query = every CDS protein translation in the H37Rv GenBank file, keyed by locus_tag.
python3 "$CODE/gbff_to_protein_fasta.py" "$GBFF" "$QUERY"

# Diamond search (bacteria scope) + eggNOG annotation, emapper defaults.
emapper.py \
  -i "$QUERY" \
  -o h37rv \
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

echo "emapper finished with exit code $? -- annotations: $OUT/h37rv.emapper.annotations"
