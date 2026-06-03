#!/usr/bin/env bash
#SBATCH --job-name=unitig_map_expand
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --error=/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/unitigs/map_work/logs/map_expand_%j.err
#SBATCH --output=/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/unitigs/map_work/logs/map_expand_%j.out

# Map every unitig to the FM211187 (ATCC 700669 / Spn23F) reference with bwa mem,
# collect ALL hit coordinates per unitig, and build a coordinate-expanded
# presence/absence matrix (one column per coordinate; unmapped -> position -10000).
# Uses the pyseer env, which provides bwa and python.

# Activate the env BEFORE enabling strict mode: sourcing ~/.bashrc / mamba activate
# run intermediate commands that return nonzero and would abort under `set -e`/`set -u`.
source ~/.bashrc
mamba activate pyseer
set -euo pipefail

THREADS=${SLURM_CPUS_PER_TASK:-16}

UNI=/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/unitigs
REF_SRC=/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/pyseer_data/pyseer_tutorial/Spn23F.fa
CODE=/nfs/research/jlees/jacqueline/thesis_code/gwas_datasets/map_expand_unitigs

WORK=$UNI/map_work
mkdir -p "$WORK/logs"

REF=$WORK/Spn23F.fa
UNITIG_FA=$WORK/unitigs.fa
SAM=$WORK/unitigs.sam

# 1. Local copy of the reference + bwa index (avoid polluting the shared tutorial dir).
if [[ ! -f "$REF.bwt" ]]; then
  cp -f "$REF_SRC" "$REF"
  bwa index "$REF"
fi

# 2. All unitigs as FASTA (id = unitig_N, seq = sequence).
if [[ ! -s "$UNITIG_FA" ]]; then
  awk -F'\t' '{print ">"$1"\n"$2}' "$UNI/spn_unitigs.unitig_ids.tsv" > "$UNITIG_FA"
fi

# 3. Map. -a: report all (incl. secondary) alignments; lowered -k/-T so short 31-mers
#    with score ~31 are reported; -c keeps repetitive seeds (multi-mapping unitigs).
bwa mem -a -k 17 -T 20 -c 10000 -t "$THREADS" "$REF" "$UNITIG_FA" > "$SAM"

# 4. Parse SAM + build coordinate-expanded matrix (MAF in [0.05,0.95]; unmapped -> -10000).
python "$CODE/parse_and_expand.py" \
  --sam          "$SAM" \
  --unitig_ids   "$UNI/spn_unitigs.unitig_ids.tsv" \
  --rtab         "$UNI/spn_unitigs.rtab" \
  --out_coords   "$UNI/spn_unitigs.coords.tsv" \
  --out_rtab     "$UNI/spn_unitigs_mapped.rtab" \
  --out_varindex "$UNI/spn_unitigs_mapped_variant_index.csv" \
  --min_af 0.05 --max_af 0.95 --min_cov 0.8 --unmapped_pos -10000

echo "Done. Outputs in $UNI:"
echo "  spn_unitigs.coords.tsv"
echo "  spn_unitigs_mapped.rtab"
echo "  spn_unitigs_mapped_variant_index.csv"
