#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=300G
#SBATCH --time=12:00:00
# Panaroo for one part of one species' split-strain dataset (analysis A0).
#   usage: sbatch run_panaroo.sh <species> <part>
# Settings match the Snakefile's panaroo_cluster rule and the original split-strain test:
# --clean-mode strict --refind-mode off --remove-invalid-genes.
set -euo pipefail

# NB: absolute path, not $(dirname $0) -- sbatch copies this script into the node's spool
# directory, so $0 does not resolve to the submit directory at run time.
source /nfs/research/jlees/jacqueline/thesis_code/pangenomerge_benchmarking/split_strain_datasets/config.sh

SPECIES=$1
PART=$2

LIST=${RESULTS_ROOT}/${SPECIES}/splits/${PART}.txt
OUTDIR=${RESULTS_ROOT}/${SPECIES}/panaroo_${PART}

if [[ ! -s "$LIST" ]]; then
  echo "missing or empty split list: $LIST (run make_splits.py first)" >&2
  exit 1
fi

mkdir -p "$OUTDIR"

source "$CONDA_SH"
conda activate "$PANAROO_ENV"

# follow whatever submit_all.sh actually requested rather than a hardcoded count
THREADS=${SLURM_CPUS_PER_TASK:-16}

echo "[panaroo] species=$SPECIES part=$PART cluster=${CLUSTER[$SPECIES]} n=$(wc -l < "$LIST") threads=$THREADS"
echo "[panaroo] outdir=$OUTDIR"

# GFF paths come from the split list, one per line
mapfile -t GFFS < "$LIST"

/usr/bin/time -v -o "${OUTDIR}/time.log" \
  panaroo -i "${GFFS[@]}" -o "$OUTDIR" \
    --clean-mode strict --refind-mode off --remove-invalid-genes \
    --threads "$THREADS"

test -f "${OUTDIR}/final_graph.gml"
test -f "${OUTDIR}/gene_data.csv"
echo "ok" > "${OUTDIR}/done.txt"
echo "[panaroo] complete: $OUTDIR"
