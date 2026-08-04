#!/usr/bin/env bash
# Submit every simulated Klebsiella inference run (both LD-pruning arms).
set -euo pipefail
for f in /nfs/research/jlees/jacqueline/thesis_code/gwas_finalruns/inference/klebsiella_homoplasic/*.sh; do
  echo "Submitting $f"
  sbatch "$f"
done
