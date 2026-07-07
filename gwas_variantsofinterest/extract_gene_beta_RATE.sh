#!/usr/bin/env bash
# Extract per-variant beta and RATE values for a gene of interest from a PPOM/POM run.
#
# For each variant in the requested gene it pulls:
#   beta  -- median effect per cutpoint from fitted_model/depruned_variant_effects.csv
#   RATE  -- variable-importance per cutpoint from cppRATE_results/RATE_values_cutpointN_depruned.txt
# and writes one tidy CSV with all cutpoints side by side.
#
# How variants are located:
#   The annotation file (from the run's submission script, e.g.
#   fields_filtered_maf05_multiallelic.txt) lists every variant by POS and gene, in the
#   same order as the phandango plot rows. The phandango row order is the variant_id order
#   used by depruned_variant_effects.csv and the _depruned RATE files (1-based), so:
#       variant_id = phandango row number (after header)
#       POS        = phandango BP column
#   We take the gene's positions from the annotation, find the matching phandango rows to
#   get variant_ids, then join beta and RATE by variant_id + cutpoint.
#
# Usage:
#   extract_gene_beta_RATE.sh <GENE> <RUN_DIR> <ANNOTATION_FILE> <N_CUTPOINTS> <OUT_CSV>
#
# Example (phpA, SPN penicillin MIC PPOM):
#   extract_gene_beta_RATE.sh phpA \
#     /nfs/research/jlees/jacqueline/thesis_results/gwas_spn_penicillin/inference/02_spn_penicillin_MIC_PPOM \
#     /nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/genotype/fields_filtered_maf05_multiallelic.txt \
#     7 \
#     /nfs/research/jlees/jacqueline/thesis_results/gwas_variantsofinterest/spn_penicillin_PPOM_phpA_beta_RATE_by_cutpoint.csv

set -euo pipefail

GENE="${1:?gene name required}"
RUN_DIR="${2:?run directory required}"
ANN="${3:?annotation file required}"
NCP="${4:?number of cutpoints required}"
OUT="${5:?output csv path required}"

EFFECTS="$RUN_DIR/fitted_model/depruned_variant_effects.csv"
PHANDANGO="$RUN_DIR/cppRATE_results/phandango_cutpoint1.plot"   # any cutpoint works: row order is shared

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

# 1) gene positions from the annotation (column 1 = POS, column 6 = GENE)
awk -F'\t' -v g="$GENE" 'NR>1 && $6==g {print $1}' "$ANN" | sort -nu > "$tmp/pos.txt"

# 2) gene positions -> variant_id via phandango row order (vid = FNR-1, BP = $3)
awk 'NR==FNR{p[$1]=1; next} FNR>1{vid=FNR-1; if($3 in p) print vid"\t"$3}' \
    "$tmp/pos.txt" "$PHANDANGO" | sort -n > "$tmp/vid.txt"

vmin=$(head -1 "$tmp/vid.txt" | cut -f1)
vmax=$(tail -1 "$tmp/vid.txt" | cut -f1)
echo "$GENE: $(wc -l < "$tmp/vid.txt") variants (variant_id $vmin-$vmax)" >&2

# 3) beta (median) per variant_id per cutpoint
awk -F, -v lo="$vmin" -v hi="$vmax" 'NR>1 && $1>=lo && $1<=hi {print $1"_"$4"\t"$2}' \
    "$EFFECTS" > "$tmp/beta.txt"

# 4) RATE per variant_id per cutpoint from each _depruned.txt
: > "$tmp/rate.txt"
for c in $(seq 1 "$NCP"); do
  awk -v c="$c" 'FNR==NR{want[$1]=1; next} /^#/{next} ($1 in want){print $1"_"c"\t"$2}' \
      "$tmp/vid.txt" "$RUN_DIR/cppRATE_results/RATE_values_cutpoint${c}_depruned.txt" >> "$tmp/rate.txt"
done

# 5) annotation detail per POS (first effect line per position)
awk -F'\t' -v g="$GENE" 'NR>1 && $6==g && !(seen[$1]++){print $1"\t"$4"|"$5"|"$8}' "$ANN" > "$tmp/ann.txt"

# 6) assemble tidy CSV: variant_id, POS, effect, impact, HGVS_P, beta_cp1..N, RATE_cp1..N
awk -F'\t' -v ncp="$NCP" '
  FILENAME~/vid\.txt$/   {bp[$1]=$2; order[++n]=$1; next}
  FILENAME~/beta\.txt$/  {beta[$1]=$2; next}
  FILENAME~/rate\.txt$/  {rate[$1]=$2; next}
  FILENAME~/ann\.txt$/   {ann[$1]=$2; next}
  END{
    printf "variant_id,POS,effect,impact,HGVS_P"
    for(c=1;c<=ncp;c++) printf ",beta_cp%d",c
    for(c=1;c<=ncp;c++) printf ",RATE_cp%d",c
    printf "\n"
    for(i=1;i<=n;i++){
      v=order[i]; p=bp[v]; split(ann[p],a,"|")
      printf "%s,%s,%s,%s,%s", v, p, a[1], a[2], a[3]
      for(c=1;c<=ncp;c++) printf ",%s", (v"_"c in beta ? beta[v"_"c] : "NA")
      for(c=1;c<=ncp;c++) printf ",%s", (v"_"c in rate ? rate[v"_"c] : "NA")
      printf "\n"
    }
  }' "$tmp/vid.txt" "$tmp/beta.txt" "$tmp/rate.txt" "$tmp/ann.txt" > "$OUT"

echo "Wrote $OUT ($(($(wc -l < "$OUT")-1)) variants)" >&2
