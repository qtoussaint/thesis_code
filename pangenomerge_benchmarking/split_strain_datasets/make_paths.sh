#!/usr/bin/env bash
# Write the --component-graphs path files for every merge configuration (analysis A0).
#
# Six configurations: 3 species x {2-way (halves), 4-way (quarters)}. Both share the same
# panaroo_all truth graph within a species, and -- because quarters nest inside halves --
# the same isolates, so the only difference between the 2-way and 4-way rows is the number
# of component graphs being merged.
set -euo pipefail
source /nfs/research/jlees/jacqueline/thesis_code/pangenomerge_benchmarking/split_strain_datasets/config.sh

for sp in "${SPECIES_LIST[@]}"; do
  base=${RESULTS_ROOT}/${sp}
  mkdir -p "${base}/paths"

  printf '%s\n' "${base}/panaroo_half_1" "${base}/panaroo_half_2" \
    > "${base}/paths/split2.txt"

  printf '%s\n' "${base}/panaroo_quarter_1" "${base}/panaroo_quarter_2" \
    "${base}/panaroo_quarter_3" "${base}/panaroo_quarter_4" \
    > "${base}/paths/split4.txt"

  for f in split2 split4; do
    n=$(wc -l < "${base}/paths/${f}.txt")
    missing=0
    while read -r d; do
      [[ -f "${d}/final_graph.gml" && -f "${d}/gene_data.csv" && -f "${d}/pan_genome_reference.fa" ]] \
        || { echo "  MISSING inputs: $d" >&2; missing=$((missing+1)); }
    done < "${base}/paths/${f}.txt"
    echo "${sp} ${f}: ${n} component graphs, ${missing} incomplete"
  done
done
