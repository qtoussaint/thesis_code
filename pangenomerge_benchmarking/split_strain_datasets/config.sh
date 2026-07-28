#!/usr/bin/env bash
# Shared configuration for the split-strain benchmark datasets (chapter 4, analysis A0).
#
# One PopPUNK strain cluster per species, split into halves and quarters, plus a Panaroo
# run over all its isolates as ground truth. Clusters chosen to be closely size-matched
# across species and to already have ggCaller GFFs, so no gene calling is repeated --
# reusing the single GFF per isolate is what lets annotation_id bridge component
# clustering_ids to truth clustering_ids.

SPECIES_LIST=(s_pneumoniae m_tuberculosis s_aureus)

# species -> PopPUNK size-balanced cluster id (n isolates)
declare -A CLUSTER=(
  [s_pneumoniae]=195      # 352 isolates
  [m_tuberculosis]=2197   # 351 isolates
  [s_aureus]=32           # 336 isolates
)

SPLIT_SEED=42

ATB=/nfs/research/jlees/jacqueline/atb_analyses/species_pangenomes
CODE_DIR=/nfs/research/jlees/jacqueline/thesis_code/pangenomerge_benchmarking/split_strain_datasets
RESULTS_ROOT=/nfs/research/jlees/jacqueline/thesis_results/pangenomerge_benchmarking/split_strain_datasets

CONDA_SH=/hps/software/users/jlees/jacqueline/etc/profile.d/conda.sh
# panaroo lives in its own env -- graph_merge carries pangenomerge's deps only.
# All 21 graphs (truth + components) must come from one panaroo build, or the
# annotation_id -> clustering_id bridge and the metrics built on it are invalid.
PANAROO_ENV=/hps/software/users/jlees/jacqueline/envs/panaroo

# every part that gets its own Panaroo graph
PARTS=(all half_1 half_2 quarter_1 quarter_2 quarter_3 quarter_4)

# ggCaller GFF directory for a species (uses that species' selected cluster)
gff_dir () {
  local sp=$1
  echo "${ATB}/${sp}/results/ggcaller/${CLUSTER[$sp]}/GFF"
}
