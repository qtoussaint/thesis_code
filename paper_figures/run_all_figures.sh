#!/bin/bash
# Remake every paper figure EXCEPT the cutpoint histograms (cutpoints_histogram_summary.R /
# cutpoints_histogram_faceted_summary.R), which are running separately.
# Reconstructed invocations for the arg-taking replot/composite scripts come from
# the constants in the manhattan_locuszoom_composite*.R scripts.
source ~/.bashrc
mamba activate gwas_pipeline
set -uo pipefail

CODE=/nfs/research/jlees/jacqueline/thesis_code/paper_figures
RES=/nfs/research/jlees/jacqueline/thesis_results
PF=$RES/paper_figures
LOG=$PF/run_all_logs
mkdir -p "$LOG"
cd "$CODE"

SPN_ANNOT=/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/genotype/fields_filtered_maf05_multiallelic.txt
TB_ANNOT=/nfs/research/jlees/jacqueline/gwas_data/tuberculosis/cryptic_regeno_snpeff/fields_filtered.txt
GOI=/nfs/research/jlees/jacqueline/thesis_code/gwas_genesofinterest
PEN_POS=$RES/gwas_datasets/inference/02_spn_penicillin_MIC/02_spn_penicillin_MIC_variant_index.csv

RUN02=$RES/gwas_spn_penicillin/inference/02_spn_penicillin_MIC_PPOM
RUN05=$RES/gwas_spn_trimethoprim/inference/05_spn_trimethoprim_MIC_PPOM
RUN07=$RES/gwas_tb_rifampicin/inference/07_tb_rifampicin_binary_logistic
RUN16=$RES/gwas_spn_penicillin/inference/16_spn_penicillin_MIC_minimabinning_PPOM

declare -a NAMES=()
declare -a STATUS=()

run() {
  local name="$1"; shift
  echo "===== START $name ====="
  if Rscript "$@" >"$LOG/$name.log" 2>&1; then
    echo "----- OK    $name"
    NAMES+=("$name"); STATUS+=("OK")
  else
    echo "##### FAIL  $name (see $LOG/$name.log)"
    NAMES+=("$name"); STATUS+=("FAIL")
  fi
}

# ---- A. self-contained summary scripts -------------------------------------
run alpha_sweep              alpha_sweep_summary.R
run tau_sweep                tau_sweep_summary.R
run regularization_prior     regularization_prior_summary.R
run dataset_split            dataset_split_summary.R
run mic_bin_histogram        mic_bin_histogram_summary.R
run heritability             heritability_summary.R
run heritability_combined    heritability_summary_combined.R
run inference_ppc            inference_ppc_summary.R
run prediction_accuracy      prediction_accuracy_summary.R
run prediction_plots         prediction_plots_summary.R
run faceted_rate_beta        faceted_rate_beta_summary.R
run violin_pbp_mray          violin_pbp_mray_summary.R
run combine_16_median        combine_16_spn_penicillin_median_manhattans.R
run combine_faceted          combine_spn_penicillin_faceted.R

# ---- B. manhattan + locus-zoom composites ----------------------------------
run composite_pen_full       manhattan_locuszoom_composite.R
run composite_pen_smaller    manhattan_locuszoom_composite.R --layout smaller
run composite_tmp_rif        manhattan_locuszoom_composite_rate_tmp_rif.R

# ---- C. per-run short overlay manhattans -----------------------------------
# 02 penicillin: unlabeled (median, exp_abs, RATE) into manhattan_short
run short_02_overlay         replot_ppom_overlay_manhattans_short.R \
    --run-dir "$RUN02" --phandango "$PEN_POS"
# 02 penicillin: labeled median + exp_abs into manhattan_short_labeled
run short_02_overlay_lab     replot_ppom_overlay_manhattans_short.R \
    --run-dir "$RUN02" --phandango "$PEN_POS" \
    --annotations "$SPN_ANNOT" --genes-of-interest "$GOI/spn_penicillin_genesofinterest.txt" \
    --label-genes "pbp1a,pbp2X,pbp2b,folA" \
    --output-dir "$PF/02_spn_penicillin_MIC_PPOM/manhattan_short_labeled"
# 02 penicillin: labeled RATE into manhattan_short_labeled
run short_02_rate_lab        replot_rate_manhattan_short.R \
    --run-dir "$RUN02" --annotations "$SPN_ANNOT" \
    --genes-of-interest "$GOI/spn_penicillin_genesofinterest.txt" \
    --label-mode gene_list --label-genes "pbp1a,pbp2X,pbp2b" \
    --output-dir "$PF/02_spn_penicillin_MIC_PPOM/manhattan_short_labeled"

# 05 trimethoprim: unlabeled RATE + labeled RATE
run short_05_rate            replot_rate_manhattan_short.R --run-dir "$RUN05"
run short_05_rate_lab        replot_rate_manhattan_short.R \
    --run-dir "$RUN05" --annotations "$SPN_ANNOT" \
    --genes-of-interest "$GOI/spn_trimethoprim_genesofinterest.txt" \
    --label-mode gene_list --label-genes "folA (dhfR),folP" \
    --output-dir "$PF/05_spn_trimethoprim_MIC_PPOM/manhattan_short_labeled"

# 07 TB rifampicin (binary): unlabeled RATE + labeled RATE
run short_07_rate            replot_rate_manhattan_short.R --run-dir "$RUN07"
run short_07_rate_lab        replot_rate_manhattan_short.R \
    --run-dir "$RUN07" --annotations "$TB_ANNOT" \
    --genes-of-interest "$GOI/tb_rifampicin_genesofinterest.txt" \
    --label-mode top_n --n-labels 10 \
    --output-dir "$PF/07_tb_rifampicin_binary_logistic/manhattan_short_labeled"

# ---- D. per-run faceted RATE cutpoints -> faceted_cutpoints/ ----------------
run faceted_02               replot_rate_faceted_cutpoints.R \
    --run-dir "$RUN02" --genes-of-interest "$GOI/spn_penicillin_genesofinterest.txt" \
    --output-dir "$PF/faceted_cutpoints"
run faceted_05               replot_rate_faceted_cutpoints.R \
    --run-dir "$RUN05" --genes-of-interest "$GOI/spn_trimethoprim_genesofinterest.txt" \
    --output-dir "$PF/faceted_cutpoints"
run faceted_16               replot_rate_faceted_cutpoints.R \
    --run-dir "$RUN16" --genes-of-interest "$GOI/spn_penicillin_genesofinterest.txt" \
    --output-dir "$PF/faceted_cutpoints"

# ---- summary ---------------------------------------------------------------
echo ""
echo "================= SUMMARY ================="
for i in "${!NAMES[@]}"; do
  printf "%-6s %s\n" "${STATUS[$i]}" "${NAMES[$i]}"
done
echo "Logs: $LOG"
