#!/usr/bin/env bash
# Example invocation for test_variant_proportionality.R.
#
# Preconditions:
#   - POM pipeline run has completed with output at $POM_DIR.
#   - PPOM pipeline run has completed with output at $PPOM_DIR.
#   - Both runs used the same dataset and (ideally) identical LD thresholds,
#     so variant indexing matches.
#
# Required per-run files (produced by gwaspipeline):
#   $POM_DIR/fitted_model/*.RDS
#   $POM_DIR/cppRATE_matrices/ld_pruning_summary.csv
#   $POM_DIR/cppRATE_matrices/bacprune_rust_results.csv
#   $PPOM_DIR/fitted_model/*.RDS
#   $PPOM_DIR/fitted_model/depruned_variant_effects.csv   (optional cross-check)
#   $PPOM_DIR/cppRATE_matrices/ld_pruning_summary.csv
#   $PPOM_DIR/cppRATE_matrices/bacprune_rust_results.csv
#
# Outputs land under $OUT_DIR.

set -euo pipefail

POM_DIR=${POM_DIR:-/path/to/pom_run_output}
PPOM_DIR=${PPOM_DIR:-/path/to/ppom_run_output}
OUT_DIR=${OUT_DIR:-/nfs/research/jlees/jacqueline/thesis_code/test_variant_proportionality/results}

Rscript "$(dirname "$0")/test_variant_proportionality.R" \
    --pom_dir     "$POM_DIR" \
    --ppom_dir    "$PPOM_DIR" \
    --out_dir     "$OUT_DIR" \
    --ci_level    0.89 \
    --top_n_plots 20
