# parameter_sweeps_prediction

Out-of-sample (prediction) versions of the four parameter sweeps that
`paper_figures/parameter_sweeps/` currently shows **in-sample** (from
`inference_ppc/`). Each sweep model is refit on a training split and scored on a
held-out **test** split (`--analysis_type prediction`), so the figures report
generalization accuracy instead of optimistic in-sample fit.

Scope: 43 prediction models / 75 runs (random split only; accuracy-only figures).

| Sweep | Models | Datasets | Runs |
|---|---|---|---|
| Regularization prior | 12 lasso/ridge variants | 02 (std K=8), 16 (minima K=5) | 24 |
| Tau sweep | 20 fixedtau×slab×λ variants | 02, 16 | 40 |
| Logistic alpha | 8 alphasd variants | 01 binary | 8 |
| Continuous alpha | 3 alphasd variants | 03 continuous | 3 |

## Layout

```
generate_prediction_stan.R         # PPOM inference->prediction stan transform (+ self-validation, stanc)
generate_alpha_prediction_stan.R   # logistic/continuous prediction stans (one-literal alpha_prior_sd edit)
generate_run_scripts.R             # emits run_scripts/ + submit_all.sh
stan_models/                       # 43 generated *_prediction.stan
run_scripts/<run_name>/run_<dataset>.sh   ;   run_scripts/submit_all.sh
figures/
  prediction_sweep_helpers.R       # shared reader + restricted-bACC resolver
  regularization_prior_prediction_summary.R
  tau_sweep_prediction_summary.R
  alpha_sweep_prediction_summary.R       # -> logistic + continuous figures
```

Results land under
`thesis_results/compare_ordinal_models/parameter_sweeps_prediction/`
(`runs/<run_name>/<dataset>/prediction_results/` and `figures/`).

## Run order

```bash
mamba activate gwas_pipeline

# 1. (re)generate Stan models + run scripts  (idempotent; safe to re-run)
Rscript generate_prediction_stan.R          # 32 PPOM models, validates + stanc-checks
Rscript generate_alpha_prediction_stan.R    # 11 logistic/continuous models, stanc-checks
Rscript generate_run_scripts.R              # 75 SLURM scripts + submit_all.sh

# 2. launch the fits on the cluster (~75 jobs, 6-12 h each)
bash run_scripts/submit_all.sh

# 3. once results land, render the figures (gaps render for any run not yet done)
Rscript figures/regularization_prior_prediction_summary.R
Rscript figures/tau_sweep_prediction_summary.R
Rscript figures/alpha_sweep_prediction_summary.R
```

Notes:
- Run scripts carry `--resume`, so re-submitting resumes from the last checkpoint.
- The Stan generator validates by regenerating the existing hand-written
  `*_tau5_prediction.stan` twin and asserting the functional blocks match, then
  `stanc`-syntax-checks every output; a residual-token check guards against any
  in-sample artifact surviving the transform.
- Figures read out-of-sample `prediction_results/prediction_accuracy_metrics.csv`;
  the dashed reference line is the production prediction PPOM / logistic /
  continuous run (random split). PPOM bACC cells that hit a missing-category
  div-by-zero are recomputed over present-only categories and starred.
