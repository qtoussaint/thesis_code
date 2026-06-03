#!/usr/bin/env python3
"""Generate wide_drift lasso and ridge PPOM models + SLURM run scripts.

These are prior comparators for the horseshoe sweep. The variant-effect prior
is copied verbatim from the existing old models:
  - lasso: to_vector(beta_variant_std) ~ double_exponential(0, 1)
  - ridge: to_vector(beta_variant_std) ~ normal(0, 5)
Both already use standardized genotype (X_std) + sd_variant back-transform and
the same y_rep / heritability generated quantities as the sweep baseline. The
ONLY change needed to drop them into the wide_drift family is widening the
cutpoint prior SD from 0.25 to 1.5 (that is all "wide_drift" means), so results
are directly comparable to the rest of the sweep.

Source stans:
  final_ordered_categorical_PPOM_free_cutpoints_{lasso,ridge}_with_centering.stan
Run scripts templated off:
  final_ordered_categorical_PPOM_free_cutpoints_wide_drift_fixedtau05_slab3
"""

import re
import stat
from pathlib import Path

CODE_ROOT = Path("/nfs/research/jlees/jacqueline/thesis_code/compare_ordinal_models")
MODELS_DIR = CODE_ROOT / "PPOM_models"
RUN_ROOT = CODE_ROOT / "run_PPOM_models"

RUN_TEMPLATE_VARIANT = (
    "final_ordered_categorical_PPOM_free_cutpoints_wide_drift_fixedtau05_slab3"
)
RUN_TEMPLATE_DIR = RUN_ROOT / RUN_TEMPLATE_VARIANT

DATASETS = [
    "02_spn_penicillin_MIC",
    "16_spn_penicillin_MIC_minimabinning",
]

CPUS = 80
MEM = "650G"
TIME = "12:00:00"

# Accurate, sweep-consistent header (the source files' headers say "tight
# Gaussian prior", which is wrong once we widen cutpoint_prior_sd to 1.5).
HEADER = {
    "lasso": (
        "// Lasso (Laplace) comparator for the wide_drift PPOM sweep. Identical to\n"
        "// final_ordered_categorical_PPOM_free_cutpoints_wide_drift.stan except the\n"
        "// regularized horseshoe on beta_variant_std is replaced by a\n"
        "// double_exponential(0, 1) prior.\n"
    ),
    "ridge": (
        "// Ridge (Gaussian) comparator for the wide_drift PPOM sweep. Identical to\n"
        "// final_ordered_categorical_PPOM_free_cutpoints_wide_drift.stan except the\n"
        "// regularized horseshoe on beta_variant_std is replaced by a normal(0, 5)\n"
        "// prior.\n"
    ),
    "lasso_nc": (
        "// Lasso (Laplace) comparator for the wide_drift PPOM sweep, on RAW genotype\n"
        "// (no centering/scaling). Same as the wide_drift_lasso model except the\n"
        "// effects act on the unstandardized variant_matrix (beta_variant_raw) rather\n"
        "// than on standardized X_std, so the double_exponential(0, 1) prior is\n"
        "// applied to raw-allele effects.\n"
    ),
    "ridge_nc": (
        "// Ridge (Gaussian) comparator for the wide_drift PPOM sweep, on RAW genotype\n"
        "// (no centering/scaling). Same as the wide_drift_ridge model except the\n"
        "// effects act on the unstandardized variant_matrix (beta_variant_raw) rather\n"
        "// than on standardized X_std, so the normal(0, 5) prior is applied to\n"
        "// raw-allele effects.\n"
    ),
}

# name -> (source stan basename, new wide_drift variant, run nickname token)
PRIORS = {
    "lasso": (
        "final_ordered_categorical_PPOM_free_cutpoints_lasso_with_centering",
        "final_ordered_categorical_PPOM_free_cutpoints_wide_drift_lasso",
        "freeCutsWDLassoPPOM",
    ),
    "ridge": (
        "final_ordered_categorical_PPOM_free_cutpoints_ridge_with_centering",
        "final_ordered_categorical_PPOM_free_cutpoints_wide_drift_ridge",
        "freeCutsWDRidgePPOM",
    ),
    "lasso_nc": (
        "final_ordered_categorical_PPOM_free_cutpoints_lasso_no_centering",
        "final_ordered_categorical_PPOM_free_cutpoints_wide_drift_lasso_no_centering",
        "freeCutsWDLassoNCPPOM",
    ),
    "ridge_nc": (
        "final_ordered_categorical_PPOM_free_cutpoints_ridge_no_centering",
        "final_ordered_categorical_PPOM_free_cutpoints_wide_drift_ridge_no_centering",
        "freeCutsWDRidgeNCPPOM",
    ),
}


def patch_stan(template, name):
    # Widen the cutpoint prior SD to match the wide_drift family (0.25 -> 1.5).
    out, n = re.subn(
        r"^(\s*real cutpoint_prior_sd\s*=\s*)0\.25(\s*;.*)$",
        lambda m: f"{m.group(1)}1.5{m.group(2)}",
        template,
        count=1,
        flags=re.M,
    )
    if n != 1:
        raise RuntimeError(f"cutpoint_prior_sd substitution failed: n={n}")
    # Replace the (now-inaccurate) source header with an accurate one.
    if "data {" not in out:
        raise RuntimeError("could not locate 'data {' to replace header")
    out = HEADER[name] + "\n" + out[out.index("data {"):]
    return out


def patch_run_script(template, new_variant, nickname):
    out = template.replace(RUN_TEMPLATE_VARIANT, new_variant)
    # Distinct analysis nickname (template uses freeCutsWDFTPPOM).
    out = out.replace("freeCutsWDFTPPOM", nickname)
    out = re.sub(
        r"^#SBATCH --cpus-per-task=\d+$",
        f"#SBATCH --cpus-per-task={CPUS}",
        out,
        count=1,
        flags=re.M,
    )
    out = re.sub(
        r"^#SBATCH --mem=\d+G$", f"#SBATCH --mem={MEM}", out, count=1, flags=re.M
    )
    out = re.sub(
        r"^#SBATCH --time=\d+:\d+:\d+$",
        f"#SBATCH --time={TIME}",
        out,
        count=1,
        flags=re.M,
    )
    out = re.sub(
        r'^THREADS="--threads \d+"$',
        f'THREADS="--threads {CPUS}"',
        out,
        count=1,
        flags=re.M,
    )
    return out


def main():
    run_templates = {
        ds: (RUN_TEMPLATE_DIR / f"run_{ds}.sh").read_text() for ds in DATASETS
    }

    n_stan = 0
    n_sh = 0
    for name, (src_base, new_variant, nickname) in PRIORS.items():
        src_stan = (MODELS_DIR / f"{src_base}.stan").read_text()
        (MODELS_DIR / f"{new_variant}.stan").write_text(patch_stan(src_stan, name))
        n_stan += 1

        run_dir = RUN_ROOT / new_variant
        run_dir.mkdir(parents=True, exist_ok=True)
        for ds in DATASETS:
            sh_path = run_dir / f"run_{ds}.sh"
            sh_path.write_text(patch_run_script(run_templates[ds], new_variant, nickname))
            sh_path.chmod(
                sh_path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH
            )
            n_sh += 1

        print(f"  generated {new_variant}")

    print(f"\nDone. Wrote {n_stan} .stan files and {n_sh} .sh run scripts.")


if __name__ == "__main__":
    main()
