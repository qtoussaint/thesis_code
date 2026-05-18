#!/usr/bin/env python3
"""Generate Stan models + SLURM run scripts for the fixed-tau x slab sweep.

Templated off the existing variant:
  final_ordered_categorical_PPOM_free_cutpoints_wide_drift_fixedtau05_slab3
For each (tau, slab) in TAUS x SLABS, writes:
  - PPOM_models/final_..._fixedtau{T}_slab{S}.stan         (2 lines patched)
  - run_PPOM_models/final_..._fixedtau{T}_slab{S}/run_{ds}.sh
      (variant slug replaced + SBATCH/threads resources overridden)

Idempotent: re-running overwrites the generated files but leaves the base
templates and the existing fixedtau05_slab3 results dir untouched.
"""

import os
import re
import stat
from itertools import product
from pathlib import Path

CODE_ROOT = Path("/nfs/research/jlees/jacqueline/thesis_code/compare_ordinal_models")
BASE_VARIANT = (
    "final_ordered_categorical_PPOM_free_cutpoints_wide_drift_fixedtau05_slab3"
)
BASE_STAN = CODE_ROOT / "PPOM_models" / f"{BASE_VARIANT}.stan"
BASE_RUN_DIR = CODE_ROOT / "run_PPOM_models" / BASE_VARIANT

TAUS = [0.001, 0.01, 0.05, 1]
SLABS = [3, 5]
DATASETS = [
    "02_spn_penicillin_MIC",
    "16_spn_penicillin_MIC_minimabinning",
]

# Resource overrides applied to every generated run script.
CPUS = 80
MEM = "650G"
TIME = "12:00:00"


def num_to_slug(x):
    """0.001 -> '0p001'; 0.05 -> '0p05'; 1 -> '1'; 50 -> '50'."""
    if float(x).is_integer():
        return str(int(x))
    return ("%g" % x).replace(".", "p")


def variant_name(tau, slab):
    return (
        "final_ordered_categorical_PPOM_free_cutpoints_wide_drift"
        f"_fixedtau{num_to_slug(tau)}_slab{num_to_slug(slab)}"
    )


def patch_stan(template, tau, slab):
    out, n_tau = re.subn(
        r"^(\s*real tau\s*=\s*)0\.05(\s*;.*)$",
        lambda m: f"{m.group(1)}{tau}{m.group(2)}",
        template,
        count=1,
        flags=re.M,
    )
    out, n_slab = re.subn(
        r"^(\s*real slab_scale\s*=\s*)3(\s*;.*)$",
        lambda m: f"{m.group(1)}{slab}{m.group(2)}",
        out,
        count=1,
        flags=re.M,
    )
    if n_tau != 1 or n_slab != 1:
        raise RuntimeError(
            f"stan template substitution failed: tau={n_tau}, slab={n_slab}"
        )
    return out


def patch_run_script(template, new_variant):
    out = template.replace(BASE_VARIANT, new_variant)
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
    stan_template = BASE_STAN.read_text()
    run_templates = {ds: (BASE_RUN_DIR / f"run_{ds}.sh").read_text() for ds in DATASETS}

    n_stan = 0
    n_sh = 0
    for tau, slab in product(TAUS, SLABS):
        new_variant = variant_name(tau, slab)

        # Stan
        stan_out = CODE_ROOT / "PPOM_models" / f"{new_variant}.stan"
        stan_out.write_text(patch_stan(stan_template, tau, slab))
        n_stan += 1

        # Run scripts
        run_dir = CODE_ROOT / "run_PPOM_models" / new_variant
        run_dir.mkdir(parents=True, exist_ok=True)
        for ds in DATASETS:
            sh_path = run_dir / f"run_{ds}.sh"
            sh_path.write_text(patch_run_script(run_templates[ds], new_variant))
            sh_path.chmod(sh_path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
            n_sh += 1

        print(f"  generated {new_variant}")

    print(f"\nDone. Wrote {n_stan} .stan files and {n_sh} .sh run scripts.")


if __name__ == "__main__":
    main()
