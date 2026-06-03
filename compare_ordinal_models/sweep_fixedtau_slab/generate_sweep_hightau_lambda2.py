#!/usr/bin/env python3
"""Generate Stan models + SLURM run scripts for the HIGH-tau sweep with the
loose lambda prior (Cauchy(0, 2)), over slab in {3, 5}.

Sibling of generate_sweep_lambda2.py. Extends the fixed-tau sweep past the
previous tau=1 endpoint to map the prediction-accuracy plateau, and fills the
slab=5 x lambda=2 cells that the earlier sweeps did not cover.

Templated off:
  final_ordered_categorical_PPOM_free_cutpoints_wide_drift_fixedtau05_slab3
For each (tau, slab), writes:
  - PPOM_models/final_..._fixedtau{T}_slab{S}_lambda2.stan        (3 lines patched)
  - run_PPOM_models/final_..._fixedtau{T}_slab{S}_lambda2/run_{ds}.sh
"""

import re
import stat
from pathlib import Path

CODE_ROOT = Path("/nfs/research/jlees/jacqueline/thesis_code/compare_ordinal_models")
BASE_VARIANT = (
    "final_ordered_categorical_PPOM_free_cutpoints_wide_drift_fixedtau05_slab3"
)
BASE_STAN = CODE_ROOT / "PPOM_models" / f"{BASE_VARIANT}.stan"
BASE_RUN_DIR = CODE_ROOT / "run_PPOM_models" / BASE_VARIANT

TAUS = [1.5, 2, 3, 5]
SLABS = [3, 5]
LAMBDA = 2
DATASETS = [
    "02_spn_penicillin_MIC",
    "16_spn_penicillin_MIC_minimabinning",
]

CPUS = 80
MEM = "650G"
TIME = "12:00:00"


def num_to_slug(x):
    if float(x).is_integer():
        return str(int(x))
    return ("%g" % x).replace(".", "p")


def variant_name(tau, slab):
    return (
        "final_ordered_categorical_PPOM_free_cutpoints_wide_drift"
        f"_fixedtau{num_to_slug(tau)}_slab{slab}_lambda{LAMBDA}"
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
    out, n_lambda = re.subn(
        r"^(\s*to_vector\(lambda_variant\)\s*~\s*cauchy\(0,\s*)1(\);.*)$",
        lambda m: f"{m.group(1)}{LAMBDA}{m.group(2)}",
        out,
        count=1,
        flags=re.M,
    )
    if n_tau != 1 or n_slab != 1 or n_lambda != 1:
        raise RuntimeError(
            f"stan substitution failed: tau={n_tau}, slab={n_slab}, lambda={n_lambda}"
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
    for tau in TAUS:
        for slab in SLABS:
            new_variant = variant_name(tau, slab)

            stan_out = CODE_ROOT / "PPOM_models" / f"{new_variant}.stan"
            stan_out.write_text(patch_stan(stan_template, tau, slab))
            n_stan += 1

            run_dir = CODE_ROOT / "run_PPOM_models" / new_variant
            run_dir.mkdir(parents=True, exist_ok=True)
            for ds in DATASETS:
                sh_path = run_dir / f"run_{ds}.sh"
                sh_path.write_text(patch_run_script(run_templates[ds], new_variant))
                sh_path.chmod(
                    sh_path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH
                )
                n_sh += 1

            print(f"  generated {new_variant}")

    print(f"\nDone. Wrote {n_stan} .stan files and {n_sh} .sh run scripts.")


if __name__ == "__main__":
    main()
