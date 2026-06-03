#!/usr/bin/env python3
"""Generate estimated-scale (adaptive) lasso & ridge PPOM comparators.

The existing wide_drift lasso/ridge use FIXED prior scales (double_exponential(0,1),
normal(0,5)), while the horseshoe estimates its global scale. These models make the
lasso/ridge scales data-adaptive too, so all three priors share the same
scale-mixture-of-normals skeleton and differ only in how the scale is generated.

For each prior family we build two parameterizations of the SAME marginal prior:
  ridge: estscale (centered, beta ~ normal(0, sigma)) and
         estscale_ncp (non-centered, beta = z * sigma)            <- adaptive ridge
  lasso: estscale (direct, beta ~ double_exponential(0, b)) and
         estscale_mixture (Park & Casella scale mixture)          <- adaptive lasso
Each in both standardization variants (X_std and raw-genotype no_centering).

Skeletons (already wide_drift: cutpoint_prior_sd=1.5, correct standardization & GQ):
  PPOM_models/final_..._wide_drift_{ridge,lasso}{,_no_centering}.stan
The four are structurally identical except the coefficient name (beta_variant_std vs
beta_variant_raw) and the single prior line; we patch the params/tparams/model blocks.

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

# Hyperpriors on the estimated scales (the key modeling choice; weakly informative,
# centered near the old fixed scales). Centered and non-centered/mixture pairs share
# the same hyperprior, so within a prior family they are the same marginal prior.
RIDGE_SCALE_PRIOR = "sigma_beta_variant ~ student_t(3, 0, 2);"
LASSO_SCALE_PRIOR = "b_lasso ~ student_t(3, 0, 1);"


def replace_once(text, old, new, what):
    n = text.count(old)
    if n != 1:
        raise RuntimeError(f"{what}: expected 1 occurrence of {old!r}, found {n}")
    return text.replace(old, new)


def insert_after_tparams(text, block):
    anchor = "transformed parameters {\n"
    return replace_once(text, anchor, anchor + block, "transformed parameters header")


def set_header(text, header):
    if "data {" not in text:
        raise RuntimeError("could not locate 'data {' to replace header")
    return header + "\n" + text[text.index("data {"):]


def patch_ridge_cp(text, coef, header):
    text = set_header(text, header)
    text = replace_once(
        text,
        f"  matrix[V, K-1] {coef};\n",
        f"  matrix[V, K-1] {coef};\n  real<lower=0> sigma_beta_variant;\n",
        "ridge cp params",
    )
    text = replace_once(
        text,
        f"  to_vector({coef}) ~ normal(0, 5);\n",
        f"  to_vector({coef}) ~ normal(0, sigma_beta_variant);\n"
        f"  {RIDGE_SCALE_PRIOR}\n",
        "ridge cp prior",
    )
    return text


def patch_ridge_ncp(text, coef, header):
    text = set_header(text, header)
    text = replace_once(
        text,
        f"  matrix[V, K-1] {coef};\n",
        "  matrix[V, K-1] z_variant;\n  real<lower=0> sigma_beta_variant;\n",
        "ridge ncp params",
    )
    text = insert_after_tparams(
        text, f"  matrix[V, K-1] {coef} = z_variant * sigma_beta_variant;\n\n"
    )
    text = replace_once(
        text,
        f"  to_vector({coef}) ~ normal(0, 5);\n",
        "  to_vector(z_variant) ~ std_normal();\n"
        f"  {RIDGE_SCALE_PRIOR}\n",
        "ridge ncp prior",
    )
    return text


def patch_lasso_direct(text, coef, header):
    text = set_header(text, header)
    text = replace_once(
        text,
        f"  matrix[V, K-1] {coef};\n",
        f"  matrix[V, K-1] {coef};\n  real<lower=0> b_lasso;\n",
        "lasso direct params",
    )
    text = replace_once(
        text,
        f"  to_vector({coef}) ~ double_exponential(0, 1);\n",
        f"  to_vector({coef}) ~ double_exponential(0, b_lasso);\n"
        f"  {LASSO_SCALE_PRIOR}\n",
        "lasso direct prior",
    )
    return text


def patch_lasso_mixture(text, coef, header):
    text = set_header(text, header)
    text = replace_once(
        text,
        f"  matrix[V, K-1] {coef};\n",
        "  matrix[V, K-1] z_variant;\n"
        "  matrix<lower=0>[V, K-1] tau2_local;\n"
        "  real<lower=0> b_lasso;\n",
        "lasso mixture params",
    )
    text = insert_after_tparams(
        text,
        f"  matrix[V, K-1] {coef};\n"
        "  for (k in 1:(K-1))\n"
        "    for (v in 1:V)\n"
        f"      {coef}[v, k] = z_variant[v, k] * sqrt(tau2_local[v, k]);\n\n",
    )
    text = replace_once(
        text,
        f"  to_vector({coef}) ~ double_exponential(0, 1);\n",
        "  to_vector(z_variant) ~ std_normal();\n"
        "  to_vector(tau2_local) ~ exponential(1 / (2 * square(b_lasso)));\n"
        f"  {LASSO_SCALE_PRIOR}\n",
        "lasso mixture prior",
    )
    return text


# Accurate headers per variant (the "_GENO" placeholder is filled per standardization).
H_RIDGE_CP = (
    "// Estimated-scale ridge comparator for the wide_drift PPOM sweep, centered\n"
    "// parameterization. Like wide_drift_ridge but the fixed normal(0, 5) is\n"
    "// replaced by normal(0, sigma_beta_variant) with sigma_beta_variant ~\n"
    "// student_t(3, 0, 2), so the L2 penalty is learned from the data._GENO\n"
)
H_RIDGE_NCP = (
    "// Adaptive (estimated-scale) ridge comparator for the wide_drift PPOM sweep,\n"
    "// non-centered parameterization (ncp): beta = z_variant * sigma_beta_variant,\n"
    "// z_variant ~ std_normal, sigma_beta_variant ~ student_t(3, 0, 2). Same marginal\n"
    "// prior as wide_drift_ridge_estscale; better ADVI/HMC geometry._GENO\n"
)
H_LASSO_DIRECT = (
    "// Estimated-scale lasso comparator for the wide_drift PPOM sweep, direct\n"
    "// parameterization. Like wide_drift_lasso but the fixed double_exponential(0, 1)\n"
    "// is replaced by double_exponential(0, b_lasso) with b_lasso ~ student_t(3, 0, 1),\n"
    "// so the L1 penalty is learned from the data._GENO\n"
)
H_LASSO_MIXTURE = (
    "// Adaptive lasso comparator for the wide_drift PPOM sweep (Park & Casella scale\n"
    "// mixture, non-centered): beta = z * sqrt(tau2), z ~ std_normal, tau2 ~\n"
    "// exponential(1 / (2 b_lasso^2)), b_lasso ~ student_t(3, 0, 1). Marginal beta ~\n"
    "// Laplace(0, b_lasso) -- same prior as wide_drift_lasso_estscale, better geometry._GENO\n"
)

GENO_NOTE = {
    "centered": "",
    "no_centering": (
        " On RAW genotype\n"
        "// (no centering/scaling): the prior acts on beta_variant_raw / variant_matrix,\n"
        "// and the estimated scale adapts to the raw-allele effect size."
    ),
}

# (prior family, standardization) -> (skeleton basename, coef name)
SKELETON = {
    ("ridge", "centered"): (
        "final_ordered_categorical_PPOM_free_cutpoints_wide_drift_ridge",
        "beta_variant_std",
    ),
    ("ridge", "no_centering"): (
        "final_ordered_categorical_PPOM_free_cutpoints_wide_drift_ridge_no_centering",
        "beta_variant_raw",
    ),
    ("lasso", "centered"): (
        "final_ordered_categorical_PPOM_free_cutpoints_wide_drift_lasso",
        "beta_variant_std",
    ),
    ("lasso", "no_centering"): (
        "final_ordered_categorical_PPOM_free_cutpoints_wide_drift_lasso_no_centering",
        "beta_variant_raw",
    ),
}

# variant token -> (family, parameterization token, patch fn, header)
PARAM = {
    "ridge_estscale": ("ridge", patch_ridge_cp, H_RIDGE_CP, "RidgeEScp"),
    "ridge_estscale_ncp": ("ridge", patch_ridge_ncp, H_RIDGE_NCP, "RidgeESncp"),
    "lasso_estscale": ("lasso", patch_lasso_direct, H_LASSO_DIRECT, "LassoESdir"),
    "lasso_estscale_mixture": ("lasso", patch_lasso_mixture, H_LASSO_MIXTURE, "LassoESmix"),
}

PREFIX = "final_ordered_categorical_PPOM_free_cutpoints_wide_drift"


def patch_run_script(template, new_variant, nickname):
    out = template.replace(RUN_TEMPLATE_VARIANT, new_variant)
    out = out.replace("freeCutsWDFTPPOM", nickname)
    out = re.sub(r"^#SBATCH --cpus-per-task=\d+$", f"#SBATCH --cpus-per-task={CPUS}", out, count=1, flags=re.M)
    out = re.sub(r"^#SBATCH --mem=\d+G$", f"#SBATCH --mem={MEM}", out, count=1, flags=re.M)
    out = re.sub(r"^#SBATCH --time=\d+:\d+:\d+$", f"#SBATCH --time={TIME}", out, count=1, flags=re.M)
    out = re.sub(r'^THREADS="--threads \d+"$', f'THREADS="--threads {CPUS}"', out, count=1, flags=re.M)
    return out


def main():
    run_templates = {ds: (RUN_TEMPLATE_DIR / f"run_{ds}.sh").read_text() for ds in DATASETS}

    n_stan = 0
    n_sh = 0
    for param_token, (family, patch_fn, header_tmpl, nick_core) in PARAM.items():
        for geno in ("centered", "no_centering"):
            skel_base, coef = SKELETON[(family, geno)]
            suffix = param_token if geno == "centered" else f"{param_token}_no_centering"
            new_variant = f"{PREFIX}_{suffix}"
            nickname = f"freeCutsWD{nick_core}{'NC' if geno == 'no_centering' else ''}PPOM"

            header = header_tmpl.replace("_GENO", GENO_NOTE[geno])
            skel = (MODELS_DIR / f"{skel_base}.stan").read_text()
            (MODELS_DIR / f"{new_variant}.stan").write_text(patch_fn(skel, coef, header))
            n_stan += 1

            run_dir = RUN_ROOT / new_variant
            run_dir.mkdir(parents=True, exist_ok=True)
            for ds in DATASETS:
                sh_path = run_dir / f"run_{ds}.sh"
                sh_path.write_text(patch_run_script(run_templates[ds], new_variant, nickname))
                sh_path.chmod(sh_path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
                n_sh += 1

            print(f"  generated {new_variant}")

    print(f"\nDone. Wrote {n_stan} .stan files and {n_sh} .sh run scripts.")


if __name__ == "__main__":
    main()
