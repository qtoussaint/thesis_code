#!/usr/bin/env python3
"""Seeded split of one strain cluster's isolates into halves and quarters (analysis A0).

Quarters are NESTED inside halves: the isolate list is shuffled once, cut into four
contiguous blocks, and half_1 = quarter_1 + quarter_2, half_2 = quarter_3 + quarter_4.
Nesting means the 2-way and 4-way merges see identical isolate content, so any difference
between them is attributable to the number of component graphs rather than to which
isolates landed where. It also gives a free consistency check: merging quarter_1 and
quarter_2 should closely reproduce the half_1 Panaroo graph.

Writes two things:
  <results>/<species>/splits/<part>.txt   one GFF path per line, fed straight to panaroo
  <code>/splits/<species>_membership.tsv  sample -> half, quarter (version-controlled record)
"""

import argparse
import random
from pathlib import Path

# species -> (cluster id, n isolates expected); mirrors config.sh
CLUSTERS = {
    "s_pneumoniae": "195",
    "m_tuberculosis": "2197",
    "s_aureus": "32",
}

ATB = Path("/nfs/research/jlees/jacqueline/atb_analyses/species_pangenomes")
CODE_DIR = Path(__file__).resolve().parent
RESULTS_ROOT = Path(
    "/nfs/research/jlees/jacqueline/thesis_results/pangenomerge_benchmarking/split_strain_datasets"
)


def sample_name(gff_path):
    """SAMEA123.fa.gff -> SAMEA123"""
    return gff_path.name.replace(".fa.gff", "")


def split_one(species, cluster, seed, results_root, code_dir):
    gff_dir = ATB / species / "results" / "ggcaller" / cluster / "GFF"
    gffs = sorted(gff_dir.glob("*.fa.gff"))
    if not gffs:
        raise SystemExit(f"no GFFs found in {gff_dir}")

    # sort first, then shuffle with a fixed seed -- filesystem order is not stable, so
    # sorting is what makes the split reproducible across machines and reruns
    rng = random.Random(seed)
    shuffled = list(gffs)
    rng.shuffle(shuffled)

    n = len(shuffled)
    # contiguous quarter boundaries; sizes differ by at most 1 when n % 4 != 0
    bounds = [round(i * n / 4) for i in range(5)]
    quarters = [shuffled[bounds[i]:bounds[i + 1]] for i in range(4)]
    halves = [quarters[0] + quarters[1], quarters[2] + quarters[3]]

    parts = {
        "all": shuffled,
        "half_1": halves[0], "half_2": halves[1],
        "quarter_1": quarters[0], "quarter_2": quarters[1],
        "quarter_3": quarters[2], "quarter_4": quarters[3],
    }

    split_dir = results_root / species / "splits"
    split_dir.mkdir(parents=True, exist_ok=True)
    for part, members in parts.items():
        (split_dir / f"{part}.txt").write_text(
            "\n".join(str(p) for p in members) + "\n")

    # version-controlled membership record
    membership = {}
    for i, q in enumerate(quarters, start=1):
        for p in q:
            membership[sample_name(p)] = (f"half_{1 if i <= 2 else 2}", f"quarter_{i}")

    out_tsv = code_dir / "splits" / f"{species}_membership.tsv"
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with open(out_tsv, "w") as fh:
        fh.write("sample\thalf\tquarter\n")
        for s in sorted(membership):
            fh.write(f"{s}\t{membership[s][0]}\t{membership[s][1]}\n")

    print(f"{species}: cluster {cluster}, {n} isolates -> "
          f"halves {len(halves[0])}/{len(halves[1])}, "
          f"quarters {'/'.join(str(len(q)) for q in quarters)}")
    print(f"  lists:      {split_dir}")
    print(f"  membership: {out_tsv}")
    return n


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--species", nargs="*", default=list(CLUSTERS),
                    help="subset of species to split (default: all three)")
    ap.add_argument("--results-root", type=Path, default=RESULTS_ROOT)
    ap.add_argument("--code-dir", type=Path, default=CODE_DIR)
    args = ap.parse_args()

    for species in args.species:
        if species not in CLUSTERS:
            raise SystemExit(f"unknown species: {species}")
        split_one(species, CLUSTERS[species], args.seed,
                  args.results_root, args.code_dir)


if __name__ == "__main__":
    main()
