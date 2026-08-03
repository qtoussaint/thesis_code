#!/usr/bin/env python3
"""Quantify run-to-run variation in the merge (analysis A3).

Scores every replicate of one configuration against the same truth, and reports:

  * mean +/- SD of pcARI / acARI / pcNMI / acNMI and the final node count
  * whether the merge *traces* agree (same number of pairs merged at each step)
  * whether the *partitions* agree (same genes grouped together)

The last two are separate questions and that distinction is the point: a trace can match
exactly while the partition differs, which is how the non-determinism was first noticed --
identical logs, different scores. Grouping by thread count tests whether the variation is
caused by parallel result ordering.

Any SD reported here is the noise floor: differences smaller than it, between methods or
configurations, are not interpretable.
"""

import argparse
import hashlib
import re
import statistics
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "graph_metrics"))
from graph_similarity import GraphScorer          # noqa: E402
from dataset_io import (merged_seqid_map, truth_seqid_map,   # noqa: E402
                        component_seqid_map, partition_signature)

SS = Path("/nfs/research/jlees/jacqueline/thesis_results/pangenomerge_benchmarking/"
          "split_strain_datasets")
DET = Path("/nfs/research/jlees/jacqueline/thesis_results/pangenomerge_benchmarking/"
           "determinism")

STEP_PAT = re.compile(
    r"Refinement merge (\d+) (identity\+context|representative|frameshift detection) step: "
    r"merging (\d+) pairs")


def trace_signature(run_log):
    """Fingerprint of the merge trace: how many pairs merged at each step, in order."""
    if not Path(run_log).exists():
        return None
    steps = STEP_PAT.findall(Path(run_log).read_text())
    return hashlib.sha256(str(steps).encode()).hexdigest()[:16], len(steps)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--species", default="s_pneumoniae")
    ap.add_argument("--split", default="split2")
    ap.add_argument("--out-csv", default=None)
    args = ap.parse_args()

    graph_all = SS / args.species / "panaroo_all"
    components = [l.strip() for l in
                  (SS / args.species / "paths" / f"{args.split}.txt").read_text().split()
                  if l.strip()]

    reps = sorted(d for d in DET.glob(f"{args.species}_{args.split}_t*_rep*")
                  if (d / "done.txt").exists())
    if not reps:
        raise SystemExit(f"no completed replicates under {DET}")
    print(f"{len(reps)} completed replicates")

    print("building baseline (truth, component, C*) once ...", flush=True)
    scorer = GraphScorer(truth_seqid_map(graph_all),
                         component_seqid_map(graph_all, components))
    print(f"  universe {len(scorer.universe)} genes; "
          f"component error {scorer.err_component:,.0f}; "
          f"reachable-floor error {scorer.err_cstar:,.0f}", flush=True)

    rows = []
    for d in reps:
        m = re.search(r"_t(\d+)_rep(\d+)$", d.name)
        threads, rep = int(m.group(1)), int(m.group(2))
        merged = merged_seqid_map(d / "final_graph.gml")
        scores = scorer.score(merged, label=d.name)
        tsig = trace_signature(d / "run.log")
        rows.append({
            "threads": threads, "rep": rep,
            "nodes": len(set(merged.values())),
            "pcARI": scores["pcARI"], "acARI": scores["acARI"],
            "pcNMI": scores["pcNMI"], "acNMI": scores["acNMI"],
            "err_merged": scores["err_merged"],
            "trace_sig": tsig[0] if tsig else None,
            "n_steps": tsig[1] if tsig else None,
            "partition_sig": partition_signature(merged),
        })
        print(f"  {d.name}: nodes={rows[-1]['nodes']} pcARI={scores['pcARI']:.6f}",
              flush=True)

    df = pd.DataFrame(rows).sort_values(["threads", "rep"])
    pd.set_option("display.width", 250); pd.set_option("display.max_columns", None)
    print("\n=== per replicate ===")
    print(df.to_string(index=False, float_format=lambda v: f"{v:.6f}"))

    print("\n=== variability by thread count ===")
    for threads, grp in df.groupby("threads"):
        n = len(grp)
        line = [f"threads={threads:<3d} n={n}"]
        for col in ("nodes", "pcARI", "acARI", "pcNMI", "acNMI"):
            vals = list(grp[col])
            sd = statistics.stdev(vals) if n > 1 else 0.0
            line.append(f"{col} {statistics.mean(vals):.6f}+/-{sd:.6f}")
        line.append(f"distinct partitions {grp['partition_sig'].nunique()}/{n}")
        line.append(f"distinct traces {grp['trace_sig'].nunique()}/{n}")
        print("  " + "  ".join(line))

    # the diagnostic that matters: identical traces but different partitions means the
    # variation is in *which* nodes merged, not how many
    same_trace = df["trace_sig"].nunique() == 1
    same_part = df["partition_sig"].nunique() == 1
    print(f"\nall traces identical: {same_trace} | all partitions identical: {same_part}")
    if same_trace and not same_part:
        print("  -> same merges performed, different genes grouped: tie-breaking differs")
    elif same_part:
        print("  -> fully deterministic across these runs")

    out = args.out_csv or str(DET / f"determinism_{args.species}_{args.split}.csv")
    df.to_csv(out, index=False)
    print(f"\nWrote {out}")


if __name__ == "__main__":
    main()
