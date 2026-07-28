#!/usr/bin/env python3
"""Calibration of the graph-similarity metrics (analysis A2d).

Feeds a set of candidate clusterings through the real production scorer and checks each
lands where theory says it must:

    random genes         standard ARI/AMI ~ 0 ; pc/ac strongly negative
    cluster-preserving   the null that matches the thesis argument -- Panaroo's clusters
      random             kept intact, only *which* clusters merge is randomised. Standard
                         ARI/AMI stay high (that is the whole problem); pc ~ 0.
    component C          pcARI = pcNMI = acARI = acNMI = 0 exactly
    best reachable C*    acARI = acNMI = 1 exactly; pcARI < 1 (the ceiling gap)
    real merge M         between 0 and 1 on both families
    truth T              pcARI = pcNMI = 1 exactly

The gap between C* and T is the point of the ac- family: it is the error no merge can fix.
"""

import argparse
import random
import sys
from pathlib import Path

import pandas as pd

REPO = "/hps/software/users/jlees/jacqueline/pangenome_merge"
sys.path.insert(0, REPO)

# metrics are thesis analysis code and live alongside this script; only graph loading and
# seqID remapping are reused from the pangenomerge package
sys.path.insert(0, str(Path(__file__).resolve().parent))
from graph_similarity import GraphScorer, best_reachable_merge  # noqa: E402

from pangenomerge.custom_functions.manipulate_seqids import (  # noqa: E402
    get_seqIDs_in_nodes, indSID_to_allSID)
from pangenomerge.panaroo_functions.load_graphs import load_graphs  # noqa: E402

REPORT_COLS = [
    "scenario", "n_seqIDs", "ARI", "AMI", "pcARI", "pcNMI", "acARI", "acNMI",
    "err_component", "err_merged", "err_cstar", "b_component", "c_component",
    "cstar_n_chained_components", "cstar_n_connected_components",
]


def load_graph(gml_path):
    """Panaroo graph, via the repo loader (parses its ';'-joined attribute strings)."""
    graphs, _, _ = load_graphs([str(gml_path)])
    return graphs[0]


def merged_seqid_map(gml_path):
    """{seqID: node} from a pangenomerge --mode test final_graph.gml.

    Read with plain nx.read_gml rather than load_graphs: pangenomerge writes 'centroid' as
    a list where panaroo writes a ';'-joined string, so the repo loader chokes on its own
    output. Mirrors node_to_seqids() in the split-strain analyse.py.
    """
    import networkx as nx
    G = nx.read_gml(str(gml_path))
    out = {}
    for n in G.nodes():
        sids = G.nodes[n].get("seqIDs", [])
        if isinstance(sids, str):
            sids = [sids] if sids else []
        elif not isinstance(sids, (list, tuple, set)):
            sids = [sids]
        for s in sids:
            out[str(s)] = n
    return out


def gid_map(graph_all_dir, component_dir):
    """component clustering_id -> truth clustering_id, bridged by annotation_id."""
    ga = pd.read_csv(Path(graph_all_dir) / "gene_data.csv",
                     usecols=["annotation_id", "clustering_id"])
    gi = pd.read_csv(Path(component_dir) / "gene_data.csv",
                     usecols=["annotation_id", "clustering_id"])
    m = ga.rename(columns={"clustering_id": "cid_all"}).merge(
        gi.rename(columns={"clustering_id": "cid_ind"}),
        on="annotation_id", how="left").dropna()
    return dict(zip(m["cid_ind"], m["cid_all"]))


def build_component_map(graph_all_dir, component_dirs):
    """C: every component cluster made globally unique by its graph index."""
    component_map = {}
    for idx, cdir in enumerate(component_dirs, start=1):
        G = indSID_to_allSID(load_graph(Path(cdir) / "final_graph.gml"),
                             gid_map(graph_all_dir, cdir))
        for node in G.nodes():
            for sid in G.nodes[node].get("seqIDs", []):
                component_map[sid] = f"{node}_g{idx}"
    return component_map


def cluster_preserving_random(component_map, merged_map, seed=0):
    """Randomise *which* component clusters merge, holding the merge pattern fixed.

    This is the null the thesis argument needs. A per-gene shuffle destroys Panaroo's
    clusters and so scores ~0 on plain ARI, which does not demonstrate anything. Here the
    component clusters stay intact and only their grouping is randomised, reproducing the
    inflated standard scores that motivated the corrected metrics.
    """
    # how many component clusters each merged node absorbed
    blocks = {}
    for sid, m_label in merged_map.items():
        c_label = component_map.get(sid)
        if c_label is not None:
            blocks.setdefault(m_label, set()).add(c_label)
    group_sizes = sorted((len(v) for v in blocks.values()), reverse=True)

    all_c = sorted({component_map[s] for s in merged_map if s in component_map})
    rng = random.Random(seed)
    rng.shuffle(all_c)

    assignment, pos = {}, 0
    for gi, size in enumerate(group_sizes):
        for c_label in all_c[pos:pos + size]:
            assignment[c_label] = f"rand_{gi}"
        pos += size
    for c_label in all_c[pos:]:            # any remainder keeps its own block
        assignment[c_label] = f"rand_solo_{c_label}"

    return {sid: assignment[component_map[sid]]
            for sid in merged_map if sid in component_map}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--graph-all", required=True,
                    help="panaroo output dir for the all-isolates truth run")
    ap.add_argument("--component-dirs", nargs="+", required=True)
    ap.add_argument("--merged-gml", required=True,
                    help="pangenomerge --mode test final_graph.gml (retains seqIDs)")
    ap.add_argument("--label", default="dataset")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out-csv", default=None)
    args = ap.parse_args()

    truth_map = get_seqIDs_in_nodes(load_graph(Path(args.graph_all) / "final_graph.gml"))
    component_map = build_component_map(args.graph_all, args.component_dirs)
    merged_map = merged_seqid_map(args.merged_gml)
    print(f"truth {len(set(truth_map.values()))} COGs / {len(truth_map)} genes; "
          f"component {len(set(component_map.values()))} clusters; "
          f"merged {len(set(merged_map.values()))} nodes")

    # C* over the genes common to all three, so it lines up with the scorer's alignment
    common = sorted((set(truth_map) & set(component_map) & set(merged_map)) - {"error"})
    cstar_labels, diag = best_reachable_merge(
        [component_map[s] for s in common], [truth_map[s] for s in common])
    cstar_map = dict(zip(common, cstar_labels))
    print(f"C*: {diag['n_cstar_clusters']} blocks from {diag['n_component_clusters']} "
          f"component clusters; {diag['n_chained_components']} of "
          f"{diag['n_connected_components']} components chain >1 truth COG")

    rng = random.Random(args.seed)
    scenarios = {
        "random_genes": {s: f"r{rng.randrange(len(set(truth_map.values())))}"
                         for s in common},
        "cluster_preserving_random": cluster_preserving_random(
            component_map, merged_map, seed=args.seed),
        "component_C": component_map,
        "best_reachable_Cstar": cstar_map,
        "real_merge_M": merged_map,
        "truth_T": truth_map,
    }

    # one scorer for all scenarios: the baseline (C* included) depends only on T and C, so
    # rebuilding it per scenario was six times the necessary work
    scorer = GraphScorer(truth_map, component_map)
    rows = []
    for name, m in scenarios.items():
        rows.append({"scenario": name, **scorer.score(m, label=name)})

    df = pd.DataFrame(rows)
    df = df[[c for c in REPORT_COLS if c in df.columns]]
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", None)
    print(f"\n=== metric calibration: {args.label} ===")
    print(df.to_string(index=False, float_format=lambda v: f"{v:.6g}"))

    out = args.out_csv or f"metric_validation_{args.label}.csv"
    df.to_csv(out, index=False)
    print(f"\nWrote {out}")


if __name__ == "__main__":
    main()
