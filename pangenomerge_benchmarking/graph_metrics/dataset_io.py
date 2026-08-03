"""Loading the three clusterings (T, C, M) out of panaroo / pangenomerge output.

Shared by the metric validation and the determinism analysis. Only graph loading and
seqID remapping are borrowed from the pangenomerge package; the metrics themselves are
thesis analysis code and live in graph_similarity.py alongside this file.
"""

import sys
from pathlib import Path

import networkx as nx
import pandas as pd

REPO = "/hps/software/users/jlees/jacqueline/pangenome_merge"
if REPO not in sys.path:
    sys.path.insert(0, REPO)

from pangenomerge.custom_functions.manipulate_seqids import (  # noqa: E402
    get_seqIDs_in_nodes, indSID_to_allSID)
from pangenomerge.panaroo_functions.load_graphs import load_graphs  # noqa: E402


def load_panaroo_graph(gml_path):
    """Panaroo graph via the repo loader, which parses its ';'-joined attribute strings."""
    graphs, _, _ = load_graphs([str(gml_path)])
    return graphs[0]


def merged_seqid_map(gml_path):
    """{seqID: node} from a pangenomerge --mode test final_graph.gml.

    Read with plain nx.read_gml rather than load_graphs: pangenomerge writes 'centroid' as
    a list where panaroo writes a ';'-joined string, so the repo loader chokes on its own
    output.
    """
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


def truth_seqid_map(graph_all_dir):
    """T: the all-isolates panaroo run."""
    return get_seqIDs_in_nodes(
        load_panaroo_graph(Path(graph_all_dir) / "final_graph.gml"))


def gid_map(graph_all_dir, component_dir):
    """component clustering_id -> truth clustering_id, bridged by annotation_id.

    This is what makes the whole comparison possible: because the GFFs were called once
    and reused, the same gene carries the same annotation_id in both runs even though
    panaroo assigns different clustering_ids.
    """
    ga = pd.read_csv(Path(graph_all_dir) / "gene_data.csv",
                     usecols=["annotation_id", "clustering_id"])
    gi = pd.read_csv(Path(component_dir) / "gene_data.csv",
                     usecols=["annotation_id", "clustering_id"])
    m = ga.rename(columns={"clustering_id": "cid_all"}).merge(
        gi.rename(columns={"clustering_id": "cid_ind"}),
        on="annotation_id", how="left").dropna()
    return dict(zip(m["cid_ind"], m["cid_all"]))


def component_seqid_map(graph_all_dir, component_dirs):
    """C: every component cluster made globally unique by its graph index."""
    out = {}
    for idx, cdir in enumerate(component_dirs, start=1):
        G = indSID_to_allSID(load_panaroo_graph(Path(cdir) / "final_graph.gml"),
                             gid_map(graph_all_dir, cdir))
        for node in G.nodes():
            for sid in G.nodes[node].get("seqIDs", []):
                out[sid] = f"{node}_g{idx}"
    return out


def partition_signature(seqid_map):
    """Order-independent fingerprint of a clustering.

    Two runs share a signature iff they group the genes identically, regardless of what
    the nodes are named. This is the definitive determinism test -- node counts can match
    while the partitions differ.
    """
    import hashlib
    blocks = {}
    for sid, node in seqid_map.items():
        blocks.setdefault(node, set()).add(sid)
    canon = sorted(tuple(sorted(b)) for b in blocks.values())
    h = hashlib.sha256()
    for block in canon:
        h.update(("|".join(block) + "\n").encode())
    return h.hexdigest()[:16]
