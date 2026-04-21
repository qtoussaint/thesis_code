"""Aggregate benchmark results from the context-search iteration sweep.

For each iter_*/ directory under RESULTS_ROOT, extract:
- runtime (wall clock) and peak memory from time.log (GNU time -v output)
- node count from final_graph.gml
- actual iterations performed + pairs-per-iteration + convergence status
  from the SLURM stdout log

Writes a single results.csv next to the iter_*/ directories.
"""

import csv
import re
from pathlib import Path

import networkx as nx

RESULTS_ROOT = Path("/nfs/research/jlees/jacqueline/thesis_results/pangenomerge_benchmarking/context_search_iterations")
LOG_DIR = RESULTS_ROOT / "logs"
OUT_CSV = RESULTS_ROOT / "results.csv"

ITER_MERGE_RE = re.compile(r"Third merge iteration (\d+): merging (\d+) pairs\.")
ITER_STOP_RE = re.compile(r"Third merge iteration (\d+): (no hits survive staleness filter|no new mergeable pairs); stopping\.")
ITER_CAP_RE = re.compile(r"Third merge: reached max iterations \((-?\d+)\); stopping\.")


def parse_time_log(path: Path):
    """Parse `/usr/bin/time -v` output. Returns (seconds, max_rss_gb) or (None, None)."""
    if not path.is_file():
        return None, None
    secs = None
    rss_gb = None
    for line in path.read_text().splitlines():
        line = line.strip()
        if line.startswith("Elapsed (wall clock) time"):
            val = line.split(":", 1)[1].strip()
            # format: h:mm:ss or m:ss.ss
            parts = val.split(":")
            parts = [float(p) for p in parts]
            if len(parts) == 3:
                secs = parts[0] * 3600 + parts[1] * 60 + parts[2]
            elif len(parts) == 2:
                secs = parts[0] * 60 + parts[1]
            else:
                secs = parts[0]
        elif line.startswith("Maximum resident set size"):
            kb = int(line.split(":", 1)[1].strip())
            rss_gb = kb / (1024 * 1024)
    return secs, rss_gb


def count_nodes(gml_path: Path):
    if not gml_path.is_file():
        return None
    G = nx.read_gml(gml_path)
    return len(G.nodes())


def parse_iteration_log(log_path: Path):
    """Return (actual_iterations, converged, pairs_per_iteration_str)."""
    if not log_path.is_file():
        return None, None, ""
    text = log_path.read_text()
    pairs_by_iter = {}
    for m in ITER_MERGE_RE.finditer(text):
        pairs_by_iter[int(m.group(1))] = int(m.group(2))
    converged = bool(ITER_STOP_RE.search(text))
    capped = bool(ITER_CAP_RE.search(text))
    if pairs_by_iter:
        actual = max(pairs_by_iter.keys())
    elif converged:
        # stopping log appears even when no merges happen (e.g. iteration 1 has no hits)
        m = ITER_STOP_RE.search(text)
        actual = int(m.group(1)) if m else 0
    else:
        actual = None
    pairs_str = ";".join(f"{k}:{pairs_by_iter[k]}" for k in sorted(pairs_by_iter))
    return actual, converged if converged or capped else None, pairs_str


def main():
    rows = []
    for iter_dir in sorted(RESULTS_ROOT.glob("iter_*")):
        tag = iter_dir.name[len("iter_"):]
        cap = -1 if tag == "unlimited" else int(tag)
        secs, rss_gb = parse_time_log(iter_dir / "time.log")
        nodes = count_nodes(iter_dir / "final_graph.gml")
        log_path = LOG_DIR / f"iter_{tag}.out"
        actual_iters, converged, pairs_str = parse_iteration_log(log_path)
        rows.append({
            "iteration_cap": cap,
            "actual_iterations": actual_iters,
            "converged": converged,
            "runtime_sec": secs,
            "max_mem_gb": rss_gb,
            "node_count": nodes,
            "pairs_per_iteration": pairs_str,
        })

    with OUT_CSV.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["iteration_cap", "actual_iterations", "converged",
                                                "runtime_sec", "max_mem_gb", "node_count",
                                                "pairs_per_iteration"])
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} rows to {OUT_CSV}")


if __name__ == "__main__":
    main()
