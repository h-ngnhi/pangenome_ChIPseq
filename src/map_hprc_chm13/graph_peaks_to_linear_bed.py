#!/usr/bin/env python3
import sys
import json
import zipfile
import numpy as np
from pathlib import Path

if len(sys.argv) != 4:
    sys.stderr.write(
        "Usage: graph_peaks_to_linear_bed.py "
        "<all_peaks.intervalcollection> <graphs_dir_with_chr*_linear_path.interval> "
        "<out.bed>\n"
    )
    sys.exit(1)

intervalcoll_path = sys.argv[1]          # multi-chrom intervalcollection
graphs_dir = Path(sys.argv[2])           # dir containing chr*_linear_path.interval
out_bed = sys.argv[3]

# Cache: chrom -> (min_node:int, node_to_distance:np.ndarray)
mapping_cache = {}

def load_mapping_for_chrom(chrom: str):
    """
    Load <chrom>_linear_path.interval and cache it.
    Expected file name: <chrom>_linear_path.interval (e.g. chr19_linear_path.interval)
    """
    if chrom in mapping_cache:
        return mapping_cache[chrom]

    interval_path = graphs_dir / f"{chrom}_linear_path.interval"
    if not interval_path.exists():
        raise FileNotFoundError(f"Missing linear mapping: {interval_path}")

    with zipfile.ZipFile(interval_path) as zf:
        with zf.open("min_node.npy") as f:
            min_node = int(np.load(f))
        with zf.open("node_to_distance.npy") as f:
            node_to_distance = np.load(f)

    mapping_cache[chrom] = (min_node, node_to_distance)
    return mapping_cache[chrom]

def node_to_dist(node_id: int, min_node: int, node_to_distance: np.ndarray):
    """Return linear distance for a node ID, or None if out of range / unmapped."""
    idx = node_id - min_node
    if idx < 0 or idx >= node_to_distance.shape[0]:
        return None
    d = node_to_distance[idx]
    if d < 0:
        return None
    return float(d)

with open(intervalcoll_path) as fin, open(out_bed, "w") as fout:
    for idx, line in enumerate(fin):
        line = line.strip()
        if not line:
            continue

        rec = json.loads(line)

        chrom = rec.get("chromosome")
        if not chrom:
            continue

        nodes = rec.get("region_paths", [])
        if not nodes:
            continue

        try:
            min_node, node_to_distance = load_mapping_for_chrom(chrom)
        except FileNotFoundError as e:
            # Skip peaks whose chromosome mapping isn't available
            sys.stderr.write(str(e) + "\n")
            continue

        dists = []
        for nid in nodes:
            d = node_to_dist(int(nid), min_node, node_to_distance)
            if d is not None:
                dists.append(d)

        if not dists:
            continue

        # --- improvement: avoid start=0 artifacts if there are real (>0) distances ---
        pos = [d for d in dists if d > 0]
        use = pos if pos else dists

        start = int(min(use))
        end = int(max(use)) + 1

        avg_q = float(rec.get("average_q_value", 0.0))
        name = f"peak_{idx}"

        fout.write(f"{chrom}\t{start}\t{end}\t{name}\t{avg_q:.2f}\n")
