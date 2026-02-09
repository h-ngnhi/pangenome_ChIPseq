#!/usr/bin/env python3
import sys
import json
import math

# Usage: vg view -aj treatment_alignments.gam \
#          | python3 match_reads_to_peaks.py all_peaks.intervalcollection

if len(sys.argv) != 2:
    sys.stderr.write(
        f"Usage:\n  vg view -aj treatment_alignments.gam | {sys.argv[0]} all_peaks.intervalcollection\n"
    )
    sys.exit(1)

peaks_file = sys.argv[1]
out_file = "results/146/vg_giraffe_vcfbub/matched_reads_to_TEpeaks.tsv"

# --- peaks of interest (indices in all_peaks.intervalcollection, 0-based) ---
# If your peak_index in the TSV is 1-based, change these to {73, 6473, ...} etc.
target_indices = {74, 6474, 6572, 7399, 9758, 17586}

# --- 1) Load region_paths for those peaks ---
peak_paths = {}  # idx -> list of node_ids
with open(peaks_file) as f:
    for idx, line in enumerate(f):
        if idx in target_indices:
            rec = json.loads(line)
            nodes = rec.get("region_paths", [])
            peak_paths[idx] = [int(n) for n in nodes]

if not peak_paths:
    sys.stderr.write("No target peaks loaded. Check indices / peaks file.\n")
    sys.exit(1)

sys.stderr.write(f"Loaded {len(peak_paths)} target peaks: {sorted(peak_paths.keys())}\n")

# --- helper: longest contiguous match between two node sequences ---
def longest_contiguous_match(read_nodes, peak_nodes):
    """
    Return length of longest contiguous alignment of read_nodes within peak_nodes.
    read_nodes and peak_nodes are lists of node_ids.
    """
    max_len = 0
    n_r = len(read_nodes)
    n_p = len(peak_nodes)
    if n_r == 0 or n_p == 0:
        return 0

    # naive but fine for short read paths / few peaks
    for start in range(n_p):
        if peak_nodes[start] != read_nodes[0]:
            continue
        j = 0
        while start + j < n_p and j < n_r and peak_nodes[start + j] == read_nodes[j]:
            j += 1
        if j > max_len:
            max_len = j
    return max_len

# --- 2) Open output file + write header (and also print header) ---
out_fh = open(out_file, "w")
header = ["read_name", "peak_index", "orientation", "match_nodes", "read_nodes"]
print("\t".join(header))
out_fh.write("\t".join(header) + "\n")

# --- 3) Stream GAM JSON from stdin ---
gam_fh = sys.stdin

for line in gam_fh:
    line = line.strip()
    if not line:
        continue
    aln = json.loads(line)

    name = aln.get("name", "")
    path = aln.get("path", {})
    mappings = path.get("mapping", [])
    if not mappings:
        continue

    # sort mappings by rank to ensure proper order
    try:
        mappings.sort(key=lambda m: int(m.get("rank", 0)))
    except Exception:
        pass

    read_nodes = []
    is_rev_flags = []
    for m in mappings:
        pos = m.get("position", {})
        nid = pos.get("node_id", None)
        if nid is None:
            continue
        read_nodes.append(int(nid))
        is_rev_flags.append(bool(pos.get("is_reverse", False)))

    if not read_nodes:
        continue

    read_len = len(read_nodes)
    required = math.ceil(0.6 * read_len)  # ≥60% of read nodes

    # orientation from GAM (majority of mappings)
    rev_count = sum(is_rev_flags)
    read_orientation = "reverse" if rev_count > (len(is_rev_flags) / 2.0) else "forward"

    # also test reverse node order against peaks
    read_nodes_rev = list(reversed(read_nodes))

    for peak_idx, peak_nodes in peak_paths.items():
        fwd_match = longest_contiguous_match(read_nodes, peak_nodes)
        rev_match = longest_contiguous_match(read_nodes_rev, peak_nodes)

        if fwd_match >= rev_match:
            match_len = fwd_match
            matched_orient = "forward"
        else:
            match_len = rev_match
            matched_orient = "reverse"

        if match_len >= required:
            out_line = [
                name,
                str(peak_idx),       # index in all_peaks.intervalcollection (0-based)
                matched_orient,      # best matching orientation vs peak path
                str(match_len),
                str(read_len),
            ]
            line_str = "\t".join(out_line)
            print(line_str)
            out_fh.write(line_str + "\n")

out_fh.close()
sys.stderr.write(f"Wrote matches to {out_file}\n")
