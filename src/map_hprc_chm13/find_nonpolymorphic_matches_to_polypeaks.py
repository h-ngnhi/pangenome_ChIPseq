#!/usr/bin/env python3
import sys
import json
import math
from collections import defaultdict

"""
Usage:

  vg view -aj results/146/vg_giraffe_chm13/treatment_alignments.gam \
    | python3 pangenome_ChIPseq/src/find_nonpolymorphic_matches_to_polypeaks.py \
        results/146/vg_giraffe_chm13/callpeaks/all_peaks.intervalcollection \
        results/146/vg_giraffe_chm13/polymorphic_reads_chm13_peak_matches.tsv
"""

if len(sys.argv) != 3:
    sys.stderr.write(
        f"Usage:\n  vg view -aj CHM13.gam | {sys.argv[0]} chm13_all_peaks.intervalcollection polymorphic_reads_chm13_peak_matches.tsv\n"
    )
    sys.exit(1)

peaks_file = sys.argv[1]
poly_match_file = sys.argv[2]
out_path = "results/146/vg_giraffe_chm13/nonpolymorphic_reads_to_polypeaks.tsv"

# ---------- 1) Read polymorphic mapping table ----------
# Expect header:
# read_name  hprc_peak_index  chm13_peak_index  orientation  match_nodes  read_nodes
poly_reads = set()
hprc_to_chm13 = defaultdict(set)
chm13_poly_peaks = set()

with open(poly_match_file) as f:
    header = f.readline().rstrip("\n").split("\t")
    if len(header) < 3:
        sys.stderr.write(f"{poly_match_file}: header has <3 columns, check format\n")
        sys.exit(1)

    for line in f:
        line = line.strip()
        if not line:
            continue
        cols = line.split("\t")
        if len(cols) < 3:
            continue
        rname = cols[0]
        try:
            hprc_idx = int(cols[1])
            chm13_idx = int(cols[2])
        except ValueError:
            continue

        poly_reads.add(rname)
        hprc_to_chm13[hprc_idx].add(chm13_idx)
        chm13_poly_peaks.add(chm13_idx)

sys.stderr.write(f"Loaded {len(poly_reads)} polymorphic reads from {poly_match_file}\n")
sys.stderr.write(f"CHM13 peaks hit by polymorphic reads (count): {len(chm13_poly_peaks)}\n")

# --- Check mapping multiplicity per HPRC peak ---
sys.stderr.write("HPRC_peak_index -> #CHM13_peaks  {list}\n")
for hprc_idx in sorted(hprc_to_chm13.keys()):
    peaks = sorted(hprc_to_chm13[hprc_idx])
    sys.stderr.write(f"  {hprc_idx}: {len(peaks)}  {peaks}\n")

# ---------- 2) Load CHM13 peaks, but focus on those used ----------
peak_paths = {}                    # peak_idx -> [node_ids]
node_to_peaks = defaultdict(set)   # node_id -> set(peak_idx)

with open(peaks_file) as f:
    for idx, line in enumerate(f):
        line = line.strip()
        if not line:
            continue
        # only store peaks that polymorphic reads mapped to
        if idx not in chm13_poly_peaks:
            continue
        rec = json.loads(line)
        nodes = rec.get("region_paths", [])
        if not nodes:
            continue
        node_ids = [int(n) for n in nodes]
        peak_paths[idx] = node_ids
        for nid in node_ids:
            node_to_peaks[nid].add(idx)

sys.stderr.write(f"Loaded {len(peak_paths)} CHM13 peaks of interest from {peaks_file}\n")

# ---------- helper: longest contiguous match ----------
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

    for start in range(n_p):
        if peak_nodes[start] != read_nodes[0]:
            continue
        j = 0
        while start + j < n_p and j < n_r and peak_nodes[start + j] == read_nodes[j]:
            j += 1
        if j > max_len:
            max_len = j
    return max_len

# ---------- 3) Open output, write header, stream GAM; find NON-polymorphic reads ----------
out_header = [
    "read_name",
    "chm13_peak_index",
    "orientation",
    "match_nodes",
    "read_nodes"
]
header_line = "\t".join(out_header)

# open output file (line-buffered)
out_fh = open(out_path, "w", buffering=1)

# write header to stdout + file
sys.stdout.write(header_line + "\n")
sys.stdout.flush()
out_fh.write(header_line + "\n")

# To avoid duplicate rows if the same read/peak combo is found multiple times
emitted = set()  # (read_name, chm13_peak_index)

for line in sys.stdin:
    line = line.strip()
    if not line:
        continue
    aln = json.loads(line)

    name = aln.get("name", "")
    # Skip polymorphic reads; we only want non-polymorphic
    if name in poly_reads:
        continue

    path = aln.get("path", {})
    mappings = path.get("mapping", [])
    if not mappings:
        continue

    # sort mappings by rank
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
    required = math.ceil(0.6 * read_len)  # >= 60% of read nodes

    rev_count = sum(is_rev_flags)
    read_orientation = "reverse" if rev_count > (len(is_rev_flags) / 2.0) else "forward"

    # Candidate peaks: any of the CHM13 peaks of interest that share a node with this read
    candidate_peaks = set()
    for nid in read_nodes:
        if nid in node_to_peaks:
            candidate_peaks.update(node_to_peaks[nid])

    if not candidate_peaks:
        continue

    read_nodes_rev = list(reversed(read_nodes))

    for peak_idx in candidate_peaks:
        # only peaks that were hit by polymorphic reads
        peak_nodes = peak_paths.get(peak_idx)
        if not peak_nodes:
            continue

        fwd_match = longest_contiguous_match(read_nodes, peak_nodes)
        rev_match = longest_contiguous_match(read_nodes_rev, peak_nodes)

        if fwd_match >= rev_match:
            match_len = fwd_match
            matched_orient = "forward"
        else:
            match_len = rev_match
            matched_orient = "reverse"

        if match_len >= required:
            key = (name, peak_idx)
            if key in emitted:
                continue
            emitted.add(key)

            out_fields = [
                name,
                str(peak_idx),
                matched_orient,
                str(match_len),
                str(read_len),
            ]
            line_str = "\t".join(out_fields)

            # write to stdout
            sys.stdout.write(line_str + "\n")
            sys.stdout.flush()

            # write to file
            out_fh.write(line_str + "\n")
            out_fh.flush()

out_fh.close()
sys.stderr.write(f"Wrote non-polymorphic matches to {out_path}\n")
