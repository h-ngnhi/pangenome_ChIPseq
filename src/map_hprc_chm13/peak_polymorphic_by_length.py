#!/usr/bin/env python3
import argparse, gzip, json, re, sys
from typing import Dict, List, Set, Tuple, Optional

AT_RE = re.compile(r"(?:^|;)AT=([^;]+)")
NH_RE = re.compile(r"(?:^|;)n_hits=([^;]+)")

def open_text(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")

def parse_at(info: str) -> Optional[Tuple[List[int], List[int]]]:
    m = AT_RE.search(info)
    if not m:
        return None
    at = m.group(1)
    parts = at.split(",")
    if len(parts) < 2:
        return None

    def to_nodes(s: str) -> List[int]:
        s = s.strip()
        if s.startswith(">"):
            s = s[1:]
        if not s:
            return []
        return [int(x) for x in s.split(">") if x]

    p0 = to_nodes(parts[0])
    p1 = to_nodes(parts[1])
    if not p0 or not p1:
        return None
    return p0, p1

def extract_core(ref: str, alt: str) -> Tuple[str, str]:
    ref = ref.upper()
    alt = alt.upper()

    # common prefix
    i = 0
    n = min(len(ref), len(alt))
    while i < n and ref[i] == alt[i]:
        i += 1
    ref2, alt2 = ref[i:], alt[i:]

    # common suffix (avoid overlap)
    j = 0
    n2 = min(len(ref2), len(alt2))
    while j < n2 and ref2 and alt2 and ref2[-(j+1)] == alt2[-(j+1)]:
        j += 1
    if j > 0:
        return ref2[:-j], alt2[:-j]
    return ref2, alt2

def parse_fasta_peaks(fa_path: str) -> List[Tuple[str, Set[int], int]]:
    """
    Returns list of (peak_id, node_set, peak_len_bp)
    Expects headers like:
    >peak0 {"start":..., "end":..., "region_paths":[...], ...}
    """
    peaks = []
    with open_text(fa_path) as f:
        header = None
        seq_chunks = []
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                # flush previous
                if header is not None:
                    peak_id, nodes = header
                    seq = "".join(seq_chunks).strip()
                    peaks.append((peak_id, nodes, len(seq)))
                # parse new
                seq_chunks = []
                # peak id is first token after >
                first = line[1:].split()[0]
                # json is after peak id (first space)
                js = line[1+len(first):].strip()
                nodes = set()
                if js:
                    try:
                        obj = json.loads(js)
                        nodes = set(int(x) for x in obj.get("region_paths", []))
                    except Exception:
                        nodes = set()
                header = (first, nodes)
            else:
                seq_chunks.append(line.strip())
        # last
        if header is not None:
            peak_id, nodes = header
            seq = "".join(seq_chunks).strip()
            peaks.append((peak_id, nodes, len(seq)))
    return peaks

def load_altonly_bp_per_node(vcf_path: str, require_nhits1: bool) -> Dict[int, float]:
    """
    Build a mapping: node_id -> bp_weight
    bp_weight is (ALT_core_len / #ALT_only_nodes) for that variant.
    If a node appears in multiple variants, we keep the MAX weight (conservative).
    """
    node_w: Dict[int, float] = {}
    with open_text(vcf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue
            _chrom, _pos, _id, ref, alt, *_ = cols[:6]
            info = cols[7]

            if require_nhits1:
                m = NH_RE.search(info)
                if not m:
                    continue
                try:
                    if int(m.group(1)) != 1:
                        continue
                except ValueError:
                    continue

            at = parse_at(info)
            if at is None:
                continue
            p0, p1 = at

            # We don't know which AT path corresponds to ALT string, but for "ALT-only nodes"
            # we only need the symmetric difference between paths.
            s0, s1 = set(p0), set(p1)
            diff = (s0 - s1) | (s1 - s0)
            if not diff:
                continue

            # core length from VCF REF/ALT strings (first ALT only)
            alt1 = alt.split(",")[0]
            _rc, ac = extract_core(ref, alt1)
            core_len = len(ac)
            if core_len <= 0:
                continue

            w = core_len / float(len(diff))
            for n in diff:
                prev = node_w.get(n)
                if prev is None or w > prev:
                    node_w[n] = w
    return node_w

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--peaks-fasta", required=True)
    ap.add_argument("--snap-vcf", required=True, help="snapped VCF (optionally already filtered)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--require-nhits1", action="store_true",
                    help="require INFO/n_hits=1 (use if you didn't pre-filter)")
    ap.add_argument("--min-frac", type=float, default=0.0,
                    help="polymorphic if (variant_bp / peak_len) >= this")
    ap.add_argument("--min-bp", type=float, default=1.0,
                    help="and variant_bp >= this (default 1bp)")
    args = ap.parse_args()

    peaks = parse_fasta_peaks(args.peaks_fasta)
    node_w = load_altonly_bp_per_node(args.snap_vcf, args.require_nhits1)

    with open(args.out, "w") as out:
        out.write("\t".join([
            "peak_id","peak_len_bp","variant_bp_est","frac_variant","is_polymorphic"
        ]) + "\n")

        for pid, pnodes, plen in peaks:
            if plen <= 0:
                out.write(f"{pid}\t0\t0\t0\t0\n")
                continue

            var_bp = 0.0
            # sum weights for overlapped variant-only nodes
            for n in pnodes:
                w = node_w.get(n)
                if w is not None:
                    var_bp += w

            frac = var_bp / float(plen)
            is_poly = 1 if (var_bp >= args.min_bp and frac >= args.min_frac) else 0
            out.write(f"{pid}\t{plen}\t{var_bp:.3f}\t{frac:.6f}\t{is_poly}\n")

if __name__ == "__main__":
    main()
