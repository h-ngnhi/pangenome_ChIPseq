import json

ALIGNMENTS_JSON = "results/146/vg_giraffe_vcfbub/treatment_alignments.filtered.json"
CHR1_NODE_IDS = "Pangenomes/vg_giraffe/vcfbub/graphs_chunks/chr1_ids.txt"
OUT = "results/146/vg_giraffe_vcfbub/chr1_mapped_read_node_ids.txt"

# Load node IDs into a set
with open(CHR1_NODE_IDS) as f:
    node_ids = set(line.strip() for line in f if line.strip())

found = set()
with open(ALIGNMENTS_JSON) as f:
    data = json.load(f)
    for aln in data["alignment"]:
        for mapping in aln.get("path", {}).get("mapping", []):
            nid = str(mapping.get("position", {}).get("node_id"))
            if nid in node_ids:
                found.add(nid)
with open(OUT, "w") as fout:
    for nid in sorted(found, key=int):
        fout.write(nid + "\n")
