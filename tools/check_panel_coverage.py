import csv, re
from pathlib import Path
P = Path("research/data/processed")
panel = {}
with open(P/"benchmark_panel.csv", newline="", encoding="utf-8") as f:
    for r in csv.DictReader(f):
        panel[r["gene"]] = (r["intogen_driver"] == "True")
drv = {g for g, v in panel.items() if v}
non = {g for g, v in panel.items() if not v}
print("panel drivers:", len(drv), "| not-listed:", len(non))

heads = [l[1:].strip() for l in open(P/"canonical_proteins.fasta", encoding="utf-8", errors="replace") if l.startswith(">")]
print("sequences in canonical_proteins.fasta:", len(heads))
print("example headers:", heads[:3])
def in_fasta(g):
    pat = r"(?<![A-Za-z0-9])" + re.escape(g) + r"(?![A-Za-z0-9])"
    return any(re.search(pat, h) for h in heads)
print("drivers in fasta:", sum(in_fasta(g) for g in drv), "| not-listed in fasta:", sum(in_fasta(g) for g in non))

def count_genes(name):
    with open(P/name, newline="", encoding="utf-8", errors="replace") as f:
        rd = csv.reader(f, delimiter="\t")
        hdr = next(rd)
        if "gene" not in hdr:
            print(name, "has no gene column:", hdr)
            return None
        i = hdr.index("gene")
        s = {}
        for row in rd:
            if len(row) > i:
                s[row[i]] = s.get(row[i], 0) + 1
    return s

for name in ("ptm_sites.tsv", "alphafold_ca_coordinates.tsv", "ptm_sites_mapped.tsv"):
    s = count_genes(name)
    if s is not None:
        print(name, "| drivers with rows:", sum(g in s for g in drv), "| not-listed with rows:", sum(g in s for g in non))
