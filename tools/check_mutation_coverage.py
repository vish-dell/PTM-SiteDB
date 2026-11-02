import csv, sys
from pathlib import Path
P = Path("research/data/processed")
drv, non = set(), set()
with open(P/"benchmark_panel.csv", newline="", encoding="utf-8") as f:
    for r in csv.DictReader(f):
        (drv if r["intogen_driver"] == "True" else non).add(r["gene"])
mapped = set()
with open(P/"ptm_sites_mapped.tsv", newline="", encoding="utf-8") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        mapped.add(r["gene"])
for path in sys.argv[1:]:
    cnt = {}
    with open(path, newline="", encoding="utf-8", errors="replace") as f:
        lines = (l for l in f if not l.startswith("#"))
        for r in csv.DictReader(lines, delimiter="\t"):
            if r.get("Variant_Classification") == "Missense_Mutation":
                g = r["Hugo_Symbol"]
                cnt[g] = cnt.get(g, 0) + 1
    print(Path(path).name)
    for label, grp in (("drivers", drv), ("not-listed", non)):
        mut = [g for g in grp if g in cnt]
        both = [g for g in mut if g in mapped]
        print(f"  {label}: genes with missense={len(mut)}, of those with mapped PTMs={len(both)}, "
              f"variants={sum(cnt[g] for g in mut)}")
