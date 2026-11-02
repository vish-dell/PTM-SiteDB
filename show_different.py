import csv
with open("verification_report.csv", newline="", encoding="utf-8-sig") as f:
    rows = [{k.lstrip("\ufeff").strip(): (v or "").strip() for k, v in r.items()} for r in csv.DictReader(f)]
hits = [r for r in rows if "DIFFERENT" in r["exp_detail"].upper()]
print("rows with DIFFERENT:", len(hits))
for r in hits:
    print(r["uniprot"], r["site"], "|", r["verdict"], "|", r["flank_match"], "|", r["exp_detail"][:110])
