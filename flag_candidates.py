import csv, re, collections
rows = []
with open("verification_report.csv", newline="", encoding="utf-8-sig") as f:
    for raw in csv.DictReader(f):
        rows.append({k.lstrip("\ufeff").strip(): (v or "").strip() for k, v in raw.items()})
print("rows in CSV:", len(rows))
print("verdicts:", dict(collections.Counter(r["verdict"].split(":")[0] for r in rows)))

def ratio(s):
    m = re.match(r"(\d+)\s*/\s*(\d+)", s)
    return int(m.group(1)) / int(m.group(2)) if m and int(m.group(2)) else None

low = [(r["uniprot"], r["site"], r["flank_match"]) for r in rows if ratio(r["flank_match"]) is not None and ratio(r["flank_match"]) < 0.8]
print("flank match below 80%:", len(low), low[:15])

tok = collections.Counter(t for r in rows for t in re.findall(r"[A-Z][A-Z_]{5,}", r["exp_detail"]))
print("experimental-check labels:", dict(tok))
odd = [(r["uniprot"], r["site"], r["exp_detail"][:70]) for r in rows
       if r["experimental"].upper() == "VERIFIED" and "PHOSPHO_PRESENT" not in r["exp_detail"].upper()]
print("verified, but no 'PHOSPHO_PRESENT' label:", len(odd), odd[:8])

lowaf = [(r["uniprot"], r["site"], r["af_plddt"]) for r in rows
         if r["verdict"].startswith("SEQUENCE_OK") and r["af_plddt"].replace(".", "").isdigit() and float(r["af_plddt"]) < 50]
print("AlphaFold-only with pLDDT below 50:", len(lowaf), lowaf[:10])
