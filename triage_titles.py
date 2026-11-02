import csv, re
G = {"P04637":["p53","tp53"],"P31749":["akt","pkb"],"P35222":["catenin"],"P00533":["egfr","epidermal growth factor"],
 "P28482":["erk","mapk","map kinase","raf"],"P27361":["erk","mapk","map kinase","raf"],"P15056":["raf"],
 "P38398":["brca1"],"O14757":["chk1","chek1","checkpoint"],"P42336":["pi3k","p110","phosphoinositide","pik3ca"]}
GEN = ["phospho","proteom","kinase","signal","mass spectrom","phosphatase","cancer","tumor","tumour","cell cycle","atlas","enrich"]
rows = list(csv.DictReader(open("pmid_check.csv", encoding="utf-8-sig")))
out = []
for r in rows:
    if r["flags"] == "OK": continue
    t = r["title"].lower()
    if not any(k in t for k in GEN + G.get(r["uniprot"], [])):
        out.append(r)
print("rows to review:", len(out))
for r in out: print(r["uniprot"], r["site"], r["pmid"], "|", r["kinase"], "|", r["title"][:80])
with open("review_titles.csv","w",newline="",encoding="utf-8-sig") as f:
    w = csv.DictWriter(f, fieldnames=rows[0].keys()); w.writeheader(); w.writerows(out)
