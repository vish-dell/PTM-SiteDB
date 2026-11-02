import csv, collections
rows = list(csv.DictReader(open("pmid_check.csv", encoding="utf-8-sig")))
cnt = collections.Counter()
for r in rows:
    for f in r["flags"].split(";"):
        cnt[f.split(":")[0]] += 1
print(dict(cnt))
for key in ("PMID_NOT_FOUND", "GENE_NOT_IN_ABSTRACT", "NOT_TAGGED_HUMAN", "VERY_OLD"):
    hits = [r for r in rows if key in r["flags"]]
    print("\n==", key, len(hits))
    for r in hits[:40]:
        print(r["uniprot"], r["site"], r["pmid"], r["kinase"], "|", r["flags"], "|", r["title"][:60])
