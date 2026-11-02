import sqlite3, json, urllib.request
from pathlib import Path

CACHE = Path("validation/uniprot_cache")
def seq(acc):
    lines = (CACHE / f"{acc}.fasta").read_text().splitlines()
    return "".join(l.strip() for l in lines if not l.startswith(">"))

c = sqlite3.connect("file:instance/ptm.db?mode=ro", uri=True)
rows = c.execute("select id,gene,uniprot,site,position,flanking,pmid from ptm_sites order by id").fetchall()

print("== flank check ==")
bad = 0
for id_, gene, acc, site, pos, fl, pmid in rows:
    s = seq(acc)
    exp = s[max(0, pos - 11):pos + 10]
    ok = (fl == exp) and s[pos - 1] == site[0]
    bad += (not ok)
    print(f"{id_:>3} {gene:6} {site:6} ok={ok}")
    if not ok:
        print("    stored ", fl)
        print("    uniprot", exp)
print("flank mismatches:", bad, "of", len(rows))

print()
print("== PMID titles (NCBI) ==")
pm = sorted({str(r[6]) for r in rows})
url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&retmode=json&id="
       + ",".join(pm))
d = json.load(urllib.request.urlopen(url, timeout=30))["result"]
for id_, gene, acc, site, pos, fl, pmid in rows:
    e = d.get(str(pmid), {})
    title = e.get("title") or "NO RECORD FOUND"
    print(f"{gene} {site} | PMID {pmid} | {title[:100]}")
