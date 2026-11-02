import sqlite3, json, re, time, urllib.request, shutil, sys
from pathlib import Path

DB = Path("instance/ptm.db")
APPLY = "--apply" in sys.argv
ALIAS = {"AKT1": ["akt", "pkb", "protein kinase b"], "CHEK1": ["chk1", "chek1", "checkpoint kinase 1"],
         "CTNNB1": ["catenin"], "EGFR": ["egfr", "epidermal growth factor receptor", "erbb1"],
         "MAPK1": ["erk2", "mapk1", "erk"], "MAPK3": ["erk1", "mapk3", "erk"], "TP53": ["p53", "tp53"]}
NAMES = {"S": ["ser", "serine", "s"], "T": ["thr", "threonine", "t"], "Y": ["tyr", "tyrosine", "y"]}
BAD = ("uniprot", "database", "mass spectrometry", "ion suppression", "phosphoproteom", "imac")

def score(gene, site, title):
    t = title.lower(); s = 0
    if any(re.search(r"\b" + re.escape(a), t) for a in ALIAS.get(gene, [gene.lower()])): s += 3
    if "phosphorylat" in t: s += 2
    num = site[1:]
    if any(re.search(rf"\b{n}\s*-?\s*{num}\b", t) for n in NAMES[site[0]]): s += 3
    if any(w in t for w in BAD): s -= 3
    return s

def get_titles(ids):
    out = {}
    ids = sorted(ids)
    for i in range(0, len(ids), 150):
        chunk = ids[i:i + 150]
        url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&retmode=json&id="
               + ",".join(chunk))
        d = json.load(urllib.request.urlopen(url, timeout=60))["result"]
        for k in chunk:
            out[k] = d.get(k, {}).get("title") or ""
        time.sleep(0.4)
    return out

con = sqlite3.connect(DB)
rows = con.execute("select id,gene,site,pmid,pmids_dbptm from ptm_sites "
                   "where pmids_dbptm is not null order by id").fetchall()
allids = {p for r in rows for p in r[4].split(";") if p}
print("fetching titles for", len(allids), "PMIDs...")
titles = get_titles(allids)

plan = []
for id_, gene, site, cur, lst in rows:
    cands = [p for p in lst.split(";") if p]
    best = max(cands, key=lambda p: (score(gene, site, titles.get(p, "")), p == str(cur)))
    cs = score(gene, site, titles.get(str(cur), ""))
    bs = score(gene, site, titles.get(best, ""))
    change = bs >= 3 and bs > cs
    print(f"{gene:6} {site:5} now={cur}({cs}) best={best}({bs}) {'CHANGE' if change else 'keep  '} | {titles.get(best,'')[:70]}")
    if change:
        plan.append((id_, best, bs))
print("rows to change:", len(plan), "of", len(rows))

if APPLY:
    bak = Path("instance/ptm.db.bak_before_pmid2")
    if bak.exists():
        sys.exit("backup exists, stopping: " + str(bak))
    con.close(); shutil.copy2(DB, bak)
    con = sqlite3.connect(DB)
    if "pmid_score" not in {r[1] for r in con.execute("pragma table_info(ptm_sites)")}:
        con.execute("alter table ptm_sites add column pmid_score INTEGER")
    for id_, best, bs in plan:
        con.execute("update ptm_sites set pmid=?, pmid_score=? where id=?", (best, bs, id_))
    con.commit()
    print("APPLIED. Backup:", bak)
else:
    print("DRY RUN only. Nothing changed. Add --apply to write.")
