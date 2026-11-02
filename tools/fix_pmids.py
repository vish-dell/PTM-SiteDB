import sqlite3, shutil, sys
from pathlib import Path

DB = Path("instance/ptm.db")
SRC = Path("research/data/processed/dbptm_phosphorylation_source.txt")
APPLY = "--apply" in sys.argv

con = sqlite3.connect(DB)
rows = con.execute("select id,gene,uniprot,site,position,pmid from ptm_sites "
                   "where ptm_type='Phosphorylation' order by id").fetchall()
want = {(r[2], r[4]): r for r in rows}
accs = {r[2] for r in rows}

found = {}
with open(SRC, encoding="utf-8", errors="replace") as f:
    for line in f:
        p = line.rstrip("\n").split("\t")
        if len(p) < 6 or p[1] not in accs or p[3] != "Phosphorylation":
            continue
        try:
            pos = int(p[2])
        except ValueError:
            continue
        if (p[1], pos) not in want:
            continue
        lst = found.setdefault((p[1], pos), [])
        for x in p[4].split(";"):
            x = x.strip()
            if x.isdigit() and x not in lst:
                lst.append(x)

plan = []
for id_, gene, acc, site, pos, old in rows:
    ids = found.get((acc, pos), [])
    if not ids:
        print(f"{gene:6} {site:6} NOT FOUND in dbPTM file - left unchanged")
        continue
    keep_old = str(old) in ids
    new = str(old) if keep_old else ids[0]
    plan.append((id_, str(old), new, ids))
    print(f"{gene:6} {site:6} old={old} in_dbptm={keep_old} -> new={new} (dbPTM lists {len(ids)})")

print("rows to update:", len(plan), "of", len(rows))

if APPLY:
    bak = Path("instance/ptm.db.bak_before_pmid")
    if bak.exists():
        sys.exit("backup already exists, stopping: " + str(bak))
    con.close()
    shutil.copy2(DB, bak)
    con = sqlite3.connect(DB)
    cols = {r[1] for r in con.execute("pragma table_info(ptm_sites)")}
    for name, typ in (("pmid_old", "TEXT"), ("pmids_dbptm", "TEXT"), ("pmid_count", "INTEGER")):
        if name not in cols:
            con.execute(f"alter table ptm_sites add column {name} {typ}")
    for id_, old, new, ids in plan:
        con.execute("update ptm_sites set pmid_old=coalesce(pmid_old,?), pmid=?, "
                    "pmids_dbptm=?, pmid_count=? where id=?",
                    (old, new, ";".join(ids), len(ids), id_))
    con.commit()
    print("APPLIED. Backup:", bak)
else:
    print("DRY RUN only. Nothing changed. Add --apply to write.")
