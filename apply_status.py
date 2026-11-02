import sqlite3, csv, sys, collections
apply = "--apply" in sys.argv
verdict = {}
with open("verification_report.csv", newline="", encoding="utf-8-sig") as f:
    for raw in csv.DictReader(f):
        r = {k.lstrip("\ufeff").strip(): v for k, v in raw.items()}
        verdict[(r["uniprot"].strip(), r["site"].strip())] = r["verdict"].strip()
c = sqlite3.connect("instance/ptm.db")
rows = c.execute("select id, uniprot, site, validation_status from ptm_sites").fetchall()
todo, missing = [], []
for i, u, s, st in rows:
    v = verdict.get((u, s))
    if v is None:
        missing.append((u, s))
        continue
    if st is None or st == "":
        todo.append((v, i))
count = collections.Counter(v for v, _ in todo)
print("CSV entries read:", len(verdict))
print("rows to update:", len(todo), dict(count))
print("rows with no entry in the CSV:", len(missing), missing[:10])
if apply:
    c.executemany("update ptm_sites set validation_status = ? where id = ?", todo)
    c.commit()
    print("updated", len(todo), "rows")
else:
    print("dry run only - add --apply to write")

sites = c.execute("select distinct uniprot, site from ptm_sites").fetchall()
cnt = collections.Counter(verdict.get(k, "NOT_IN_CSV") for k in sites)
print("rows in database:", len(rows), "| unique sites:", len(sites))
print("unique sites by verification result:")
for k, v in cnt.most_common():
    print("  ", k, ":", v)
