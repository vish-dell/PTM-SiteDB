import sqlite3, re
c = sqlite3.connect("instance/ptm.db")
rows = c.execute("select id, gene, site, pmid, pmids_dbptm from ptm_sites where validation_status is null").fetchall()
nums, per = set(), {}
for r in rows:
    txt = " ".join(str(x) for x in (r[3], r[4]) if x)
    for m in re.findall(r"\d{5,9}", txt):
        n = int(m); nums.add(n); per.setdefault(n, []).append((r[1], r[2]))
s = sorted(nums)
runs, cur = [], ([s[0]] if s else [])
for a in s[1:]:
    if a == cur[-1] + 1:
        cur.append(a)
    else:
        if len(cur) >= 3: runs.append(cur)
        cur = [a]
if len(cur) >= 3: runs.append(cur)
print("distinct PMIDs among unchecked rows:", len(s))
print("runs of 3+ consecutive PMIDs:", len(runs))
for r in runs[:10]:
    print(r[0], "...", r[-1], "(", len(r), "PMIDs ) e.g.", sorted({x for n in r for x in per[n]})[:4])
print("rows touched by runs:", len({x for r in runs for n in r for x in per[n]}))
