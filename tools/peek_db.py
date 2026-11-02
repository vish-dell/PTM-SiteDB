import sqlite3
c = sqlite3.connect("file:instance/ptm.db?mode=ro", uri=True)
cur = c.execute("select * from ptm_sites")
print([d[0] for d in cur.description])
for r in cur:
    print(" | ".join(str(x) for x in r))
