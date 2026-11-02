import sqlite3
c = sqlite3.connect("instance/ptm.db")
print("quarantine columns:", [r[1] for r in c.execute("pragma table_info(ptm_sites_quarantine)")])
print("ptm_sites columns:", [r[1] for r in c.execute("pragma table_info(ptm_sites)")])
print("rows:", c.execute("select count(*) from ptm_sites_quarantine").fetchone()[0])
for r in c.execute("select * from ptm_sites_quarantine limit 3"): print(r)
