import sqlite3
c = sqlite3.connect("instance/ptm.db")
rows = c.execute("select kinase, count(*) from ptm_sites where kinase is not null and kinase not in ('', '?') group by kinase order by 2 desc").fetchall()
print("distinct kinases:", len(rows))
for r in rows: print(r)
print("proteins:", c.execute("select count(distinct uniprot) from ptm_sites").fetchone()[0])
