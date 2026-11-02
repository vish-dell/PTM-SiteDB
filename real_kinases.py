import sqlite3
c = sqlite3.connect("instance/ptm.db")
rows = c.execute("select uniprot, site, kinase, pmid from ptm_sites where kinase is not null and trim(kinase) != '' and length(trim(kinase)) > 1 order by kinase").fetchall()
print("records with real kinase:", len(rows))
for r in rows: print(r)
