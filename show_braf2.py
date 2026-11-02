import sqlite3
c = sqlite3.connect("instance/ptm.db")
print("BRAF rows left:", c.execute("select count(*) from ptm_sites where uniprot='P15056'").fetchone()[0])
for r in c.execute("select id, site, validation_status from ptm_sites where uniprot='P15056' order by id"):
    print(r)
