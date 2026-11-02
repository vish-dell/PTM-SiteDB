import sqlite3
c = sqlite3.connect("instance/ptm.db")
for r in c.execute("select id, site, kinase, pmid, validation_status from ptm_sites where uniprot='P04637' and site='S37'"):
    print(r)
