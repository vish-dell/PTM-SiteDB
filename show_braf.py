import sqlite3
c = sqlite3.connect("instance/ptm.db")
for r in c.execute("select id, site, pmid, validation_status from ptm_sites where uniprot='P15056' and site in ('S449','S486')"):
    print(r)
