import sqlite3
c = sqlite3.connect("instance/ptm.db")
for s in ("S2","T312","Y437","S475","S477","T479"):
    for r in c.execute("select id, site, pmid, kinase, validation_status from ptm_sites where uniprot='P31749' and site=?", (s,)):
        print(r)
