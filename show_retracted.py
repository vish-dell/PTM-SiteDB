import sqlite3
c = sqlite3.connect("instance/ptm.db")
for r in c.execute("select id, site, pmid, validation_status from ptm_sites where uniprot='P35222' and site in ('Y86','Y654')"): print(r)
