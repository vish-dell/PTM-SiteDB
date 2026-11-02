import sqlite3
c = sqlite3.connect("instance/ptm.db")
for r in c.execute("select gene, count(*), sum(disease_context is not null and trim(disease_context)!='') from ptm_sites group by gene"): print(r)
for r in c.execute("select gene, site, disease_context, confidence, pmid from ptm_sites where disease_context like '%Diabetes%'"): print(r)
