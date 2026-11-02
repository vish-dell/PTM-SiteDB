import sqlite3
c = sqlite3.connect("instance/ptm.db")
print("records:", c.execute("select count(*) from ptm_sites").fetchone()[0])
print("unique sites:", c.execute("select count(*) from (select distinct uniprot, site from ptm_sites)").fetchone()[0])
print("proteins (uniprot):", c.execute("select count(distinct uniprot) from ptm_sites").fetchone()[0])
print("genes:", c.execute("select count(distinct gene) from ptm_sites").fetchone()[0])
print("kinases:", c.execute("select count(distinct kinase) from ptm_sites where kinase is not null and length(trim(kinase))>1").fetchone()[0])
import re
p = set()
for (x,) in c.execute("select pmid from ptm_sites"): p |= set(re.findall(r"\d+", x or ""))
print("distinct PMIDs:", len(p))
