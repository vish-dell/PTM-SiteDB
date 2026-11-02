import sqlite3
c = sqlite3.connect("instance/ptm.db")
for up, site in [("P35222","T393"),("P31749","Y176")]:
    n = c.execute("update ptm_sites set kinase='', kinase_pmid=NULL, kinase_evidence=NULL where uniprot=? and site=?", (up, site)).rowcount
    print(up, site, "rows reset:", n)
c.commit()
print("records with kinase:", c.execute("select count(*) from ptm_sites where kinase is not null and length(trim(kinase))>1").fetchone()[0])
