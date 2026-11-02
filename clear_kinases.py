import sqlite3
c = sqlite3.connect("instance/ptm.db")
pairs = [("P04637","T18"),("P04637","S215"),("P04637","S183"),("P04637","S371"),("P35222","Y489"),("P31749","T450")]
for up, site in pairs:
    c.execute("update ptm_sites set kinase='' where uniprot=? and site=?", (up, site))
c.commit()
print("records with kinase:", c.execute("select count(*) from ptm_sites where kinase is not null and length(trim(kinase))>1").fetchone()[0])
print("distinct kinases:", c.execute("select count(distinct kinase) from ptm_sites where kinase is not null and length(trim(kinase))>1").fetchone()[0])
