import sqlite3
c = sqlite3.connect("instance/ptm.db")
print("total sites:", c.execute("select count(*) from ptm_sites").fetchone()[0])
print("by status:", c.execute("select validation_status, count(*) from ptm_sites group by validation_status").fetchall())
print("by gene:", c.execute("select gene, count(*) from ptm_sites group by gene order by 2 desc").fetchall())
print("with kinase:", c.execute("select count(*) from ptm_sites where kinase is not null and kinase <> '' and kinase <> char(8212)").fetchone()[0])
