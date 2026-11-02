import sqlite3
c = sqlite3.connect("instance/ptm.db")
n = c.execute("update ptm_sites set kinase = char(8212) where uniprot = 'P04637' and site = 'S37' and kinase = 'CDK5'").rowcount
c.commit()
print("rows updated:", n)
