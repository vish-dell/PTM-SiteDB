import sqlite3
c = sqlite3.connect("instance/ptm.db")
bad = [("P15056","S313"),("P15056","S358"),("P15056","S412"),("P15056","S424"),("P15056","S449"),("P15056","S486")]
n = 0
for u, s in bad:
    n += c.execute("delete from ptm_sites where uniprot=? and site=?", (u, s)).rowcount
c.commit()
print("rows removed:", n)
