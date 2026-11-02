import sqlite3
con = sqlite3.connect("instance/ptm.db")
q = """select ptm_type, count(*) as records, count(distinct position) as residues
       from ptm_sites where uniprot = 'P04637' group by ptm_type"""
for row in con.execute(q):
    print(row)