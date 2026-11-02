import sqlite3
c = sqlite3.connect("instance/ptm.db")
print([r[0] for r in c.execute("select name from sqlite_master where type='table'")])
