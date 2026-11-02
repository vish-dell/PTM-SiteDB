import sqlite3
for f in ("instance/ptm_experimental.db", "instance/ptm.db"):
    c = sqlite3.connect(f)
    print(f, "records:", c.execute("select count(*) from ptm_sites").fetchone()[0],
          "| unique sites:", c.execute("select count(*) from (select distinct uniprot, site from ptm_sites)").fetchone()[0],
          "| proteins:", c.execute("select count(distinct uniprot) from ptm_sites").fetchone()[0])
