import sqlite3
c = sqlite3.connect("instance/ptm.db")
shared = [r[1] for r in c.execute("pragma table_info(ptm_sites_quarantine)") if r[1] != "quarantine_reason"]
live = [r[1] for r in c.execute("pragma table_info(ptm_sites)")]
shared = [x for x in shared if x in live]
cols = ",".join(shared)
reason = "Only citation PMID 11279024 is retracted on PubMed"
before = c.execute("select count(*) from ptm_sites").fetchone()[0]
try:
    c.execute(f"insert into ptm_sites_quarantine ({cols}, quarantine_reason) select {cols}, ? from ptm_sites where id=632 and uniprot='P35222' and site='Y654'", (reason,))
except sqlite3.IntegrityError:
    other = ",".join(x for x in shared if x != "id")
    c.execute(f"insert into ptm_sites_quarantine ({other}, quarantine_reason) select {other}, ? from ptm_sites where id=632 and uniprot='P35222' and site='Y654'", (reason + " (original id 632)",))
c.execute("delete from ptm_sites where id=632 and uniprot='P35222' and site='Y654'")
c.commit()
print("records:", before, "->", c.execute("select count(*) from ptm_sites").fetchone()[0])
print("unique sites:", c.execute("select count(*) from (select distinct uniprot, site from ptm_sites)").fetchone()[0])
print("proteins:", c.execute("select count(distinct uniprot) from ptm_sites").fetchone()[0])
print("quarantine rows:", c.execute("select count(*) from ptm_sites_quarantine").fetchone()[0])
