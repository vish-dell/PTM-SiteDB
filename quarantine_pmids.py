import sqlite3, shutil
shutil.copy("instance/ptm.db", "instance/ptm_backup_before_pmid_cleanup.db")
c = sqlite3.connect("instance/ptm.db")
shared = [r[1] for r in c.execute("pragma table_info(ptm_sites_quarantine)") if r[1] != "quarantine_reason"]
cols = ",".join(shared)
bad = {336:"PMID 101545721 does not exist in PubMed",
       349:"PMID 11725543 is about perlecan, unrelated to AKT1",
       354:"PMID 119577 is a 1979 rat cardiac study, unrelated to AKT1",
       361:"PMID 16173467 is about a rice protein, unrelated to AKT1",
       362:"PMID 16173479 is an acupuncture study, unrelated to AKT1",
       363:"PMID 16173491 is a drug pharmacokinetics study, unrelated to AKT1"}
q_before = c.execute("select count(*) from ptm_sites_quarantine").fetchone()[0]
r_before = c.execute("select count(*) from ptm_sites").fetchone()[0]
for id_, why in bad.items():
    if not c.execute("select 1 from ptm_sites where id=? and uniprot='P31749'", (id_,)).fetchone():
        print("skip, not found:", id_); continue
    reason = "PMID check: " + why
    try:
        c.execute(f"insert into ptm_sites_quarantine ({cols}, quarantine_reason) select {cols}, ? from ptm_sites where id=?", (reason, id_))
    except sqlite3.IntegrityError:
        other = ",".join(x for x in shared if x != "id")
        c.execute(f"insert into ptm_sites_quarantine ({other}, quarantine_reason) select {other}, ? from ptm_sites where id=?",
                  (reason + f" (original id {id_})", id_))
    c.execute("delete from ptm_sites where id=? and uniprot='P31749'", (id_,))
# clear kinase labels the cited papers do not support (sites are kept)
for site in ("S315", "S6", "S33"):
    c.execute("update ptm_sites set kinase='' where uniprot='P04637' and site=? and kinase in ('HIPK2','CDK5')", (site,))
c.commit()
print("records:", r_before, "->", c.execute("select count(*) from ptm_sites").fetchone()[0])
print("quarantine rows:", q_before, "->", c.execute("select count(*) from ptm_sites_quarantine").fetchone()[0])
print("unique sites:", c.execute("select count(*) from (select distinct uniprot, site from ptm_sites)").fetchone()[0])
print("with real kinase:", c.execute("select count(*) from ptm_sites where kinase is not null and length(trim(kinase))>1").fetchone()[0])
