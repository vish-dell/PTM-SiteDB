import sqlite3, shutil, csv, collections
shutil.copy("instance/ptm.db", "instance/ptm_experimental.db")

rows = list(csv.DictReader(open("verification_report.csv", encoding="utf-8-sig")))
keep = {(r["uniprot"].strip(), r["site"].strip())
        for r in rows if r["verdict"].strip().startswith("VERIFIED_EXPERIMENTAL_STRUCTURE")}
print("sites passing in CSV:", len(keep))

c = sqlite3.connect("instance/ptm_experimental.db")
before = c.execute("select count(*) from ptm_sites").fetchone()[0]
for id_, up, site in c.execute("select id, uniprot, site from ptm_sites").fetchall():
    if (up, site) not in keep:
        c.execute("delete from ptm_sites where id=?", (id_,))
c.execute("update ptm_sites set validation_status='VERIFIED_EXPERIMENTAL_STRUCTURE'")
c.commit()

print("records:", before, "->", c.execute("select count(*) from ptm_sites").fetchone()[0])
print("unique sites:", c.execute("select count(*) from (select distinct uniprot, site from ptm_sites)").fetchone()[0])
print("proteins:", c.execute("select count(distinct uniprot) from ptm_sites").fetchone()[0])
print("with real kinase:", c.execute("select count(*) from ptm_sites where kinase is not null and length(trim(kinase))>1").fetchone()[0])
for g, n in c.execute("select gene, count(distinct site) from ptm_sites group by gene order by 2 desc"):
    print("  ", g, n)
