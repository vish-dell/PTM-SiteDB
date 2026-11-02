import sqlite3
OPTIONAL = False   # set True to also add PTK6 (AKT1 Y315/Y326) and ERK (EGFR T693)
c = sqlite3.connect("instance/ptm.db")
cols = [r[1] for r in c.execute("pragma table_info(ptm_sites)")]
for col in ("kinase_pmid", "kinase_evidence"):
    if col not in cols: c.execute(f"alter table ptm_sites add column {col} TEXT")
new = [("P04637","S46","HIPK2","11780126","direct phosphorylation; antisense knockdown, human cells", None),
       ("P35222","S45","CK1","12051714","in vitro assay; CK1 yes, CK2 and GSK3B no; human cells", None),
       ("P35222","T393","CK2","12700239","major CK2 site; T393 mutant reduces phosphorylation", None),
       ("P31749","Y176","ACK1","20333297","direct phosphorylation; ACK1 siRNA lowers Y176 in human cancer cells", "YYAM")]
opt = [("P31749","Y315","PTK6","20606012","direct phosphorylation; Src-independent", None),
       ("P31749","Y326","PTK6","20606012","direct phosphorylation; Src-independent", None),
       ("P00533","T693","ERK","32782566","MEK inhibitor and T669A mutant; ERK inferred (T693 = T669 mature numbering)", "PSGE")]
for up, site, kin, pm, ev, guard in new + (opt if OPTIONAL else []):
    row = c.execute("select flanking from ptm_sites where uniprot=? and site=? limit 1", (up, site)).fetchone()
    if row is None: print(up, site, "not in database, skipped"); continue
    if guard and guard not in (row[0] or "").upper():
        print(up, site, "flank", row[0], "does not contain", guard, "- skipped, check numbering"); continue
    n = c.execute("update ptm_sites set kinase=?, kinase_pmid=?, kinase_evidence=? "
                  "where uniprot=? and site=? and (kinase is null or length(trim(kinase))<=1)",
                  (kin, pm, ev, up, site)).rowcount
    print(up, site, kin, "rows updated:", n)
c.commit()
print("records with kinase:", c.execute("select count(*) from ptm_sites where kinase is not null and length(trim(kinase))>1").fetchone()[0])
print("distinct kinases:", c.execute("select count(distinct kinase) from ptm_sites where kinase is not null and length(trim(kinase))>1").fetchone()[0])
