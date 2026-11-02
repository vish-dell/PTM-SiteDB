import sqlite3
REASON = "RECORD_ERROR: UniProt has a different residue here"
bad = [("P15056","S313"),("P15056","S358"),("P15056","S412"),("P15056","S424"),("P15056","S449"),("P15056","S486")]
c = sqlite3.connect("instance/ptm.db")
c.execute("attach database 'instance/ptm.db.bak_465' as bak")
have = {(u, s) for u, s in c.execute("select uniprot, site from ptm_sites_quarantine where validation_status like 'RECORD_ERROR%'")}
n = 0
for u, s in bad:
    if (u, s) in have:
        continue
    n += c.execute("""insert into ptm_sites_quarantine (gene, uniprot, ptm_type, site, position, flanking, kinase,
        disease_context, evidence, confidence, pmid, source_url, protein_id, validation_status)
        select gene, uniprot, ptm_type, site, position, flanking, kinase, disease_context, evidence, confidence,
        pmid, source_url, protein_id, ? from bak.ptm_sites where uniprot = ? and site = ?""", (REASON, u, s)).rowcount
c.commit()
print("rows added to quarantine:", n)
print("quarantine total:", c.execute("select count(*) from ptm_sites_quarantine").fetchone()[0])
