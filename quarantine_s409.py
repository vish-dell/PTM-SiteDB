import sqlite3
REASON = "FLANK_MISMATCH: 0/20 flank match against UniProt (numbering offset suspected)"
c = sqlite3.connect("instance/ptm.db")
cols = "gene, uniprot, ptm_type, site, position, flanking, kinase, disease_context, evidence, confidence, pmid, source_url, protein_id"
n = c.execute("insert into ptm_sites_quarantine (" + cols + ", validation_status) select " + cols + ", ? from ptm_sites where uniprot = 'P15056' and site = 'S409'", (REASON,)).rowcount
d = c.execute("delete from ptm_sites where uniprot = 'P15056' and site = 'S409'").rowcount
c.commit()
print("copied to quarantine:", n, "| removed from ptm_sites:", d)
print("rows now:", c.execute("select count(*) from ptm_sites").fetchone()[0])
