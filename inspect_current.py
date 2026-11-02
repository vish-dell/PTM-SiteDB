import sqlite3

db = r"instance\ptm.db"

conn = sqlite3.connect(db)
cur = conn.cursor()

print("\n===== ACTIVE DATABASE =====")

print("Records:",
      cur.execute("SELECT COUNT(*) FROM ptm_sites").fetchone()[0])

print("Unique sites:",
      cur.execute("""
      SELECT COUNT(*)
      FROM (
          SELECT DISTINCT uniprot, ptm_type, site
          FROM ptm_sites
      )
      """).fetchone()[0])

print("Genes:")

for row in cur.execute("""
    SELECT gene, COUNT(*), COUNT(DISTINCT site)
    FROM ptm_sites
    GROUP BY gene
    ORDER BY gene
"""):
    print(row)

print("\n===== CURRENT KINASE ANNOTATIONS =====")

for row in cur.execute("""
    SELECT gene, site, kinase, kinase_pmid,
           kinase_evidence
    FROM ptm_sites
    WHERE kinase IS NOT NULL
    AND TRIM(kinase) NOT IN ('', '—')
    ORDER BY gene, site
"""):
    print(row)

print("\n===== QUARANTINE =====")

print("Quarantined records:",
      cur.execute("""
      SELECT COUNT(*) FROM ptm_sites_quarantine
      """).fetchone()[0])

print("\n===== RETRACTED PMID CHECK =====")

for row in cur.execute("""
    SELECT id, gene, site, pmid
    FROM ptm_sites
    WHERE pmid LIKE '%11279024%'
"""):
    print("STILL ACTIVE:", row)

conn.close()

print("\nREAD ONLY AUDIT COMPLETE")