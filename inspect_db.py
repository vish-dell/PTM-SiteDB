import sqlite3, csv
c = sqlite3.connect("instance/ptm.db")
total = c.execute("select count(*) from ptm_sites").fetchone()[0]
uniq = c.execute("select count(*) from (select distinct uniprot, site from ptm_sites)").fetchone()[0]
print("rows:", total, "| unique sites (protein + residue):", uniq, "| extra rows for the same site:", total - uniq)
print("unique sites per gene:", c.execute("select gene, count(*) from (select distinct gene, uniprot, site from ptm_sites) group by gene order by 2 desc").fetchall())
dups = c.execute("select gene, site, count(*) from ptm_sites group by uniprot, site having count(*) > 1 order by gene, site").fetchall()
print("sites with more than one row:", len(dups), dups[:12])
with open("verification_report.csv", newline="", encoding="utf-8") as f:
    r = csv.reader(f)
    print("CSV header:", next(r))
    print("first row :", next(r, None))
