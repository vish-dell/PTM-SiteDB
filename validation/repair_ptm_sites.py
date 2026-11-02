"""Repair the legacy 67-row ptm_sites table using the dbPTM table you already downloaded.

Why: diagnostics showed that many legacy rows cannot map to any structure because the
record itself is wrong (position beyond the protein, wrong residue letter, or a flanking
sequence that belongs to a different site). No offset or viewer change can fix a wrong record.

What it does (DRY RUN by default; nothing is written unless you pass --apply):
  1. Backs up the database (instance/ptm.db.bak_before_repair).
  2. Adds a validation_status column (nothing is deleted).
  3. For every PHOSPHORYLATION row, looks up gene+position in dbPTM and checks the
     residue letter (centre of dbPTM's 21-mer flank).
        corroborated  -> kept; flanking replaced by dbPTM's real flank.
        not found     -> moved to table ptm_sites_quarantine (reversible).
        letter differs-> moved to ptm_sites_quarantine.
  4. Non-phospho rows (ubiquitin / acetyl / ...) are kept and marked
     'unverified_not_in_phospho_source' because the dbPTM file is phosphorylation-only.

Run from C:\\ptm:
  .\\.venv312\\Scripts\\python.exe .\\validation\\repair_ptm_sites.py            (dry run)
  .\\.venv312\\Scripts\\python.exe .\\validation\\repair_ptm_sites.py --apply    (write)
"""
import argparse
import csv
import shutil
import sqlite3
import sys
from pathlib import Path

ap = argparse.ArgumentParser()
ap.add_argument("--db", default="instance/ptm.db")
ap.add_argument("--dbptm", default="research/data/processed/ptm_sites.tsv")
ap.add_argument("--apply", action="store_true")
a = ap.parse_args()

db_path, ref_path = Path(a.db), Path(a.dbptm)
for p in (db_path, ref_path):
    if not p.exists():
        sys.exit(f"Missing file: {p}")

# gene+position -> (residue letter, 21-mer flank) from dbPTM
genes = set()
con = sqlite3.connect(db_path)
con.row_factory = sqlite3.Row
rows = con.execute("SELECT * FROM ptm_sites ORDER BY gene, position").fetchall()
genes = {r["gene"] for r in rows}

ref = {}
with open(ref_path, newline="", encoding="utf-8") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        if r["gene"] in genes:
            fl = (r["flanking_sequence"] or "").strip().upper()
            if len(fl) >= 11:
                ref[(r["gene"], int(r["position"]))] = (fl[10], fl)

plan = []
for r in rows:
    site = (r["site"] or "").strip().upper()
    letter, pos = site[:1], r["position"]
    if (r["ptm_type"] or "") != "Phosphorylation":
        plan.append((r, "keep", "unverified_not_in_phospho_source", None))
        continue
    hit = ref.get((r["gene"], pos))
    if hit is None:
        plan.append((r, "quarantine", "position_not_in_dbPTM", None))
    elif hit[0] != letter:
        plan.append((r, "quarantine", f"letter_mismatch_dbPTM_has_{hit[0]}", None))
    else:
        plan.append((r, "keep", "corroborated_by_dbPTM", hit[1]))

print(f"{'ID':>3} {'site':<8} {'action':<11} reason")
for r, act, why, _ in plan:
    print(f"{r['id']:>3} {r['site']:<8} {act:<11} {why}")
keep = sum(1 for p in plan if p[1] == "keep")
print(f"\nkeep {keep}, quarantine {len(plan) - keep}  (total {len(plan)})")

if not a.apply:
    print("\nDRY RUN - nothing written. Re-run with --apply to write.")
    sys.exit(0)

bak = db_path.with_name(db_path.name + ".bak_before_repair")
if not bak.exists():
    shutil.copy2(db_path, bak)
    print("Backup written:", bak)

cols = [c[1] for c in con.execute("PRAGMA table_info(ptm_sites)")]
if "validation_status" not in cols:
    con.execute("ALTER TABLE ptm_sites ADD COLUMN validation_status TEXT")
con.execute("CREATE TABLE IF NOT EXISTS ptm_sites_quarantine AS SELECT * FROM ptm_sites WHERE 0")
qcols = [c[1] for c in con.execute("PRAGMA table_info(ptm_sites_quarantine)")]
if "validation_status" not in qcols:
    con.execute("ALTER TABLE ptm_sites_quarantine ADD COLUMN validation_status TEXT")

for r, act, why, flank in plan:
    if act == "keep":
        if flank:
            con.execute("UPDATE ptm_sites SET validation_status=?, flanking=? WHERE id=?",
                        (why, flank, r["id"]))
        else:
            con.execute("UPDATE ptm_sites SET validation_status=? WHERE id=?", (why, r["id"]))
    else:
        con.execute("UPDATE ptm_sites SET validation_status=? WHERE id=?", (why, r["id"]))
        names = ",".join(c for c in [c[1] for c in con.execute("PRAGMA table_info(ptm_sites)")])
        con.execute(f"INSERT INTO ptm_sites_quarantine ({names}) SELECT {names} FROM ptm_sites WHERE id=?",
                    (r["id"],))
        con.execute("DELETE FROM ptm_sites WHERE id=?", (r["id"],))
con.commit()
con.close()
print("Applied. Quarantined rows are in table ptm_sites_quarantine; original file is", bak.name)