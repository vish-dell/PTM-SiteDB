"""FINAL: keep only fully verified PTM rows.

A row stays in ptm_sites ONLY if ALL of these pass:
  1. phosphorylation row (the only type a source here can validate)
  2. UniProt (the row's own accession) has the recorded residue letter at the recorded position
  3. dbPTM has the same gene + same accession + same position, and its residue letter agrees
  4. dbPTM's full 21-residue flank matches the UniProt sequence EXACTLY around that position
Everything else is moved to ptm_sites_quarantine (reversible; nothing is deleted).
Rows that stay get dbPTM's verified flank (padding '-' removed).

If UniProt cannot be reached for an accession, the rows cannot be verified and --apply REFUSES to run.
A per-row report is always written to validation/logs/final_verification.csv.

Run from C:\\ptm:
  .\\.venv312\\Scripts\\python.exe .\\validation\\finalize_ptm_sites.py            (check only)
  .\\.venv312\\Scripts\\python.exe .\\validation\\finalize_ptm_sites.py --apply    (write)
Needs internet (UniProt). Sequences are cached in validation/uniprot_cache/.
"""
import argparse, csv, json, shutil, sqlite3, sys, urllib.request
from pathlib import Path

ap = argparse.ArgumentParser()
ap.add_argument("--db", default="instance/ptm.db")
ap.add_argument("--dbptm", default="research/data/processed/ptm_sites.tsv")
ap.add_argument("--cache", default="validation/uniprot_cache")
ap.add_argument("--report", default="validation/logs/final_verification.csv")
ap.add_argument("--apply", action="store_true")
a = ap.parse_args()

for p in (a.db, a.dbptm):
    if not Path(p).exists():
        sys.exit(f"Missing file: {p}")
Path(a.cache).mkdir(parents=True, exist_ok=True)
Path(a.report).parent.mkdir(parents=True, exist_ok=True)


def uniprot_seq(acc):
    f = Path(a.cache) / f"{acc}.fasta"
    if not f.exists():
        req = urllib.request.Request(f"https://rest.uniprot.org/uniprotkb/{acc}.fasta",
                                     headers={"User-Agent": "PTMSiteDB-final-verify/1.0"})
        with urllib.request.urlopen(req, timeout=30) as r:
            f.write_text(r.read().decode())
    return "".join(l.strip() for l in f.read_text().splitlines() if not l.startswith(">")).upper()


con = sqlite3.connect(a.db)
con.row_factory = sqlite3.Row
rows = con.execute("SELECT * FROM ptm_sites ORDER BY gene, position").fetchall()
genes = {r["gene"] for r in rows}

ref = {}
with open(a.dbptm, newline="", encoding="utf-8") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        if r["gene"] in genes:
            ref[(r["gene"], r["uniprot"], int(r["position"]))] = (r["flanking_sequence"] or "").strip().upper()

seqs, fetch_errors = {}, {}
for acc in sorted({r["uniprot"] for r in rows}):
    try:
        seqs[acc] = uniprot_seq(acc)
    except Exception as e:  # network or accession problem: never guessed around
        fetch_errors[acc] = str(e)

report, keep_ids, unverifiable = [], {}, 0
for r in rows:
    site, pos, acc = (r["site"] or "").strip().upper(), r["position"], r["uniprot"]
    letter, why = site[:1], None
    seq = seqs.get(acc)
    fl = ref.get((r["gene"], acc, pos))
    if (r["ptm_type"] or "") != "Phosphorylation":
        why = "not_phosphorylation_no_validation_source"
    elif seq is None:
        why, unverifiable = "UNVERIFIABLE_uniprot_not_reachable", unverifiable + 1
    elif not (1 <= pos <= len(seq)):
        why = f"position_beyond_protein_length_{len(seq)}"
    elif seq[pos - 1] != letter:
        why = f"uniprot_has_{seq[pos - 1]}_not_{letter}"
    elif fl is None:
        why = "not_in_dbPTM_for_this_accession"
    elif len(fl) != 21 or fl[10] != letter:
        why = "dbPTM_letter_or_flank_length_disagrees"
    else:
        window = "".join(seq[i - 1] if 1 <= i <= len(seq) else "-" for i in range(pos - 10, pos + 11))
        if window != fl:
            why = "dbPTM_flank_does_not_match_uniprot_window"
    ok = why is None
    if ok:
        keep_ids[r["id"]] = fl.replace("-", "")
    report.append({"id": r["id"], "gene": r["gene"], "accession": acc, "site": site, "position": pos,
                   "verdict": "VERIFIED" if ok else "REJECTED", "reason": "all_checks_passed" if ok else why})

with open(a.report, "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=list(report[0]))
    w.writeheader(); w.writerows(report)

print(f"{'ID':>3} {'gene':<7} {'site':<7} verdict  reason")
for x in report:
    print(f"{x['id']:>3} {x['gene']:<7} {x['site']:<7} {x['verdict']:<8} {x['reason']}")
print(f"\nVERIFIED {len(keep_ids)} | REJECTED {len(rows) - len(keep_ids)} | total {len(rows)}")
print("Report:", a.report)
if fetch_errors:
    print("UniProt fetch failed for:", fetch_errors)

if not a.apply:
    print("\nCHECK ONLY - nothing written. Re-run with --apply to write.")
    sys.exit(0)
if unverifiable:
    sys.exit(f"\nREFUSING to apply: {unverifiable} row(s) could not be checked against UniProt. "
             "Fix your internet connection and run again. Nothing was written.")

bak = Path(a.db + ".bak_before_final")
if not bak.exists():
    shutil.copy2(a.db, bak)
    print("Backup:", bak)
cols = [c[1] for c in con.execute("PRAGMA table_info(ptm_sites)")]
if "validation_status" not in cols:
    con.execute("ALTER TABLE ptm_sites ADD COLUMN validation_status TEXT")
    cols.append("validation_status")
con.execute("CREATE TABLE IF NOT EXISTS ptm_sites_quarantine AS SELECT * FROM ptm_sites WHERE 0")
if "validation_status" not in [c[1] for c in con.execute("PRAGMA table_info(ptm_sites_quarantine)")]:
    con.execute("ALTER TABLE ptm_sites_quarantine ADD COLUMN validation_status TEXT")
names = ",".join(cols)
by_id = {x["id"]: x for x in report}
for r in rows:
    x = by_id[r["id"]]
    if r["id"] in keep_ids:
        con.execute("UPDATE ptm_sites SET flanking=?, validation_status='verified_uniprot_and_dbPTM' WHERE id=?",
                    (keep_ids[r["id"]], r["id"]))
    else:
        con.execute("UPDATE ptm_sites SET validation_status=? WHERE id=?", (x["reason"], r["id"]))
        con.execute(f"INSERT INTO ptm_sites_quarantine ({names}) SELECT {names} FROM ptm_sites WHERE id=?", (r["id"],))
        con.execute("DELETE FROM ptm_sites WHERE id=?", (r["id"],))
con.commit(); con.close()
print("Applied. ptm_sites now holds only verified rows; the rest are in ptm_sites_quarantine.")