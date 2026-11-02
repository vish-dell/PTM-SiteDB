"""Import a TSV exported from the ORIGINAL PTM-SiteDB into the local database.

    python import_tsv.py export.tsv

Required columns (header names, case-insensitive): gene, uniprot, ptm_type, site.
Optional: position, flanking, kinase, disease_context, evidence, confidence, pmid.
Rows missing a required field are skipped and reported. Nothing is filled in or guessed.
"""
import csv, datetime, re, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from services.dbutil import db, ensure_schema

REQ = ["gene", "uniprot", "ptm_type", "site"]
OPT = ["position", "flanking", "kinase", "disease_context", "evidence", "confidence", "pmid"]

def main(path):
    ensure_schema()
    conn = db()
    before = conn.execute("SELECT COUNT(*) FROM ptm_sites").fetchone()[0]
    ok = bad = 0
    with open(path, encoding="utf-8-sig", newline="") as f:
        rd = csv.DictReader(f, delimiter="\t")
        rd.fieldnames = [h.strip().lower() for h in rd.fieldnames or []]
        missing = [c for c in REQ if c not in rd.fieldnames]
        if missing:
            sys.exit(f"Missing required column(s): {missing}. Found: {rd.fieldnames}")
        for r in rd:
            if not all((r.get(c) or "").strip() for c in REQ):
                bad += 1; continue
            pos = (r.get("position") or "").strip()
            if not pos.isdigit():
                m = re.search(r"(\d+)$", r["site"].strip()); pos = m.group(1) if m else ""
            vals = [r[c].strip() for c in REQ] + [int(pos) if pos else None] + [(r.get(c) or "").strip() for c in OPT[1:]]
            conn.execute("""INSERT OR IGNORE INTO ptm_sites
                (gene,uniprot,ptm_type,site,position,flanking,kinase,disease_context,evidence,confidence,pmid)
                VALUES (?,?,?,?,?,?,?,?,?,?,?)""", vals)
            ok += 1
    after = conn.execute("SELECT COUNT(*) FROM ptm_sites").fetchone()[0]
    conn.execute("INSERT INTO sync_log(synced_at,source,pages,records_parsed,records_new) VALUES (?,?,?,?,?)",
                 (datetime.datetime.now(datetime.timezone.utc).replace(tzinfo=None).isoformat(timespec="seconds") + "Z", f"file:{Path(path).name}", None, ok, after - before))
    conn.commit(); conn.close()
    print(f"Read {ok} valid rows ({bad} skipped); {after - before} new records added.")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
