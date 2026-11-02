"""Local PTM-SiteDB connection diagnostics. No network calls and no invented data."""
from pathlib import Path
import sqlite3, sys

ROOT = Path(__file__).resolve().parent
DB = ROOT / "instance" / "ptm.db"
print("Project:", ROOT)
print("Database:", DB)
print("Database file exists:", DB.exists())
if not DB.exists():
    print("RESULT: database file not created yet; run the site once or sync_data.py")
    sys.exit(2)
try:
    con = sqlite3.connect(str(DB))
    tables = [r[0] for r in con.execute("SELECT name FROM sqlite_master WHERE type='table' ORDER BY name")]
    print("Tables:", ", ".join(tables) if tables else "none")
    if "ptm_sites" not in tables:
        print("RESULT: ptm_sites table is missing. Check services/dbutil.py and restart app.py.")
        sys.exit(2)
    n = con.execute("SELECT COUNT(*) FROM ptm_sites").fetchone()[0]
    genes = con.execute("SELECT COUNT(DISTINCT gene) FROM ptm_sites").fetchone()[0]
    print("PTM records:", n)
    print("Distinct genes:", genes)
    if "sync_log" in tables:
        last = con.execute("SELECT synced_at, source, pages, records_parsed, records_new FROM sync_log ORDER BY id DESC LIMIT 1").fetchone()
        print("Last sync:", last if last else "never")
    con.close()
    print("RESULT:", "database connected and populated" if n else "database connected but EMPTY; run sync_data.py or import a genuine TSV export")
except Exception as exc:
    print("DATABASE ERROR:", repr(exc))
    sys.exit(2)

proc = ROOT / "research" / "data" / "processed"
print("\nResearch source files:")
for name in ("dbptm_phosphorylation_source.txt", "intogen_drivers_source.tsv", "ptm_sites_mapped.tsv", "benchmark_panel.csv", "source_manifest.json"):
    p = proc / name
    print(f"  {name}: " + (f"{p.stat().st_size} bytes" if p.exists() else "MISSING"))
