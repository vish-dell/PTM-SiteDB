"""Database helpers shared by the app and the analysis blueprint."""
import os, sqlite3
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent
# PTM_DB_PATH lets you point the app at another SQLite file (used for testing only).
DB = Path(os.environ.get("PTM_DB_PATH", BASE / "instance" / "ptm.db"))

SCHEMA = """
CREATE TABLE IF NOT EXISTS ptm_sites (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    gene TEXT NOT NULL, uniprot TEXT NOT NULL, ptm_type TEXT NOT NULL,
    site TEXT NOT NULL, position INTEGER, flanking TEXT, kinase TEXT,
    disease_context TEXT, evidence TEXT, confidence TEXT, pmid TEXT,
    source_url TEXT, protein_id TEXT,
    UNIQUE(uniprot, ptm_type, site, pmid)
);
CREATE TABLE IF NOT EXISTS sync_log (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    synced_at TEXT NOT NULL, source TEXT NOT NULL,
    pages INTEGER, records_parsed INTEGER, records_new INTEGER
);
"""

def db():
    conn = sqlite3.connect(DB)
    conn.row_factory = sqlite3.Row
    return conn

def ensure_schema():
    DB.parent.mkdir(parents=True, exist_ok=True)
    conn = db()
    conn.executescript(SCHEMA)
    conn.commit()
    conn.close()
