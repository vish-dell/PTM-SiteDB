import csv, os, sqlite3, sys
import pytest
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from test_proximity_v2 import FakeCbio, model, rec     # noqa: E402


def test_batch_columns_are_populated_and_failures_are_labelled(tmp_path, monkeypatch):
    from services import dbutil, external as ext
    monkeypatch.setattr(dbutil, "DB", tmp_path / "b.db"); dbutil.ensure_schema()
    conn = sqlite3.connect(dbutil.DB)
    for i, site in enumerate(("S15", "S16", "S17")):
        conn.execute("INSERT INTO ptm_sites(gene,uniprot,ptm_type,site,position,pmid) VALUES('TP53','P04637','Phosphorylation',?,?,?)",
                     (site, 15 + i, str(i)))
    conn.commit(); conn.close()
    st = {15: ("SER", (0, 0, 0), 92, []), 16: ("SER", (3, 0, 0), 92, []), 17: ("SER", (6, 0, 0), 92, []),
          248: ("ARG", (5, 0, 0), 91, []), 175: ("ARG", (60, 0, 0), 88, []), 282: ("ARG", (90, 0, 0), 85, [])}
    monkeypatch.setattr(ext, "alphafold_model", lambda acc: ("AF-X-F1", model(st)))
    import batch_proximity as B
    monkeypatch.setattr(ext, "fetch", FakeCbio(mutations=[rec(1, "R248Q"), rec(2, "R175H")]))
    rows, _ = B.run("s1", 20, 10.0, 3, 0, 100, 1)
    r = rows[0]
    assert r["n_eligible_residues"] == 3 and r["n_qualifying_residues"] == 1          # eligible != qualifying
    assert r["fold_enrichment"] is not None and r["retrieval_status"] == "ok"
    monkeypatch.setattr(ext, "fetch", FakeCbio(profiles=[]))
    rows, _ = B.run("s1", 20, 10.0, 3, 0, 100, 1)
    assert rows[0]["status"].startswith("not retrieved (no_mutation_profile)") and rows[0].get("p_value") is None
