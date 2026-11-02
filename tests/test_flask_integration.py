"""End-to-end Flask tests through the real routes. Only external network sources are mocked."""
import os
import sqlite3
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from test_proximity_v2 import FakeCbio, model, rec     # noqa: E402  (shared fixtures)

UNIPROT = "P04637"
# R248 sits 5 A from S15 (233 residues away in sequence) -> the only qualifying eligible residue.
STRUCT = {15: ("SER", (0, 0, 0), 92, []), 248: ("ARG", (5, 0, 0), 91, []),
          175: ("ARG", (60, 0, 0), 88, []), 282: ("ARG", (90, 0, 0), 85, []), 300: ("ALA", (150, 0, 0), 80, [])}


@pytest.fixture()
def client(tmp_path, monkeypatch):
    from services import dbutil, external as ext
    dbfile = tmp_path / "t.db"
    monkeypatch.setattr(dbutil, "DB", dbfile)
    dbutil.ensure_schema()
    conn = sqlite3.connect(dbfile)
    conn.execute("INSERT INTO ptm_sites(gene,uniprot,ptm_type,site,position,flanking,evidence,confidence,pmid) "
                 "VALUES('TP53',?,'Phosphorylation','S15',15,'','in vivo','high','1')", (UNIPROT,))
    conn.commit(); conn.close()
    monkeypatch.setattr(ext, "alphafold_model", lambda acc: ("AF-P04637-F1", model(STRUCT)))
    import app as appmod
    appmod.app.config["TESTING"] = True
    return appmod.app.test_client(), ext, monkeypatch


def cb(**kw):
    kw.setdefault("genes", {"TP53": 7157})
    kw.setdefault("mutations", [rec(1, "R248Q"), rec(2, "R175H"), rec(3, "R282W"),
                                rec(4, "S15S", "Silent"), rec(5, "R248*", "Nonsense_Mutation")])
    return FakeCbio(**kw)


def test_vendored_3dmol_is_served_and_pages_use_it(client):
    c, _, _ = client
    r = c.get("/static/vendor/3Dmol-min.js")
    assert r.status_code == 200 and len(r.data) > 100_000
    for url in (f"/viewer/{UNIPROT}", "/proximity", "/spatial"):
        html = c.get(url).get_data(as_text=True)
        assert "vendor/3Dmol-min.js" in html and "js/ptmviewer.js" in html, url
        assert html.index("vendor/3Dmol-min.js") < html.index("cdn.jsdelivr.net")      # local copy first


def test_proximity_page_exposes_configurable_controls(client):
    c, _, _ = client
    html = c.get(f"/proximity?uniprot={UNIPROT}").get_data(as_text=True)
    for ident in ("p-threshold", "p-minsep", "p-distmode", "p-plddt", "p-nperm", "p-seed", "p-study"):
        assert f'id="{ident}"' in html


def test_proximity_api_end_to_end(client):
    c, ext, mp = client
    mp.setattr(ext, "fetch", cb())
    j = c.get(f"/api/f/proximity/{UNIPROT}?study=s1&n_perm=400&seed=3").get_json()
    assert j["available"] and j["structure_source"] == "alphafold" and "PREDICTED" in j["structure_label"]
    assert j["retrieval"]["profile_id"] == "s1_mutations" and j["retrieval"]["complete"]
    assert j["retrieval"]["n_sequenced"] == 321 and j["cbio_check"]["n_records_raw"] == 5
    # eligible = {248,175,282,300}; qualifying = {248}; mutated eligible = 3 -> P(X>=1) = 1 - 1/4
    assert (j["n_eligible_residues"], j["n_qualifying_residues"], j["n_mutated_positions"]) == (4, 1, 3)
    assert j["p_value"] == pytest.approx(0.75) and j["fold_enrichment"] == pytest.approx(4 / 3, abs=1e-3)
    assert j["n_permutations"] == 400 and j["seed"] == 3 and 0 < j["p_perm"] <= 1
    assert [h["pos"] for h in j["hits"]] == [248] and j["hits"][0]["dist"] == pytest.approx(5.0)
    assert j["n_synonymous"] == 1                                                # "Silent" recognised
    assert any("Nonsense_Mutation" in k for k in j["exclusions"]["missense_mutations"])
    assert j["highlight"]["ptm"][0]["resi"] == 15 and j["highlight"]["mutations"][0]["resi"] == 248


def test_zero_mutation_study_is_not_estimable_and_says_so(client):
    c, ext, mp = client
    mp.setattr(ext, "fetch", cb(mutations=[]))
    j = c.get(f"/api/f/proximity/{UNIPROT}?study=s1").get_json()
    assert j["available"] and j["cbio_check"]["n_records_raw"] == 0
    assert j["p_value"] is None and j["p_perm"] is None and j["fold_enrichment"] is None
    assert "not estimable" in j["p_value_reason"] and j["retrieval"]["status"] == "ok_no_mutations"


@pytest.mark.parametrize("kw,status", [
    ({"studies": set()}, "study_not_found"),
    ({"profiles": []}, "no_mutation_profile"),
    ({"lists": []}, "no_sample_list"),
    ({"mutation_response": {"err": 1}}, "malformed_response"),
])
def test_retrieval_failures_are_reported_not_turned_into_zero(client, kw, status):
    c, ext, mp = client
    mp.setattr(ext, "fetch", cb(**kw))
    j = c.get(f"/api/f/proximity/{UNIPROT}?study=s1").get_json()
    assert j["available"] is False and j["retrieval_status"] == status
    assert "p_value" not in j and "n_missense" not in j
    assert "no mutation data from any other study" in " ".join(j["notes"]).lower()


def test_bad_parameters_and_missing_study(client):
    c, ext, mp = client
    mp.setattr(ext, "fetch", cb())
    assert c.get(f"/api/f/proximity/{UNIPROT}").get_json()["available"] is False
    assert c.get(f"/api/f/proximity/{UNIPROT}?study=s1&dist_mode=bogus").get_json()["available"] is False
    j = c.get(f"/api/f/proximity/{UNIPROT}?study=s1&min_sep=9999&n_perm=999999").get_json()   # clamped, not crashed
    assert j["min_sep"] == 200 and j["n_perm_requested"] == 20000


def test_model_endpoint_serves_alphafold_coordinates(client):
    c, _, _ = client
    r = c.get(f"/api/f/spatial/{UNIPROT}/model?source=alphafold&id=AF-P04637-F1")
    assert r.status_code == 200 and r.get_data(as_text=True).lstrip().startswith("ATOM")
    assert c.get(f"/api/f/spatial/{UNIPROT}/model?source=nope&id=x").status_code == 400
