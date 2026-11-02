"""Regression tests for the corrected 3D-proximal mutation-PTM pipeline.

Everything here uses small hand-built fixtures (never app data). Expected values are derived by hand or by an
independent implementation (scipy), not by running the code under test.
"""
import math
import os
import random
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from services import external as ext          # noqa: E402
from services import proximity as P           # noqa: E402


# ------------------------------------------------------------------ fixture helpers
def atom(serial, name, res, resi, xyz, b, chain="A"):
    x, y, z = xyz
    return (f"ATOM  {serial:5d} {name:<4s} {res:>3s} {chain}{resi:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00{b:6.2f}           {name.strip()[0]}")


def model(residues):
    """residues: {resi: (resname, ca_xyz, plddt, [extra (name, xyz)...])} -> PDB text."""
    lines, n = [], 0
    for resi in sorted(residues):
        resn, ca, b, extra = residues[resi]
        n += 1
        lines.append(atom(n, " CA ", resn, resi, ca, b))
        for name, xyz in extra:
            n += 1
            lines.append(atom(n, f" {name:<3s}", resn, resi, xyz, b))
    return "\n".join(lines)


def use_model(monkeypatch, text):
    monkeypatch.setattr(P.ext, "alphafold_model", lambda acc: ("AF-TEST-F1", text))


def mut(pos, wt, alt, typ="Missense_Mutation"):
    return {"pos": pos, "pos_end": pos, "wt": wt, "alt": alt, "protein_change": f"{wt}{pos}{alt}", "type": typ}


PTM = [{"pos": 1, "site": "S1", "ptm_type": "Phosphorylation", "aa": "S"}]


# ------------------------------------------------------------------ mathematics (independent checks)
def test_hypergeometric_matches_scipy():
    from scipy.stats import hypergeom
    rng = random.Random(7)
    for _ in range(200):
        N = rng.randint(5, 400)
        K = rng.randint(0, N)
        n = rng.randint(0, N)
        k = rng.randint(0, min(K, n) + 1)
        assert P.hyper_sf(k, N, K, n) == pytest.approx(float(hypergeom.sf(k - 1, N, K, n)), abs=1e-9)


def test_bh_matches_statsmodels_style_reference():
    ps = [0.001, 0.008, 0.039, 0.041, 0.042, 0.06, 0.074, 0.205, 0.212, 0.216]
    n = len(ps)
    ref, prev = [0] * n, 1.0                      # textbook step-up definition
    for rank in range(n, 0, -1):
        prev = min(prev, ps[rank - 1] * n / rank)
        ref[rank - 1] = prev
    assert P.bh_fdr(ps) == pytest.approx(ref)
    assert all(0 <= q <= 1 for q in P.bh_fdr(ps))


def test_permutation_null_agrees_with_exact_when_composition_is_uniform():
    # one amino-acid letter => composition matching changes nothing => permutation p ~ hypergeometric p
    L, K, n, k = 200, 30, 12, 4
    seq = "A" * L
    universe = list(range(1, L + 1))
    target = set(range(1, K + 1))
    mutated = list(range(1, k + 1)) + list(range(100, 100 + n - k))
    hits, why = P._perm_hits(seq, universe, mutated, target, 20000, seed=11)
    assert why is None and len(hits) == 20000
    assert P._empirical_p(hits, k) == pytest.approx(P.hyper_sf(k, L, K, n), abs=0.02)


def test_permutation_is_reproducible_and_seed_dependent():
    seq = "SAGTKLE" * 30
    uni, mutated, tgt = list(range(1, 211)), [3, 9, 15, 40, 77], set(range(1, 40))
    a, _ = P._perm_hits(seq, uni, mutated, tgt, 500, seed=5)
    b, _ = P._perm_hits(seq, uni, mutated, tgt, 500, seed=5)
    c, _ = P._perm_hits(seq, uni, mutated, tgt, 500, seed=6)
    assert (a == b).all() and not (a == c).all()


def test_permutation_keeps_reference_residue_composition():
    seq = ("S" * 10) + ("A" * 90)
    # 2 of the mutated positions are Ser; every permutation must place exactly 2 draws on Ser residues
    hits, _ = P._perm_hits(seq, range(1, 101), [1, 2, 50], set(range(1, 11)), 300, seed=1)
    assert (hits >= 2).all() and (hits <= 2).all()    # target = the Ser block -> exactly the 2 Ser draws


def test_composition_null_impossible_reports_reason_not_a_number():
    seq = "S" * 3 + "A" * 50
    hits, why = P._perm_hits(seq, [1, 2, 10, 11, 12], [1, 2, 3], {1}, 100, seed=1)   # 3 mutated Ser, 2 eligible Ser
    assert hits is None and "too few eligible residues" in why and "S: 3 mutated vs 2 eligible" in why


# ------------------------------------------------------------------ null universe
def structure_for_universe():
    far = (200.0, 0, 0)
    return {1: ("SER", (0, 0, 0), 90, []), 3: ("ALA", (3, 0, 0), 90, []),       # 3: close in 3D, 2 apart in sequence
            50: ("ALA", (5, 0, 0), 90, []), 60: ("ALA", (40, 0, 0), 90, []),
            70: ("ALA", (80, 0, 0), 90, []), 80: ("ALA", far, 90, [])}


def test_null_universe_is_eligible_residues_only(monkeypatch):
    use_model(monkeypatch, model(structure_for_universe()))
    cb = {"mutations": [mut(50, "A", "T"), mut(3, "A", "V")], "check": {}, "n_sequenced": 100}
    r = P.proximity_analysis("TEST", PTM, cb, threshold=10, min_sep=20, n_perm=200)
    # eligible = {50, 60, 70, 80} (3 is within 20 of the PTM, 1 is the PTM): N=4, K=1 (residue 50)
    assert (r["n_eligible_residues"], r["n_qualifying_residues"]) == (4, 1)
    # residue 3 is sequence-proximal by definition: reported, not tested
    assert r["n_mutated_sequence_proximal"] == 1 and r["n_mutated_positions"] == 1
    # n=1, k=1 -> P = K/N = 1/4  (old code: N=6, n=2 -> 1/3)
    assert r["p_value"] == pytest.approx(0.25)
    assert r["expected"] == pytest.approx(0.25)
    assert r["fold_enrichment"] == pytest.approx(4.0)


# ------------------------------------------------------------------ mutation mapping checks
def test_reference_residue_mismatch_is_excluded_with_reason(monkeypatch):
    st = structure_for_universe()
    st[50] = ("SER", (5, 0, 0), 90, [])                    # model has Ser at 50
    use_model(monkeypatch, model(st))
    cb = {"mutations": [mut(50, "A", "T"),                 # record says Ala at 50 -> numbering/isoform mismatch
                        mut(60, "A", "T"),
                        mut(60, "A", "V", "Nonsense_Mutation"),
                        mut(500, "A", "T")],               # outside the model
          "check": {}, "n_sequenced": 10}
    r = P.proximity_analysis("TEST", PTM, cb, threshold=10, min_sep=20, n_perm=100)
    ex = r["exclusions"]["missense_mutations"]
    assert any("disagrees with the model" in k for k in ex)
    assert any("outside the structure model" in k for k in ex)
    assert any("Nonsense_Mutation" in k for k in ex)
    assert r["n_missense"] == 3 and r["n_missense_mapped"] == 1      # only A60T is verified
    assert [h["pos"] for h in r["hits"]] == []                        # 60 is 35 A away


def test_low_confidence_residues_and_ptms_are_excluded(monkeypatch):
    st = structure_for_universe()
    st[50] = ("ALA", (5, 0, 0), 55, [])                    # pLDDT 55
    use_model(monkeypatch, model(st))
    r = P.proximity_analysis("TEST", PTM + [{"pos": 70, "site": "S70", "ptm_type": "Phosphorylation", "aa": "S"}],
                             {"mutations": [mut(50, "A", "T")], "check": {}}, n_perm=50)
    assert any("low confidence" in k for k in r["exclusions"]["missense_mutations"])
    # S70 record says Ser but the model has Ala -> PTM excluded, with the reason
    assert any(e["site"] == "S70" and "model has A" in e["reason"] for e in r["ptm_exclusions"])
    assert r["n_ptm_confident"] == 1


def test_silent_records_are_recognised_as_synonymous_control(monkeypatch):
    use_model(monkeypatch, model(structure_for_universe()))
    cb = {"mutations": [mut(50, "A", "A", "Silent"), mut(60, "A", "T")], "check": {}}
    r = P.proximity_analysis("TEST", PTM, cb, n_perm=50)
    assert r["n_synonymous"] == 1 and r["negative_control"]["n_observed"] == 1
    r2 = P.proximity_analysis("TEST", PTM, {"mutations": [mut(60, "A", "T")], "check": {}}, n_perm=50)
    assert r2["n_synonymous"] == 0 and "not estimable" in r2["negative_control_note"]


# ------------------------------------------------------------------ never turn "missing" into a number
def test_zero_mutations_is_not_estimable_not_zero_or_one(monkeypatch):
    use_model(monkeypatch, model(structure_for_universe()))
    r = P.proximity_analysis("TEST", PTM, {"mutations": [], "check": {"n_records_raw": 0}}, n_perm=50)
    assert r["p_value"] is None and r["p_perm"] is None
    assert r["fold_enrichment"] is None and r["fold_ci95"] is None
    assert "zero mutation records" in r["p_value_reason"] and "not estimable" in r["p_value_reason"]
    assert "not evidence of no effect" in r["note"]


def test_no_eligible_residue_reports_reason(monkeypatch):
    tight = {1: ("SER", (0, 0, 0), 90, []), 5: ("ALA", (4, 0, 0), 90, []), 9: ("ALA", (6, 0, 0), 90, [])}
    use_model(monkeypatch, model(tight))
    r = P.proximity_analysis("TEST", PTM, {"mutations": [mut(5, "A", "T")], "check": {}}, min_sep=20, n_perm=50)
    assert r["p_value"] is None and "no eligible residue" in r["note"]
    assert r["n_eligible_residues"] == 0 and r["n_mutated_sequence_proximal"] == 1


def test_untestable_permutation_has_specific_reason():
    cls = {1: {"eligible": False, "qualifies": False, "dist": 0, "nearest_ptm_pos": 1, "min_seq_sep": 0},
           60: {"eligible": True, "qualifies": True, "dist": 4.0, "nearest_ptm_pos": 1, "min_seq_sep": 59},
           61: {"eligible": True, "qualifies": False, "dist": 50.0, "nearest_ptm_pos": 1, "min_seq_sep": 60}}
    r = P.proximity_test(cls, [60], sequence=None)
    assert r["p_value"] == pytest.approx(0.5) and r["p_perm"] is None
    assert "no usable protein sequence" in r["p_perm_note"]


# ------------------------------------------------------------------ distance definition
def test_distance_modes_differ_and_are_configurable(monkeypatch):
    # C-alpha 14 A apart, but side-chain atoms 3 A apart
    st = {1: ("SER", (0, 0, 0), 90, [("OG", (7, 0, 0))]),
          50: ("ALA", (14, 0, 0), 90, [("CB", (10, 0, 0))]),
          90: ("ALA", (300, 0, 0), 90, [])}
    use_model(monkeypatch, model(st))
    cb = {"mutations": [mut(50, "A", "T")], "check": {}}
    ca = P.proximity_analysis("TEST", PTM, cb, threshold=10, min_sep=20, n_perm=50, dist_mode="ca")
    ma = P.proximity_analysis("TEST", PTM, cb, threshold=10, min_sep=20, n_perm=50, dist_mode="min_atom")
    assert ca["n_qualifying_residues"] == 0 and ca["dist_mode"] == "ca"
    assert ma["n_qualifying_residues"] == 1 and ma["hits"][0]["dist"] == pytest.approx(3.0, abs=0.01)
    assert "C-alpha" in ca["distance_definition"] and "any heavy-atom pair" in ma["distance_definition"]
    bad = P.proximity_analysis("TEST", PTM, cb, dist_mode="nonsense")
    assert bad["available"] is False


def test_reported_distance_equals_geometry(monkeypatch):
    st = {1: ("SER", (0, 0, 0), 90, []), 50: ("ALA", (3, 4, 0), 90, []), 90: ("ALA", (300, 0, 0), 90, [])}
    use_model(monkeypatch, model(st))
    r = P.proximity_analysis("TEST", PTM, {"mutations": [mut(50, "A", "T")], "check": {}}, n_perm=50)
    assert r["hits"][0]["dist"] == pytest.approx(5.0)          # 3-4-5 triangle
    assert r["highlight"]["mutations"][0]["resi"] == 50 and r["highlight"]["ptm"][0]["resi"] == 1


def test_results_are_reproducible_for_a_seed(monkeypatch):
    use_model(monkeypatch, model(structure_for_universe()))
    cb = {"mutations": [mut(50, "A", "T"), mut(60, "A", "T")], "check": {}}
    a = P.proximity_analysis("TEST", PTM, cb, n_perm=300, perm_seed=4)
    b = P.proximity_analysis("TEST", PTM, cb, n_perm=300, perm_seed=4)
    assert a["p_perm"] == b["p_perm"] and a["n_permutations"] == 300 and 0 < a["p_perm"] <= 1


# ------------------------------------------------------------------ cBioPortal retrieval (mock server)
class FakeCbio:
    """Routes ext.fetch calls to canned cBioPortal answers."""

    def __init__(self, **kw):
        self.studies = kw.get("studies", {"s1"})
        self.profiles = kw.get("profiles", [{"molecularProfileId": "s1_mutations", "molecularAlterationType": "MUTATION_EXTENDED"}])
        self.lists = kw.get("lists", [{"sampleListId": "s1_sequenced", "category": "all_cases_with_mutation_data", "sampleCount": 321}])
        self.genes = kw.get("genes", {"TP53": 7157})
        self.mutations = kw.get("mutations", [])
        self.mutation_response = kw.get("mutation_response")        # override raw payload
        self.ignore_page = kw.get("ignore_page", False)
        self.calls = []

    def __call__(self, url, params=None, **kw):
        self.calls.append((url, params))
        tail = url.split("/api/")[-1]
        if tail.startswith("studies/") and tail.count("/") == 1:
            if tail.split("/")[1] in self.studies:
                return {"studyId": tail.split("/")[1]}
            raise ext.ExternalError("www.cbioportal.org: not found (404)")
        if tail.endswith("/molecular-profiles") and tail.startswith("studies/"):
            return self.profiles
        if tail.endswith("/sample-lists"):
            return self.lists
        if tail.startswith("genes/"):
            sym = tail.split("/")[1]
            if sym in self.genes:
                return {"hugoGeneSymbol": sym, "entrezGeneId": self.genes[sym]}
            raise ext.ExternalError("www.cbioportal.org: not found (404)")
        if tail.endswith("/mutations"):
            if self.mutation_response is not None:
                return self.mutation_response
            size, page = params["pageSize"], 0 if self.ignore_page else params["pageNumber"]
            return self.mutations[page * size:(page + 1) * size]
        raise AssertionError("unexpected URL " + url)


def rec(i, change="R248Q", typ="Missense_Mutation", start=None):
    import re
    found = re.search(r"\d+", change)
    pos = int(found.group()) if found else 0
    return {"sampleId": f"S{i}", "patientId": f"P{i}", "proteinChange": change, "mutationType": typ,
            "proteinPosStart": pos if start is None else start, "proteinPosEnd": pos if start is None else start}


def test_retrieval_success_reports_audit_trail(monkeypatch):
    f = FakeCbio(mutations=[rec(1), rec(2, "R175H"), rec(3, "R282W", "Nonsense_Mutation"), rec(4, "", "Splice_Site", 0)])
    monkeypatch.setattr(ext, "fetch", f)
    out = ext.cbio_gene_mutations("s1", "TP53")
    ret = out["retrieval"]
    assert ret["status"] == "ok" and ret["complete"] and ret["pages"] == 1
    assert ret["profile_id"] == "s1_mutations" and ret["sample_list_id"] == "s1_sequenced"
    assert out["n_sequenced"] == 321 and out["check"]["n_records_raw"] == 4
    assert out["check"]["n_skipped_no_position"] == 1                       # splice record without protein position
    assert {m["pos"] for m in out["mutations"]} == {248, 175, 282}
    assert out["mutations"][0]["wt"] == "R" and out["mutations"][0]["alt"] == "Q"


def test_zero_mutations_is_a_valid_result_not_an_error(monkeypatch):
    monkeypatch.setattr(ext, "fetch", FakeCbio(mutations=[]))
    out = ext.cbio_gene_mutations("s1", "TP53")
    assert out["mutations"] == [] and out["retrieval"]["status"] == "ok_no_mutations" and out["retrieval"]["complete"]


def test_discovers_nonstandard_profile_id(monkeypatch):
    f = FakeCbio(profiles=[{"molecularProfileId": "s1_mutations_custom", "molecularAlterationType": "MUTATION_EXTENDED"},
                           {"molecularProfileId": "s1_gistic", "molecularAlterationType": "COPY_NUMBER_ALTERATION"}],
                 mutations=[rec(1)])
    monkeypatch.setattr(ext, "fetch", f)
    out = ext.cbio_gene_mutations("s1", "TP53")
    assert out["retrieval"]["profile_id"] == "s1_mutations_custom"
    assert any("s1_mutations_custom/mutations" in c[0] for c in f.calls)


@pytest.mark.parametrize("kwargs,gene,status", [
    ({"studies": set()}, "TP53", "study_not_found"),
    ({"profiles": [{"molecularProfileId": "s1_gistic", "molecularAlterationType": "COPY_NUMBER_ALTERATION"}]}, "TP53", "no_mutation_profile"),
    ({"lists": [{"sampleListId": "s1_all", "category": "all_cases_in_study"}]}, "TP53", "no_sample_list"),
    ({}, "NOTAGENE", "gene_not_found"),
    ({"mutation_response": {"message": "oops"}}, "TP53", "malformed_response"),
])
def test_failures_have_distinct_statuses_and_are_not_zero_mutations(monkeypatch, kwargs, gene, status):
    monkeypatch.setattr(ext, "fetch", FakeCbio(**kwargs))
    with pytest.raises(ext.MutationRetrievalError) as e:
        ext.cbio_gene_mutations("s1", gene)
    assert e.value.status == status and e.value.retrieval["status"] == status


def test_pagination_collects_all_pages(monkeypatch):
    f = FakeCbio(mutations=[rec(i, f"R{100 + i}Q") for i in range(25)])
    monkeypatch.setattr(ext, "fetch", f)
    out = ext.cbio_gene_mutations("s1", "TP53", page_size=10)
    assert out["check"]["n_records_raw"] == 25 and out["retrieval"]["pages"] == 3 and out["retrieval"]["complete"]


def test_pagination_that_does_not_advance_is_reported_incomplete(monkeypatch):
    f = FakeCbio(mutations=[rec(i, f"R{100 + i}Q") for i in range(10)], ignore_page=True)
    monkeypatch.setattr(ext, "fetch", f)
    with pytest.raises(ext.MutationRetrievalError) as e:
        ext.cbio_gene_mutations("s1", "TP53", page_size=10)
    assert e.value.status == "incomplete"


def test_alphafold_entry_selection_uses_exact_accession(monkeypatch):
    rows = [{"entryId": "AF-Q99999-F1", "uniprotAccession": "Q99999", "pdbUrl": "https://x/other.pdb"},
            {"entryId": "AF-P04637-F2", "uniprotAccession": "P04637", "pdbUrl": "https://x/f2.pdb"},
            {"entryId": "AF-P04637-F1", "uniprotAccession": "P04637", "pdbUrl": "https://x/f1.pdb"}]
    monkeypatch.setattr(ext, "fetch", lambda url, **k: rows if "prediction" in url else "ATOM  1")
    assert ext.alphafold_entry("P04637")["entryId"] == "AF-P04637-F1"
    monkeypatch.setattr(ext, "fetch", lambda url, **k: rows[:1])
    with pytest.raises(ext.ExternalError):
        ext.alphafold_entry("P04637")
