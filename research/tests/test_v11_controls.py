"""Tests for the v1.1 controls: recurrence collapse, decoy-anchor control, burden-stratified test,
family-wise grid correction, and an end-to-end CLI run.

IMPORTANT: every dataset below is a tiny hand-built UNIT-TEST FIXTURE used only to check that code paths
behave as specified. None of it is, or may be reported as, biological data or a project result.
"""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("ptm_analysis_v11", ROOT / "scripts" / "annotate_and_analyze_maf.py")
analysis = importlib.util.module_from_spec(spec)
spec.loader.exec_module(analysis)


def _line_gene(n=60, cluster=range(45, 56), anchor=20, gene="G"):
    """Residues on a 30 A grid except a cluster placed within 5 A of the anchor."""
    coords = {}
    for p in range(1, n + 1):
        coords[(gene, p)] = (30.0 * p, 0.0, 0.0)
    coords[(gene, anchor)] = (0.0, 0.0, 0.0)
    for i, p in enumerate(cluster):
        coords[(gene, p)] = (3.0 + 0.1 * i, 0.0, 0.0)
    return coords


def _opp(gene="G", n=60, aa_pattern=True):
    rows = []
    for p in range(1, n + 1):
        for alt in "ACG":
            rows.append({"gene": gene, "protein_pos": p, "context96": "X", "alt": alt,
                         "genomic_position": 1000 + p * 3 + "ACG".index(alt),
                         "ref_aa": "S" if (p % 2 == 0 and aa_pattern) else "A"})
    return pd.DataFrame(rows)


def test_recurrent_hotspot_counts_once_under_primary_mode():
    gene = "G"
    coords = _line_gene(gene=gene)
    opp = _opp(gene)
    # one residue (50) next to the PTM mutated ten times
    mut = pd.DataFrame({"gene": gene, "protein_pos": [50] * 10, "context96": "X",
                        "genomic_position": 1150, "alt": "A"})
    r = analysis.setting_gene_result(gene, [20], mut, opp, coords, sep=8, dist=10.0, n_perm=2000, seed=1,
                                     structure_mapping_ok=True)
    assert r["recurrence_mode"] == "distinct_residue"
    assert r["observed_proximal_mutation_events"] == 1
    assert r["observed_proximal_events_all_variants"] == 10
    assert r["p_context_null_all_variants"] < 0.01          # legacy counting is wildly "significant"
    assert r["p_context_null"] > 0.05                        # one residue is weak evidence
    legacy = analysis.setting_gene_result(gene, [20], mut, opp, coords, sep=8, dist=10.0, n_perm=2000, seed=1,
                                          structure_mapping_ok=True, recurrence="all_events")
    assert legacy["observed_proximal_mutation_events"] == 10


def test_distinct_residues_still_accumulate_evidence():
    gene = "G"
    coords = _line_gene(gene=gene)
    opp = _opp(gene)
    mut = pd.DataFrame({"gene": gene, "protein_pos": [46, 48, 50, 52, 54], "context96": "X",
                        "genomic_position": [1138, 1144, 1150, 1156, 1162], "alt": "A"})
    r = analysis.setting_gene_result(gene, [20], mut, opp, coords, sep=8, dist=10.0, n_perm=2000, seed=1,
                                     structure_mapping_ok=True)
    assert r["observed_proximal_mutation_events"] == 5
    assert r["p_context_null"] < 0.01


def test_decoy_control_flags_ptm_specific_excess():
    gene = "G"
    coords = _line_gene(gene=gene, cluster=[60, 62, 64], n=70)       # even residues are 'S'
    opp = _opp(gene, n=70)
    mut = pd.DataFrame({"gene": gene, "protein_pos": [60, 62, 64, 10, 12], "context96": "X",
                        "genomic_position": [1180, 1186, 1192, 1030, 1036], "alt": "A"})
    r = analysis.setting_gene_result(gene, [20], mut, opp, coords, sep=8, dist=10.0, n_perm=500, seed=1,
                                     structure_mapping_ok=True, n_decoy=200)
    assert r["decoy_status"] == "ok"
    assert r["observed_proximal_mutation_events"] == 3
    assert r["p_decoy"] < 0.02
    assert r["decoy_z"] > 3 or not np.isfinite(r["decoy_z"])


def test_decoy_control_is_not_significant_for_unrelated_mutations():
    gene = "G"
    coords = _line_gene(gene=gene, cluster=[60, 62, 64], n=70)
    opp = _opp(gene, n=70)
    mut = pd.DataFrame({"gene": gene, "protein_pos": [10, 12, 14, 40], "context96": "X",
                        "genomic_position": [1030, 1036, 1042, 1120], "alt": "A"})
    r = analysis.setting_gene_result(gene, [20], mut, opp, coords, sep=8, dist=10.0, n_perm=500, seed=1,
                                     structure_mapping_ok=True, n_decoy=200)
    assert r["decoy_status"] == "ok"
    assert r["observed_proximal_mutation_events"] == 0
    assert r["p_decoy"] > 0.2


def test_decoy_control_reports_unavailable_when_pool_too_small():
    gene = "G"
    coords = _line_gene(gene=gene, n=60)
    opp = _opp(gene, n=60)
    anchors = list(range(2, 60, 2))                                  # every S residue is an anchor
    mut = pd.DataFrame({"gene": gene, "protein_pos": [50], "context96": "X", "genomic_position": 1150, "alt": "A"})
    r = analysis.setting_gene_result(gene, anchors, mut, opp, coords, sep=8, dist=10.0, n_perm=100, seed=1,
                                     structure_mapping_ok=True, n_decoy=50)
    assert r["decoy_status"] == "unavailable_decoy_pool_too_small"
    assert "decoy_z" not in r


def test_stratified_test_removes_burden_confounding():
    rng = np.random.default_rng(3)
    high = rng.normal(3.0, 0.4, 24)
    low = rng.normal(0.5, 0.4, 24)
    scores = np.concatenate([high, low])
    strata = np.array([1] * 24 + [0] * 24)
    # drivers are enriched among high-burden genes, but labels are independent of score WITHIN a stratum
    is_driver = np.array([True] * 18 + [False] * 6 + [True] * 6 + [False] * 18)
    from scipy.stats import mannwhitneyu
    naive = mannwhitneyu(scores[is_driver], scores[~is_driver], alternative="greater").pvalue
    adj = analysis.stratified_label_permutation_test(scores, is_driver, strata, n_perm=4000, seed=1)
    assert naive < 0.01
    assert adj["status"] == "ok" and adj["one_sided_p_value"] > 0.05


def test_stratified_test_detects_genuine_within_stratum_shift():
    rng = np.random.default_rng(5)
    scores = np.concatenate([rng.normal(3, 0.4, 24), rng.normal(0.5, 0.4, 24)])
    strata = np.array([1] * 24 + [0] * 24)
    is_driver = np.array([True] * 12 + [False] * 12 + [True] * 12 + [False] * 12)
    scores = scores + 1.5 * is_driver
    adj = analysis.stratified_label_permutation_test(scores, is_driver, strata, n_perm=4000, seed=1)
    assert adj["one_sided_p_value"] < 0.01


def test_group_comparison_reports_v11_primary_and_burden_test():
    n = 30
    rng = np.random.default_rng(1)
    result = pd.DataFrame({
        "gene": [f"G{i}" for i in range(n)], "status": "ok", "primary_test_eligible": True,
        "p_context_null": rng.uniform(0.01, 1, n), "q_context_null": 0.5,
        "n_distinct_mutated_residues": rng.integers(1, 40, n),
        "ptm_anchors_with_coordinates": rng.integers(1, 20, n),
        "decoy_z": rng.normal(0, 1, n),
    })
    panel = pd.DataFrame({"gene": result["gene"], "intogen_driver": [i % 2 == 0 for i in range(n)]})
    s = analysis.primary_group_comparison(result, panel, n_strat_perm=500)
    assert s["status"] == "ok"
    assert s["burden_stratified_raw_score_test"]["status"] == "ok"
    v11 = s["primary_v1_1_decoy_calibrated_stratified_test"]
    assert v11["status"] == "ok" and 0 < v11["one_sided_p_value"] <= 1


def test_cli_end_to_end_smoke(tmp_path, monkeypatch):
    """Runs main() with --robustness and --decoys on a tiny hand-built genome (code-path check only)."""
    sys.path.insert(0, str(ROOT / "scripts"))
    bspec = importlib.util.spec_from_file_location("builder_smoke", ROOT / "scripts" / "build_coding_opportunities.py")
    builder = importlib.util.module_from_spec(bspec)
    bspec.loader.exec_module(builder)
    rng = np.random.default_rng(11)
    codons = [c for c, a in builder.CODON_TABLE.items() if a != "*"]
    n_cod = 150
    s_codons = [c for c, a in builder.CODON_TABLE.items() if a == "S"]
    cds = "".join(rng.choice(s_codons) if rng.random() < 0.3 else rng.choice(codons) for _ in range(n_cod))
    flank = lambda k: "".join(rng.choice(list("ACGT"), k))
    genome = flank(100) + cds + flank(100)
    fa = tmp_path / "g.fa"
    fa.write_text(">1\n" + "\n".join(genome[i:i + 60] for i in range(0, len(genome), 60)) + "\n")
    start = 101  # 1-based genomic position of first CDS base
    opp_rows = []
    genes = [f"GENE{i}" for i in range(8)]
    for gene in genes:
        for ci in range(n_cod):
            codon = cds[3 * ci:3 * ci + 3]
            for k in range(3):
                gpos = start + 3 * ci + k
                for alt in "ACGT":
                    if alt == codon[k]:
                        continue
                    new = codon[:k] + alt + codon[k + 1:]
                    aa_new, aa_old = builder.CODON_TABLE[new], builder.CODON_TABLE[codon]
                    if aa_new == "*" or aa_new == aa_old:
                        continue
                    left, right = genome[gpos - 2], genome[gpos]
                    opp_rows.append({"gene": gene, "chromosome": "1", "genomic_position": gpos, "ref": codon[k],
                                     "alt": alt, "protein_pos": ci + 1, "context96": analysis.context96(left, codon[k], alt, right),
                                     "ref_aa": aa_old, "alt_aa": aa_new, "transcript_id": "T1"})
    opp = pd.DataFrame(opp_rows)
    opp_path = tmp_path / "opp.tsv"
    opp.to_csv(opp_path, sep="\t", index=False)
    aa_of = {ci + 1: builder.CODON_TABLE[cds[3 * ci:3 * ci + 3]] for ci in range(n_cod)}
    sty = [p for p, a in aa_of.items() if a in "STY"]
    coord_rows, ptm_rows, maf_rows = [], [], []
    for gi, gene in enumerate(genes):
        for p in range(1, n_cod + 1):
            coord_rows.append({"gene": gene, "position": p, "x": 4.0 * p + gi, "y": (p * 7) % 13, "z": (p * 3) % 5, "plddt": 90})
        for p in [p for p in sty if aa_of[p] == 'S'][:2]:
            ptm_rows.append({"gene": gene, "position": p, "ptm_type": "Phosphorylation"})
        sub = opp[opp["gene"] == gene].sample(25, random_state=gi, replace=False)
        sub = pd.concat([sub, sub.iloc[:5]])   # recurrent variants on purpose
        for r in sub.itertuples():
            maf_rows.append({"Hugo_Symbol": gene, "Chromosome": "1", "Start_Position": r.genomic_position,
                             "Reference_Allele": r.ref, "Tumor_Seq_Allele2": r.alt,
                             "Variant_Classification": "Missense_Mutation", "Tumor_Sample_Barcode": f"S{gi}"})
    pd.DataFrame(coord_rows).to_csv(tmp_path / "coords.tsv", sep="\t", index=False)
    pd.DataFrame(ptm_rows).to_csv(tmp_path / "ptm.tsv", sep="\t", index=False)
    pd.DataFrame(maf_rows).to_csv(tmp_path / "m.maf", sep="\t", index=False)
    pd.DataFrame({"gene": genes, "intogen_driver": [i < 4 for i in range(8)]}).to_csv(tmp_path / "panel.csv", index=False)
    res = tmp_path / "res"
    res.mkdir()
    monkeypatch.setattr(analysis, "RESULTS", res)
    monkeypatch.setattr(sys, "argv", ["x", str(tmp_path / "m.maf"), str(fa), "--cohort", "smoke",
                                      "--opportunities", str(opp_path), "--structure", str(tmp_path / "coords.tsv"),
                                      "--ptm", str(tmp_path / "ptm.tsv"), "--panel", str(tmp_path / "panel.csv"),
                                      "--permutations", "200", "--decoys", "50", "--robustness"])
    assert analysis.main() == 0
    summ = json.loads((res / "smoke_analysis_summary.json").read_text())
    assert summ["primary_settings"]["recurrence_mode"] == "distinct_residue"
    assert summ["genes_tested_primary"] > 0
    prim = pd.read_csv(res / "smoke_ptm_proximity.tsv", sep="\t")
    assert {"q_decoy", "decoy_z", "p_context_null_all_variants"}.issubset(prim.columns)
    assert (prim["decoy_status"] == "ok").any()
    rob = pd.read_csv(res / "smoke_robustness_summary.tsv", sep="\t")
    assert len(rob) == 12
    assert "genes_q_lt_005_family_corrected" in rob.columns
    assert (rob["genes_q_lt_005_family_corrected"] <= rob["genes_q_lt_005_per_setting"]).all()
