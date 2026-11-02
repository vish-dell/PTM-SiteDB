from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "annotate_and_analyze_maf.py"
spec = importlib.util.spec_from_file_location("ptm_analysis", SCRIPT)
analysis = importlib.util.module_from_spec(spec)
assert spec and spec.loader
spec.loader.exec_module(analysis)


def test_sbs96_reverse_complement_is_canonicalized():
    # A[C>A]T is reverse-complement-equivalent to A[G>T]T? Compare generated
    # examples where the latter is the reverse complement of the former.
    assert analysis.context96("A", "C", "A", "G") == analysis.context96("C", "G", "T", "T")


def test_invalid_context_is_rejected():
    assert analysis.context96("N", "C", "A", "G") is None
    assert analysis.context96("A", "C", "C", "G") is None


def test_benjamini_hochberg_monotonic_and_nan_safe():
    p = pd.Series([0.01, 0.04, np.nan, 0.02], index=["a", "b", "c", "d"])
    q = analysis.benjamini_hochberg(p)
    assert np.isnan(q["c"])
    assert abs(q["a"] - 0.03) < 1e-12
    assert abs(q["d"] - 0.03) < 1e-12
    assert abs(q["b"] - 0.04) < 1e-12


def test_sample_null_preserves_context_event_counts_and_is_deterministic():
    counts = np.array([2, 3])
    probs = np.array([0.0, 1.0])
    a = analysis.sample_null(25, counts, probs, np.random.default_rng(17))
    b = analysis.sample_null(25, counts, probs, np.random.default_rng(17))
    assert np.array_equal(a, b)
    assert np.all(a == 3)


def test_primary_gene_with_zero_missense_events_remains_testable_if_qc_passes():
    gene = "GENE1"
    ptm_positions = [1]
    empty_mutations = pd.DataFrame(columns=["gene", "protein_pos", "context96"])
    opportunities = pd.DataFrame({
        "gene": [gene, gene],
        "protein_pos": [1, 10],
        "context96": ["A[C>A]A", "A[C>A]A"],
    })
    coords = {(gene, 1): (0.0, 0.0, 0.0), (gene, 10): (1.0, 0.0, 0.0)}
    result = analysis.setting_gene_result(
        gene, ptm_positions, empty_mutations, opportunities, coords,
        sep=8, dist=10.0, n_perm=100, seed=1, structure_mapping_ok=True,
    )
    assert result["status"] == "ok_no_observed_missense_events"
    assert result["primary_test_eligible"] is True
    assert result["observed_proximal_mutation_events"] == 0
    assert result["p_context_null"] == 1.0


def test_sequence_separation_is_applied_before_structural_distance():
    gene = "GENE2"
    coords = {
        (gene, 10): (0.0, 0.0, 0.0),
        (gene, 15): (0.1, 0.0, 0.0),
        (gene, 30): (0.2, 0.0, 0.0),
    }
    positions = np.array([10, 15, 30])
    mask = analysis.near_mask(positions, np.array([10]), coords, gene, sep=8, dist=10.0)
    assert mask.tolist() == [False, False, True]


def test_coding_opportunity_builder_uses_same_sbs96_canonicalization():
    builder_path = Path(__file__).resolve().parents[1] / "scripts" / "build_coding_opportunities.py"
    builder_spec = importlib.util.spec_from_file_location("ptm_coding_opportunities", builder_path)
    builder = importlib.util.module_from_spec(builder_spec)
    assert builder_spec and builder_spec.loader
    builder_spec.loader.exec_module(builder)
    assert builder.context96("A", "C", "A", "G") == analysis.context96("A", "C", "A", "G")
    assert builder.context96("C", "G", "T", "T") == analysis.context96("C", "G", "T", "T")


def test_ptm_panel_ingestion_accepts_site_labels_like_s127(tmp_path, monkeypatch):
    builder_path = Path(__file__).resolve().parents[1] / "scripts" / "build_200gene_panel.py"
    builder_spec = importlib.util.spec_from_file_location("ptm_panel_builder", builder_path)
    builder = importlib.util.module_from_spec(builder_spec)
    assert builder_spec and builder_spec.loader
    builder_spec.loader.exec_module(builder)
    src = tmp_path / "dbptm.tsv"
    src.write_text(
        "Gene Symbol\tUniProt Accession\tModified Position\tPTM Type\tFlanking Sequence\n"
        "TP53\tP04637\tS127\tPhosphorylation\tABCDE\n",
        encoding="utf-8",
    )
    monkeypatch.setattr(builder, "PTM_SOURCE", src)
    out = builder.load_ptm()
    assert len(out) == 1
    assert out.iloc[0]["gene"] == "TP53"
    assert out.iloc[0]["position"] == 127
    assert out.iloc[0]["uniprot"] == "P04637"


def test_primary_group_comparison_records_driver_background_test():
    result = pd.DataFrame({
        "gene": ["A", "B", "C", "D"],
        "status": ["ok", "ok", "ok", "ok_no_observed_missense_events"],
        "primary_test_eligible": [True, True, True, True],
        "p_context_null": [0.01, 0.02, 0.5, 1.0],
        "q_context_null": [0.04, 0.04, 0.5, 1.0],
    })
    panel = pd.DataFrame({"gene": ["A", "B", "C", "D"], "intogen_driver": [True, True, False, False]})
    summary = analysis.primary_group_comparison(result, panel)
    assert summary["status"] == "ok"
    assert summary["driver_genes_tested"] == 2
    assert summary["not_listed_genes_tested"] == 2
    assert summary["one_sided_p_value"] is not None
    assert summary["secondary_driver_q_lt_005_hits"] == 2
    assert summary["secondary_not_listed_q_lt_005_hits"] == 0
