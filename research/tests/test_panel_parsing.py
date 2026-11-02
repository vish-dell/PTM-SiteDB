"""Panel-builder tests. All rows below are hand-written UNIT-TEST FIXTURES in the real dbPTM file layout
(ENTRY_SPECIES, accession, position, type, PubMed IDs, flank); they are not biological data."""
import pandas as pd
import pytest

import scripts.build_200gene_panel as bp

REAL_LAYOUT = (
    "AAA1_HUMAN\tP00001\t10\tPhosphorylation\t111;222\tAAAAASAAAAAAAAAAAAAAA\n"
    "AAA1_HUMAN\tP00001\t55\tPhosphorylation\t111\tAAAAATAAAAAAAAAAAAAAA\n"
    "BBB2_HUMAN\tP00002\t7\tPhosphorylation\t333\tAAAAAYAAAAAAAAAAAAAAA\n"
    "CCC3_HUMAN\tP00003\t99\tPhosphorylation\t444\tAAAAASAAAAAAAAAAAAAAA\n"
    "DDD4_HUMAN\tP00004\t5\tPhosphorylation\t555\tAAAAASAAAAAAAAAAAAAAA\n"
    "AAA1_MOUSE\tQ00001\t10\tPhosphorylation\t666\tAAAAASAAAAAAAAAAAAAAA\n"
    "ZZZ9_YEAST\tQ00009\t3\tPhosphorylation\t777\tAAAAASAAAAAAAAAAAAAAA\n"
)


@pytest.fixture()
def paths(tmp_path, monkeypatch):
    (tmp_path / "dbptm.txt").write_text(REAL_LAYOUT, encoding="utf-8")
    (tmp_path / "intogen.tsv").write_text("SYMBOL\tTRANSCRIPT\nGENEA\tENST1\nGENEB\tENST2\n", encoding="utf-8")
    for name, value in {"PROC": tmp_path, "PTM_SOURCE": tmp_path / "dbptm.txt", "INTOGEN_SOURCE": tmp_path / "intogen.tsv",
                        "OUT_PTM": tmp_path / "ptm_sites.tsv", "OUT_DRIVERS": tmp_path / "intogen_drivers.tsv",
                        "OUT_PANEL": tmp_path / "benchmark_panel.csv", "GENE_MAP_CACHE": tmp_path / "gene_map.tsv"}.items():
        monkeypatch.setattr(bp, name, value)
    monkeypatch.setattr(bp, "resolve_gene_symbols",
                        lambda accs: {"P00001": "GENEA", "P00002": "GENEB", "P00003": "GENEC", "P00004": "GENED"})
    return tmp_path


def test_real_dbptm_layout_is_parsed_and_filtered_to_human(paths):
    df = bp.load_ptm()
    assert set(df["uniprot"]) == {"P00001", "P00002", "P00003", "P00004"}      # mouse/yeast rows excluded
    assert set(df["gene"]) == {"GENEA", "GENEB", "GENEC", "GENED"}
    assert df["position"].dtype.kind == "i" and sorted(df["position"]) == [5, 7, 10, 55, 99]
    assert (df["ptm_type"] == "Phosphorylation").all()


def test_panel_has_driver_and_not_listed_groups(paths, capsys):
    assert bp.main() == 0
    panel = pd.read_csv(paths / "benchmark_panel.csv")
    assert set(panel.loc[panel["intogen_driver"], "gene"]) == {"GENEA", "GENEB"}
    assert set(panel.loc[~panel["intogen_driver"], "gene"]) == {"GENEC", "GENED"}
    assert "do not call this a 200-gene panel" in capsys.readouterr().out      # honest about the small fixture


def test_gene_cache_avoids_network_and_unmapped_are_excluded(tmp_path, monkeypatch):
    monkeypatch.setattr(bp, "GENE_MAP_CACHE", tmp_path / "gene_map.tsv")
    pd.DataFrame({"accession": ["P1", "P2"], "gene": ["G1", "G2"]}).to_csv(bp.GENE_MAP_CACHE, sep="\t", index=False)
    monkeypatch.setattr(bp, "fetch_reviewed_human_gene_map", lambda: pytest.fail("network must not be used"))
    monkeypatch.setattr(bp, "map_uniprot_to_gene", lambda a: pytest.fail("network must not be used"))
    assert bp.resolve_gene_symbols(["P1", "P2"]) == {"P1": "G1", "P2": "G2"}


def test_bulk_lookup_used_for_many_accessions_and_missing_are_not_guessed(tmp_path, monkeypatch):
    monkeypatch.setattr(bp, "GENE_MAP_CACHE", tmp_path / "gene_map.tsv")
    accs = [f"P{i:05d}" for i in range(600)]
    monkeypatch.setattr(bp, "fetch_reviewed_human_gene_map", lambda: {a: f"G{a[1:]}" for a in accs[:550]})
    monkeypatch.setattr(bp, "map_uniprot_to_gene", lambda a: {})
    out = bp.resolve_gene_symbols(accs)
    assert len(out) == 550 and accs[599] not in out
    assert len(pd.read_csv(bp.GENE_MAP_CACHE, sep="\t")) == 550              # cached for the next run
