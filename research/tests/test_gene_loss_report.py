"""Gene-loss accounting. The GTF/FASTA below are hand-built UNIT-TEST FIXTURES, not biological data."""
import json
import sys

import pandas as pd

import scripts.build_coding_opportunities as bco
from scripts.fasta_utils import IndexedFasta

CDS = "ATGGCTAAACCCGGGAAA"                         # M A K P G K (18 nt)
SEQ = "ACGTACGTAC" + CDS + "TGCATGCATGCA" * 5


def _gtf(lines):
    return "\n".join(lines) + "\n"


def test_every_panel_gene_is_accounted_for(tmp_path, monkeypatch, capsys):
    fa = tmp_path / "g.fa"
    fa.write_text(">chr1\n" + "\n".join(SEQ[i:i + 60] for i in range(0, len(SEQ), 60)) + "\n")
    IndexedFasta(fa).close()
    attrs = 'gene_id "{g}"; transcript_id "{t}"; gene_name "{n}"; transcript_type "{ty}";'
    gtf = tmp_path / "g.gtf"
    gtf.write_text(_gtf([
        "chr1\tH\tCDS\t11\t28\t.\t+\t0\t" + attrs.format(g="G1", t="T1", n="GENEA", ty="protein_coding"),
        "chr1\tH\texon\t11\t28\t.\t+\t.\t" + attrs.format(g="G2", t="T2", n="GENEB", ty="lncRNA"),          # present, not coding
        "chr1\tH\tCDS\t11\t27\t.\t+\t0\t" + attrs.format(g="G4", t="T4", n="GENED", ty="protein_coding"),   # 17 nt, not a codon multiple
    ]))
    panel = tmp_path / "panel.csv"
    pd.DataFrame({"gene": ["GENEA", "GENEB", "GENEC", "GENED"], "intogen_driver": [True, True, False, False]}).to_csv(panel, index=False)
    out = tmp_path / "opp.tsv"
    monkeypatch.setattr(sys, "argv", ["x", "--gtf", str(gtf), "--fasta", str(fa), "--panel", str(panel), "--out", str(out),
                                      "--proteins-out", str(tmp_path / "p.fa"), "--manifest-out", str(tmp_path / "m.json")])
    assert bco.main() == 0

    report = pd.read_csv(tmp_path / "gene_loss_report.tsv", sep="\t").set_index("gene")["reason"].to_dict()
    assert report == {"GENEB": "no_protein_coding_transcript_with_cds",
                      "GENEC": "gene_symbol_not_found_in_gtf",
                      "GENED": "non_triplet_complete_cds"}
    assert set(pd.read_csv(out, sep="\t")["gene"]) == {"GENEA"}
    acc = json.loads((tmp_path / "m.json").read_text())["gene_accounting"]
    assert (acc["panel_genes_requested"], acc["panel_genes_retained"], acc["panel_genes_lost"]) == (4, 1, 3)
    assert acc["driver_genes_retained"] == "1/2" and acc["not_listed_genes_retained"] == "0/2"
    text = capsys.readouterr().out
    assert "Gene accounting: 1/4 panel genes retained" in text and "WARNING: more than 20%" in text


def test_report_written_even_when_nothing_is_selected(tmp_path, monkeypatch):
    fa = tmp_path / "g.fa"
    fa.write_text(">chr1\n" + SEQ + "\n")
    IndexedFasta(fa).close()
    gtf = tmp_path / "g.gtf"
    gtf.write_text('chr1\tH\tCDS\t11\t28\t.\t+\t0\tgene_id "G1"; transcript_id "T1"; gene_name "OTHER"; transcript_type "protein_coding";\n')
    panel = tmp_path / "panel.csv"
    pd.DataFrame({"gene": ["GENEA"], "intogen_driver": [True]}).to_csv(panel, index=False)
    monkeypatch.setattr(sys, "argv", ["x", "--gtf", str(gtf), "--fasta", str(fa), "--panel", str(panel), "--out", str(tmp_path / "o.tsv"),
                                      "--proteins-out", str(tmp_path / "p.fa"), "--manifest-out", str(tmp_path / "m.json")])
    try:
        bco.main()
    except SystemExit as exc:
        assert "gene_loss_report" in str(exc) or "per-gene reasons" in str(exc)
    else:
        raise AssertionError("expected SystemExit")
    assert pd.read_csv(tmp_path / "gene_loss_report.tsv", sep="\t")["reason"].tolist() == ["gene_symbol_not_found_in_gtf"]
