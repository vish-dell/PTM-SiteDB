from pathlib import Path

from scripts.fasta_utils import IndexedFasta


def test_fasta_index_build_and_fetch_across_wrapped_lines(tmp_path):
    fasta = tmp_path / "tiny.fa"
    # CRLF handling matters for Windows-created FASTA files.
    fasta.write_bytes(b">chr1 chromosome one\r\nACGT\r\nTGCA\r\nAAAA\r\n>chr2\r\nttac\r\nccgg\r\n")
    with IndexedFasta(fasta) as ref:
        assert ref.references == ["chr1", "chr2"]
        assert ref.fetch("chr1", 0, 4) == "ACGT"
        assert ref.fetch("chr1", 2, 9) == "GTTGCAA"
        assert ref.fetch("chr2", 1, 7) == "TACCCG"
    assert Path(str(fasta) + ".fai").exists()


def test_fasta_bounds_and_unknown_contig_raise(tmp_path):
    fasta = tmp_path / "tiny.fa"
    fasta.write_text(">chr1\nACGTAC\n", encoding="ascii")
    with IndexedFasta(fasta) as ref:
        try:
            ref.fetch("chrX", 0, 1)
        except KeyError:
            pass
        else:
            raise AssertionError("Unknown contig must fail")
        try:
            ref.fetch("chr1", 0, 7)
        except ValueError:
            pass
        else:
            raise AssertionError("Out-of-range FASTA slice must fail")


import pytest


def test_index_rejects_inconsistent_line_lengths(tmp_path):
    fasta = tmp_path / "bad.fa"
    fasta.write_text(">chr1\nACGT\nAC\nACGT\n", encoding="ascii")        # short line in the middle
    with pytest.raises(ValueError, match="inconsistent"):
        IndexedFasta(fasta)


def test_index_rejects_duplicate_ids_and_orphan_sequence(tmp_path):
    dup = tmp_path / "dup.fa"
    dup.write_text(">chr1\nACGT\n>chr1\nACGT\n", encoding="ascii")
    with pytest.raises(ValueError, match="Duplicate"):
        IndexedFasta(dup)
    orphan = tmp_path / "orphan.fa"
    orphan.write_text("ACGT\n>chr1\nACGT\n", encoding="ascii")
    with pytest.raises(ValueError, match="before the first"):
        IndexedFasta(orphan)


def test_index_handles_no_trailing_newline_blank_lines_and_matches_samtools_columns(tmp_path):
    fasta = tmp_path / "ok.fa"
    fasta.write_bytes(b">chr1 desc\nACGTAC\nGTACGT\nAC")                 # last line short, no final newline
    with IndexedFasta(fasta) as ref:
        assert ref.fetch("chr1", 0, 14) == "ACGTACGTACGTAC"
    # standard samtools .fai: name, length, offset, bases per line, bytes per line
    assert Path(str(fasta) + ".fai").read_text().strip() == "chr1\t14\t11\t6\t7"
