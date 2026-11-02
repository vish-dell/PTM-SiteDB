"""Regression: FASTA indexing from the web app failed with "No module named 'scripts'" when the site was
started from the project root (research/ was not on sys.path)."""
import subprocess, sys, textwrap
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_scripts_package_importable_after_mounting_research_app(tmp_path):
    fa = tmp_path / "t.fa"
    fa.write_text(">chr1\nACGTACGTAC\nGTACGTACGT\n")
    code = textwrap.dedent(f"""
        import sys
        sys.path[:] = [p for p in sys.path if 'research' not in p]       # as when launched from the project root
        import app                                                       # the real startup entry point
        from scripts.fasta_utils import IndexedFasta                     # the exact import the upload handler uses
        with IndexedFasta({str(fa)!r}):
            pass
        print("OK")
    """)
    r = subprocess.run([sys.executable, "-c", code], cwd=ROOT, capture_output=True, text=True, timeout=120)
    assert r.returncode == 0 and "OK" in r.stdout, r.stdout + r.stderr
    assert (tmp_path / "t.fa.fai").exists()          # the index was really built
