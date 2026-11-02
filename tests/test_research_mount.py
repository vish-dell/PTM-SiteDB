"""The research app is mounted under /research of the root app; uploads must work from that launch mode."""
import io
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def test_fasta_upload_indexes_when_mounted_under_root():
    import app as root_app
    client = root_app.app.test_client()
    upload = ROOT / "research" / "instance" / "uploads" / "MOUNTTEST_reference.fa"
    try:
        r = client.post("/research/upload-analysis", content_type="multipart/form-data", data={
            "cohort": "MOUNTTEST",
            "maf": (io.BytesIO(b"Hugo_Symbol\tChromosome\n"), "t.maf"),
            "fasta": (io.BytesIO(b">1\nACGTACGTACGTACGT\n"), "t.fa")})
        text = r.get_data(as_text=True)
        assert r.status_code == 200
        assert "No module named" not in text and "Could not index" not in text
    finally:
        for p in (upload, Path(str(upload) + ".fai"), upload.with_name("MOUNTTEST.maf")):
            if p.exists():
                p.unlink()
