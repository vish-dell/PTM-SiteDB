"""Guards against the same mutation file being analysed under two different cohort labels.

Why: two labelled cohorts whose MAF files have an identical SHA-256 are one dataset, not two. Any
"replication" between them is meaningless (they differ only by permutation noise). This was observed in
the shipped cohort_results (TCGA_PAAD and TCGA_PRAD). These helpers make that impossible to miss.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path


def file_sha256(path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def recorded_maf_hashes(results_dir) -> dict:
    """{cohort: maf_sha256} from existing <cohort>_analysis_summary.json provenance blocks."""
    out = {}
    for f in sorted(Path(results_dir).glob("*_analysis_summary.json")):
        try:
            s = json.loads(f.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            continue
        h = (s.get("input_provenance") or {}).get("maf_sha256")
        if h:
            out[s.get("cohort") or f.name.replace("_analysis_summary.json", "")] = h
    return out


def duplicate_groups(hashes: dict) -> list:
    """Groups of cohort labels sharing one MAF hash (only groups with >= 2 labels)."""
    by = {}
    for cohort, h in hashes.items():
        by.setdefault(h, []).append(cohort)
    return [sorted(v) for v in by.values() if len(v) > 1]


def find_duplicate_input(maf_path, cohort: str, results_dir):
    """Cohort label (other than `cohort`) already analysed with a byte-identical MAF, else None."""
    h = file_sha256(maf_path)
    for other, oh in recorded_maf_hashes(results_dir).items():
        if oh == h and other != cohort:
            return other
    return None
