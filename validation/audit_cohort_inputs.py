"""Audit cohort inputs for independence. Usage: python validation/audit_cohort_inputs.py
Reports MAF SHA-256 per cohort, duplicate groups, row/sample counts and TCGA tissue-source-site prefixes."""
import sys, collections
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT / "research"))
import pandas as pd
from scripts import cohort_guard as G
RES = ROOT / "research/data/processed/cohort_results"; UP = ROOT / "research/instance/uploads"
h = G.recorded_maf_hashes(RES)
print("Recorded MAF hashes:"); [print(f"  {c}: {x}") for c, x in h.items()]
print("Duplicate groups (NOT independent cohorts):", G.duplicate_groups(h) or "none")
for maf in sorted(UP.glob("*.maf")):
    d = pd.read_csv(maf, sep="\t", comment="#", low_memory=False, usecols=["Tumor_Sample_Barcode"])
    pre = collections.Counter(d["Tumor_Sample_Barcode"].str.split("-").str[1])
    print(f"{maf.name}: rows={len(d)} samples={d['Tumor_Sample_Barcode'].nunique()} sha256={G.file_sha256(maf)[:16]}… top TSS codes={pre.most_common(5)}")
