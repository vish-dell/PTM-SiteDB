import hashlib
import sys
import pandas as pd

PAAD = r"C:\ptm\research\instance\uploads\TCGA_PAAD.maf"
PRAD = r"C:\ptm\research\instance\uploads\TCGA_PRAD.maf"
KEY = ["Tumor_Sample_Barcode", "Chromosome", "Start_Position",
       "Reference_Allele", "Tumor_Seq_Allele2"]


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def load(path):
    m = pd.read_csv(path, sep="\t", comment="#", low_memory=False)
    m["patient"] = m["Tumor_Sample_Barcode"].str[:12]
    return m


print("Checking file identity (read-only)...")
h1, h2 = sha256(PAAD), sha256(PRAD)
print("PAAD sha256:", h1)
print("PRAD sha256:", h2)
print("IDENTICAL FILES" if h1 == h2 else "Files differ")
print()

patients = {}
for name, path in (("PAAD", PAAD), ("PRAD", PRAD)):
    m = load(path)
    patients[name] = set(m["patient"])
    print(f"== {name} ==")
    print("rows:", len(m))
    print("exact duplicate variant rows:", int(m.duplicated(KEY).sum()))
    print("tumor samples:", m["Tumor_Sample_Barcode"].nunique(),
          "| patients:", m["patient"].nunique())
    print("example barcodes:", list(m["Tumor_Sample_Barcode"].drop_duplicates().head(3)))
    print()

shared = patients["PAAD"] & patients["PRAD"]
print("patients appearing in BOTH files:", len(shared))
print("(a clean pair of cohorts should share 0)")
