import sys, pandas as pd
a, b = sys.argv[1], sys.argv[2]
sa = set(pd.read_csv(a, sep="\t", usecols=["Tumor_Sample_Barcode"], comment="#")["Tumor_Sample_Barcode"])
sb = set(pd.read_csv(b, sep="\t", usecols=["Tumor_Sample_Barcode"], comment="#")["Tumor_Sample_Barcode"])
ov = len(sa & sb)
print(f"A: {len(sa)} samples | B: {len(sb)} samples | overlap: {ov}")
print("OK: independent" if ov == 0 else "WARNING: cohorts share samples, not independent")
