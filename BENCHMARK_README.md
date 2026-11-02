# Independent benchmark workflow

This project uses the **IntOGen driver catalog release 2024.09.20** as an independent gene-level reference. The official IntOGen download page identifies the release date and provides the driver catalog; the benchmark script never uses the project's own hit list to define positives.

## One command on Windows

From the project folder:

```powershell
.\run_benchmark.ps1
```

The script will:
1. verify the BRCA result exists;
2. download the official IntOGen 2024.09.20 catalog automatically;
3. build `benchmark.csv` from the genes actually present in the BRCA result;
4. run `benchmark_report.py` using permutation q-values at q <= 0.05;
5. write `results/benchmark_brca.csv` and `results/benchmark_brca.summary.csv`.

## If the IntOGen website is unreachable

Download `IntOGen-Drivers-20240920.zip` from the official IntOGen download page and then run:

```powershell
python build_intogen_benchmark.py results/proximity_brca_tcga_sep8_dist10.csv "C:\full\path\IntOGen-Drivers-20240920.zip" --out benchmark.csv
python benchmark_report.py results/proximity_brca_tcga_sep8_dist10.csv benchmark.csv --out results/benchmark_brca.csv --metric perm --q-cutoff 0.05
```

## Interpretation

`positive` means the gene is listed as a driver in IntOGen 2024.09.20. `negative` means only that the gene is **not listed** in that release; it does not prove that the gene is biologically non-driver.

Rows for which the project has no usable q-value are reported as `not_comparable` and are never counted as negatives.

## Reproducibility

Record the IntOGen release (2024.09.20), benchmark metric (`perm`), q-value cutoff (0.05), and the exact project result filename in the report.
