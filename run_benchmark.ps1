$ErrorActionPreference = 'Stop'
$project = 'results/proximity_brca_tcga_sep8_dist10.csv'
if (-not (Test-Path $project)) { throw "Missing $project. Run batch_proximity.py first." }
python build_intogen_benchmark.py $project --out benchmark.csv
python benchmark_report.py $project benchmark.csv --out results/benchmark_brca.csv --metric perm --q-cutoff 0.05
Write-Host "`nBenchmark complete. Inspect:" -ForegroundColor Green
Write-Host "results/benchmark_brca.summary.csv"
