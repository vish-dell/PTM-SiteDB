# Run commands (Windows PowerShell, project in C:\ptm)

## 1. One-time setup + sync + start (simplest)
```powershell
cd C:\ptm
.\INSTALL_AND_RUN.ps1
```
That script creates `.venv312`, installs `requirements.txt`, runs `sync_data.py`, downloads/builds the research panel, then starts the site.

## 2. Same steps manually
```powershell
cd C:\ptm
py -3.12 -m venv .venv312
.\.venv312\Scripts\python.exe -m pip install -r requirements.txt -r research\requirements-dev.txt
.\.venv312\Scripts\python.exe diagnose_connections.py            # checks cBioPortal/AlphaFold/RCSB reachability first
.\.venv312\Scripts\python.exe sync_data.py                       # SYNC the PTM-SiteDB records into instance\ptm.db
.\.venv312\Scripts\python.exe research\scripts\bootstrap_sources.py
.\.venv312\Scripts\python.exe research\scripts\build_200gene_panel.py
.\.venv312\Scripts\python.exe app.py                             # RUN the website
```
Open http://127.0.0.1:5000/ (main), http://127.0.0.1:5000/research/ (research), http://127.0.0.1:5000/connection-check

## 3. Tests and simulation (return the output files listed in section 6)
```powershell
.\.venv312\Scripts\python.exe -m pytest tests research\tests -v -rs -p no:cacheprovider *> validation\logs\pytest_windows.txt
.\.venv312\Scripts\python.exe validation\simulation\run_simulation.py *> validation\logs\simulation_windows.txt
.\.venv312\Scripts\python.exe validation\audit_cohort_inputs.py *> validation\logs\cohort_input_audit_windows.txt
```
Optional browser test: `npm install jsdom` then re-run pytest (the skipped test).

## 4. Correct cohorts (PAAD, BRCA, LUAD) — you must supply a DIFFERENT real MAF per cohort
Download each cohort's somatic MAF (e.g. GDC/TCGA open-access MAF, build GRCh37 to match `hg19.fa`), record the source URL/version, and confirm sample counts:
```powershell
$py = ".\.venv312\Scripts\python.exe"; $fa = "C:\ptm-data\hg19.fa"
foreach ($c in @(@("TCGA_PAAD","C:\ptm-data\PAAD.maf"), @("TCGA_BRCA","C:\ptm-data\BRCA.maf"), @("TCGA_LUAD","C:\ptm-data\LUAD.maf"))) {
  & $py research\scripts\annotate_and_analyze_maf.py $c[1] $fa --cohort $c[0] --permutations 5000 --decoys 1000 --recurrence distinct_residue --robustness *> "validation\logs\cohort_$($c[0]).txt"
}
& $py research\scripts\make_results_report.py TCGA_PAAD TCGA_BRCA TCGA_LUAD
& $py validation\audit_cohort_inputs.py
```
The report script now refuses to run if two cohorts share a byte-identical MAF. Delete the old invalid `TCGA_PAAD_*`/`TCGA_PRAD_*` outputs first (move them to a backup folder).

## 5. Plagiarism review
Upload the similarity report (PDF/HTML) and your document; I will complete the table in FINAL_REPORT.md section F.

## 6. Files to send back
validation\logs\pytest_windows.txt, simulation_windows.txt, cohort_input_audit_windows.txt, cohort_TCGA_*.txt, research\data\processed\cohort_results\*_analysis_summary.json, RESULTS.md, and the similarity report.
