$ErrorActionPreference = "Stop"
Set-Location -LiteralPath $PSScriptRoot

Write-Host "PTM-SiteDB — original interface + Research extension" -ForegroundColor Cyan
Write-Host "Original styles and Analysis Suite are preserved." -ForegroundColor Cyan

# --- 1. Folder integrity check -------------------------------------------------------------
# An interrupted or over-long-path extraction (Windows limit: 260 characters) silently drops folders.
$required = @("app.py", "requirements.txt", "templates\base.html", "templates\protein.html",
              "services\proximity.py", "static\js\viewer.js",
              "research\app.py", "research\templates\home.html", "research\scripts\annotate_and_analyze_maf.py")
$missing = @($required | Where-Object { -not (Test-Path -LiteralPath (Join-Path $PSScriptRoot $_)) })
if ($missing.Count -gt 0) {
    Write-Host "`nThis copy of the project is INCOMPLETE. Missing:" -ForegroundColor Red
    $missing | ForEach-Object { Write-Host "  $_" -ForegroundColor Red }
    Write-Host "Extract the zip again into a SHORT folder, e.g.  tar -xf <zipfile> -C C:\  (or use 7-Zip)." -ForegroundColor Yellow
    throw "Incomplete extraction - see the list above."
}
if ($PSScriptRoot.Length -gt 80) {
    Write-Warning "This folder path is long ($($PSScriptRoot.Length) characters). Windows can fail on paths over 260 characters. A short path such as C:\ptm is safer."
}

# --- 2. Python virtual environment (created, or repaired if broken) ---------------------------
function Test-Venv {
    if (!(Test-Path ".\.venv312\Scripts\python.exe")) { return $false }
    try { & ".\.venv312\Scripts\python.exe" -m pip --version *> $null; return ($LASTEXITCODE -eq 0) } catch { return $false }
}
if (-not (Test-Venv)) {
    if (Test-Path ".\.venv312") {
        $stamp = Get-Date -Format "yyyyMMddHHmmss"
        Write-Warning "The existing .venv312 is broken (it cannot run pip). Renaming it to .venv312_broken_$stamp and creating a fresh one."
        Rename-Item ".\.venv312" ".venv312_broken_$stamp"
    }
    $launchers = @({ py -3.12 -m venv .venv312 }, { py -3 -m venv .venv312 }, { python -m venv .venv312 })
    foreach ($launcher in $launchers) {
        try { & $launcher } catch { }
        if (Test-Venv) { break }
        if (Test-Path ".\.venv312") { Remove-Item ".\.venv312" -Recurse -Force -ErrorAction SilentlyContinue }
    }
    if (-not (Test-Venv)) {
        throw "Could not create a Python virtual environment. Install Python 3.12 from python.org (tick 'Add python.exe to PATH') and run this script again."
    }
}
$py = (Resolve-Path ".\.venv312\Scripts\python.exe").Path

Write-Host "`n=== INSTALL DEPENDENCIES ===" -ForegroundColor Cyan
& $py -m pip install --upgrade pip
if ($LASTEXITCODE -ne 0) { throw "pip upgrade failed." }
& $py -m pip install -r requirements.txt
if ($LASTEXITCODE -ne 0) { throw "Dependency installation failed. Check the pip error above." }

Write-Host "`n=== DATABASE ===" -ForegroundColor Cyan
if (Test-Path .\instance\ptm.db) {
    Write-Host "Using the supplied curated database instance\ptm.db. Skipping sync_data.py so the cleaned data is not overwritten." -ForegroundColor Green
} else {
    throw "instance\ptm.db is missing. Re-extract the submitted zip. Do not run sync_data.py: it reloads uncleaned records."
}

Write-Host "`n=== PREPARE EXPANDED RESEARCH SOURCES ===" -ForegroundColor Cyan
& $py research\scripts\bootstrap_sources.py
if ($LASTEXITCODE -ne 0) { Write-Warning "Research-source download failed; check source manifest and internet access." }

$ptm = ".\research\data\processed\dbptm_phosphorylation_source.txt"
$intogen = ".\research\data\processed\intogen_drivers_source.tsv"
if ((Test-Path $ptm) -and (Test-Path $intogen)) {
    & $py research\scripts\build_200gene_panel.py
    if ($LASTEXITCODE -ne 0) { Write-Warning "Research panel build failed; inspect the schema/QC messages." }
} else {
    Write-Warning "Expanded research data are incomplete. Counts will remain unavailable until real sources load."
}

Write-Host "`nStarting website at http://127.0.0.1:5000/" -ForegroundColor Green
Write-Host "Connection Check: http://127.0.0.1:5000/connection-check" -ForegroundColor Green
Write-Host "Research extension: http://127.0.0.1:5000/research/" -ForegroundColor Green
& $py app.py
