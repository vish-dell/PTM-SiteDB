$ErrorActionPreference = "Stop"
Set-Location $PSScriptRoot

if (!(Test-Path ".\.venv312\Scripts\python.exe")) {
    py -3.12 -m venv .venv312
}

& .\.venv312\Scripts\python.exe -m pip install --upgrade pip
& .\.venv312\Scripts\python.exe -m pip install -r requirements.txt
& .\.venv312\Scripts\python.exe scripts\bootstrap_sources.py
& .\.venv312\Scripts\python.exe scripts\build_200gene_panel.py
Write-Host "`nBootstrap complete. Start the website with .\start.bat" -ForegroundColor Green
