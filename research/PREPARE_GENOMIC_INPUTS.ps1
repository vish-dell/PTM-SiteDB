param(
    [Parameter(Mandatory=$true)][string]$ReferenceFasta,
    [Parameter(Mandatory=$true)][string]$Gtf,
    [int]$MaxGenes = 0
)
$ErrorActionPreference = "Stop"

# Resolve user-provided paths before changing the current directory.
$callerDirectory = (Get-Location).Path
if (![System.IO.Path]::IsPathRooted($ReferenceFasta)) {
    $ReferenceFasta = [System.IO.Path]::GetFullPath((Join-Path $callerDirectory $ReferenceFasta))
}
if (![System.IO.Path]::IsPathRooted($Gtf)) {
    $Gtf = [System.IO.Path]::GetFullPath((Join-Path $callerDirectory $Gtf))
}
$projectRoot = Split-Path -Parent $PSScriptRoot
$pythonPath = Join-Path $projectRoot ".venv312\Scripts\python.exe"
if (!(Test-Path $pythonPath)) {
    throw "Environment not found at '$pythonPath'. Run .\INSTALL_AND_RUN.ps1 from the top-level project folder first."
}
$py = (Resolve-Path $pythonPath).Path
if (!(Test-Path $ReferenceFasta)) { throw "Reference FASTA not found: $ReferenceFasta" }
if (!(Test-Path $Gtf)) { throw "GTF/GTF.GZ not found: $Gtf" }
$ref = (Resolve-Path $ReferenceFasta).Path
$gtfPath = (Resolve-Path $Gtf).Path
Set-Location -LiteralPath $PSScriptRoot

# Build the benchmark panel first (human dbPTM sites x IntOGen drivers); it needs internet once for UniProt gene symbols.
if (!(Test-Path ".\data\processed\benchmark_panel.csv")) {
    Write-Host "Building the benchmark panel from the real dbPTM + IntOGen sources (internet needed once for UniProt gene symbols)..." -ForegroundColor Cyan
    & $py scripts\build_200gene_panel.py
    if ($LASTEXITCODE -ne 0 -or !(Test-Path ".\data\processed\benchmark_panel.csv")) {
        throw "Panel build failed. Check your internet connection and that data\processed\dbptm_phosphorylation_source.txt and intogen_drivers_source.tsv exist."
    }
}

if (!(Test-Path ($ref + ".fai"))) {
    Write-Host "Creating FASTA index (.fai). A whole genome takes a few minutes; progress is shown - do NOT press Ctrl+C." -ForegroundColor Cyan
    & $py -c "from scripts.fasta_utils import IndexedFasta; import sys; f=IndexedFasta(sys.argv[1]); f.close()" $ref
    if ($LASTEXITCODE -ne 0) { throw "Could not index FASTA. Check FASTA file formatting and write permissions." }
}

Write-Host "`n[1/3] Building SBS96 coding missense opportunity map..." -ForegroundColor Cyan
$cmd = @("scripts\build_coding_opportunities.py", "--gtf", $gtfPath, "--fasta", $ref)
if ($MaxGenes -gt 0) { $cmd += @("--max-genes", [string]$MaxGenes) }
& $py @cmd
if ($LASTEXITCODE -ne 0) { throw "Opportunity map generation failed. Confirm the GTF and FASTA use the same genome assembly." }

Write-Host "`n[2/3] Fetching AlphaFold models and mapping coordinates/PTM anchors..." -ForegroundColor Cyan
$fetchArgs = @("scripts\fetch_alphafold_coords.py")
if ($MaxGenes -gt 0) { $fetchArgs += @("--max-genes", [string]$MaxGenes) }
& $py @fetchArgs
if ($LASTEXITCODE -ne 0) { throw "AlphaFold coordinate preparation failed. Inspect the log above and mapping manifest." }

Write-Host "`n[3/3] Recording exact reference path and preparation QC..." -ForegroundColor Cyan
& $py -c "import json,sys; from pathlib import Path; p=Path('data/processed/reference_config.json'); p.write_text(json.dumps({'reference_fasta':str(Path(sys.argv[1]).resolve()),'gtf':str(Path(sys.argv[2]).resolve()),'genome_assembly_note':'User-declared matched FASTA/GTF; verify assembly/build manually'},indent=2),encoding='utf-8')" $ref $gtfPath
if ($LASTEXITCODE -ne 0) { throw "Could not record reference configuration." }
& $py scripts\preflight_check.py
if ($LASTEXITCODE -ne 0) { throw "Preflight reported a fatal missing-input error. Review the result." }
Write-Host "`nGenomic input preparation finished." -ForegroundColor Green
Write-Host "Review data\processed\opportunity_manifest.json, structure_manifest.json, and preparation_qc.json before interpreting cohort results." -ForegroundColor Yellow
Write-Host "AlphaFold coordinates are predictions, not experimentally solved structures." -ForegroundColor Yellow
