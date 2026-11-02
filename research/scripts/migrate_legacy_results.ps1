param(
  [Parameter(Mandatory=$true)] [string]$OldProject
)

$ErrorActionPreference = "Stop"
$new = Split-Path -Parent (Split-Path -Parent $PSScriptRoot)
$src = Join-Path $OldProject "results"
$dst = Join-Path $new "data\processed\legacy_results"

if (!(Test-Path $src)) { throw "Old results folder not found: $src" }
New-Item -ItemType Directory -Force -Path $dst | Out-Null

Get-ChildItem $src -File | Where-Object {
    $_.Name -match '^(proximity_|sites_|independent_|hotspot_)'
} | ForEach-Object {
    Copy-Item $_.FullName (Join-Path $dst $_.Name) -Force
}

Write-Host "Legacy result files copied to $dst" -ForegroundColor Green
