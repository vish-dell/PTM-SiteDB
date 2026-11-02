# Runs the pre-specified v1.1 analysis on TWO independent same-cancer cohorts and builds RESULTS.md.
# You must supply real somatic MAFs (same genome build as the FASTA) and check sample independence yourself.
param(
  [Parameter(Mandatory=$true)][string]$MafDiscovery,
  [Parameter(Mandatory=$true)][string]$MafReplication,
  [Parameter(Mandatory=$true)][string]$Fasta,
  [string]$CohortDiscovery = "discovery",
  [string]$CohortReplication = "replication"
)
$ErrorActionPreference = "Stop"
$py = Join-Path (Split-Path $PSScriptRoot -Parent) ".venv312\Scripts\python.exe"
if (-not (Test-Path $py)) { $py = "python" }
foreach ($pair in @(@($MafDiscovery,$CohortDiscovery), @($MafReplication,$CohortReplication))) {
  & $py (Join-Path $PSScriptRoot "scripts\annotate_and_analyze_maf.py") $pair[0] $Fasta `
      --cohort $pair[1] --permutations 5000 --decoys 1000 --recurrence distinct_residue --robustness
  if ($LASTEXITCODE -ne 0) { throw "Analysis failed for cohort $($pair[1])" }
}
& $py (Join-Path $PSScriptRoot "scripts\make_results_report.py") $CohortDiscovery $CohortReplication
Write-Host "Done. Open research\RESULTS.md. Report the primary result even if it is null."
