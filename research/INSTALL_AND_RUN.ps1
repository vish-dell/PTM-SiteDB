$ErrorActionPreference = "Stop"
Write-Host "This folder is the Research sub-application. To launch the complete original-styled website, run the parent project installer:" -ForegroundColor Cyan
Write-Host (Join-Path (Split-Path -Parent $PSScriptRoot) "INSTALL_AND_RUN.ps1") -ForegroundColor Yellow
& (Join-Path (Split-Path -Parent $PSScriptRoot) "INSTALL_AND_RUN.ps1")
