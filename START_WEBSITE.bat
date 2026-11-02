@echo off
setlocal
cd /d "%~dp0"
if not exist ".venv312\Scripts\python.exe" goto setup
".venv312\Scripts\python.exe" -c "import flask,requests,pandas,numpy,scipy,Bio" >nul 2>&1
if errorlevel 1 goto setup
echo Starting the complete PTM-SiteDB site (original interface + research extension)...
echo Main site: http://127.0.0.1:5000/
echo Research extension: http://127.0.0.1:5000/research/
".venv312\Scripts\python.exe" app.py
goto done
:setup
echo First launch or missing dependency: installing/checking dependencies and real data sources...
powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0INSTALL_AND_RUN.ps1"
:done
pause
