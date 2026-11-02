# PTM-SiteDB integrated repair notes — 10 October 2026

## What was repaired

- Rebuilt from the clean `PTM-SiteDB_OriginalStyle_Research_Integrated_v2` source archive instead of carrying forward the locally corrupted `app.py` / front-end files.
- Kept the original website at `/` and mounts the cohort/research application at `/research/`; `START_WEBSITE.bat` launches the project-root `app.py` and does not run the research sub-app as the root page.
- Fixed the broken genomic-preparation environment path: it now looks for `.venv312` at the project root and resolves FASTA/GTF paths before changing working directory.
- Made research status pages tolerate missing, empty, or malformed result CSVs without crashing or turning invalid files into results.
- Made source ZIP discovery handle nested ZIPs, gzip-compressed tables, extensionless text tables, CSV/TSV/TXT/DAT/TAB files, and spreadsheet tables; failures now report an archive inventory to help diagnose the actual official archive contents.
- Restored defensive hypothesis response handling and added cache-busting for the browser feature script to avoid stale JavaScript.
- Made every protein detail page and every feature page offer the 3D structure viewer, including curated proteins with no PTM records currently loaded.
- Reworked the viewer to try mapped experimental PDB coordinates first and fall back to a clearly labelled AlphaFold predicted model. A PTM marker is only highlighted if the residue exists and its residue identity matches the site label; unmatched sites remain unhighlighted.
- Added regression tests for startup wiring, empty-row 3D routes, hypothesis JavaScript safety, and genomic input paths.

## Validation performed in this build environment

- `pytest`: 35 tests passed.
- Python compilation: all 38 Python files passed.
- Jinja parsing: all 37 HTML templates parsed.
- Node syntax checks passed for `viewer.js`, `features.js`, `circmap.js`, and `app.js`.

A live Flask HTTP smoke test was not possible in this build container because Flask is not installed in its Python environment. The included Windows installer installs the packages in `requirements.txt`; run it once when setting up the project.

## Scientific/data limitation (important)

This ZIP does not bundle the official dbPTM experimental table, IntOGen driver archive/table, a reference FASTA/GTF, a completed panel, or user cohort data. The installer attempts official downloads, but network/source layout/format can prevent ingestion. In that case the site must show unavailable or zero source-backed results until ingestion succeeds. No placeholder PTM records, invented cohort results, fabricated structures, or guessed site mappings are added. AlphaFold models are predictions, not experimental structures, and neither model availability nor residue coverage can be guaranteed for every protein/site.
