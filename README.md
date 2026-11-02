# Decode the Phosphoproteome of Human Cancer

An evidence-first analysis platform built on the records of PTM-SiteDB.
Stack: HTML5 · CSS3 · JavaScript · Flask · SQLite · Cytoscape.js · 3Dmol.js.

## Quick start (Windows PowerShell)
```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe sync_data.py          # imports the REAL records from the original PTM-SiteDB
.\.venv\Scripts\python.exe app.py                # http://127.0.0.1:5000
```
If `sync_data.py` cannot read the original site, export a TSV from it yourself and run
`python import_tsv.py export.tsv` (columns: gene, uniprot, ptm_type, site, + optional position, flanking,
kinase, disease_context, evidence, confidence, pmid).

## Data policy — no synthetic, simulated or assumed data
* **No PTM records ship with this project.** The old `seed_phosphorylation_page1.tsv` was deleted: its
  flanking sequences did not look like the proteins they were attached to, and it was silently loaded
  whenever the database was empty. The database is now filled only by `sync_data.py` / `import_tsv.py`.
* Every analysis is either a deterministic computation over the imported records or uses a live public source
  (UniProt, Ensembl Compara, EBI Proteins API, RCSB/PDBe-SIFTS, AlphaFold DB, cBioPortal), cached in `instance/cache`.
* If a source cannot answer, the page says **"Not available"** and shows why. Nothing is estimated or filled in.
* Run **Data Integrity Audit** after importing: it checks every record's residue letter and flanking window
  against the UniProt sequence and lists any that disagree.

## The nine modules (each has its own page under *Analysis Suite*)
| Page | What it computes | Real data used |
|---|---|---|
| PTM Hotspots | Windows holding more modified residues than a uniform null predicts (binomial tail, Bonferroni) | imported records |
| 3D Spatial Distance | Cα–Cα and closest-atom distances between PTM residues; flags pairs far in sequence but close in 3D; interactive 3D view | PDB coordinates via PDBe-SIFTS (each residue letter verified); AlphaFold only as a labelled fallback |
| Kinase–PTM Network | kinase → PTM site → protein → cancer context (Cytoscape.js), kinase ranking, cascade candidates | imported records |
| Conservation Explorer | Is each PTM residue identical in Ensembl Compara orthologues (11 model species) | Ensembl + UniProt |
| Mutation–PTM–Cancer | (A) variants at/near PTM residues with disease/COSMIC/ClinVar cross-refs; (B) per-site mutation counts in a **user-selected** cBioPortal study with an enrichment test | EBI Proteins API, cBioPortal |
| PTM Crosstalk | Different PTM types on the same/nearby residues, shared kinase/disease | imported records |
| Hypothesis Generator | Ranked, testable research *questions*; every card lists the criteria met/unmet, suggested experiments, and PMIDs | all layers above (optional) |
| PhosphoRank | Transparent, adjustable checklist score for choosing sites to follow up | imported records |
| Data Integrity Audit | Record-by-record verification against UniProt | UniProt |

Every page has a *Method & assumptions* section. Hotspots, proximity and crosstalk are labelled **candidate**
relationships; hypotheses are questions, not findings.

## Known limits
* Hotspot/mutation enrichment tests use a uniform-placement null; PTM databases are curation-biased and known
  driver hotspots violate it, so treat results as screens.
* Gene-symbol lookups (Ensembl, cBioPortal) depend on the symbol matching the current HGNC symbol.
* Structure mapping ignores PDB insertion codes; residues that cannot be verified stay unmapped.
* The first run of a page that contacts external services can take a minute; results are cached afterwards.

## v4 additions

- **3D-Proximal Mutations** (`/proximity`): cancer missense mutations (cBioPortal) that are >= 20 residues from every PTM site in sequence but <= 10 A from one in the AlphaFold model (pLDDT >= 70 only). Exact hypergeometric test against uniform placement. This is a screen: real mutations are not uniform (signatures, driver hotspots) and distance alone does not show functional effect.
- **Hotspots**: Benjamini-Hochberg `fdr_q` plus a family-wise permutation p (`perm_p_fwer`) that re-places sites on residues of the same amino-acid letter, removing the S/T/Y-composition bias of the uniform null.
- **Fix**: conserved-site hypotheses now compare identity against orthologues aligned at that site (gaps excluded, >= 3 required).
- **Tests**: `python -m unittest discover -s tests` (FDR, exact test, classification, composition null, mocked pipeline). These use tiny hand-built fixtures; no fixture is ever loaded into the app database.

Run order: `pip install -r requirements.txt`, `python sync_data.py` (or `python import_tsv.py export.tsv`), `python app.py`.

- **Batch screen**: `python batch_proximity.py brca_tcga luad_tcga --sep 8 12 20` runs the 3D-proximal test on every protein with >= 3 PTM units, BH-corrects across proteins, and writes `results/*.csv`. Use one study to discover and an independent one to replicate; sensitivity over `--sep` shows whether hits depend on the arbitrary cutoff.

- **v4.1**: PTM-centric 3D-neighbourhood test per site (`sites_*.csv`, and a table on `/proximity`) and Fisher pooling of independent cohorts (`pooled_sep*.csv`). All inputs are real public data fetched live; nothing is simulated.

## Research-grade 3D PTM–mutation workflow

The primary contribution is an integrated, reproducible, structure-guided PTM–mutation proximity screen. The production pipeline uses real PTM records, real cBioPortal cohort mutations, residue-verified experimental PDB/SIFTS mappings where available, and an explicitly labelled AlphaFold fallback. No biological results are hard-coded.

### Recommended analysis

```bash
python sync_data.py
python batch_proximity.py brca_tcga luad_tcga --sep 8 12 20 --dist 10 --n-perm 5000 --perm-seed 1
python replication_report.py results/proximity_brca_tcga_sep8_dist10.csv results/proximity_luad_tcga_sep8_dist10.csv --out results/replication_sep8_dist10.csv
```

The permutation null preserves the reference-amino-acid composition of distinct mutated positions. The current cBioPortal protein-level endpoint does **not** provide enough genomic context to construct a trinucleotide-context-matched null, so outputs explicitly label that analysis as unavailable rather than substituting a fabricated context.

Permutation p-values are reported separately from the exact hypergeometric baseline. Missing permutation p-values carry an explanatory note. Synonymous mutations, when returned by the same real cohort endpoint, are shown only as a positional calibration control and are not assumed to be biologically neutral.

Spatial analysis reports per-site structural provenance and distinguishes experimental PDB/SIFTS mappings from AlphaFold predicted mappings and unmapped residues. Distances are calculated only from actual coordinates in one selected structure/chain; no placeholder coordinates are created.

### Scientific interpretation

A statistically enriched PTM neighbourhood is a candidate association, not proof of functional causality. Discovery/replication labels should only be assigned from actual cohort result files, and benchmark agreements should only be reported after computing them from verified reference data.

## 27+ validation workflow

The research-facing contribution is the integrated, reproducible, structure-guided and composition-aware PTM–mutation proximity analysis with cohort replication and explicit structural provenance.

### 1. Load the real PTM records

```powershell
python sync_data.py
```

The importer reads the original PTM-SiteDB search pages and writes real records to `instance/ptm.db`; it does not create replacement records.

### 2. Run discovery/replication analysis

```powershell
python batch_proximity.py brca_tcga luad_tcga --sep 8 12 20 --dist 10 --n-perm 5000 --perm-seed 1
```

Use one cohort as discovery and the other as replication only after checking that the cohorts are appropriate and non-overlapping for the intended comparison. A result is not called replicated unless it passes the same pre-specified criterion in both cohorts.

### 3. Build the replication summary

```powershell
python replication_report.py results/proximity_brca_tcga_sep8_dist10.csv results/proximity_luad_tcga_sep8_dist10.csv --out results/replication_sep8_dist10.csv --require-permutation
```

### 4. Sensitivity analysis

Run the batch with several distance cutoffs in one command (`--dist 6 8 10 12`); each (study, sep, dist) is written to its own file, e.g. `proximity_brca_tcga_sep8_dist10.csv`. Do not pick the threshold that gives significance. Summarize with:

```powershell
python sensitivity_report.py --study brca_tcga --seps 8 12 20 --dists 6 8 10 12
```

This is a sensitivity analysis, not a substitute for a trinucleotide-aware cancer mutation model. The current cBioPortal protein-level records do not contain enough genomic-context information to claim a trinucleotide-matched null.

### 5. Benchmarking

Benchmark calls must be supplied from a real, cited reference resource. Do not put hand-entered positive examples into production results. Benchmark CSV columns: `gene, uniprot (optional), site (optional), reference_result (positive|negative), reference_source`. `benchmark_report.py` computes TP/FP/FN/TN, precision, recall and specificity (written to `<out>.summary.csv`); rows our run could not test are reported as `not_comparable`, never as negatives.

```powershell
python benchmark_report.py results/sites_brca_tcga_sep8_dist10.csv benchmark.csv --out results/benchmark_brca_sep8_dist10.csv
```

### 6. Interpretation rule

A candidate is a statistical candidate, not a proven functional PTM–mutation interaction. Report observed count, expected count, fold enrichment, exact p, permutation p where valid, BH q, confidence interval where estimable, structural provenance, sequence separation, distance threshold and replication status.


### Related work
See `RELATED_WORK.md` for a sourced comparison with ActiveDriverDB, HotMAPS and CLUMPS, and for the claims that may and may not be made.


## Expanded Research extension (keeps the original website design)

The original Phosphoproteome Atlas UI, navigation, Search, Browse, Statistics, protein pages, and Analysis Suite remain served from the root routes. The new cohort/benchmark workflow is mounted separately under `/research/`; it has its own `research/data/raw`, `research/data/processed`, and cohort-results files to avoid collisions with the legacy PTM catalogue database. The Research pages reuse the original teal-and-white `main.css` look rather than replacing the entire site's theme.

### Windows setup

From this project root, run `INSTALL_AND_RUN.ps1` or double-click `START_WEBSITE.bat`. The installer creates `.venv312`, installs `requirements.txt` (no `pysam` or `pyfaidx` compiler dependency), attempts to download dbPTM and IntOGen, and builds the real-source panel if both tables are available. Open `http://127.0.0.1:5000/` for the original site and `http://127.0.0.1:5000/research/` for the extended cohort workflow.

No demo biological observations are bundled. If sources or genomic reference inputs are unavailable, the Research pages show a not-loaded/not-ready status; missing values are never replaced with synthetic counts. IntOGen fallback archive location: `research/data/raw/IntOGen-Drivers-20240920.zip`.

## Correct startup path and evidence rules (October 2026)

- Start the complete site from the **project root** with `START_WEBSITE.bat`, or run `.venv312\Scripts\python.exe app.py`. The root `app.py` serves the original PTM-SiteDB interface and mounts the research extension below `/research/`. Do **not** run `python research\app.py` to launch the complete site; that starts the research sub-application by itself and replaces the root page at port 5000.
- The protein detail pages and every feature module link to the 3D viewer. It prefers structures from RCSB PDB when residue mappings are verified and falls back to an explicitly labelled AlphaFold prediction when a compatible model is available. A model or site mapping is never invented; unresolved/mismatched PTM positions remain unhighlighted.
- This source bundle does not include a generated database, the official dbPTM table, an IntOGen driver table, a reference FASTA/GTF, or cohort files. Setup attempts to acquire source data; if that fails, results remain unavailable/zero and the source manifest reports why. Do not interpret missing inputs as negative biological findings.
- Run checks from the project root with `.venv312\Scripts\python.exe -m pytest -q` (after the environment/dependencies are installed).

## Update 2026-10-11 — research analysis controls v1.1
The `research/` module now counts recurrent mutations once per residue, adds a same-gene decoy-anchor control, tests drivers vs not-listed genes within mutation-burden × PTM-anchor strata, and corrects the 12 threshold settings as one family. No biological result is claimed until two real cohorts have been run (`research/RUN_TWO_COHORTS.ps1`) and reported, null results included. Details: `research/CHANGELOG_v1_1.md`.
