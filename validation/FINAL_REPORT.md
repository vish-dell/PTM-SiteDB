# PTM-SiteDB — Final Research, Validation and Integrity Report
Generated 2026-10-25. Every number below was produced by a command whose raw output is in `validation/logs/`.
Labels: **[VERIFIED]** executed here · **[SYNTHETIC]** simulated data · **[SOURCE]** read from a cited publication · **[PENDING]** not done.

## A. Baseline audit
* Structure: root Flask site (`app.py`, `features.py`, `services/`, `templates/`, `static/`) + research Flask app mounted at `/research`
  (`research/`: scripts, tests, 200-gene panel, cohort outputs). 11 analysis modules exist (hotspots, spatial, network, conservation,
  mutations, proximity, crosstalk, hypotheses, priority, integrity, features index). **No new module was added** (Phase 3 below).
* Baseline reproduced [VERIFIED]: **5 failed, 90 passed, 1 skipped**. After fixes [VERIFIED]: **101 passed, 1 skipped, 0 failed, 0 warnings** (`logs/pytest_final.txt`).
  The skip needs Node + `npm install jsdom` (browser rebinding test) — not run.
* Production database [VERIFIED]: `instance/ptm.db` now has **30 verified phosphorylation sites** (cleaned 25 Oct 2026; originally 67 records: 51 phospho, 11 acetyl, 5 ubiquitin; 37 rows kept in ptm_sites_quarantine) from the original site; the research module
  separately holds 463,956 dbPTM records. Any claim about "database scale" must say which.
* Defects found and fixed: see section D. Defects found, NOT fixable offline: cohort data (section E.3).
* Housekeeping [VERIFIED]: stray 0-byte file named `-ErrorAction`, a UTF-16 `flask_error.log`, and a bundled `.venv312` (1.2 GB total ZIP) — delete before submission.

## B. Research-gap statement (documented capabilities only)
* HotPho [SOURCE: Huang et al., Nat Commun 2021;12:2313, doi:10.1038/s41467-021-22481-w]: takes mutation and phosphosite lists, maps them to **PDB** structures, computes pairwise
  distances, clusters, and reports prioritized hybrid clusters; built on the authors' HotSpot3D; handles complexes; applied to 225,151 phosphosites (PTMcosmos/CPTAC) with MC3 mutations.
* ActiveDriverDB [SOURCE: Krassowski et al., NAR 2018;46:D901, doi:10.1093/nar/gkx973; 2021 update, Front Cell Dev Biol 9:626821, doi:10.3389/fcell.2021.626821]:
  database annotating TCGA/ClinVar/population variants against PTM sites (±7 aa windows; "direct" = substitution of the modified residue), MIMP kinase-motif rewiring, networks, upload + API; hg19.
* Gap PTM-SiteDB can defensibly target: a **small, fully inspectable, reproducible** per-gene test that (i) works on predicted (AlphaFold, pLDDT-filtered) structures,
  (ii) uses an explicit eligible-residue null with exact + composition-matched permutation p-values, (iii) reports "not estimable" instead of 0/1, and (iv) ships a
  calibration/power study of its own test. I could **not** verify from retrieved text how HotPho/ActiveDriverDB handle their nulls or multiple testing, so no superiority claim is made.
* Unresolved: no head-to-head benchmark; no validated cohort evidence (E.3).

## C. Methodology actually implemented (read from source)
Per gene, AlphaFold Cα (or closest-atom) distances; residues with pLDDT ≥ cutoff (default 70). **Eligible** = not a PTM residue and ≥ `min_sep` (default 20) from every PTM.
**Qualifying** = eligible and ≤ `threshold` Å (default 10) from the nearest PTM. Test: P(X ≥ k), X ~ Hypergeometric(N eligible, K qualifying, n distinct mutated eligible positions) — exact,
equals `scipy.stats.hypergeom.sf(k−1,N,K,n)` [VERIFIED on every N ≤ 12 grid case, `tests/test_proximity.py`]. Plus a seeded composition-matched permutation test; Wilson CI on fold enrichment (approximate).
Cohort-level (research module): SBS96-context null, decoy-anchor control, burden-stratified permutation test, BH correction, driver (IntOGen 2024.09.20 release) vs not-listed panel genes.
Limitations: AlphaFold only for the statistic; no regional mutation-rate model in the single-gene test (see E.2-C); "not listed" ≠ non-driver.

## D. Contributions, separated
**Pre-existing:** all analysis modules, cohort pipeline, 3Dmol viewer, hypergeometric test, SBS96 null.
**Done in this work [VERIFIED]:**
1. Fixed the 5 failing tests (3Dmol assets; 4 missing proximity controls; stale null-universe expectations; `scripts` import path).
2. Independent verification of the hypergeometric tail + exhaustive SciPy comparison tests.
3. Phase-4 null/planted-signal/confounding simulation of the production `proximity_test` (E.2).
4. **Cohort-independence guard** (`research/scripts/cohort_guard.py`, wired into the upload route and `make_results_report.py`, 3 tests) after discovering duplicated cohort inputs (E.3).
5. Reproducible audit script `validation/audit_cohort_inputs.py`.

## E. Validation results
### E.1 Tests
`validation/logs/pytest_final.txt` — 101 passed, 1 skipped, 0 failed. Python 3.12.3, Flask 3.1.3, SciPy 1.17.1, NumPy 2.4.4, pandas 3.0.2 (sandbox Linux, not your Windows machine).
Two tests were **corrected, not bypassed**: `test_site_neighbourhood` 0.4→2/3 and `test_hit_found_and_reported` 0.5→1.0. Hand derivations plus SciPy show the old values counted the PTM residue and sequence-near residues
(which can never qualify) in N. Under the project's own stated eligible-residue null the new values are exact. If your supervisor intended the all-residue universe, that is a different (flawed) null — decide deliberately.

### E.2 Simulation [SYNTHETIC] (`validation/simulation/`, seed 20261025, `logs/simulation_run.txt`)
Tests the production `proximity_test` (exact test), not a separate demo. Synthetic data never touch the DB.
* **A. Null** — 24 scenarios (N 100–1000, 5%/15% qualifying, n 5–100), 4,000 replicates each: max FPR 0.0112 @α=0.01, 0.049 @0.05, 0.091 @0.10;
  **no scenario's 95 % Wilson lower bound exceeds α** (`null_fpr.csv`, `null_calibration.png`). The test is **conservative** (discrete p-values; often well below α), not miscalibrated upward.
* **B. Planted signal** (N=300, 10 % qualifying, 2,000 reps), power at α=0.05:

| n mutated | RR 1 (null) | 1.5 | 2 | 3 | 5 |
|---|---|---|---|---|---|
| 10 | 0.012 | 0.038 | 0.076 | 0.220 | 0.492 |
| 20 | 0.036 | 0.130 | 0.262 | 0.560 | 0.876 |
| 50 | 0.048 | 0.218 | 0.480 | 0.873 | 0.994 |
| 100 | 0.030 | 0.299 | 0.696 | 0.979 | 1.000 |

  Interpretation: adequate power only for strong effects (RR ≥ 3) with ≥ 20 mutated positions; **a per-gene null result at typical TCGA counts is weak evidence of absence** (`power.csv`, `power_curves.png`).
* **C. Alternative setting — regional mutation-rate heterogeneity** (no true proximity effect; qualifying residues in one contiguous patch; n=30; 3,000 reps):
  FPR@0.05 = 0.017 (σ=0), 0.060 (σ=0.5), **0.124 (σ=1.0), 0.143 (σ=1.5)** (`confounding_fpr.csv`). **The exact test is not robust to regional rate variation.** This is a real limitation; it is the reason the
  decoy/context controls matter and must not be dropped. The method is therefore validated for a *uniform-propensity null only*.
* **D. Permutation vs exact** (300 pairs, uniform composition): correlation 0.90, mean |Δp| 0.094; permutation FPR@0.05 = 0.020 vs exact 0.047. They agree in rank but **not** closely in value — the repair notes' claim that they "agree" is only partly supported. Needs investigation (the +1 smoothing and the small permutation count (499) are candidate causes; not tested).

### E.3 Cohorts — **the shipped PAAD/PRAD results are invalid as labelled** [VERIFIED]
`logs/cohort_input_audit.txt`: `TCGA_PAAD.maf` and `TCGA_PRAD.maf` have **identical SHA-256** (`9cb009b1…`), identical rows (126,252), 1,009 samples and identical top tissue-source-site codes (BH, D8, AN, AC, A8; these look like breast sites — verify against the GDC TSS table).
1,009 samples fits TCGA-BRCA, not PAAD (~185) or PRAD (~500). Both runs report 581 candidate missense SNVs; their p-values differ only through permutation count (1000 vs 5000: 0.1455 vs 0.1434).
Consequences: no PAAD result, no PRAD result and no two-cohort "replication" exists. At most one BRCA-like result (pending identity confirmation) exists, null (decoy-calibrated stratified p = 0.80; 0 genes q<0.05).
`make_results_report.py` now **refuses** to report these. **BRCA/LUAD/PAAD with correct MAFs: [PENDING]** — needs downloads I cannot make from the sandbox. Commands are in `RUN_COMMANDS.md`.
Do not pool, do not cite PAAD/PRAD results.

### E.4 Comparison (documented-feature comparison only — no benchmark, no performance claims)
| Aspect | PTM-SiteDB (implemented) | HotPho | ActiveDriverDB |
|---|---|---|---|
| Objective | per-gene test: are missense mutations 3D-near PTM sites more than expected; cohort driver-vs-non-driver benchmark | discover 3D co-clusters of phosphosites and mutations [SOURCE] | annotate/interpret variants via PTM sites, kinase-motif rewiring [SOURCE] |
| Structures | AlphaFold, pLDDT-filtered (experimental structures shown, not used in the statistic) | PDB structures, incl. complexes [SOURCE] | not a structural-distance method (sequence windows ±7 aa) [SOURCE] |
| PTMs | dbPTM (research), 30 verified phosphosites (DB, after cleanup) | phosphosites only (225,151) [SOURCE] | multiple PTM types, >385k (2018) / >260k experimental (2021 update) [SOURCE] |
| Mutations | uploaded MAF / cBioPortal | MC3 [SOURCE] | TCGA, ClinVar, population; direct vs indirect effects [SOURCE] |
| Null / multiple testing | exact hypergeometric + permutation, SBS96 context null, BH | **not verified** | **not verified** |
| Validation | own null/power simulation (this report); no external benchmark | not assessed here | not assessed here |
| Known limits | uniform-propensity null (E.2-C), tiny DB, no valid cohort yet | — | hg19 [SOURCE] |

### E.5 Phase 3 audit
Hotspots, spatial, network, conservation, mutations, crosstalk, priority and integrity pages already exist; I did **not** add features or edit their statistics (no offline way to validate sources/alignments, and
adding more would raise feature count, not rigor). Recommended next improvement: add multiple-testing correction and a stated background to the hotspot module, and re-run its tests against analytical examples.

## F. Academic-integrity review — **[PENDING: no similarity report supplied]**
I cannot do a match-by-match review without the report. Upload the PDF/HTML report and the document. Until then: 6 % is unreviewed; I will not rewrite text blindly, invent citations or claim "plagiarism-free".
| Matched passage | Source | Reason | Needs fix? | Revision | Citation | Action |
|---|---|---|---|---|---|---|
| (awaiting report) | | | | | | none |
Verified references ready to cite: Huang 2021 (HotPho); Krassowski 2018 and 2021 (ActiveDriverDB). Reimand & Bader 2013 (ActiveDriver) is from memory — verify before citing. 3Dmol.js (BSD-3) licence is bundled in `static/vendor/`.
Also: your repo's `RELATED_WORK.md`/`NOVELTY.md` were not rewritten.

## G. Scores (honest, current evidence)
| | /10 | Evidence | To reach 9 |
|---|---|---|---|
| Novelty | 5 | transparent eligible-residue null, not-estimable reporting, calibration study, independence guard | a validated cohort result + a head-to-head benchmark against HotPho/ActiveDriver on shared data; fix regional-rate robustness |
| Research gap | 6 | gap is real and sourced for HotPho/ADDB features | verify comparators' nulls from Methods; show a case where the method finds what they miss (or honestly doesn't) |
| Methodology | 6 | exact test verified vs SciPy; null FPR controlled; power quantified; honest confounding failure | model regional rate (e.g., stratified/segment-matched null) and show FPR restored; resolve perm-vs-exact gap; real independent cohorts |
| **Total** | **17/30** | | **27+/30 is not supported by present evidence**; it needs the pending items |

## H. Changed / new files
Changed: `templates/viewer.html`, `templates/feature_spatial.html`, `templates/feature_proximity.html`, `research/app.py`, `research/scripts/make_results_report.py`, `tests/test_proximity.py`, `tests/test_proximity_pipeline.py`.
New: `research/scripts/cohort_guard.py`, `research/tests/test_cohort_guard.py`, `validation/` (simulation, audit script, logs, this report), `RUN_COMMANDS.md`.
Remaining: Phase 5 cohorts, Phase 7 similarity review, real-data comparator benchmark, jsdom test, regional-rate-robust null.

## Addendum (25 Oct 2026): cohort inputs corrected
Section E.3 describes the earlier duplicated MAF inputs. Those were replaced. audit_cohort_inputs.py on the current files reports distinct hashes, TCGA_PAAD = 176 samples, TCGA_PRAD = 493 samples, and no duplicate groups. Results were regenerated from the current files.

