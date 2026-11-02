# Repair notes — 17 October 2026 (3D viewer + 3D-proximal mutation–PTM analysis)

All changes are in the existing files; no second Flask app, no duplicate routes, original layout/templates kept.
**83 tests pass** (`pytest tests research/tests`), including real-route Flask tests and a jsdom browser test.

## 1. Blank 3D viewer — root cause and fix
* **Cause (`static/js/features.js`, `loadSpatialOverview`)**: the 3Dmol viewer was cached in a module variable and created
  once. Every analysis run rebuilds the result area with `innerHTML`, so the cached viewer stayed bound to the *old,
  detached* `<div id="v3d">`. It drew into nothing while the status line still said "structure loaded, N sites marked".
  Reproduced in jsdom on the original file (1 viewer created over 2 runs, live container has no canvas) and fixed.
* **Fix**: new `static/js/ptmviewer.js` keys viewers by the live DOM element, checks the model parsed to >0 atoms, and
  after drawing verifies a non-zero-size canvas exists; otherwise it says so instead of claiming success.
* 3Dmol was loaded from one host (`3Dmol.csb.pitt.edu`). It is now vendored (`static/vendor/3Dmol-min.js`, v2.5.3,
  BSD-3) with CDN fallbacks (`templates/_3dmol.html`).
* Protein viewer (`viewer.js`): uses the shared helper, shows pLDDT legend + "Reset view", and if no structure verifies any
  PTM it still shows the real structure with **no** highlights (it previously reported "unavailable").
* Highlighting is only ever applied after the loaded model itself shows the expected amino acid at that residue.
  Predicted models are labelled PREDICTED everywhere; experimental and AlphaFold remain separate.
* `alphafold_entry()` now selects the API row for *this accession* (not blindly row 0).

## 2. Mutation retrieval (`services/external.py::cbio_gene_mutations`)
Previously assumed profile `<study>_mutations` and sample list `<study>_sequenced`, no pagination handling, and any
failure looked the same. Now, with a full audit trail (`retrieval`) shown in the UI:
study exists → real MUTATION_EXTENDED profile discovered → sequenced sample list discovered (sample count read) →
gene resolved → paginated fetch (full last page / repeated page / non-list ⇒ reported). Distinct statuses:
`ok`, `ok_no_mutations` (valid zero), `study_not_found`, `gene_not_found`, `no_mutation_profile`, `no_sample_list`,
`request_failed`, `malformed_response`, `incomplete`. A failure never becomes "0 mutations" and no other study is substituted.

## 3. Statistics (`services/proximity.py`) — real defects fixed
1. **Wrong null universe.** N included PTM residues and residues within `min_sep` of a PTM, which can never qualify, and the
   draws `n` included mutations on them. Now N = *eligible* residues (confident, not PTM, ≥ min_sep from every verified
   PTM); mutations on PTM/sequence-proximal residues are counted and reported, not tested. Site-level tests use the
   per-site universe. Two old tests encoded the flawed numbers (0.5, 0.4) and were corrected with hand derivations.
2. **`hyper_sf` crashed** (math domain error) when `n − i > N − K`; support bounds now respected (cross-checked against scipy
   on 200 random cases).
3. **Reference-residue check** (previously parsed but never used): mutations whose reference amino acid disagrees with the
   model are excluded with a reason; >10 % mismatch raises a transcript/isoform warning. PTM records are verified the same way.
4. **Distance definition** configurable: `ca` (default, unchanged) or `min_atom` (closest atom pair); thresholds, pLDDT
   cut-off, permutations (≤20 000) and seed are request parameters and UI controls. Defaults unchanged (10 Å, 20, 70).
5. **Synonymous control** looked for the word "synonymous", but MAF/cBioPortal writes "Silent"; both are now recognised.
6. **Not estimable**: every p-value/CI/fold/q that cannot be validly computed is `null` with a specific `*_reason`; the UI
   shows "Not estimable ⓘ" + reason, never 0, 1 or "—". Permutation null is numpy-based, seeded, reports valid permutation
   count, and states exactly which amino-acid stratum made the null impossible. Permutation vs exact agree under uniform composition (tested).
7. `batch_proximity.py`: `n_eligible_residues` held the *qualifying* count; `fold_enrichment`/`fold_ci95` columns were never
   filled; retrieval failures now labelled `not retrieved (<status>)`. New `--dist-mode`, `--plddt-min`.

## 4. UI (`proximity` page)
Controls for distance mode, pLDDT, permutations, seed; retrieval audit; exclusion-reason table; 3D view highlighting verified
PTM residues (teal) and qualifying mutated residues (red) with click-to-focus and reset.

## Not verified / known limitations (be honest in the report)
* **No live-network check was possible in the build sandbox** (cBioPortal, AlphaFold, RCSB, PDBe unreachable). Retrieval is
  tested against a mock that follows the documented cBioPortal v3 API; run `python diagnose_connections.py` and one real
  study once on your machine. A browser WebGL smoke test also needs your machine (jsdom uses a fake 3Dmol).
* The 3D-proximal statistic uses the **AlphaFold** model only. Experimental structures remain in the viewer and the Spatial
  module (SIFTS-verified) but are not used for this test (partial, multi-chain coverage).
* Null model does not include trinucleotide/signature context; panel-coverage (targeted sequencing) does not bias the test
  (it uses mutation *positions*, not frequencies) but does affect any frequency shown elsewhere.
* Wilson-based fold CI is approximate (labelled so). Exact/permutation p-values are screening evidence, not causal.
* The earlier research-plan items (200+ gene scale-up, IntOGen benchmark, MAF trinucleotide null) live under `research/` and
  need your downloaded data; they were not re-run here. `research/data/` is not included in this zip — copy your existing folder back.

## Run
`pip install -r requirements.txt` → `python app.py`. Optional JS test: `npm install jsdom` then `pytest`.

## Addendum — research upload: "Could not index the uploaded FASTA: No module named 'scripts'"
`research/app.py` imports `scripts.fasta_utils`, but `research/` was not on `sys.path` when the site is started from the
project root (`python app.py`). `research/app.py` now adds its own folder to `sys.path`; regression test
`tests/test_research_import_path.py` (fresh process, real startup import, real FASTA index). The separate banner
"Analysis inputs are not fully prepared" is expected until `PREPARE_GENOMIC_INPUTS.ps1` has been run with a genome
FASTA + matching GTF; it is a data-preparation step, not a bug.
