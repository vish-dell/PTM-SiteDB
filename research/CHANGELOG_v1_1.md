# v1.1 — review fixes (11 Oct 2026)

| Review finding | Fix | Where |
|---|---|---|
| Recurrent hotspots counted as independent events (toy check: 10 hits on one residue gave p≈0.0002) | Each mutated residue counts once (`--recurrence distinct_residue`); legacy counting reported as `p_context_null_all_variants` | `collapse_recurrence`, `setting_gene_result`, `site_results` |
| Driver-vs-background comparison confounded by mutation burden / functional regions | (a) decoy-anchor control: same-gene, same-amino-acid non-PTM residues, `decoy_z`/`p_decoy`/`q_decoy`; (b) label permutation within mutation-burden × PTM-anchor strata | `decoy_calibration`, `stratified_label_permutation_test`, `primary_group_comparison` |
| 12 threshold settings not corrected as a family | One BH family across all gene×setting tests (`q_family_all_12_settings`, `genes_q_lt_005_family_corrected`) | `main()` robustness block |
| No results section / threshold figure | `scripts/make_results_report.py` → `RESULTS.md` + `threshold_sensitivity.png` from real outputs only; refuses if none exist | new script |
| Two-cohort run not scripted | `RUN_TWO_COHORTS.ps1` | new script |
| **Bug found while testing:** `load_coordinates` passed `low_memory=False` to pandas' Python parser, which pandas rejects | Removed the argument | `load_coordinates` |
| Speed | Per-position loops replaced with vectorised/cached distance matrices (3,000-residue gene, 100 anchors, 1,000 decoys, 5,000 permutations ≈ 1 s) | `near_mask`, `near_matrix`, `CoordMap` |

Tests: 21 pass (9 new in `tests/test_v11_controls.py`, incl. an end-to-end CLI run on a hand-built fixture genome).
The fixtures are unit-test scaffolding only and are not data or results.

## NOT fixed (cannot be fixed without real data or new inputs)
- **No real-data result exists yet.** Run `RUN_TWO_COHORTS.ps1` on two independent real cohorts.
- CA–CA distance only (side-chain geometry would need Cβ/heavy-atom coordinates from the AlphaFold files).
- No regional mutation-rate / replication-timing model; decoys do not remove that confounder.
- "Not listed in IntOGen" is still not a validated non-driver label.
- Novelty is unchanged by code: HotPho, ActiveDriverDB, HotMAPS, CLUMPS and others already cover PTM–mutation and 3D-clustering ideas. The defensible claim is the
  context-matched, decoy-controlled, burden-stratified, pre-specified driver-vs-background test, supported by real replicated results.
