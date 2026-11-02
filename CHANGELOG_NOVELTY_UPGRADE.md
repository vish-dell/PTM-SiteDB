# Method upgrade — 3D PTM Mutation Analysis (file formerly titled "Novelty Upgrade")

## What changed

The original exact hypergeometric screen is preserved as the baseline. A second, composition-matched permutation null has been added to reduce bias caused by uneven reference-amino-acid composition.

### New method
- Keeps the observed number of distinct missense-mutated positions for each reference amino-acid class.
- Randomly samples positions from confident AlphaFold residues (pLDDT >= 70).
- Computes one-sided empirical p-values using a +1 correction.
- Applies Benjamini-Hochberg correction to permutation p-values.
- Uses a reproducible seed.
- Defaults to 5,000 permutations.

### New command options
```powershell
python batch_proximity.py brca_tcga luad_tcga --sep 8 12 20 --dist 10 --n-perm 5000 --perm-seed 1
```

### New output columns
Protein-level CSV:
- `p_perm`
- `q_perm`
- `n_permutations`

Site-level CSV:
- `p_perm`
- `q_perm_protein`
- `q_perm_all_sites`
- `n_permutations`

The existing `p_value`, `q_bh`, and `q_bh_all_sites` remain unchanged.

## Scientific limitation
This null controls amino-acid composition only. It does not model trinucleotide mutation context, mutational signatures, replication timing, chromatin state, cancer-specific hotspot processes, or patient-level clonality. The permutation analysis therefore improves the baseline null but does not establish causality.

## Validation
The project test suite passes after the upgrade.


## Targeted execution upgrade — structural coverage, effect sizes, and control

### 3D structural mapping
- The spatial module now compares the best verified experimental PDB chain with the best real AlphaFold chain.
- AlphaFold is selected only when it provides better verified PTM coverage or when an experimental structure has fewer than two usable PTM sites.
- AlphaFold chain selection is based on actual residue coverage rather than assuming chain A.
- The API now reports experimental and AlphaFold PTM coverage separately.
- The spatial UI shows per-site mapping provenance and explains which real structure was selected.
- Distances are never combined across different structures/chains.

### Statistics
- Protein- and site-level outputs now include fold enrichment.
- An approximate Wilson 95% CI is reported for the observed/expected enrichment ratio when estimable.
- The CI is explicitly an effect-size summary and does not replace the exact or permutation hypothesis tests.

### Negative/calibration control
- Real synonymous mutations returned by the selected cBioPortal study are reported as an explicitly labelled positional calibration control when present.
- They are not assumed to be biologically neutral.
- No synonymous records are invented when the cohort contains none.

### Validation status
- Automated test suite: 14 tests passing.
- Python compilation: passing.
- Real PTM database sync could not be executed in this sandbox because the original PTM-SiteDB host refused the network connection. No substitute records were inserted.
- Therefore no real cohort analysis or benchmark result is claimed from this execution.

## Fix pass (post-evaluation)
- `batch_proximity.py`: `--dist` now accepts several values; output files are named `..._sep<N>_dist<D>.csv`, so a distance sweep no longer overwrites itself and `sensitivity_report.py` finds the files it expects.
- `sensitivity_report.py`: reads only the dist-tagged files (removed the ambiguous fallback to dist-less names).
- `benchmark_report.py`: now computes TP/FP/FN/TN, precision, recall and specificity from a cited positive/negative benchmark; untestable rows are `not_comparable`, not negatives.
- `tests/test_reports.py`: 3 new tests (17 total passing).
- `RELATED_WORK.md`: sourced comparison with ActiveDriverDB, HotMAPS and CLUMPS, with the defensible and non-defensible claims stated.
- Not fixed (needs network/data not available in this sandbox): real BRCA/LUAD results, a real benchmark file, and a trinucleotide-aware null.
