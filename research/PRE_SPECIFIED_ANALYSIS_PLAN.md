# Pre-specified analysis plan — v1.0

**Version date:** 9 October 2026  
**Purpose:** freeze the expanded-panel primary analysis before running a new cohort through this version of the website. This file is a local, versioned plan; it is **not** an independently time-stamped registry preregistration.

## Primary hypothesis

Among genes with experimental PTM evidence and sufficient transcript, genomic-context and structural QC, the distribution of gene-specific PTM-proximity evidence scores `-log10(p_context_null)` is shifted higher in genes listed as cancer drivers by the selected IntOGen source than in PTM-bearing genes **not listed** in that release.

## Primary comparison and settings

- **Reference groups:** IntOGen-listed PTM-bearing driver genes versus PTM-bearing genes not listed in the selected IntOGen release. “Not listed” is not asserted to mean biologically non-driver.
- **Primary cohort-level test:** one-sided Mann–Whitney U comparing each eligible gene's `-log10(p_context_null)` score between the two groups.
- **Primary gene-level proximity test:** one declared setting, sequence separation ≥8 amino acids and Cα–Cα distance ≤10 Å.
- **Null:** within each gene, draw possible coding missense SNV opportunities preserving the observed variant counts in each of the 96 trinucleotide substitution contexts (SBS96); count draws landing within the pre-specified proximity zone. Empirical p-value uses the add-one rule.
- **Monte Carlo settings:** 5,000 null replicates; seed 1 (gene-specific derived seed is logged in the output code); report the exact implementation and random seed in the generated result manifest.
- **Gene-level multiplicity:** Benjamini–Hochberg q-values across the genes that pass all primary QC gates. Gene-level q<0.05 hit fractions and a one-sided Fisher exact comparison are secondary summaries.

## Eligibility and failure rules

A gene is excluded from inferential analysis, and reported as untestable/partial, if its PTM anchors cannot be mapped to the selected sequence/coordinates; relevant Cα coordinates are missing for observed events; genomic variants do not map to the selected coding transcript; or an observed SBS96 context lacks coding missense opportunities. Missing data is never converted into a negative biological result. A gene with zero observed missense events but complete structural/transcript/context eligibility remains in the comparison with observed proximity count 0 and p=1.

## Historical result kept separate from the expanded-panel hypothesis

The earlier six-gene analysis reported no significant protein-level PTM-proximity enrichment for TP53, EGFR, AKT1 or CTNNB1 under its primary setting. BRAF and BRCA1 failed structural eligibility/QC and were untestable, not negatives. The TP53 S127 site-level sensitivity pattern is hypothesis-generating only; it did not independently replicate at the protein level. These are results from the earlier project version, not outcomes from this expanded panel.

## Secondary analyses

- Threshold matrix: sequence separation 8/12/20 × Cα distance 6/8/10/12 Å. Show every cell, including zero hits and QC-limited settings. These 12 settings per cohort are descriptive robustness checks; significance from a non-primary threshold is exploratory and is **not** a confirmatory discovery.
- Site-level PTM-anchor tests are secondary and use across-site BH correction.
- Independent same-cancer replication should use a separately sourced cohort, the same genome assembly/transcript rules, and unchanged primary thresholds. Sample overlap must be checked externally; this website does not automatically prove independence.
- TP53 S127, if present in the real mapped data, is hypothesis-generating only. No claim of biological importance or replication is authorized by a database record alone.

## Limits

SBS96 matching controls the observed trinucleotide substitution spectrum conditional on the gene's candidate coding missense opportunities. It is not COSMIC signature deconvolution, does not model gene-specific callable territory or local regional mutation rate, and does not address indels/multi-nucleotide variants. AlphaFold coordinates are predictions, not experimentally determined structures. One transcript per gene may not represent all isoforms. IntOGen catalog absence is not a gold-standard non-driver label.

---

## v1.1 amendments — 11 October 2026 (made BEFORE any run of the expanded 200-gene panel)

v1.0 above is preserved unchanged. The amendments below respond to specific weaknesses found on review. They are
dated, local and **not** an independently time-stamped registry entry; every v1.0 statistic is still computed and reported.

1. **Recurrence.** v1.0 counted every variant as an independent event, so one hotspot residue mutated in many tumours
   could dominate a gene. v1.1 counts each mutated *residue* once (`--recurrence distinct_residue`, default). The v1.0
   counting (`all_events`) is reported alongside as `p_context_null_all_variants` and may be selected with `--recurrence all_events`.
2. **PTM-specificity control (decoy anchors).** For each gene, the real PTM anchors are compared with `--decoys` (default
   1000) random sets of the same number of non-PTM residues of the same amino-acid identity (S→S, T→T, Y→Y) from the same
   gene. Statistic per set: observed proximal mutated residues minus the exact context-null expectation. Reported:
   `decoy_z`, `p_decoy`, BH `q_decoy`. Genes with fewer than 2× as many eligible decoy residues as anchors are
   `unavailable_decoy_pool_too_small` (untestable for this endpoint, never counted as negative).
3. **Primary endpoint (replaces the v1.0 primary).** One-sided stratified permutation test of `decoy_z`, driver vs
   not-listed, permuting labels only within strata of mutation burden (distinct mutated residues) × PTM-anchor count
   (tertiles when ≥60 genes). The v1.0 Mann–Whitney on `-log10(p_context_null)`, a burden-stratified test of that raw
   score, and the q<0.05 Fisher comparison are all retained as **secondary**.
4. **Grid multiplicity.** Within each cohort, the 12 threshold settings are corrected as one Benjamini–Hochberg family
   (`q_family_all_12_settings`; summary column `genes_q_lt_005_family_corrected`). Only family-corrected hits may be called
   discoveries; the primary setting stays 8 / 10 Å.
5. **Replication rule.** A positive claim requires a significant primary result in the discovery cohort and a same-direction
   one-sided p < 0.05 in an independent same-cancer cohort at unchanged thresholds. A null primary result is reported as null.
6. **Unchanged limits.** CA–CA distance (side-chain geometry not modelled), AlphaFold predictions, one transcript per gene,
   no regional mutation-rate model, "not listed" ≠ validated non-driver, SNV missense only. Decoys match residue identity and
   anchor count but cannot absorb unmeasured regional mutation-rate variation.
