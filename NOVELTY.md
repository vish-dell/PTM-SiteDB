# Contribution statement (for the report / viva)

**Project:** Decode the Phosphoproteome of Human Cancer

No numeric novelty or quality score is claimed anywhere in this project: that is a judgement for evaluators and depends on
what they compare against. What follows is the evidence you can present, and the checks to run first.
Nothing here is described as first-of-kind, novel as a concept, or unique.

## What is different from existing resources
Reference sites (PTM-SiteDB, PhosphoSitePlus, dbPTM, UniProt) are mainly search/browse catalogues.
This project adds an **analysis layer that is traceable to evidence**:

| # | Contribution | Why it is more than a catalogue |
|---|---|---|
| 1 | Hotspot detection with a stated statistical null and multiple-testing correction | reports significance, not just counts |
| 2 | Residue-verified 3D distances, incl. "far in sequence, close in 3D" pairs | distances only where the residue letter is confirmed in real coordinates |
| 3 | Kinase → site → protein → cancer network with cascade candidates | connects records across proteins |
| 4 | Site-level conservation from Compara orthologue alignments | per-site, per-species evidence |
| 5 | PTM-site × cancer-cohort mutation overlay with enrichment test (user-chosen cBioPortal study) | links PTM sites to patient mutation counts |
| 6 | Crosstalk candidates with shared kinase/disease | integrates records |
| 7 | Evidence-linked hypothesis generator showing every criterion, assay idea and PMID | reasoning is inspectable |
| 8 | PhosphoRank — adjustable, transparent checklist | not a black-box score |
| 9 | Data Integrity Audit against UniProt | quantifies curation errors in the records |

These are integrations of established ideas in one interactive, no-synthetic-data platform. The scope of the contribution is
integration, transparency and reproducibility; it is not a claim that any single analysis is new.

## Before you claim novelty (from your own guide)
Search the Nucleic Acids Research Database Issue, the journal *Database* (Oxford) and bio.tools for a resource
that already offers the same combination, and cite what you find in a "related work" table.

## 3D PTM–mutation proximity (research module)

Sequence-only mutation-PTM tools cannot flag mutations that are distant in sequence but adjacent in 3D to a modified residue. `/proximity` does this with an explicit exact test and confidence-filtered AlphaFold coordinates. No claim is made that no existing tool does this: HotPho, PRISMOID, ClusTar (2026), ActiveDriverDB, HotMAPS, CLUMPS and others already work on closely related questions (see `RELATED_WORK.md`). Report results from real runs only, including null results.

## Implemented methodological upgrade: composition-matched null

The 3D PTM–mutation screen now reports two complementary null models:

1. **Exact hypergeometric baseline:** the original uniform-placement test over confident AlphaFold residues.
2. **Amino-acid-composition-matched permutation test:** the same number of distinct missense-mutated positions is randomly re-sampled within each reference-amino-acid class from confident residues. The default is 5,000 reproducible permutations (`--n-perm 5000`, `--perm-seed 1`).

Permutation p-values are reported as `p_perm`; BH-adjusted values are reported as `q_perm` for protein-level results and `q_perm_all_sites` for site-level results. The original exact p/q values are retained, so the new method does not replace or hide the baseline.

This is a methodological refinement of existing approaches, not a claim that the concept is new or unique in the literature. This original null does not model trinucleotide context, mutational signatures, replication timing, clonality, or patient-level dependence; the research module (`research/`) adds a trinucleotide-matched null.

## Research module v1.1 controls (see `research/CHANGELOG_v1_1.md`)
Recurrent mutations count once per residue; a decoy-anchor control compares PTM anchors with same-gene, same-amino-acid non-PTM residues; the driver-vs-not-listed test is stratified by mutation burden and PTM-anchor count; the 12 threshold settings are corrected as one family. No biological result is claimed until real cohorts have been run and reported (null results included). TP53 S127 remains a hypothesis-generating database/site candidate only.

## Related work
A sourced comparison is in `RELATED_WORK.md`. Do not use "no tool does this", "first" or "novel" wording.
