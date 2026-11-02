# PTM-SiteDB Integrated Explorer v3 — scalable, source-audited rebuild

A local Flask website for testing whether cancer-cohort missense mutations occur unusually close to experimentally recorded PTM residues in protein 3D space. This is a **new analysis workflow**; it does not overwrite or reinterpret the previous six-gene findings as results from the expanded panel.

## What v3 adds

1. **A scalable PTM source:** bootstraps the official dbPTM experimental phosphorylation download and records local source hashes.
2. **Balanced 200-gene target:** selects up to 100 genes from the real intersection of PTM-bearing genes and the downloaded IntOGen driver catalog, plus up to 100 PTM-bearing genes not listed in that release. Selection is deterministic (seed 1), with the actual available and selected counts recorded. The app warns if the real sources cannot support 200 genes.
3. **A real driver/background comparison:** the primary cohort-level statistic is a one-sided Mann–Whitney U test of per-gene `-log10(p_context_null)` scores for eligible IntOGen-listed drivers versus PTM-bearing genes not listed in that release. The q<0.05 hit-fraction Fisher test is secondary.
4. **An SBS96 context-matched coding-opportunity null:** a matched genome FASTA/GTF pair is used to enumerate coding missense SNV opportunities. The null preserves the observed 96-channel trinucleotide substitution counts within each gene and samples from eligible coding opportunities.
5. **Strict mapping QC:** mismatched assembly/context, missing structural positions, incomplete transcript mapping and untestable PTM anchors are explicit statuses, not negative biological results. Genes with zero observed missense events but complete QC remain tested with observed proximity count 0 and p=1.
6. **A compact threshold matrix:** all 12 settings (separation 8/12/20 × distance 6/8/10/12 Å) are shown per cohort and can be viewed together for multiple cohorts. The single primary setting remains separation ≥8 and distance ≤10 Å.
7. **An expanded related-work review:** PTMcode, dbPTM, PhosphoSitePlus terms, Mutation3D, SpacePAC, HotSpot3D, ActiveDriverDB, HotPho, PRISMOID and ClusTar (2026) are acknowledged directly. No first-of-kind claim is made.
8. **No synthetic biological results:** the package ships code and blank data directories, not artificial mutations/PTM rows. Cohort results are generated only from source downloads and user-provided MAF/FASTA/GTF inputs.

## Historical six-gene result

The previous small-panel primary analysis found no significant protein-level PTM-proximity enrichment for TP53, EGFR, AKT1 or CTNNB1. BRAF and BRCA1 were untestable due to structural eligibility/QC. TP53 S127 remains hypothesis-generating only and did not independently replicate at the protein level. Those are prior results, not results from the newly expanded panel.

## Defensible contribution

> This project implements a reproducible, auditable hypothesis-testing workflow that maps experimental PTM residues and cohort mutations onto a declared transcript/structure coordinate system, tests proximity under a gene-specific SBS96-matched coding missense opportunity null, compares gene-level evidence between an external driver catalog and a PTM-bearing not-listed background, and exposes mapping/QC failures plus a complete pre-specified sensitivity matrix.

This is a methodological integration/transparency contribution—not a claim that PTM–mutation spatial co-clustering is novel. See `RELATED_WORK.md` and `PRE_SPECIFIED_ANALYSIS_PLAN.md`.

## Windows setup — integrated original-style site

This repository is an add-on to the original PTM-SiteDB web site. The original Home/Search/Browse/Statistics/protein pages and Analysis Suite remain at the root path; this research module is mounted under `/research/` and reuses the same teal-and-white stylesheet. Do **not** run this folder's `app.py` on its own when you want the complete combined website.

From the **top-level project folder** (the folder containing the original `app.py` and `INSTALL_AND_RUN.ps1`), run:

```powershell
Set-ExecutionPolicy -Scope Process Bypass
.\INSTALL_AND_RUN.ps1
```

The installer creates the root `.venv312`, installs the combined root `requirements.txt` (no `pysam`/`pyfaidx` compilation dependency), attempts source downloads, and builds a real-source panel only when both source tables are present. Then it launches the combined app at `http://127.0.0.1:5000/`; the research pages are at `http://127.0.0.1:5000/research/`.

If the IntOGen download is blocked, place the official `IntOGen-Drivers-20240920.zip` at `research\data\raw\IntOGen-Drivers-20240920.zip`, then run from the project root:

```powershell
.\.venv312\Scripts\python.exe research\scripts\bootstrap_sources.py
.\.venv312\Scripts\python.exe research\scripts\build_200gene_panel.py
```

The archive's internal declared release is logged where possible. A ZIP filename date is not authoritative; check `research/data/processed/source_manifest.json`, especially if its README reports a release different from the download-page label.

## Genome opportunity and structure preparation

The installer cannot safely guess which genome assembly your cohort uses. Obtain a reference FASTA and gene annotation GTF/GTF.GZ from the same build (GRCh38 recommended), and confirm that the MAF uses that same build. Then run:

```powershell
.\research\PREPARE_GENOMIC_INPUTS.ps1 -ReferenceFasta "C:\data\GRCh38.fa" -Gtf "C:\data\gencode.annotation.gtf.gz"
```

This enumerates all possible missense single-nucleotide coding opportunities for the selected transcript per panel gene, calculates SBS96 context, writes canonical protein sequences, retrieves AlphaFold structures and maps coordinates/PTM residues. It can take substantial time and disk space for 200+ genes. Inspect `opportunity_manifest.json`, `structure_manifest.json`, and `preparation_qc.json` before running a cohort. If only a quick development check is needed, use `-MaxGenes 10`; results from a capped set are not 200-gene study results.

## Run a cohort

Use **Run Analysis** to upload a real somatic MAF (`.maf` or gzip-compressed MAF). Required common columns are gene symbol, chromosome, start position, reference allele and alternate allele; consequence classification is used when present to filter to missense SNVs. The genome FASTA is required to independently verify the reference base and calculate flanking trinucleotide context.

Primary defaults: sequence separation ≥8 residues, Cα–Cα distance ≤10 Å, 5,000 Monte Carlo replicates, seed 1, BH q cutoff 0.05. Tick robustness to write all 12 sensitivity settings. The analysis uses only sites whose PTM/source-to-sequence mapping is supported. PTM sites with absent structure coordinates keep the gene untestable rather than being treated as distant.

## Statistics and limitations

- The SBS96 null is a **trinucleotide spectrum-matched empirical null**, not fitted COSMIC signature deconvolution. It does not model regional replication timing, gene-specific callable territory, local mutation-rate variation or all selection effects.
- Only missense SNVs are analyzed in the present implementation; indels and multi-nucleotide variants are excluded.
- AlphaFold structures are predicted models, not experimental structures; mapping identity/coverage and pLDDT are retained.
- The primary group comparison excludes genes that fail prespecified QC. It is a hypothesis test, not proof of causality. “Not listed by IntOGen” is not a gold-standard biological non-driver label.
- The 12-setting grid is a secondary descriptive sensitivity analysis; non-primary threshold hits are exploratory and are not corrected as a family by this version.
- Separate cohort files can be analyzed, but the website does not automatically verify sample overlap. Independent validation must be designed and documented by the researcher.
- TP53 S127 is hypothesis-generating only and may be shown only if supported by current real inputs/results. See `S127_EVIDENCE.md`; do not claim specific biological validation without it.
- PhosphoSitePlus is not downloaded by this package. Its current terms must be checked before bulk use/redistribution.

## Verification

Run checks after installing dependencies:

```powershell
.\.venv312\Scripts\python.exe research\scripts\preflight_check.py
.\.venv312\Scripts\python.exe -m pytest -q research
```

Review all warning/error fields in `data/processed/preparation_qc.json`. The application binds to `127.0.0.1` by default so uploaded cohort files remain local unless you deliberately change deployment settings.

## v1.1 changes (see `CHANGELOG_v1_1.md` and the amendment at the end of `PRE_SPECIFIED_ANALYSIS_PLAN.md`)

- Primary endpoint is now the **decoy-calibrated PTM-specific score** compared between IntOGen drivers and not-listed genes with a **mutation-burden × PTM-anchor stratified permutation test**.
- Recurrent mutations count once per residue; the 12-setting grid is BH-corrected as one family.
- Run two independent cohorts with `.\research\RUN_TWO_COHORTS.ps1 -MafDiscovery a.maf.gz -MafReplication b.maf.gz -Fasta GRCh38.fa`, then read `research\RESULTS.md`.
- Until real cohorts have been run, this project has **no biological result**; do not describe it as having one.
