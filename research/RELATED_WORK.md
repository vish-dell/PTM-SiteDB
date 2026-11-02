# Related-work review and contribution scope

**Web/PubMed check date:** 9 October 2026. Selected primary papers and resource pages were opened directly. The searches were targeted rather than a registered systematic review; completeness across all databases/languages is not claimed.

## Comparison of established work

| Resource / method | Existing contribution | Consequence for this project |
|---|---|---|
| **PTMcode** — https://pmc.ncbi.nlm.nih.gov/articles/PMC3531129/ | PTM functional associations supported by multiple evidence channels, including structural proximity. | Structural proximity around PTMs is not itself a new concept. |
| **dbPTM 2025** — https://biomics.lab.nycu.edu.tw/dbPTM/ and https://biomics.lab.nycu.edu.tw/dbPTM/download.php | Large resource of experimentally reported and other PTM evidence with UniProt-linked annotations. | Used here as a scalable source; the project does not claim to originate PTM data. |
| **PhosphoSitePlus** — https://www.phosphosite.org/psp/legal/terms-and-conditions | Major curated PTM resource. Its data-use terms must be checked before any bulk download or redistribution. | Not automatically downloaded by this package; use only if the applicable agreement allows it. |
| **Mutation3D** — https://pmc.ncbi.nlm.nih.gov/articles/PMC4833594/ | Structural clustering of cancer-associated variants. | Existing mutation-to-structure clustering prior art. |
| **SpacePAC** — https://pmc.ncbi.nlm.nih.gov/articles/PMC4227039/ | Statistical analysis of spatial somatic-mutation clusters in protein structures. | Existing spatial clustering and null-model prior art. |
| **HotSpot3D** — https://www.nature.com/articles/ng.3586 | Protein-structure-guided mutation hotspot detection, including intramolecular/intermolecular neighborhoods. | Existing 3D hotspot framework. |
| **ActiveDriverDB** — https://pubmed.ncbi.nlm.nih.gov/29126202/ | PTM-aware interpretation of cancer mutations and disease-associated variation. | Direct conceptual neighbor; historical prior work, not a dependency of this website. |
| **HotPho** — https://pubmed.ncbi.nlm.nih.gov/33875650/ and https://pmc.ncbi.nlm.nih.gov/articles/PMC8055881/ | Directly identifies spatial co-clustering of phosphorylation sites and cancer mutations, including hundreds of hybrid clusters. | PTM–mutation 3D co-clustering predates this project. |
| **PRISMOID** — https://pubmed.ncbi.nlm.nih.gov/31161204/ and https://pmc.ncbi.nlm.nih.gov/articles/PMC7299293/ | 3D resource for PTM sites, nearby/on-site mutations, and structural/functional context. | Existing integrated 3D PTM–mutation resource. |
| **ClusTar (2026)** — https://pubs.acs.org/doi/10.1021/acs.jcim.6c01767 | Pan-cancer framework for spatial PTM–mutation clusters in protein–protein interaction context, published online July 2026. | Further invalidates an unqualified “first of its kind” claim. |

## Defensible contribution

> This project implements a reproducible, auditable hypothesis-testing workflow that maps experimentally recorded PTM residues and cohort mutations onto a declared protein-transcript/structure coordinate system, tests PTM-proximal mutation enrichment using a gene-specific SBS96-matched coding missense opportunity null, compares the resulting gene-level evidence scores between IntOGen-listed drivers and PTM-bearing genes not listed in that release, and exposes QC failures and a complete predeclared threshold-sensitivity matrix.

This is a **methodological integration and transparency claim**, not a claim that PTM–mutation spatial analysis, 3D clustering, or PTM-aware cancer mutation interpretation is unprecedented. The contribution must be judged on the correctness of mapping, null-model assumptions, panel size, QC pass rate, and results from real independent cohorts—not on the number of UI features.

## PubMed query log

Queries opened on 9 October 2026. These URLs allow the searches to be reproduced and updated:

1. PTM/phosphosite + 3D/structural + mutation/clustering: https://pubmed.ncbi.nlm.nih.gov/?term=%28PTM+OR+phosphosite%29+AND+%283D+OR+structural%29+AND+%28mutation+OR+mutation+clustering%29
2. PTM + structure-guided + cancer mutation: https://pubmed.ncbi.nlm.nih.gov/?term=%28post-translational+modification%29+AND+%28structure-guided+OR+protein+structure%29+AND+cancer+mutation
3. TP53 Ser127 phosphorylation: https://pubmed.ncbi.nlm.nih.gov/?term=%22TP53%22+AND+%22Ser127%22+AND+phosphorylation
4. PTMcode original paper: https://pmc.ncbi.nlm.nih.gov/articles/PMC3531129/
5. HotPho original paper/full text: https://pmc.ncbi.nlm.nih.gov/articles/PMC8055881/
6. PRISMOID original paper/full text: https://pmc.ncbi.nlm.nih.gov/articles/PMC7299293/
7. ClusTar 2026 article: https://pubs.acs.org/doi/10.1021/acs.jcim.6c01767

## Source and license rules

- dbPTM data downloads should be cited and used under their published terms.
- PhosphoSitePlus is not bundled or silently scraped; its current terms may require a separate agreement for downloads.
- IntOGen release identity is recorded from the archive contents when available. A ZIP filename/date alone is not treated as authoritative.
- Search results and site-level hints do not establish causality or experimental validation.
