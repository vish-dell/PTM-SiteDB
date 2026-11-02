# Related work (sourced; verify before citing)

Entries below come from the resources' own abstracts/pages, read on 2026-10-06. This is **not** an exhaustive
literature review. PTMcode, dbPTM, PhosphoSitePlus and Mutation3D/SpacePAC/Hotspot3D were not read in detail here;
check them directly. No uniqueness is claimed. The fuller 9 Oct 2026 review is `research/RELATED_WORK.md`.

| Resource | What it does (per its own description) | How this project differs |
|---|---|---|
| ActiveDriverDB (Krassowski et al., NAR 2018; Front Cell Dev Biol 2021) | Annotates TCGA/ClinVar/population substitutions at PTM sites. Impact classes are defined **in sequence**: direct (at the PTM residue), proximal (1-2 positions), distal (3-7 positions), plus kinase-motif network rewiring. | Looks only along the sequence. This project asks about mutations **>= `--sep` residues away in sequence but within `--dist` A in 3D** of a PTM site, which a sequence-window definition cannot express. |
| HotMAPS (Tokheim et al., Cancer Res 2016) | Finds 3D hotspot regions of missense mutations in PDB structures across TCGA; hotspots linked to conservation, driver scores and interfaces. | Clusters **mutations**; it is not PTM-site-centred. This project fixes the PTM site as the anchor and tests enrichment in its 3D neighbourhood. |
| CLUMPS (Kamburov et al., PNAS 2015) | Tests 3D clustering of missense mutations in PDB structures and at interaction interfaces (4,742 tumours). | Same difference: cluster-of-mutations vs PTM-anchored neighbourhood test. |
| HotPho (2021) | 3D co-clustering of phosphosites and cancer mutations; reports hundreds of hybrid spatial clusters. | **Direct prior art** for PTM–mutation 3D co-clustering. This project differs only in method: a context-matched exact null, decoy-anchor control and a pre-specified driver-vs-not-listed test. |
| PRISMOID (2019) | 3D PTM database with mutation context (mutations on and near PTM sites). | Direct prior art for 3D PTM–mutation context; a database, not a hypothesis test. |
| ClusTar (2026) | Pan-cancer framework for spatial PTM–mutation clusters in a protein–protein-interaction context. | Closely related recent work; read before comparing. |
| Other 3D-clustering tools (SpacePAC, Mutation3D, Hotspot3D) | Named in the literature as 3D mutation-cluster detectors. | Not reviewed in detail here. |

## Statement this review supports
"ActiveDriverDB defines PTM-proximal mutations in sequence (within 7 residues); structure-based tools such as
HotMAPS and CLUMPS cluster mutations without anchoring on PTM sites; HotPho and PRISMOID already relate PTM sites and
mutations in 3D. This project is a reproducible hypothesis-testing workflow for sequence-distant mutations in the 3D
neighbourhood of PTM sites, with context-matched and decoy-controlled nulls, AlphaFold/PDB provenance and cohort
replication, and it reports null results."

## Claims you must NOT make
"No existing tool does PTM-anchored 3D mutation proximity" (HotPho and PRISMOID do closely related work), "first", "novel",
"unique", or "validated S127 discovery".

## Still-open limitation
The null models amino-acid composition only. It does not model trinucleotide context, signatures, replication timing or
clonality (cBioPortal protein-level data lacks genomic context). Treat hits as screening candidates.

## Sources
- Krassowski M, et al. Nucleic Acids Res 2018;46(D1):D901-D910. doi:10.1093/nar/gkx973
- Krassowski M, et al. Front Cell Dev Biol 2021;9:626821. doi:10.3389/fcell.2021.626821
- Tokheim C, et al. Cancer Res 2016;76(13):3719-3731. doi:10.1158/0008-5472.CAN-15-3190
- Kamburov A, et al. Proc Natl Acad Sci U S A 2015 (CLUMPS)
