# PTM-SiteDB audit (2 Nov 2026)
Copy the files in this folder over the same paths in your project (instance/ptm.db is already corrected; backup first).

## Fixed
| Area | Problem | Fix |
|---|---|---|
| Code | services/analysis.py: sites keyed without protein, so ALL-network merged TP53 S33/S37 into CTNNB1 | uniprot added to key |
| Code | batch_verify.py: edge-truncated flanks scored 0/12 (TP53 S6) | now 15/15 |
| Code | app.py: "Proteins 15" on home/statistics/API = catalogue size; only 7 have records | now 7, `catalogue_proteins`=15 added |
| CSS | global `table{min-width:1000px}` | scoped to `.table-wrap>table`; small/JS tables exempt; connection table compact |
| Data | 4 impossible PMIDs inside pmids_dbptm (92065884, 91184134 x2, 976) | tokens removed, pmid_count fixed, logged in table data_audit_log; no row deleted |

## Passed (0 errors, audit_ptm_db.py)
30/30 sites: residue letter = position = UniProt residue, none beyond protein length, flank centre OK (4 flanks are shorter only because the site is at a protein end), no duplicates, all Phosphorylation on S/T/Y, accessions valid and unique, catalogue lengths = UniProt for the 7 proteins with records, no PMID duplicated within a cell, primary PMID always in its list, no tyrosine-only kinase on S/T sites.
All pages and API routes return 200 (404/400 only for ALL on per-protein APIs and spatial/model without parameters, by design). Counts for sites (30), kinases (16), publications (26) equal the database.

## Not verifiable offline (run these with internet)
- `python pubmed_check.py` (then `--all`) : PubMed existence, retraction, gene mention. There was no pubmed_check.py in the zip, so this is new.
- `python batch_verify.py instance/ptm.db` after importing expanded_sites_min2pmids.tsv.
- Catalogue lengths of AKT3, BRAF, BRCA1, MAPK8, PIK3CA, PIM1, RB1, SMAD4 (no cached UniProt file); the values match my recollection of UniProt but are unchecked.
- /api/f/integrity, conservation, variants, cbio, and live AlphaFold spatial calls (network).

## Needs your decision
- 8 kinase annotations (TP53 S6, T18, S37, S183, S215, S315, S371; AKT1 T450) look doubtful from my own recollection of the literature; they are kept in the DB but listed under kinase_pending_review in ptm_sitedb_clean_sites.csv. Check them against the cited papers.
- 7 sites have no kinase ('—'). Only 15 of 30 sites have an accepted kinase.
- Crosstalk is empty (all records are Phosphorylation). The 4 test modules that need pytest were not run here (20 other tests pass; your own log shows 101 passed).
- 8 of 15 catalogue proteins have 0 sites; importing expanded_sites_min2pmids.tsv fills BRAF, BRCA1, PIK3CA.
- Novelty/plagiarism: I cannot compute a novelty or plagiarism score. ptm_sitedb_clean_sites.csv is built only from your own verified database; use a similarity checker on your report text.

## AlphaFold / mapping (research panel, 100 genes, 79,267 residues)
Accession-gene match with the UniProt map: 100/100. pLDDT in range: yes. Duplicate positions: none. Models with sequence_identity < 1: 3 (>=0.995). 36 flags are intentional model-vs-canonical numbering offsets (e.g. AFF3, DSP, BDNF), handled by the sequence-mapping step; 737 of 4,151 sites carry a shifted position for that reason and 0 fail the S/T/Y centre check. 6 genes failed mapping QC and are excluded (CHD7, FAT4, MROH9, PDLIM2, RPP14, VPS13B). Only 32% of sites have pLDDT >= 70, so 3D-proximity results use the confident minority. Of your 7 database genes only EGFR is in the panel (TP53, AKT1 etc. are fetched live via /api/f/spatial).
