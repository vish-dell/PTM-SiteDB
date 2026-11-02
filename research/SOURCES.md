# Source and citation register

Checked 9 October 2026. URLs are listed for citation and reproducibility; not all datasets are redistributed with this source package.

## Primary data

1. **dbPTM** — official home: https://biomics.lab.nycu.edu.tw/dbPTM/ ; download page: https://biomics.lab.nycu.edu.tw/dbPTM/download.php ; experimental phosphorylation download directory: https://biomics.lab.nycu.edu.tw/dbPTM/download/experiment/ . The download is retrieved at setup time. Data source release/retrieval metadata and local file hashes are logged.
2. **IntOGen Cancer Drivers** — official download page: https://www.intogen.org/download ; release notes: https://www.intogen.org/releasednotes . The website records the release declared inside an archive README when one is present. If a bundled archive filename suggests one release but its internal metadata says another, cite the internal declaration and investigate the discrepancy.
3. **GDC / TCGA MAF** — data format: https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/ and https://docs.gdc.cancer.gov/Encyclopedia/pages/Mutation_Annotation_Format/ . The user must supply a legally accessible cohort MAF.
4. **Genome FASTA and GTF** — user-supplied matched assembly (GRCh38 recommended). Record source/version and confirm MAF, FASTA and GTF builds agree.
5. **UniProtKB** — API used for gene-symbol mapping and primary sequence comparison: https://rest.uniprot.org/ ; human TP53 reviewed entry: https://www.uniprot.org/uniprot/P04637 .
6. **AlphaFold Protein Structure Database** — https://alphafold.ebi.ac.uk/ ; API/coordinates are predicted structural models and are explicitly labelled as such.

## Related work

7. PTMcode: https://pmc.ncbi.nlm.nih.gov/articles/PMC3531129/
8. Mutation3D: https://pmc.ncbi.nlm.nih.gov/articles/PMC4833594/
9. SpacePAC: https://pmc.ncbi.nlm.nih.gov/articles/PMC4227039/
10. HotSpot3D: https://www.nature.com/articles/ng.3586
11. ActiveDriverDB: https://pubmed.ncbi.nlm.nih.gov/29126202/
12. HotPho: https://pubmed.ncbi.nlm.nih.gov/33875650/ ; https://pmc.ncbi.nlm.nih.gov/articles/PMC8055881/
13. PRISMOID: https://pubmed.ncbi.nlm.nih.gov/31161204/ ; https://pmc.ncbi.nlm.nih.gov/articles/PMC7299293/
14. ClusTar (2026): https://pubs.acs.org/doi/10.1021/acs.jcim.6c01767

## Important restrictions and study limits

PhosphoSitePlus terms: https://www.phosphosite.org/psp/legal/terms-and-conditions . The package does not automatically retrieve PhosphoSitePlus data. SBS96 matching preserves observed trinucleotide substitution counts within gene-specific candidate coding missense opportunities; it is not a fitted COSMIC mutational signature mixture and does not model regional callability or all mutation processes. The website ships no synthetic biological observations.
