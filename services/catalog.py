"""Protein catalogue: standard public UniProt metadata (gene, accession, name, length).

This is NOT PTM data. It only lets Browse list every curated protein.
The integrity checker (services/external.py + /integrity) re-verifies the
lengths against live UniProt sequences.
"""
PROTEIN_CATALOG = [
    ("AKT1",   "P31749", "RAC-alpha serine/threonine-protein kinase", 480),
    ("AKT3",   "Q9Y243", "RAC-gamma serine/threonine-protein kinase", 479),
    ("BRAF",   "P15056", "Serine/threonine-protein kinase B-raf", 766),
    ("BRCA1",  "P38398", "Breast cancer type 1 susceptibility protein", 1863),
    ("CHEK1",  "O14757", "Serine/threonine-protein kinase Chk1", 476),
    ("CTNNB1", "P35222", "Catenin beta-1", 781),
    ("EGFR",   "P00533", "Epidermal growth factor receptor", 1210),
    ("MAPK1",  "P28482", "Mitogen-activated protein kinase 1 (ERK2)", 360),
    ("MAPK3",  "P27361", "Mitogen-activated protein kinase 3 (ERK1)", 379),
    ("MAPK8",  "P45983", "Mitogen-activated protein kinase 8 (JNK1)", 427),
    ("PIK3CA", "P42336", "Phosphatidylinositol 4,5-bisphosphate 3-kinase catalytic subunit alpha", 1068),
    ("PIM1",   "P11309", "Serine/threonine-protein kinase pim-1", 313),
    ("RB1",    "P06400", "Retinoblastoma-associated protein", 928),
    ("SMAD4",  "Q13485", "Mothers against decapentaplegic homolog 4", 552),
    ("TP53",   "P04637", "Cellular tumor antigen p53", 393),
]
