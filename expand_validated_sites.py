#!/usr/bin/env python3
"""Build an import file of MORE experimentally supported phospho-sites from the dbPTM download
you already have (research/data/processed/dbptm_phosphorylation_source.txt).

    python expand_validated_sites.py --genes TP53 AKT1 EGFR CTNNB1 MAPK1 MAPK3 CHEK1 BRAF BRCA1 PIK3CA \
        --min-pmids 2 --out expanded_sites.tsv
    python import_tsv.py expanded_sites.tsv
    python batch_verify.py instance/ptm.db          # UniProt + PDB + AlphaFold check of every site

Nothing is invented: gene/UniProt/position/flank/PMIDs are copied from dbPTM; the residue letter is read from
the centre of dbPTM's own 21-mer and must be S/T/Y. Kinase / disease / confidence are left blank
(dbPTM's download does not carry them). Keep only what batch_verify marks VERIFIED/SEQUENCE_OK.
"""
import argparse, csv, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SRC = HERE / "research/data/processed/dbptm_phosphorylation_source.txt"
MAP = HERE / "research/data/processed/uniprot_gene_map.tsv"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genes", nargs="+", required=True)
    ap.add_argument("--source", default=str(SRC)); ap.add_argument("--gene-map", default=str(MAP))
    ap.add_argument("--min-pmids", type=int, default=2, help="require >= N distinct PMIDs (default 2)")
    ap.add_argument("--out", default="expanded_sites.tsv")
    a = ap.parse_args()
    want = {g.upper() for g in a.genes}
    acc2gene = {}
    with open(a.gene_map, encoding="utf-8-sig") as f:
        rd = csv.reader(f, delimiter="\t"); hdr = [h.lower() for h in next(rd)]
        ia = next(i for i, h in enumerate(hdr) if "uniprot" in h or "accession" in h)
        ig = next(i for i, h in enumerate(hdr) if "gene" in h)
        for r in rd:
            if len(r) > max(ia, ig) and r[ig].upper() in want:
                acc2gene[r[ia]] = r[ig]
    seen, rows = set(), []
    with open(a.source, encoding="utf-8", errors="replace") as f:
        for line in f:
            p = line.rstrip("\r\n").split("\t")
            if len(p) < 6 or p[1] not in acc2gene or p[3] != "Phosphorylation" or not p[2].isdigit():
                continue
            fl = p[5].strip().upper()
            pm = sorted({x for x in p[4].split(";") if x.strip().isdigit()})
            if len(fl) != 21 or fl[10] not in "STY" or len(pm) < a.min_pmids:
                continue
            key = (p[1], int(p[2]))
            if key in seen:
                continue
            seen.add(key)
            rows.append([acc2gene[p[1]], p[1], "Phosphorylation", f"{fl[10]}{p[2]}", p[2], fl, "", "",
                         "dbPTM experimental", "", pm[0]])
    rows.sort(key=lambda r: (r[0], int(r[4])))
    with open(a.out, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["gene", "uniprot", "ptm_type", "site", "position", "flanking", "kinase",
                    "disease_context", "evidence", "confidence", "pmid"])
        w.writerows(rows)
    from collections import Counter
    print(f"{len(rows)} candidate sites ->", a.out, dict(Counter(r[0] for r in rows)))

if __name__ == "__main__":
    main()
