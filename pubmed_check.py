#!/usr/bin/env python3
"""Check every PMID in instance/ptm.db against PubMed (needs internet).

    python pubmed_check.py            # primary PMID of each site  (26 unique)
    python pubmed_check.py --all      # every PMID in pmids_dbptm   (~1,900 unique, ~2 min)

Writes pmid_check.csv with one status per (PMID, gene): OK | NOT_FOUND | RETRACTED | NO_GENE_MENTION.
NO_GENE_MENTION = title/abstract mentions neither the gene symbol nor its protein name (needs a human look;
it is NOT proof of a wrong citation, MS papers often name proteins only in supplements).
"""
import argparse, csv, json, re, sqlite3, sys, time, urllib.parse, urllib.request
import xml.etree.ElementTree as ET
EU = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

def get(url, tries=3):
    for i in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "ptm-audit/1.0"}), timeout=60) as r:
                return r.read().decode("utf-8", "replace")
        except Exception:
            time.sleep(2 * (i + 1))
    return None

def parse_articles(xml_text):
    """{pmid: {'text': title+abstract, 'retracted': bool}} from an efetch XML string."""
    out = {}
    for art in ET.fromstring(xml_text).iter("PubmedArticle"):
        pmid = art.findtext(".//MedlineCitation/PMID")
        title = "".join(art.find(".//ArticleTitle").itertext()) if art.find(".//ArticleTitle") is not None else ""
        abst = " ".join("".join(a.itertext()) for a in art.iter("AbstractText"))
        pubtypes = [p.text or "" for p in art.iter("PublicationType")]
        out[pmid] = {"text": f"{title} {abst}", "retracted": any("Retracted Publication" in p for p in pubtypes)}
    return out

def classify(rec, gene, names):
    if rec is None:
        return "NOT_FOUND"
    if rec["retracted"]:
        return "RETRACTED"
    pat = "|".join(re.escape(x) for x in [gene] + names if x)
    return "OK" if re.search(rf"(?<![A-Za-z0-9])({pat})(?![A-Za-z0-9])", rec["text"], re.I) else "NO_GENE_MENTION"

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--all", action="store_true"); ap.add_argument("--db", default="instance/ptm.db")
    ap.add_argument("--out", default="pmid_check.csv"); a = ap.parse_args()
    sys.path.insert(0, ".")
    from services.catalog import PROTEIN_CATALOG
    names = {g: [nm.split(",")[0].split(" (")[0]] for g, u, nm, l in PROTEIN_CATALOG}
    ALIAS = {"TP53": ["p53"], "CTNNB1": ["beta-catenin", "β-catenin", "catenin"], "MAPK1": ["ERK2", "ERK", "MAPK"],
             "MAPK3": ["ERK1", "ERK", "MAPK"], "AKT1": ["Akt", "PKB", "protein kinase B"], "EGFR": ["ErbB1", "HER1", "EGF receptor"],
             "CHEK1": ["Chk1", "CHK1"], "BRAF": ["B-Raf", "B-RAF"], "BRCA1": ["BRCA-1"], "PIK3CA": ["PI3K", "p110"]}
    for g, al in ALIAS.items():
        names.setdefault(g, []).extend(al)
    con = sqlite3.connect(a.db); con.row_factory = sqlite3.Row
    pairs = set()
    for r in con.execute("select gene,pmid,pmids_dbptm from ptm_sites"):
        toks = (r["pmids_dbptm"] or "").split(";") if a.all else [r["pmid"]]
        pairs |= {(t.strip(), r["gene"]) for t in toks if t and t.strip().isdigit()}
    ids = sorted({p for p, _ in pairs}); recs = {}
    for i in range(0, len(ids), 150):
        xml = get(EU + "efetch.fcgi?db=pubmed&retmode=xml&id=" + ",".join(ids[i:i + 150]))
        if xml is None:
            sys.exit("PubMed unreachable - run again with internet access.")
        recs.update(parse_articles(xml)); time.sleep(0.4)
    with open(a.out, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f); w.writerow(["pmid", "gene", "status", "link"])
        for p, g in sorted(pairs):
            w.writerow([p, g, classify(recs.get(p), g, names.get(g, [])), f"https://pubmed.ncbi.nlm.nih.gov/{p}/"])
    print("wrote", a.out, len(pairs), "PMID-gene pairs")

if __name__ == "__main__":
    main()
