#!/usr/bin/env python3
"""Offline re-runnable audit of instance/ptm.db (sites, proteins, PMID structure, kinase format).
    python audit_ptm_db.py       # exit code 1 if any hard error is found
Needs verify_cache/<ACC>.fasta (created by batch_verify.py) for the residue/length checks."""
import os, re, sqlite3, sys, collections
sys.path.insert(0, ".")
from services.catalog import PROTEIN_CATALOG as PC
VALID = {"Phosphorylation", "Acetylation", "Ubiquitination", "Methylation", "SUMOylation", "Glycosylation"}
ACC = re.compile(r"^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$")
con = sqlite3.connect("instance/ptm.db"); con.row_factory = sqlite3.Row
err, warn = [], []
seqs = {f[:-6]: "".join(open("verify_cache/" + f).read().splitlines()[1:]) for f in os.listdir("verify_cache") if f.endswith(".fasta")}
seen = collections.Counter()
for r in con.execute("select * from ptm_sites"):
    m = re.fullmatch(r"([A-Z])(\d+)", r["site"] or ""); tag = f"{r['gene']} {r['site']}"
    if not m: err.append((tag, "bad site format")); continue
    if int(m.group(2)) != r["position"]: err.append((tag, "site number != position"))
    s = seqs.get(r["uniprot"])
    if s is None: warn.append((tag, "no cached UniProt sequence"))
    elif r["position"] > len(s): err.append((tag, "position > protein length"))
    elif s[r["position"] - 1] != m.group(1): err.append((tag, "residue != UniProt"))
    fl = (r["flanking"] or "").upper()
    if len(fl) == 21 and fl[10] != m.group(1): err.append((tag, "flank centre != residue"))
    if r["ptm_type"] not in VALID: err.append((tag, "invalid ptm_type"))
    if r["ptm_type"] == "Phosphorylation" and m.group(1) not in "STY": err.append((tag, "phospho site not S/T/Y"))
    if not ACC.match(r["uniprot"] or ""): err.append((tag, "bad accession"))
    seen[(r["uniprot"], r["site"], r["ptm_type"])] += 1
    toks = (r["pmids_dbptm"] or "").split(";")
    if len(toks) != len(set(toks)): err.append((tag, "duplicate PMID in cell"))
    for t in set(toks) | {r["pmid"] or ""}:
        if not t.isdigit() or not 1000 <= int(t) <= 45_000_000: err.append((tag, f"implausible PMID {t!r}"))
    if r["pmid"] not in toks: err.append((tag, "primary PMID not in dbPTM list"))
    if r["kinase"] is None or str(r["kinase"]).strip() in ("", "None", "nan"): err.append((tag, "blank kinase value"))
err += [(k, "duplicate site row") for k, v in seen.items() if v > 1]
used = {r[0] for r in con.execute("select uniprot from ptm_sites")}
err += [(u, "site uniprot not in protein catalogue") for u in used - {p[1] for p in PC}]
for g, u, n, l in PC:
    if not ACC.match(u): err.append((g, "bad catalogue accession"))
    if u in seqs and len(seqs[u]) != l: err.append((g, f"catalogue length {l} != UniProt {len(seqs[u])}"))
    if u not in used: warn.append((g, "catalogue protein has 0 site records"))
print("ERRORS:", len(err)); [print("  ", e) for e in err]
print("WARNINGS:", len(warn)); [print("  ", w) for w in warn]
sys.exit(1 if err else 0)
