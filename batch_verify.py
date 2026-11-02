#!/usr/bin/env python3
"""
batch_verify.py - verify EVERY PTM site in the SQLite database.

For each site it checks, using only real public data (no assumptions):
  1. UniProt: is the residue at that position the amino acid the record says,
     and does the record's flanking sequence match UniProt?
  2. PDB (via RCSB search + PDBe SIFTS numbering): is there an experimental structure
     whose residue at that position is the right amino acid (not a mutant, not disordered)?
  3. AlphaFold: residue identity and pLDDT confidence at that position.

Usage (PowerShell, from the project folder):
  python batch_verify.py --inspect                 # find DB and show its tables/columns
  python batch_verify.py path\\to\\your.db           # run on every site
Options if auto-detection fails: --table --acc-col --site-col --flank-col
Output: verification_report.csv  (downloads are cached in .\\verify_cache, so reruns are fast)
"""
import argparse, csv, glob, json, os, re, sqlite3, sys, time
import urllib.request, urllib.error
from concurrent.futures import ThreadPoolExecutor

T2O = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
       "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
       "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}
MODS = {"SEP": "S", "TPO": "T", "PTR": "Y", "MSE": "M", "ALY": "K", "M3L": "K", "MLY": "K", "MLZ": "K",
        "M2L": "K", "KCX": "K", "CSO": "C", "HYP": "P", "SNN": "N"}
THREE_TO_ONE = {**T2O, **MODS}
PTM_RES = {"SEP", "TPO", "PTR", "ALY", "M3L", "MLY", "MLZ", "M2L", "KCX", "CSO", "HYP"}
ONE = set("ACDEFGHIKLMNPQRSTVWY")
CACHE = "verify_cache"
os.makedirs(CACHE, exist_ok=True)


# ----------------------------------------------------------------- network + cache
def http(url, data=None, tries=3):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, data=data, headers={
                "Content-Type": "application/json", "User-Agent": "ptm-verify/1.0"})
            with urllib.request.urlopen(req, timeout=60) as r:
                return "" if r.status == 204 else r.read().decode("utf-8", "replace")
        except urllib.error.HTTPError as e:
            if e.code in (400, 404):
                return None
            time.sleep(2 * (i + 1))
        except Exception:
            time.sleep(2 * (i + 1))
    return None


def cached(key, url, data=None):
    path = os.path.join(CACHE, re.sub(r"[^\w.-]", "_", key))
    if os.path.exists(path):
        with open(path, encoding="utf-8") as f:
            return f.read()
    txt = http(url, data)
    if txt:
        with open(path, "w", encoding="utf-8") as f:
            f.write(txt)
    return txt


# ----------------------------------------------------------------- database
def find_db():
    found = [p for ext in ("*.db", "*.sqlite", "*.sqlite3")
             for p in glob.glob(os.path.join("**", ext), recursive=True)
             if "venv" not in p and "node_modules" not in p]
    return found


def schema(con):
    cur = con.cursor()
    tabs = [r[0] for r in cur.execute("select name from sqlite_master where type='table'")]
    return {t: [r[1] for r in cur.execute(f'pragma table_info("{t}")')] for t in tabs}


def pick(cols, *needles):
    for n in needles:
        for c in cols:
            if n in c.lower():
                return c
    return None


def load_sites(con, a):
    sc = schema(con)
    table = a.table or next((t for t, c in sc.items() if pick(c, "site") and pick(c, "flank")), None) \
        or next((t for t, c in sc.items() if pick(c, "site")), None)
    if not table:
        sys.exit("Could not find a table with a 'site' column. Run --inspect and pass --table/--site-col.")
    cols = sc[table]
    site = a.site_col or pick(cols, "site")
    flank = a.flank_col or pick(cols, "flank")
    acc = a.acc_col or pick(cols, "uniprot", "accession")
    if acc:
        sql = f'select "{acc}", "{site}", ' + (f'"{flank}"' if flank else "''") + f' from "{table}"'
    else:  # accession lives in a protein table joined by protein_id
        fk = pick(cols, "protein_id")
        ptab = next((t for t, c in sc.items() if t != table and pick(c, "uniprot", "accession") and "id" in [x.lower() for x in c]), None)
        if not (fk and ptab):
            sys.exit(f"No UniProt column in '{table}' (columns: {cols}). Pass --acc-col, or send me the --inspect output.")
        pacc = pick(sc[ptab], "uniprot", "accession")
        sql = (f'select p."{pacc}", t."{site}", ' + (f't."{flank}"' if flank else "''") +
               f' from "{table}" t join "{ptab}" p on p.id = t."{fk}"')
    print(f"Using table '{table}': {sql}")
    return con.execute(sql).fetchall()


# ----------------------------------------------------------------- parsing helpers
def parse_site(s):
    s = str(s or "").strip()
    m = re.match(r"^(Ala|Arg|Asn|Asp|Cys|Gln|Glu|Gly|His|Ile|Leu|Lys|Met|Phe|Pro|Ser|Thr|Trp|Tyr|Val)\s*[-_:]?\s*(\d+)", s, re.I)
    if m:
        return T2O[m.group(1).upper()], int(m.group(2))
    m = re.match(r"^(?:ph|p)?([ACDEFGHIKLMNPQRSTVWY])\s*[-_:]?\s*(\d+)", s)
    return (m.group(1), int(m.group(2))) if m else (None, None)


def flank_score(seq, pos, aa, flank):
    f = re.sub(r"\s+", "", str(flank or "")).upper()
    if len(f) < 5:
        return None
    best = None
    for c in range(len(f)):
        if f[c] != aa:
            continue
        cmp = hit = 0
        for i, ch in enumerate(f):
            if i == c or ch not in ONE:
                continue
            j = pos - 1 + i - c
            if 0 <= j < len(seq):
                cmp += 1
                hit += seq[j] == ch
        if cmp and (best is None or hit / cmp > best[0] / best[1]):
            best = (hit, cmp)
    return best


# ----------------------------------------------------------------- PDB probing
def probe(txt, chain, resi):
    """Look at one residue in a PDB file: residue name, B-factor, SEQADV note, REMARK 465 (missing)."""
    out = {"resn": None, "b": None, "seqadv": None, "missing": False}
    for l in txt.splitlines():
        rec = l[:6]
        if rec == "ENDMDL":
            break
        if rec in ("ATOM  ", "HETATM"):
            if l[21] == chain and l[22:26].strip() == str(resi) and l[26] == " " and l[16] in " A1":
                rn = l[17:20].strip()
                if rec == "HETATM" and rn not in MODS:
                    continue
                out["resn"] = rn
                if l[12:16].strip() == "CA":
                    try:
                        out["b"] = float(l[60:66])
                    except ValueError:
                        pass
        elif rec == "SEQADV" and l[16] == chain and l[18:22].strip() == str(resi):
            out["seqadv"] = l[49:].strip()
        elif rec == "REMARK":
            m = re.match(r"^REMARK 465\s+\w{3}\s+(\S)\s+(-?\d+)\s*$", l)
            if m and m.group(1) == chain and int(m.group(2)) == resi:
                out["missing"] = True
    return out


def judge(p, aa):
    if p["resn"] is None:
        return "DISORDERED" if p["missing"] else "NO_COORDS"
    one = THREE_TO_ONE.get(p["resn"])
    if one != aa:
        return f"DIFFERENT({p['resn']})"
    if p["seqadv"] and re.search(r"MUTATION|CONFLICT", p["seqadv"], re.I):
        return "ENGINEERED_MUTANT"
    return "OK_PHOSPHO_PRESENT" if p["resn"] in PTM_RES else "OK"


class Protein:
    """Everything fetched once per UniProt accession."""
    def __init__(self, acc):
        self.acc = acc
        self.base = acc.split("-")[0]
        fa = cached(f"{acc}.fasta", f"https://rest.uniprot.org/uniprotkb/{acc}.fasta")
        self.seq = "".join(fa.splitlines()[1:]).strip().upper() if fa else ""
        self.blocks = {}      # pdb id -> SIFTS mapping blocks
        self.ids = self._search()
        with ThreadPoolExecutor(8) as ex:
            list(ex.map(self._sifts, self.ids))
        self.af = self._alphafold()

    def _search(self):
        t = lambda a, v: {"type": "terminal", "service": "text",
                          "parameters": {"attribute": a, "operator": "exact_match", "value": v}}
        q = {"query": {"type": "group", "logical_operator": "and", "nodes": [
            t("rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_name", "UniProt"),
            t("rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession", self.base)]},
            "return_type": "entry",
            "request_options": {"paginate": {"start": 0, "rows": 1000}, "results_content_type": ["experimental"]}}
        txt = cached(f"search_{self.base}.json", "https://search.rcsb.org/rcsbsearch/v2/query", json.dumps(q).encode())
        try:
            return [r["identifier"] for r in json.loads(txt).get("result_set", [])] if txt else []
        except Exception:
            return []

    def _sifts(self, pdb):
        txt = cached(f"sifts_{pdb}.json", f"https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb.lower()}")
        try:
            d = json.loads(txt)[pdb.lower()]
            self.blocks[pdb] = d.get("UniProt", {})[self.base]["mappings"]
        except Exception:
            self.blocks[pdb] = []

    def _alphafold(self):
        txt = cached(f"af_{self.base}.json", f"https://alphafold.ebi.ac.uk/api/prediction/{self.base}")
        try:
            url = json.loads(txt)[0]["pdbUrl"]
            return cached(f"af_{self.base}.pdb", url)
        except Exception:
            return None

    def covering(self, pos):
        out = []
        for pdb, blocks in self.blocks.items():
            for b in blocks:
                try:
                    if b["unp_start"] <= pos <= b["unp_end"]:
                        out.append((pdb, str(b["chain_id"]), b["start"]["author_residue_number"] + (pos - b["unp_start"])))
                except (KeyError, TypeError):
                    pass
        return out


# ----------------------------------------------------------------- per-site verification
def verify(prot, aa, pos, flank, max_pdb):
    r = {"uniprot_residue": "", "sequence_ok": "", "flank_match": "", "n_entries_checked_for_search": len(prot.ids),
         "n_covering": 0, "experimental": "", "exp_detail": "", "alphafold": "", "af_plddt": ""}
    if not prot.seq:
        r["sequence_ok"] = "UNIPROT_UNAVAILABLE"
    else:
        have = prot.seq[pos - 1] if 0 < pos <= len(prot.seq) else "(beyond end)"
        r["uniprot_residue"] = have
        r["sequence_ok"] = "YES" if have == aa else "NO"
        fs = flank_score(prot.seq, pos, aa, flank)
        r["flank_match"] = "no flank in record" if fs is None else f"{fs[0]}/{fs[1]}"

    cov = prot.covering(pos)
    seen, tried, notes = set(), 0, []
    r["n_covering"] = len({c[0] for c in cov})
    if not prot.ids:
        r["experimental"] = "NO_PDB_ENTRIES"
    elif not cov:
        r["experimental"] = "OUTSIDE_ALL_SOLVED_SEGMENTS"
    else:
        r["experimental"] = "COVERED_BUT_NOT_VERIFIED"
        for pdb, chain, resi in cov:
            if pdb in seen and tried >= max_pdb:
                continue
            if pdb not in seen:
                if tried >= max_pdb:
                    break
                seen.add(pdb)
                tried += 1
            txt = cached(f"{pdb}.pdb", f"https://files.rcsb.org/download/{pdb}.pdb")
            if not txt:
                notes.append(f"{pdb}:no .pdb file")
                continue
            st = judge(probe(txt, chain, resi), aa)
            if st.startswith("OK"):
                r["experimental"] = "VERIFIED"
                r["exp_detail"] = f"{pdb} chain {chain} res {resi} {st}" + (f" (+{r['n_covering']-1} other covering entries)" if r["n_covering"] > 1 else "")
                break
            notes.append(f"{pdb}/{chain}:{st}")
        if r["experimental"] != "VERIFIED":
            r["exp_detail"] = f"checked {tried} of {r['n_covering']} covering entries: " + "; ".join(notes[:8])

    if prot.af:
        p = probe(prot.af, "A", pos)
        if p["resn"]:
            r["alphafold"] = "OK" if THREE_TO_ONE.get(p["resn"]) == aa else f"DIFFERENT({p['resn']})"
            r["af_plddt"] = "" if p["b"] is None else round(p["b"])
        else:
            r["alphafold"] = "NOT_IN_MODEL"
    else:
        r["alphafold"] = "NO_MODEL"

    if r["sequence_ok"] == "NO":
        r["verdict"] = "RECORD_ERROR: UniProt has a different residue here"
    elif r["experimental"] == "VERIFIED":
        r["verdict"] = "VERIFIED_EXPERIMENTAL_STRUCTURE"
    elif r["alphafold"] == "OK":
        r["verdict"] = "SEQUENCE_OK + ALPHAFOLD_ONLY" if r["sequence_ok"] == "YES" else "ALPHAFOLD_ONLY"
    else:
        r["verdict"] = "UNVERIFIED"
    return r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("db", nargs="?")
    ap.add_argument("--inspect", action="store_true")
    ap.add_argument("--table"); ap.add_argument("--acc-col"); ap.add_argument("--site-col"); ap.add_argument("--flank-col")
    ap.add_argument("--max-pdb", type=int, default=15, help="max covering PDB entries downloaded per site")
    ap.add_argument("--out", default="verification_report.csv")
    a = ap.parse_args()

    if not a.db:
        dbs = find_db()
        if a.inspect or len(dbs) != 1:
            print("Database files found:", *dbs, sep="\n  ")
            if not a.inspect:
                sys.exit("Pass the one you want: python batch_verify.py <path>")
        a.db = dbs[0] if dbs else sys.exit("No .db/.sqlite file found.")
    con = sqlite3.connect(a.db)
    if a.inspect:
        for t, c in schema(con).items():
            print(f"{t}: {c}")
        return

    rows = load_sites(con, a)
    print(f"{len(rows)} site records")
    prots, results = {}, []
    for n, (acc, site, flank) in enumerate(rows, 1):
        acc = str(acc or "").strip()
        aa, pos = parse_site(site)
        base = {"uniprot": acc, "site": site}
        if not acc or not aa:
            results.append({**base, "verdict": "UNREADABLE_SITE_OR_ACCESSION"})
            print(f"[{n}/{len(rows)}] {acc} {site}: unreadable")
            continue
        if acc not in prots:
            print(f"   fetching data for {acc} ...")
            prots[acc] = Protein(acc)
        res = verify(prots[acc], aa, pos, flank, a.max_pdb)
        results.append({**base, **res})
        print(f"[{n}/{len(rows)}] {acc} {site}: {res['verdict']}")

    fields = ["uniprot", "site", "verdict", "uniprot_residue", "sequence_ok", "flank_match", "experimental",
              "exp_detail", "n_covering", "n_entries_checked_for_search", "alphafold", "af_plddt"]
    with open(a.out, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(results)
    print("\nSummary:")
    for v in sorted({r["verdict"] for r in results}):
        print(f"  {sum(r['verdict'] == v for r in results):3d}  {v}")
    print(f"\nFull table written to {a.out}")


if __name__ == "__main__":
    main()