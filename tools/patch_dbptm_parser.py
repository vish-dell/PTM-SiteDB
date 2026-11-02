"""Fix the dbPTM parser in research/scripts/build_200gene_panel.py.

Run from the project root (the folder that contains the 'research' folder):

    .\\.venv312\\Scripts\\python.exe patch_dbptm_parser.py

What it fixes (all three were reproduced on a file with your real layout):
  1. dbPTM rows are  entry_name, accession, position, ptm_type, pmids, flanking_sequence  (6 columns).
     The old parser assumed  accession, position, type, flank, pmid  and so read the entry name as the
     UniProt ID and took digits from the accession as the position.
  2. The file mixes species. Only rows whose UniProt entry name ends in _HUMAN (UniProt's organism code
     for Homo sapiens) are kept. Species is never guessed from gene names.
  3. `out["gene"] = ""` on an empty DataFrame creates zero rows, so the UniProt gene-symbol lookup never
     ran and every row was dropped without a message.

Every dropped row is counted and printed with its reason. Nothing is dropped silently.
The layout is accepted only if at least 95% of the first 5,000 lines match the 6-column pattern
(entry name, valid UniProt accession, numeric position). Any other file falls back to the original
parser, which keeps its old behaviour except for fix 3.

A backup is saved next to the file before anything is changed, and the new code is syntax-checked first.
Running it twice is safe.
"""
import sys
from pathlib import Path

TARGET = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("research/scripts/build_200gene_panel.py")
MARK = "# --- dbPTM entry-name layout (patched) ---"

NEW_CODE = r'''# --- dbPTM entry-name layout (patched) ---
_ACC_RE = re.compile(r"^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})(-[0-9]+)?$")
_ENTRY_RE = re.compile(r"^[A-Z0-9]+_[A-Z0-9]+$")


def _sniff_dbptm_layout(path: Path, n: int = 5000) -> str:
    """'entry6' only if >=95% of the first n lines are: entry_name, UniProt accession, integer position, ..."""
    hits = total = 0
    with path.open(encoding="utf-8", errors="replace", newline="") as fh:
        for i, line in enumerate(fh):
            if i >= n:
                break
            f = line.rstrip("\r\n").split("\t")
            total += 1
            if len(f) == 6 and _ENTRY_RE.match(f[0]) and _ACC_RE.match(f[1]) and f[2].isdigit():
                hits += 1
    return "entry6" if total and hits / total >= 0.95 else "unknown"


def _load_ptm_entry6() -> pd.DataFrame:
    """Real dbPTM layout: entry_name, accession, position, ptm_type, pmids, flanking_sequence (21 aa window)."""
    from collections import Counter
    n = Counter()
    rows = []
    with PTM_SOURCE.open(encoding="utf-8", errors="replace", newline="") as fh:
        for line in fh:
            n["lines read"] += 1
            f = line.rstrip("\r\n").split("\t")
            if len(f) != 6:
                n["dropped: not 6 tab-separated fields"] += 1
                continue
            entry, acc, pos, ptm, pmid, flank = (x.strip() for x in f)
            if not entry.endswith("_HUMAN"):
                n["dropped: not Homo sapiens (entry name does not end in _HUMAN)"] += 1
                continue
            if not _ACC_RE.match(acc):
                n["dropped: invalid UniProt accession"] += 1
                continue
            if not pos.isdigit() or int(pos) < 1:
                n["dropped: invalid position"] += 1
                continue
            rows.append((acc, int(pos), ptm, pmid, flank))

    out = pd.DataFrame(rows, columns=["uniprot", "position", "ptm_type", "pubmed_id", "flanking_sequence"])
    n["human rows kept so far"] = len(out)
    n["  of which isoform-suffixed accessions (e.g. P12345-2)"] = int(out["uniprot"].str.contains("-").sum()) if len(out) else 0

    base = out["uniprot"].str.split("-").str[0]
    accessions = sorted(set(base))
    print(f"Resolving gene symbols for {len(accessions)} distinct human UniProt accessions via UniProt REST (75 per request)...")
    mapping = map_uniprot_to_gene(accessions)
    out["gene"] = base.map(mapping).fillna("").astype(str).str.upper()
    no_gene = out["gene"] == ""
    n["dropped: no gene symbol returned by UniProt"] = int(no_gene.sum())
    out = out[~no_gene]
    bad_fmt = ~out["gene"].str.match(r"^[A-Z0-9][A-Z0-9-]*$", na=False)
    n["dropped: gene symbol format"] = int(bad_fmt.sum())
    out = out[~bad_fmt].copy()

    # Internal consistency check, computed from the file itself: the flanking window is 21 residues with the
    # modified residue in the middle (index 10); for phosphorylation it should be S, T or Y.
    fl = out["flanking_sequence"]
    has21 = (fl.str.len() == 21) & (out["ptm_type"] == "Phosphorylation")
    out["flank_center_ok"] = pd.Series(pd.NA, index=out.index, dtype="boolean")
    out.loc[has21, "flank_center_ok"] = fl[has21].str[10].isin(list("STY"))
    n["flank-window check: rows with a 21-residue phospho window"] = int(has21.sum())
    n["flank-window check: centre residue is S/T/Y"] = int((out["flank_center_ok"] == True).sum())
    n["flank-window check: centre residue is NOT S/T/Y (kept, flagged flank_center_ok=False)"] = int((out["flank_center_ok"] == False).sum())

    def _join(series):
        return ";".join(sorted({p for v in series for p in str(v).split(";") if p}))

    out = (out.groupby(["gene", "uniprot", "position", "ptm_type"], as_index=False, sort=True)
              .agg(flanking_sequence=("flanking_sequence", "first"),
                   pubmed_id=("pubmed_id", _join),
                   flank_center_ok=("flank_center_ok", "first")))
    out["source"] = "dbPTM experimental PTM download"
    n["final unique (gene, accession, position, type) records"] = len(out)
    n["final distinct genes"] = int(out["gene"].nunique())

    print("dbPTM parsing report")
    for k, v in n.items():
        print(f"  {k}: {v:,}")
    if out.empty:
        raise RuntimeError(f"No usable human PTM records after parsing {PTM_SOURCE}. See the counts above.")
    return out[["gene", "uniprot", "position", "ptm_type", "flanking_sequence", "source", "pubmed_id", "flank_center_ok"]]


def load_ptm() -> pd.DataFrame:
    if _sniff_dbptm_layout(PTM_SOURCE) == "entry6":
        return _load_ptm_entry6()
    print("dbPTM layout not recognised as the 6-column entry-name format; using the original parser.")
    return _load_ptm_legacy()
'''


def main() -> int:
    if not TARGET.exists():
        print(f"Cannot find {TARGET}. Run this from the project root, or pass the path to build_200gene_panel.py.")
        return 2
    raw = TARGET.read_bytes().decode("utf-8")
    crlf = "\r\n" in raw
    text = raw.replace("\r\n", "\n")
    if MARK in text:
        print("Already patched. Nothing changed.")
        return 0
    a, b = "def load_ptm() -> pd.DataFrame:", "def load_drivers() -> pd.DataFrame:"
    if text.count(a) != 1 or text.count(b) != 1:
        print("The file does not look like the expected version (could not find load_ptm / load_drivers once each). Nothing changed.")
        return 3
    i, j = text.index(a), text.index(b)
    old = text[i:j]
    if old.count("    out = pd.DataFrame()\n") != 1:
        print("Unexpected load_ptm body (could not find the single 'out = pd.DataFrame()' line). Nothing changed.")
        return 3
    old = old.replace(a, "def _load_ptm_legacy() -> pd.DataFrame:")
    old = old.replace("    out = pd.DataFrame()\n", "    out = pd.DataFrame(index=df.index)\n")
    new = text[:i] + NEW_CODE + "\n\n" + old + text[j:]
    try:
        compile(new, str(TARGET), "exec")
    except SyntaxError as e:
        print(f"Patched code failed the syntax check ({e}). Nothing changed.")
        return 4
    backup = TARGET.with_name(TARGET.name + ".bak_before_dbptm_patch")
    if not backup.exists():
        backup.write_bytes(TARGET.read_bytes())
    TARGET.write_bytes((new.replace("\n", "\r\n") if crlf else new).encode("utf-8"))
    print(f"Patched {TARGET}\nBackup: {backup}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())