"""Build an independent gene-level benchmark from the IntOGen 2024.09.20 driver catalog.

The script accepts either the official IntOGen ZIP already downloaded locally or,
by default, downloads the official CC0 release automatically. It never creates a
benchmark from the project's own hits.

Usage:
  python build_intogen_benchmark.py results/proximity_brca_tcga_sep8_dist10.csv
  python build_intogen_benchmark.py results/proximity_brca_tcga_sep8_dist10.csv path/to/IntOGen-Drivers-20240920.zip

The output benchmark is gene-level. Genes absent from IntOGen are labelled
"negative" only in the narrow sense "not listed as a driver in this release";
they are NOT claimed to be biologically non-drivers.
"""
import argparse, csv, io, re, sys, zipfile
from pathlib import Path
from urllib.request import Request, urlopen

OFFICIAL_URL = "https://www.intogen.org/download?file=IntOGen-Drivers-20240920.zip"
DEFAULT_ZIP = "IntOGen-Drivers-20240920.zip"


def read_project_genes(path):
    with open(path, newline="", encoding="utf-8-sig") as f:
        rows = list(csv.DictReader(f))
    if not rows or "gene" not in rows[0]:
        raise ValueError("Project CSV must contain a 'gene' column.")
    return sorted({(r.get("gene") or "").strip() for r in rows if (r.get("gene") or "").strip()})


def load_zip(path=None):
    if path:
        p = Path(path)
        if not p.exists():
            raise FileNotFoundError(f"IntOGen ZIP not found: {p}")
        return p.read_bytes()
    print("Downloading the official IntOGen 2024.09.20 driver catalog...", flush=True)
    req = Request(OFFICIAL_URL, headers={"User-Agent": "PTM-SiteDB-benchmark/1.0"})
    try:
        with urlopen(req, timeout=60) as r:
            data = r.read()
    except Exception as e:
        raise RuntimeError(
            "Automatic IntOGen download failed. The official site may be unreachable "
            "from this network. Download IntOGen-Drivers-20240920.zip manually from "
            "https://www.intogen.org/download and pass its full path as the second argument."
        ) from e
    if not data.startswith(b"PK"):
        preview = data[:120].decode("utf-8", errors="replace")
        raise RuntimeError("IntOGen returned a non-ZIP response (possibly a login/404 page): " + preview)
    return data


def pick_member(names):
    preferred = [
        n for n in names
        if n.lower().endswith("compendium_cancer_genes.tsv")
    ]
    if preferred:
        return preferred[0]
    candidates = [n for n in names if n.lower().endswith(".tsv") and "driver" in n.lower()]
    if candidates:
        return candidates[0]
    raise ValueError("Could not find Compendium_Cancer_Genes.tsv or a driver TSV in the IntOGen ZIP.\n" + "\n".join(names))


def read_drivers(zip_bytes):
    with zipfile.ZipFile(io.BytesIO(zip_bytes)) as z:
        member = pick_member(z.namelist())
        raw = z.read(member).decode("utf-8-sig")
    reader = csv.DictReader(io.StringIO(raw), delimiter="\t")
    fields = reader.fieldnames or []
    symbol_col = next((c for c in fields if c.strip().upper() in {"SYMBOL", "GENE", "GENE_SYMBOL"}), None)
    if not symbol_col:
        symbol_col = fields[0] if fields else None
    if not symbol_col:
        raise ValueError("Driver TSV has no columns.")
    genes = set()
    for row in reader:
        g = (row.get(symbol_col) or "").strip()
        if g:
            genes.add(g)
    return genes, member


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("project_csv")
    ap.add_argument("intogen_zip", nargs="?", help="Optional local IntOGen-Drivers-20240920.zip; otherwise downloaded from the official URL.")
    ap.add_argument("--out", default="benchmark.csv")
    a = ap.parse_args()

    project = read_project_genes(a.project_csv)
    data = load_zip(a.intogen_zip)
    drivers, member = read_drivers(data)

    rows = []
    for gene in project:
        positive = gene in drivers
        rows.append({
            "gene": gene,
            "uniprot": "",
            "site": "",
            "reference_result": "positive" if positive else "negative",
            "reference_source": (
                "IntOGen driver catalog release 2024.09.20; " +
                ("listed as a driver" if positive else "not listed as a driver in this release; not evidence of non-driver status")
            )
        })

    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["gene","uniprot","site","reference_result","reference_source"])
        w.writeheader(); w.writerows(rows)

    print(f"Created {out}")
    print(f"Project genes tested in benchmark: {len(rows)}")
    print(f"IntOGen positive: {sum(r['reference_result']=='positive' for r in rows)}")
    print(f"Not listed in IntOGen: {sum(r['reference_result']=='negative' for r in rows)}")
    print(f"Source member: {member}")
    print("Untestable project rows are handled by benchmark_report.py and are never converted to negatives.")

if __name__ == "__main__":
    main()
