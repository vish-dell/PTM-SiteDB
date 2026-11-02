from __future__ import annotations

import hashlib
import json
import os
import re
import time
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "data" / "processed"
PTM_SOURCE = PROC / "dbptm_phosphorylation_source.txt"
INTOGEN_SOURCE = PROC / "intogen_drivers_source.tsv"
OUT_PTM = PROC / "ptm_sites.tsv"
OUT_DRIVERS = PROC / "intogen_drivers.tsv"
OUT_PANEL = PROC / "benchmark_panel.csv"
GENE_MAP_CACHE = PROC / "uniprot_gene_map.tsv"
# The dbPTM download covers many species; the genome-based analysis is human, so only *_HUMAN entries are used.
SPECIES_SUFFIX = os.environ.get("PTM_SPECIES_SUFFIX", "_HUMAN").upper()
TARGET_TOTAL_GENES = 200
TARGET_PER_GROUP = 100
SEED = 1


def read_flexible(path: Path) -> pd.DataFrame:
    errors = []
    for sep in ("\t", ",", None):
        try:
            df = pd.read_csv(path, sep=sep, engine="python" if sep is None else "c", low_memory=False)
            if df.shape[1] >= 2:
                return df
        except Exception as exc:
            errors.append(str(exc))
    raise RuntimeError(f"Unable to parse {path}: {'; '.join(errors[-2:])}")


def norm(s: object) -> str:
    return re.sub(r"[^a-z0-9]", "", str(s).lower())


def choose_col(df: pd.DataFrame, names: list[str]) -> str | None:
    lookup = {norm(c): c for c in df.columns}
    for n in names:
        if norm(n) in lookup:
            return lookup[norm(n)]
    for c in df.columns:
        nc = norm(c)
        if any(norm(n) in nc for n in names):
            return c
    return None


def map_uniprot_to_gene(accessions: list[str]) -> dict[str, str]:
    """Resolve missing gene symbols through the official UniProtKB REST API."""
    if not accessions:
        return {}
    try:
        import requests
    except Exception:
        return {}
    out: dict[str, str] = {}
    session = requests.Session()
    session.headers.update({"User-Agent": "PTM-SiteDB-Integrated-Explorer/3.0"})
    for i in range(0, len(accessions), 75):
        batch = accessions[i:i+75]
        query = " OR ".join(f"accession:{a}" for a in batch)
        try:
            r = session.get(
                "https://rest.uniprot.org/uniprotkb/search",
                params={"query": query, "format": "tsv", "fields": "accession,gene_names", "size": 500},
                timeout=60,
            )
            r.raise_for_status()
            lines = [x.split("\t") for x in r.text.splitlines() if x.strip()]
            for row in lines[1:]:
                if len(row) >= 2 and row[1].strip():
                    accession = row[0].split(";")[0].strip()
                    gene = row[1].strip().split()[0].upper()
                    if gene:
                        out[accession] = gene
        except Exception as exc:
            print(f"WARNING: UniProt gene-symbol mapping batch {i//75+1} failed: {exc}")
    return out


def _load_gene_cache() -> dict[str, str]:
    if not GENE_MAP_CACHE.exists():
        return {}
    try:
        d = pd.read_csv(GENE_MAP_CACHE, sep="\t", dtype=str).dropna()
        return dict(zip(d["accession"], d["gene"]))
    except Exception:
        return {}


def _save_gene_cache(mapping: dict[str, str]) -> None:
    pd.DataFrame(sorted(mapping.items()), columns=["accession", "gene"]).to_csv(GENE_MAP_CACHE, sep="\t", index=False)


def fetch_reviewed_human_gene_map() -> dict[str, str]:
    """One request for all reviewed human UniProt entries (accession -> primary gene symbol)."""
    import requests
    last = None
    for attempt in range(3):
        try:
            r = requests.get("https://rest.uniprot.org/uniprotkb/stream",
                             params={"query": "(organism_id:9606) AND (reviewed:true)",
                                     "fields": "accession,gene_primary", "format": "tsv"},
                             headers={"User-Agent": "PTM-SiteDB-Integrated-Explorer/3.0"}, timeout=300)
            r.raise_for_status()
            out: dict[str, str] = {}
            for line in r.text.splitlines()[1:]:
                row = line.split("\t")
                if len(row) >= 2 and row[1].strip():
                    gene = row[1].strip().split(";")[0].split()[0].upper()
                    if gene:
                        out[row[0].split(";")[0].strip()] = gene
            return out
        except Exception as exc:  # network problems are reported, never silently turned into empty data
            last = exc
            time.sleep(2 * (attempt + 1))
    print(f"WARNING: bulk UniProt gene-symbol download failed: {last}")
    return {}


def resolve_gene_symbols(accessions: list[str]) -> dict[str, str]:
    """accession -> gene symbol, using a local cache, then one bulk request, then batched lookups."""
    mapping = _load_gene_cache()
    missing = [a for a in accessions if a not in mapping]
    if missing:
        print(f"Resolving gene symbols for {len(missing):,} UniProt accessions (cached: {len(mapping):,})...")
        if len(missing) > 500:
            mapping.update(fetch_reviewed_human_gene_map())
            missing = [a for a in accessions if a not in mapping]
        if missing:
            base = sorted({a.split("-")[0] for a in missing})
            batch = map_uniprot_to_gene(base)
            for a in missing:
                g = batch.get(a.split("-")[0])
                if g:
                    mapping[a] = g
        _save_gene_cache(mapping)
        still = [a for a in accessions if a not in mapping]
        if still:
            print(f"NOTE: {len(still):,} accessions have no gene symbol and are excluded (not guessed).")
    return mapping


def load_ptm() -> pd.DataFrame:
    df = read_flexible(PTM_SOURCE)
    gene_c = choose_col(df, ["gene", "genesymbol", "symbol", "hgnc"])
    uniprot_c = choose_col(df, ["uniprot", "uniprotid", "accession", "protein"])
    pos_c = choose_col(df, ["position", "modifiedposition", "site", "residueposition"])
    type_c = choose_col(df, ["ptmtype", "modification", "ptm"])
    flank_c = choose_col(df, ["sequence", "flankingsequence", "flanking"])

    # dbPTM describes four core tab-delimited fields (UniProt ID, modified
    # position, PTM type and a flanking sequence). Some download builds expose
    # a header; if a release is headerless, re-read by documented column order.
    if uniprot_c is None or pos_c is None:
        raw = None
        for sep in ("\t", None, ","):
            try:
                candidate = pd.read_csv(PTM_SOURCE, sep=sep, engine="python" if sep is None else "c",
                                        header=None, low_memory=False, dtype=str)
                if candidate.shape[1] >= 3:
                    raw = candidate
                    break
            except Exception:
                continue
        if raw is None:
            raise RuntimeError("Could not parse dbPTM table as a headered or headerless tabular file.")
        # Remove a header row if one is actually present.
        first = " ".join(raw.iloc[0].astype(str).str.lower().tolist())
        if "uniprot" in first and ("position" in first or "modified" in first):
            raw = raw.iloc[1:].reset_index(drop=True)
        # Real dbPTM download layout (6 tab-separated columns, no header):
        #   ENTRY_NAME_SPECIES | UniProt accession | position | PTM type | PubMed IDs | flanking sequence
        first_col = raw.iloc[:, 0].astype(str).head(2000)
        entry_like = bool(first_col.str.match(r"^[A-Za-z0-9]+_[A-Za-z0-9]+$").mean() > 0.9)
        if raw.shape[1] >= 6 and entry_like:
            fallback_names = ["entry_name", "uniprot", "position", "ptm_type", "pubmed_id", "flanking_sequence"]
        else:
            fallback_names = ["uniprot", "position", "ptm_type", "flanking_sequence", "pubmed_id"]
        if raw.shape[1] > len(fallback_names):
            fallback_names.extend([f"extra_{i}" for i in range(raw.shape[1]-len(fallback_names))])
        raw.columns = fallback_names[:raw.shape[1]]
        df = raw
        if "entry_name" in df.columns and SPECIES_SUFFIX:
            keep = df["entry_name"].astype(str).str.upper().str.endswith(SPECIES_SUFFIX)
            print(f"dbPTM rows: {len(df):,} total; {int(keep.sum()):,} with species suffix {SPECIES_SUFFIX}.")
            df = df[keep].reset_index(drop=True)
            if df.empty:
                raise RuntimeError(f"No dbPTM rows with species suffix {SPECIES_SUFFIX}.")
        gene_c = None
        uniprot_c = "uniprot"
        pos_c = "position"
        type_c = "ptm_type" if "ptm_type" in df.columns else None
        flank_c = "flanking_sequence" if "flanking_sequence" in df.columns else None

    out = pd.DataFrame(index=df.index)   # index matters: scalar columns assigned to an empty frame stay empty
    out["gene"] = df[gene_c].astype(str).str.strip().str.upper() if gene_c else ""
    out["uniprot"] = df[uniprot_c].astype(str).str.strip()
    # Accept positions encoded either as integer text or site labels like S127.
    out["position"] = pd.to_numeric(df[pos_c].astype(str).str.extract(r"(\d+)", expand=False), errors="coerce")
    # Normalize FASTA-like accession strings while preserving isoform suffixes.
    out["uniprot"] = out["uniprot"].str.replace(r"^(?:sp|tr)\|([^|]+)\|.*$", r"\1", regex=True)
    out["ptm_type"] = df[type_c].astype(str).str.strip() if type_c else "Phosphorylation"
    out["flanking_sequence"] = df[flank_c].astype(str) if flank_c else ""
    out["source"] = "dbPTM experimental PTM download"
    out = out.dropna(subset=["position"])
    out["position"] = out["position"].astype(int)
    missing_gene = out["gene"].isin(["", "nan", "NA", "None"])
    if missing_gene.any():
        accessions = sorted(set(out.loc[missing_gene, "uniprot"]))
        mapping = resolve_gene_symbols(accessions)
        out.loc[missing_gene, "gene"] = out.loc[missing_gene, "uniprot"].map(mapping).fillna("")
    out = out[(out["gene"].str.match(r"^[A-Z0-9][A-Z0-9-]*$", na=False)) & (out["uniprot"] != "") & (out["uniprot"] != "nan") & (out["position"] > 0)]
    out = out.drop_duplicates(subset=["gene", "uniprot", "position", "ptm_type"]).reset_index(drop=True)
    if out.empty:
        raise RuntimeError(f"No usable PTM records after parsing {PTM_SOURCE}. Inspect the source file header/column layout before analysis.")
    return out


def load_drivers() -> pd.DataFrame:
    df = read_flexible(INTOGEN_SOURCE)
    gene_c = choose_col(df, ["symbol", "gene", "hugo_symbol", "genesymbol"])
    if gene_c is None:
        raise RuntimeError(f"Could not identify gene symbol column in IntOGen table. Columns: {list(df.columns)}")
    out = pd.DataFrame({"gene": df[gene_c].astype(str).str.upper().str.strip()})
    out = out[out["gene"].str.match(r"^[A-Z0-9][A-Z0-9-]+$", na=False)].drop_duplicates()
    return out.sort_values("gene").reset_index(drop=True)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    if not PTM_SOURCE.exists() or not INTOGEN_SOURCE.exists():
        raise SystemExit("Run scripts/bootstrap_sources.py first; required real source files are missing.")

    ptm = load_ptm()
    drivers = load_drivers()
    ptm.to_csv(OUT_PTM, sep="\t", index=False)
    drivers.to_csv(OUT_DRIVERS, sep="\t", index=False)
    # Precompute small indexes so the website doesn't parse millions of PTM rows
    # on every page request.
    gene_index = (ptm.groupby("gene", dropna=False)
                  .agg(site_count=("position", "size"), accession_count=("uniprot", "nunique"),
                       ptm_types=("ptm_type", lambda x: ";".join(sorted(set(x.astype(str))))))
                  .reset_index())
    gene_index.to_csv(PROC / "ptm_gene_index.tsv", sep="\t", index=False)
    ptm_by_type = ptm["ptm_type"].value_counts(dropna=False).to_dict()
    (PROC / "ptm_summary.json").write_text(json.dumps({
        "records": int(len(ptm)), "ptm_bearing_genes": int(ptm["gene"].nunique()),
        "by_type": {str(k): int(v) for k, v in ptm_by_type.items()},
        "normalization_source": "dbPTM experimental phosphorylation download",
    }, indent=2), encoding="utf-8")

    ptm_genes = set(ptm["gene"].astype(str).str.upper())
    driver_genes = set(drivers["gene"])
    driver_candidates = sorted(driver_genes & ptm_genes)
    not_listed_candidates = sorted(ptm_genes - driver_genes)

    # Balanced, reproducible panel: target 100 genes/group (200 total) when
    # sources support it. If one group has fewer candidates, use available genes
    # and fill the total from the other group without inventing any records.
    rng_seed = SEED
    n_driver = min(TARGET_PER_GROUP, len(driver_candidates))
    n_not_listed = min(TARGET_PER_GROUP, len(not_listed_candidates))
    if n_driver + n_not_listed < TARGET_TOTAL_GENES:
        extra = TARGET_TOTAL_GENES - (n_driver + n_not_listed)
        extra_driver = min(extra, max(0, len(driver_candidates) - n_driver))
        n_driver += extra_driver
        extra -= extra_driver
        n_not_listed += min(extra, max(0, len(not_listed_candidates) - n_not_listed))
    sampled_drivers = sorted(pd.Series(driver_candidates).sample(n=n_driver, random_state=rng_seed).tolist()) if n_driver else []
    sampled_not_listed = sorted(pd.Series(not_listed_candidates).sample(n=n_not_listed, random_state=rng_seed + 1).tolist()) if n_not_listed else []

    panel = pd.DataFrame([
        {"gene": g, "intogen_driver": True, "panel_group": "intogen_driver_sampled" if len(driver_candidates) > n_driver else "intogen_driver_all_available"}
        for g in sampled_drivers
    ] + [
        {"gene": g, "intogen_driver": False, "panel_group": "ptm_bearing_not_listed_in_intogen_sampled" if len(not_listed_candidates) > n_not_listed else "ptm_bearing_not_listed_in_intogen_all_available"}
        for g in sampled_not_listed
    ])
    panel = panel.drop_duplicates("gene").sort_values(["intogen_driver", "gene"], ascending=[False, True]).reset_index(drop=True)
    panel.to_csv(OUT_PANEL, index=False)

    manifest_path = PROC / "source_manifest.json"
    try:
        source_manifest = json.loads(manifest_path.read_text(encoding="utf-8")) if manifest_path.exists() else {}
    except Exception:
        source_manifest = {}
    intogen_release = source_manifest.get("sources", {}).get("IntOGen", {}).get("release_declared_inside_archive")
    summary = {
        "ptm_source": "dbPTM experimental PTM-site download",
        "ptm_source_sha256": sha256(PTM_SOURCE),
        "intogen_source_sha256": sha256(INTOGEN_SOURCE),
        "panel_random_seed": SEED,
        "target_total_genes": TARGET_TOTAL_GENES,
        "target_genes_per_group": TARGET_PER_GROUP,
        "dbptm_unique_ptm_records": int(len(ptm)),
        "ptm_bearing_genes": int(ptm["gene"].nunique()),
        "intogen_catalog_genes": int(drivers["gene"].nunique()),
        "driver_genes_with_ptm_available": int(len(driver_candidates)),
        "driver_genes_selected": int(len(sampled_drivers)),
        "not_listed_ptm_bearing_genes_available": int(len(not_listed_candidates)),
        "not_listed_in_intogen_panel_genes": int(len(sampled_not_listed)),
        "total_panel_genes": int(panel["gene"].nunique()),
        "intogen_archive_declared_release": intogen_release,
        "panel_is_200_plus": bool(panel["gene"].nunique() >= TARGET_TOTAL_GENES),
        "interpretation": "The panel is a reproducible seeded sample of PTM-bearing IntOGen-listed drivers and PTM-bearing genes not listed in that release, targeting 100 per group and 200 total. Genes absent from IntOGen are not proven biological non-drivers. The actual panel can be smaller if the real source intersection does not contain enough genes.",
    }
    (PROC / "panel_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    if len(panel) < TARGET_TOTAL_GENES:
        print(f"WARNING: only {len(panel)} genes can be assembled from these source files; do not call this a 200-gene panel.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
