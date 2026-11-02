from __future__ import annotations
import json
from datetime import datetime, timezone
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "data" / "processed"


def read_table(path: Path, sep="\t") -> pd.DataFrame:
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    try:
        return pd.read_csv(path, sep=sep, low_memory=False)
    except Exception:
        return pd.DataFrame()


def main() -> int:
    checks: list[dict] = []
    warnings: list[str] = []
    errors: list[str] = []
    paths = {
        "source_manifest": PROC / "source_manifest.json",
        "ptm_source": PROC / "ptm_sites.tsv",
        "intogen_drivers": PROC / "intogen_drivers.tsv",
        "panel": PROC / "benchmark_panel.csv",
        "opportunities": PROC / "context_opportunities.tsv",
        "coordinates": PROC / "alphafold_ca_coordinates.tsv",
        "mapped_ptm": PROC / "ptm_sites_mapped.tsv",
        "reference_config": PROC / "reference_config.json",
        "opportunity_manifest": PROC / "opportunity_manifest.json",
        "structure_manifest": PROC / "structure_manifest.json",
    }
    for name, path in paths.items():
        checks.append({"item": name, "exists": path.exists(), "bytes": path.stat().st_size if path.exists() else 0})

    ptm = read_table(paths["ptm_source"])
    drivers = read_table(paths["intogen_drivers"])
    panel = pd.read_csv(paths["panel"]) if paths["panel"].exists() else pd.DataFrame()
    opp = read_table(paths["opportunities"])
    coords = read_table(paths["coordinates"])
    mapped = read_table(paths["mapped_ptm"])

    if ptm.empty:
        errors.append("Normalized dbPTM experimental PTM source is missing/empty. Rerun source bootstrap.")
    if drivers.empty:
        errors.append("Normalized IntOGen driver table is missing/empty. Download the official ZIP into data/raw and rerun source bootstrap.")
    if panel.empty:
        errors.append("Benchmark panel is missing/empty. Run scripts/build_200gene_panel.py after source bootstrap.")
    if opp.empty:
        warnings.append("Coding opportunity map is not prepared yet. Run PREPARE_GENOMIC_INPUTS.ps1 with a matched FASTA/GTF before cohort analysis.")
    if coords.empty:
        warnings.append("AlphaFold CA coordinate map is not prepared yet. Run PREPARE_GENOMIC_INPUTS.ps1 and inspect structure mapping QC.")
    if mapped.empty:
        warnings.append("Sequence-mapped PTM anchors are not prepared yet. Structural inference must not use raw UniProt positions as a fallback.")
    if not opp.empty and not {"gene", "protein_pos", "context96", "chromosome", "genomic_position", "ref", "alt"}.issubset(opp.columns):
        errors.append("Opportunity map has an invalid schema.")
    if not coords.empty and not {"gene", "position", "x", "y", "z"}.issubset(coords.columns):
        errors.append("Coordinate map has an invalid schema.")

    panel_summary = {}
    panel_summary_path = PROC / "panel_summary.json"
    if panel_summary_path.exists():
        try: panel_summary = json.loads(panel_summary_path.read_text(encoding="utf-8"))
        except Exception: pass
    panel_n = int(panel["gene"].astype(str).str.upper().nunique()) if not panel.empty and "gene" in panel.columns else 0
    drivers_n = int(panel["intogen_driver"].astype(str).str.lower().isin({"true", "1", "yes"}).sum()) if not panel.empty and "intogen_driver" in panel.columns else 0
    not_listed_n = max(0, panel_n - drivers_n)
    if panel_n < 200:
        warnings.append(f"The observed source intersection yielded {panel_n} panel genes (<200). Do not label this as a completed 200-gene study.")
    if not opp.empty and not coords.empty:
        opp_genes = set(opp["gene"].astype(str).str.upper()) if "gene" in opp.columns else set()
        coord_genes = set(coords["gene"].astype(str).str.upper()) if "gene" in coords.columns else set()
        if len(opp_genes & coord_genes) < 100:
            warnings.append(f"Only {len(opp_genes & coord_genes)} genes overlap the opportunity map and structure coordinates. Structural eligibility is below the aspirational target; inspect manifests.")
    if not mapped.empty and "gene" in mapped.columns and not panel.empty:
        mapped_outside = set(mapped["gene"].astype(str).str.upper()) - set(panel["gene"].astype(str).str.upper())
        if mapped_outside:
            warnings.append(f"{len(mapped_outside)} mapped PTM genes are outside the current panel; verify the panel/mapping source version before analysis.")

    summary = {
        "checked_utc": datetime.now(timezone.utc).isoformat(),
        "source_counts": {
            "normalized_ptm_records": int(len(ptm)),
            "ptm_bearing_genes": int(ptm["gene"].astype(str).str.upper().nunique()) if not ptm.empty and "gene" in ptm.columns else 0,
            "intogen_genes": int(drivers["gene"].astype(str).str.upper().nunique()) if not drivers.empty and "gene" in drivers.columns else 0,
            "panel_genes": panel_n,
            "panel_driver_genes": drivers_n,
            "panel_not_listed_genes": not_listed_n,
            "coding_missense_opportunities": int(len(opp)),
            "opportunity_genes": int(opp["gene"].astype(str).str.upper().nunique()) if not opp.empty and "gene" in opp.columns else 0,
            "coordinate_residues": int(len(coords)),
            "coordinate_genes": int(coords["gene"].astype(str).str.upper().nunique()) if not coords.empty and "gene" in coords.columns else 0,
            "mapped_experimental_ptm_records": int(len(mapped)),
        },
        "panel_summary": panel_summary,
        "checks": checks,
        "warnings": warnings,
        "errors": errors,
        "status": "ERROR" if errors else "WARNING" if warnings else "PASS",
    }
    out = PROC / "preparation_qc.json"
    out.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 2 if errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
