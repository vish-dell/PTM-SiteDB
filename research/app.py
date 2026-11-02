from __future__ import annotations

import json
import os
import re
import subprocess
import sys
from pathlib import Path

import pandas as pd
import numpy as np
from flask import Flask, jsonify, redirect, render_template, request, send_file, url_for

try:
    from scipy.stats import fisher_exact, mannwhitneyu
except Exception:
    fisher_exact = None
    mannwhitneyu = None

ROOT = Path(__file__).resolve().parent
# Research scripts (and the root-site upload handler) import `scripts.*`. Make that package resolvable whether this
# app runs standalone or is mounted under /research by the root site. Appended (never inserted) so it cannot shadow
# root-level modules such as `app` or `services`.
if str(ROOT) not in sys.path:
    sys.path.append(str(ROOT))
PROC = ROOT / "data" / "processed"
RAW = ROOT / "data" / "raw"
RESULTS = PROC / "cohort_results"
APP_UPLOADS = ROOT / "instance" / "uploads"
for folder in (PROC, RAW, RESULTS, APP_UPLOADS):
    folder.mkdir(parents=True, exist_ok=True)

PTM_MAPPED = PROC / "ptm_sites_mapped.tsv"
PTM_RAW = PROC / "ptm_sites.tsv"
PANEL_PATH = PROC / "benchmark_panel.csv"
DRIVERS_PATH = PROC / "intogen_drivers.tsv"
OPPORTUNITIES_PATH = PROC / "context_opportunities.tsv"
COORDS_PATH = PROC / "alphafold_ca_coordinates.tsv"
REFERENCE_CONFIG = PROC / "reference_config.json"

app = Flask(__name__)
app.secret_key = os.environ.get("PTMSITEDB_SECRET_KEY", "local-research-app-change-this-key")
app.config["MAX_CONTENT_LENGTH"] = 8 * 1024 * 1024 * 1024


def read_tsv(path: Path) -> pd.DataFrame:
    try:
        return pd.read_csv(path, sep="\t", low_memory=False) if path.exists() else pd.DataFrame()
    except Exception:
        return pd.DataFrame()


def read_json(path: Path) -> dict:
    try:
        return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {}
    except Exception:
        return {}


def read_csv_table(path: Path) -> pd.DataFrame:
    """Read a generated CSV defensively; absent/empty/corrupt files are unavailable, not results."""
    try:
        return pd.read_csv(path, low_memory=False) if path.exists() and path.stat().st_size > 0 else pd.DataFrame()
    except Exception:
        return pd.DataFrame()


def ptm_path() -> Path:
    """Full normalized experimental PTM table for browse/search pages."""
    return PTM_RAW if PTM_RAW.exists() and PTM_RAW.stat().st_size > 0 else PTM_MAPPED


def analysis_ptm_path() -> Path:
    """Structural inference must use sequence-mapped PTM anchors only; never fall back to raw isoform positions."""
    return PTM_MAPPED


def reference_path() -> Path | None:
    cfg = read_json(REFERENCE_CONFIG)
    val = cfg.get("reference_fasta")
    if val and Path(val).exists():
        return Path(val)
    for p in (RAW / "reference.fa", RAW / "reference.fasta", RAW / "GRCh38.fa", RAW / "GRCh38.fasta"):
        if p.exists():
            return p
    return None


def data_status() -> dict:
    # Read precomputed manifests/indexes rather than loading the full million-row PTM
    # table for every page request.
    drivers = read_tsv(DRIVERS_PATH)
    panel = read_csv_table(PANEL_PATH)
    summary = read_json(PROC / "panel_summary.json")
    ptm_summary = read_json(PROC / "ptm_summary.json")
    structure_manifest = read_json(PROC / "structure_manifest.json")
    gene_results = structure_manifest.get("gene_results", {}) if isinstance(structure_manifest.get("gene_results", {}), dict) else {}
    mapped_gene_count = sum(int(v.get("sequence_mapped_ptm_records", 0) or 0) > 0 for v in gene_results.values() if isinstance(v, dict))
    return {
        **summary,
        "ptm_source_file": ptm_path().name if ptm_path().exists() else "not loaded",
        "ptm_records": int(ptm_summary.get("records", summary.get("dbptm_unique_ptm_records", 0)) or 0),
        "ptm_bearing_genes": int(ptm_summary.get("ptm_bearing_genes", summary.get("ptm_bearing_genes", 0)) or 0),
        "ptm_summary_by_type": ptm_summary.get("by_type", {}),
        "mapped_ptm_records": int(structure_manifest.get("mapped_experimental_ptm_records", 0) or 0),
        "mapped_ptm_genes": int(mapped_gene_count),
        "analysis_ptm_file": analysis_ptm_path().name if analysis_ptm_path().exists() else "not prepared",
        "intogen_genes": int(drivers["gene"].nunique()) if not drivers.empty and "gene" in drivers.columns else 0,
        "panel_genes": int(panel["gene"].astype(str).str.upper().nunique()) if not panel.empty and "gene" in panel.columns else 0,
        "driver_panel": int(truth_series(panel["intogen_driver"]).sum()) if not panel.empty and "intogen_driver" in panel.columns else 0,
        "non_driver_panel_genes": int((~truth_series(panel["intogen_driver"])).sum()) if not panel.empty and "intogen_driver" in panel.columns else 0,
        "panel_is_200_plus": bool(int(panel["gene"].astype(str).str.upper().nunique()) >= 200) if not panel.empty and "gene" in panel.columns else False,
        "intogen_release_declared": source_manifest().get("sources", {}).get("IntOGen", {}).get("release_declared_inside_archive"),
        "opportunity_map_ready": OPPORTUNITIES_PATH.exists(),
        "coordinates_ready": COORDS_PATH.exists(),
        "canonical_proteins_ready": (PROC / "canonical_proteins.fasta").exists(),
        "reference_fasta_ready": bool(reference_path() and Path(str(reference_path()) + ".fai").exists()),
        "source_manifest": source_manifest(),
    }


def read_tsv_head(path: Path, nrows: int = 300) -> pd.DataFrame:
    try:
        return pd.read_csv(path, sep="\t", nrows=nrows, low_memory=False) if path.exists() else pd.DataFrame()
    except Exception:
        return pd.DataFrame()


def search_ptm_records(term: str, *, exact_gene: bool = False, limit: int = 300) -> list[dict]:
    path = ptm_path()
    if not path.exists() or not term:
        return []
    matches = []
    try:
        for chunk in pd.read_csv(path, sep="\t", chunksize=100000, low_memory=False):
            gene = chunk.get("gene", pd.Series("", index=chunk.index)).astype(str)
            if exact_gene:
                mask = gene.str.upper().eq(term.upper())
            else:
                mask = gene.str.contains(term, case=False, na=False, regex=False)
                for col in ("uniprot", "ptm_type"):
                    if col in chunk.columns:
                        mask |= chunk[col].astype(str).str.contains(term, case=False, na=False, regex=False)
            if mask.any():
                matches.extend(chunk.loc[mask].head(max(0, limit-len(matches))).to_dict(orient="records"))
            if len(matches) >= limit:
                break
    except Exception:
        return matches[:limit]
    return matches[:limit]


def source_manifest() -> dict:
    return read_json(PROC / "source_manifest.json")


def available_cohorts() -> list[str]:
    names = [p.name.removesuffix("_ptm_proximity.tsv") for p in RESULTS.glob("*_ptm_proximity.tsv")]
    return sorted(names, key=lambda n: (-(RESULTS / f"{n}_ptm_proximity.tsv").stat().st_mtime, n.lower()))


def truth_series(s: pd.Series) -> pd.Series:
    return s.astype(str).str.lower().isin({"true", "1", "yes"})


def cohort_result(cohort: str) -> Path | None:
    if not re.fullmatch(r"[A-Za-z0-9_-]+", cohort or ""):
        return None
    p = RESULTS / f"{cohort}_ptm_proximity.tsv"
    return p if p.exists() else None


def load_indexed_fasta():
    """Import IndexedFasta by file path so it works whether the app runs on its own or mounted under the root site."""
    import importlib.util
    name = "ptm_research_fasta_utils"
    module = sys.modules.get(name)
    if module is None:
        path = ROOT / "scripts" / "fasta_utils.py"
        if not path.exists():
            raise RuntimeError(f"Included FASTA indexer is missing: {path}")
        spec = importlib.util.spec_from_file_location(name, path)
        module = importlib.util.module_from_spec(spec)
        sys.modules[name] = module          # dataclasses need the module registered before it executes
        spec.loader.exec_module(module)
    return module.IndexedFasta


def benchmark_metrics(cohort: str) -> tuple[dict, list[dict]]:
    p = cohort_result(cohort)
    if not p or not PANEL_PATH.exists():
        return {}, []
    result = read_tsv(p)
    panel = read_csv_table(PANEL_PATH)
    if result.empty or panel.empty or "gene" not in result.columns:
        return {}, []
    panel["gene"] = panel["gene"].astype(str).str.upper()
    result["gene"] = result["gene"].astype(str).str.upper()
    merged = panel.merge(result, on="gene", how="left", suffixes=("_panel", ""))
    if "intogen_driver_panel" in merged.columns:
        driver_col = "intogen_driver_panel"
    elif "intogen_driver" in merged.columns:
        driver_col = "intogen_driver"
    else:
        driver_col = None
    if driver_col is None:
        merged["intogen_driver"] = False
        driver_col = "intogen_driver"
    eligible = truth_series(merged.get("primary_test_eligible", pd.Series(False, index=merged.index)))
    if "status" in merged.columns:
        eligible &= merged["status"].astype(str).isin({"ok", "ok_no_observed_missense_events"})
    p_context = pd.to_numeric(merged.get("p_context_null", pd.Series(np.nan, index=merged.index)), errors="coerce")
    q = pd.to_numeric(merged.get("q_context_null", pd.Series(np.nan, index=merged.index)), errors="coerce")
    eligible &= p_context.notna() & q.notna()
    tested = merged[eligible].copy()
    tested["hit"] = (pd.to_numeric(tested.get("q_context_null", pd.Series(np.nan, index=tested.index)), errors="coerce") < 0.05)
    drivers = tested[truth_series(tested[driver_col])]
    not_listed = tested[~truth_series(tested[driver_col])]
    d_hit = int(drivers["hit"].sum()) if len(drivers) else 0
    n_hit = int(not_listed["hit"].sum()) if len(not_listed) else 0
    odds = fisher_p = None
    if fisher_exact is not None and len(drivers) and len(not_listed):
        table = [[d_hit, len(drivers)-d_hit], [n_hit, len(not_listed)-n_hit]]
        try:
            result_fisher = fisher_exact(table, alternative="greater")
            odds, fisher_p = float(result_fisher.statistic), float(result_fisher.pvalue)
        except Exception:
            pass

    # Primary group-level benchmark: do driver genes show stronger calibrated evidence?
    group_p = group_u = median_difference = cliffs_delta = None
    d_scores = n_scores = []
    if "p_context_null" in tested.columns:
        pvals = pd.to_numeric(tested["p_context_null"], errors="coerce")
        tested = tested.loc[pvals.notna()].copy()
        pvals = pd.to_numeric(tested["p_context_null"], errors="coerce").clip(lower=1e-300, upper=1.0)
        tested["evidence_score"] = -np.log10(pvals)
        is_driver = truth_series(tested[driver_col])
        d_scores = tested.loc[is_driver, "evidence_score"].to_numpy(float)
        n_scores = tested.loc[~is_driver, "evidence_score"].to_numpy(float)
        if mannwhitneyu is not None and len(d_scores) and len(n_scores):
            try:
                test_u = mannwhitneyu(d_scores, n_scores, alternative="greater", method="auto")
                group_u, group_p = float(test_u.statistic), float(test_u.pvalue)
                median_difference = float(np.median(d_scores) - np.median(n_scores))
                pairwise = d_scores[:, None] - n_scores[None, :]
                cliffs_delta = float((np.sum(pairwise > 0) - np.sum(pairwise < 0)) / pairwise.size)
            except Exception:
                pass
    summary_json = read_json(RESULTS / f"{cohort}_analysis_summary.json") or {}
    group_cmp = summary_json.get("primary_group_comparison", {}) if isinstance(summary_json, dict) else {}
    v11 = group_cmp.get("primary_v1_1_decoy_calibrated_stratified_test", {}) or {}
    burden = group_cmp.get("burden_stratified_raw_score_test", {}) or {}
    def _fmt(x):
        return "—" if x is None or (isinstance(x, float) and pd.isna(x)) else (f"{x:.4g}" if isinstance(x, float) else str(x))
    metrics = {
        "v11_status": v11.get("status", "not_run"),
        "v11_p": _fmt(v11.get("one_sided_p_value")),
        "v11_driver_median": _fmt(v11.get("driver_median_decoy_z")),
        "v11_not_listed_median": _fmt(v11.get("not_listed_median_decoy_z")),
        "v11_driver_n": v11.get("driver_genes_with_decoy_score", "—"),
        "v11_not_listed_n": v11.get("not_listed_genes_with_decoy_score", "—"),
        "v11_strata": v11.get("n_informative_strata", "—"),
        "burden_p": _fmt(burden.get("one_sided_p_value")),
        "recurrence_mode": (summary_json.get("primary_settings", {}) or {}).get("recurrence_mode", "—") if isinstance(summary_json, dict) else "—",
        "cohort": cohort,
        "panel_genes": int(len(panel)),
        "testable_genes": int(len(tested)),
        "untestable_or_partial_genes": int(len(merged)-len(tested)),
        "driver_genes_tested": int(len(drivers)),
        "not_listed_genes_tested": int(len(not_listed)),
        "driver_hits": d_hit,
        "not_listed_hits": n_hit,
        "driver_hit_rate": f"{100*d_hit/len(drivers):.1f}%" if len(drivers) else "—",
        "not_listed_hit_rate": f"{100*n_hit/len(not_listed):.1f}%" if len(not_listed) else "—",
        "primary_test_name": "One-sided Mann–Whitney U on -log10(gene-specific context-null p-values)",
        "primary_u": "—" if group_u is None or pd.isna(group_u) else f"{group_u:.4g}",
        "primary_p": "—" if group_p is None or pd.isna(group_p) else f"{group_p:.4g}",
        "primary_p_numeric": group_p,
        "median_evidence_score_difference": "—" if median_difference is None or pd.isna(median_difference) else f"{median_difference:.4g}",
        "cliffs_delta": "—" if cliffs_delta is None or pd.isna(cliffs_delta) else f"{cliffs_delta:.4g}",
        "fisher_odds_ratio": "—" if odds is None or pd.isna(odds) else f"{odds:.4g}",
        "fisher_p_greater": "—" if fisher_p is None or pd.isna(fisher_p) else f"{fisher_p:.4g}",
        "fisher_p_numeric": fisher_p,
        "interpretation": "Primary cohort-level benchmark compares per-gene -log10(p) scores between IntOGen-listed drivers and PTM-bearing genes not listed in that release among genes passing prespecified QC. Fisher's exact test on q<0.05 hit fractions is a secondary descriptive benchmark. Neither comparison establishes causality; IntOGen absence is not a confirmed biological non-driver label.",
    }
    display_cols = [c for c in ["gene", driver_col, "status", "primary_test_eligible", "missense_mutations_mapped", "observed_proximal_mutation_events", "expected_proximal_events_context_null", "fold_enrichment", "p_context_null", "q_context_null", "decoy_z", "p_decoy", "q_decoy", "is_primary_hit"] if c in merged.columns]
    display_table = merged[display_cols].copy()
    if "q_context_null" in display_table.columns:
        display_table = display_table.sort_values("q_context_null", na_position="last")
    rows = display_table.head(100).to_dict(orient="records")
    return metrics, rows


@app.get("/")
def home():
    return render_template("home.html", status=data_status(), manifest=source_manifest(), cohorts=available_cohorts())


@app.get("/search")
def search():
    q = request.args.get("q", "").strip()
    rows = search_ptm_records(q, limit=300) if q else []
    return render_template("search.html", q=q, rows=rows, count=len(rows), status=data_status())


@app.get("/gene/<gene>")
def gene_page(gene: str):
    g = gene.upper()
    panel = read_csv_table(PANEL_PATH)
    rows = search_ptm_records(g, exact_gene=True, limit=1000)
    benchmark = {}
    if not panel.empty and g in set(panel["gene"].astype(str).str.upper()):
        benchmark = panel.loc[panel["gene"].astype(str).str.upper() == g].iloc[0].to_dict()
    gene_results = []
    for p in RESULTS.glob("*_ptm_proximity.tsv"):
        df = read_tsv(p)
        if not df.empty and "gene" in df.columns:
            hits = df[df["gene"].astype(str).str.upper() == g]
            if not hits.empty:
                row = hits.iloc[0].to_dict()
                row["cohort"] = p.name.removesuffix("_ptm_proximity.tsv")
                gene_results.append(row)
    return render_template("gene.html", gene=g, rows=rows, benchmark=benchmark, gene_results=gene_results)


@app.get("/browse")
def browse():
    rows = read_tsv_head(ptm_path(), 300)
    index_path = PROC / "ptm_gene_index.tsv"
    if index_path.exists():
        try:
            gene_idx = pd.read_csv(index_path, sep="\t", usecols=["gene"])
            genes = sorted(gene_idx["gene"].astype(str).str.upper().replace("nan", "").dropna().unique())
        except Exception:
            genes = []
    else:
        genes = []
        if not rows.empty and "gene" in rows.columns:
            genes = sorted(rows["gene"].astype(str).replace("nan", "").dropna().unique())
    return render_template("browse.html", rows=rows.to_dict(orient="records"), genes=genes[:2500], status=data_status())


@app.get("/statistics")
def statistics():
    ptm_summary = read_json(PROC / "ptm_summary.json")
    by_type = ptm_summary.get("by_type", {})
    panel = read_csv_table(PANEL_PATH)
    driver_pct = float((truth_series(panel["intogen_driver"]).mean() * 100)) if not panel.empty and "intogen_driver" in panel.columns else 0.0
    return render_template("statistics.html", status=data_status(), by_type=by_type, driver_pct=driver_pct)


@app.get("/benchmark")
def benchmark():
    cohorts = available_cohorts()
    selected = request.args.get("cohort", "")
    if selected not in cohorts:
        selected = cohorts[0] if cohorts else ""
    metrics, rows = benchmark_metrics(selected) if selected else ({}, [])
    return render_template("benchmark.html", metrics=metrics, rows=rows, cohorts=cohorts, selected=selected)


@app.get("/robustness")
def robustness():
    cohorts = available_cohorts()
    selected = request.args.get("cohort", "all")
    if selected not in cohorts and selected != "all":
        selected = "all" if len(cohorts) > 1 else (cohorts[0] if cohorts else "")
    target_cohorts = cohorts if selected == "all" else ([selected] if selected else [])
    matrices = []
    all_rows = []
    for cohort in target_cohorts:
        summary_path = RESULTS / f"{cohort}_robustness_summary.tsv"
        rows = read_tsv(summary_path) if summary_path.exists() else pd.DataFrame()
        if rows.empty:
            continue
        all_rows.extend(rows.to_dict(orient="records"))
        cohort_rows = []
        for sep in (8, 12, 20):
            cells = {int(float(r.distance_max_A)): r._asdict() for r in rows[rows["sequence_separation_min"] == sep].itertuples(index=False)}
            row = {"sep": sep, "cells": []}
            for dist in (6, 8, 10, 12):
                cell = cells.get(dist, {})
                row["cells"].append({
                    "dist": dist,
                    "hits": int(cell.get("genes_q_lt_005_per_setting", cell.get("genes_q_lt_005", 0)) or 0),
                    "hits_family": int(cell.get("genes_q_lt_005_family_corrected", 0) or 0),
                    "has_family": "genes_q_lt_005_family_corrected" in cell,
                    "tested": int(cell.get("genes_tested", 0) or 0),
                    "primary": bool(cell.get("is_primary_setting", False)),
                    "file": cell.get("result_file", ""),
                })
            cohort_rows.append(row)
        matrices.append({"cohort": cohort, "rows": cohort_rows})
    return render_template("robustness.html", cohorts=cohorts, selected=selected, matrices=matrices, rows=all_rows)


@app.get("/analysis")
def analysis():
    tables = []
    for p in sorted(RESULTS.glob("*_ptm_proximity.tsv"), key=lambda x: -x.stat().st_mtime):
        df = read_tsv(p)
        summary = read_json(RESULTS / p.name.replace("_ptm_proximity.tsv", "_analysis_summary.json"))
        if df.empty:
            tables.append({"cohort": p.stem.replace("_ptm_proximity", ""), "error": "Result file is empty", "summary": summary, "top": []})
            continue
        q = pd.to_numeric(df.get("q_context_null", pd.Series(np.nan, index=df.index)), errors="coerce")
        eligible = truth_series(df.get("primary_test_eligible", pd.Series(False, index=df.index)))
        df_display = df.sort_values(["primary_test_eligible", "q_context_null"], ascending=[False, True], na_position="last") if "primary_test_eligible" in df else df
        rows = df_display.head(100).to_dict(orient="records")
        sites_path = RESULTS / p.name.replace("_ptm_proximity.tsv", "_ptm_sites.tsv")
        sites = read_tsv(sites_path)
        s127 = sites[(sites.get("gene", pd.Series(dtype=str)).astype(str).str.upper() == "TP53") & (pd.to_numeric(sites.get("position", pd.Series(dtype=float)), errors="coerce") == 127)] if not sites.empty else pd.DataFrame()
        tables.append({
            "cohort": p.stem.replace("_ptm_proximity", ""), "rows": len(df),
            "q_lt_005": int(((q < 0.05) & eligible).sum()), "tested": int(eligible.sum()),
            "top": rows, "summary": summary,
            "s127": s127.to_dict(orient="records"), "site_hits": int((pd.to_numeric(sites.get("q_context_null", pd.Series(dtype=float)), errors="coerce") < 0.05).sum()) if not sites.empty else 0,
        })
    return render_template("analysis.html", tables=tables, status=data_status())


@app.route("/upload-analysis", methods=["GET", "POST"])
def upload_analysis():
    message = None
    error = None
    status = data_status()
    if request.method == "POST":
        maf = request.files.get("maf")
        if not maf or not maf.filename:
            error = "Choose a real MAF file first."
        else:
            cohort = request.form.get("cohort", "cohort").strip() or "cohort"
            safe = re.sub(r"[^A-Za-z0-9_-]", "_", cohort).strip("_") or "cohort"
            maf_name = maf.filename.lower()
            if not maf_name.endswith((".maf", ".txt", ".tsv", ".maf.gz", ".tsv.gz", ".txt.gz")):
                error = "MAF file must end in .maf, .txt, .tsv or .gz variant."
            else:
                maf_path = APP_UPLOADS / f"{safe}{'.maf.gz' if maf_name.endswith('.gz') else '.maf'}"
                maf.save(maf_path)
                fasta = request.files.get("fasta")
                opp = request.files.get("opportunities")
                coords_file = request.files.get("structure")
                fasta_path = reference_path()
                opp_path = OPPORTUNITIES_PATH
                coords_path = COORDS_PATH
                try:
                    from scripts.cohort_guard import find_duplicate_input
                    dup = find_duplicate_input(maf_path, safe, RESULTS)
                    if dup:
                        raise RuntimeError(f"This MAF is byte-identical to the one already analysed as '{dup}'. "
                                           f"Two cohort labels on one file are not two cohorts; upload the correct, independent MAF for '{safe}'.")
                    if fasta and fasta.filename:
                        fasta_path = APP_UPLOADS / f"{safe}_reference.fa"
                        fasta.save(fasta_path)
                        try:
                            IndexedFasta = load_indexed_fasta()
                            # IndexedFasta creates the standard .fai index when needed.
                            with IndexedFasta(fasta_path):
                                pass
                        except Exception as exc:
                            raise RuntimeError(f"Could not index the uploaded FASTA: {exc}") from exc
                    if opp and opp.filename:
                        opp_path = APP_UPLOADS / f"{safe}_context_opportunities.tsv"
                        opp.save(opp_path)
                    if coords_file and coords_file.filename:
                        coords_path = APP_UPLOADS / f"{safe}_ca_coordinates.tsv"
                        coords_file.save(coords_path)
                    if not fasta_path or not Path(fasta_path).exists():
                        raise RuntimeError("No reference FASTA configured. Use PREPARE_GENOMIC_INPUTS.ps1 first, or upload the same FASTA build used for the opportunity map.")
                    for label, path in (("context opportunity map", opp_path), ("CA coordinate map", coords_path), ("sequence-mapped PTM file", analysis_ptm_path())):
                        if not Path(path).exists():
                            raise RuntimeError(f"Missing {label}: {path}. Follow the data preparation guide and rerun bootstrap/preparation.")
                    mapped_check = read_tsv(analysis_ptm_path())
                    if mapped_check.empty or not {"gene", "position", "ptm_type"}.issubset(mapped_check.columns):
                        raise RuntimeError("The sequence-mapped PTM table is empty or has an invalid schema. Re-run PREPARE_GENOMIC_INPUTS.ps1; the app will not substitute unprojected raw UniProt coordinates into the structural analysis.")
                    cmd = [
                        sys.executable, str(ROOT / "scripts" / "annotate_and_analyze_maf.py"), str(maf_path), str(fasta_path),
                        "--cohort", safe, "--opportunities", str(opp_path), "--structure", str(coords_path),
                        "--ptm", str(analysis_ptm_path()), "--panel", str(PANEL_PATH),
                    ]
                    nperm = request.form.get("permutations", "5000")
                    if nperm not in {"1000", "5000"}:
                        nperm = "5000"
                    cmd += ["--permutations", nperm]
                    if request.form.get("robustness") == "on":
                        cmd.append("--robustness")
                    result = subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=60 * 60 * 6)
                    summary = read_json(RESULTS / f"{safe}_analysis_summary.json")
                    message = (f"Analysis finished for {safe}: {summary.get('genes_tested_primary', 0)} genes passed primary QC; "
                               f"{summary.get('genes_q_lt_005_primary', 0)} protein-level q-values were below 0.05. "
                               "Review all untestable/partial QC rows before interpretation. Results are not a biological claim until independently validated.")
                    if result.stdout:
                        (RESULTS / f"{safe}_run_log.txt").write_text(result.stdout[-50000:], encoding="utf-8")
                except subprocess.TimeoutExpired:
                    error = "Analysis exceeded the 6-hour local execution limit. Review panel size/permutations or run the CLI directly."
                except subprocess.CalledProcessError as exc:
                    error = (exc.stderr or exc.stdout or "Analysis failed; see local console")[-7000:]
                except Exception as exc:
                    error = str(exc)
    return render_template("upload.html", message=message, error=error, status=status)


@app.get("/related-work")
def related_work():
    return render_template("related_work.html")


@app.get("/s127-evidence")
def s127_evidence():
    return render_template("s127_evidence.html")


@app.get("/methodology")
def methodology():
    return render_template("methodology.html", status=data_status())


@app.get("/download")
def download_page():
    downloadable = []
    fixed_files = [PANEL_PATH, DRIVERS_PATH, PTM_RAW, PTM_MAPPED, OPPORTUNITIES_PATH, COORDS_PATH,
                   PROC / "source_manifest.json", PROC / "panel_summary.json", PROC / "opportunity_manifest.json", PROC / "structure_manifest.json", PROC / "ptm_site_mapping_qc.tsv", PROC / "preparation_qc.json", PROC / "ptm_summary.json", PROC / "ptm_gene_index.tsv"]
    seen = set()
    for p in fixed_files:
        if p.exists() and p.name not in seen:
            downloadable.append({"name": p.name, "size": p.stat().st_size})
            seen.add(p.name)
    for p in sorted(RESULTS.glob("*")):
        if p.is_file() and p.suffix.lower() in {".tsv", ".csv", ".json", ".txt"}:
            downloadable.append({"name": p.name, "size": p.stat().st_size})
    return render_template("download.html", downloadable=downloadable)


@app.get("/download-file/<name>")
def download_file(name: str):
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", name or ""):
        return ("Invalid filename", 400)
    candidates = [PANEL_PATH, DRIVERS_PATH, PTM_MAPPED, PTM_RAW, OPPORTUNITIES_PATH, COORDS_PATH,
                  PROC / "source_manifest.json", PROC / "panel_summary.json", PROC / "opportunity_manifest.json", PROC / "structure_manifest.json", PROC / "ptm_site_mapping_qc.tsv", PROC / "preparation_qc.json", PROC / "ptm_summary.json", PROC / "ptm_gene_index.tsv"]
    candidates += [p for p in RESULTS.iterdir() if p.is_file()]
    for p in candidates:
        if p.exists() and p.name == name:
            return send_file(p, as_attachment=True)
    return ("Not found", 404)


@app.get("/help")
def help_page():
    return render_template("help.html", status=data_status())


@app.get("/api/status")
def api_status():
    return jsonify({"status": data_status(), "manifest": source_manifest(), "cohorts": available_cohorts()})


@app.get("/api/results/<cohort>")
def api_results(cohort: str):
    p = cohort_result(cohort)
    if not p:
        return jsonify({"error": "cohort_not_found"}), 404
    df = read_tsv(p)
    return jsonify({"cohort": cohort, "rows": df.to_dict(orient="records")})


if __name__ == "__main__":
    port = int(os.environ.get("PORT", "5000"))
    # Local-only binding: uploaded genomic data remain on the user's machine.
    app.run(host="127.0.0.1", port=port, debug=False)
