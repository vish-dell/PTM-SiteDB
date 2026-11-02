"""Build RESULTS.md (+ threshold-sensitivity figure) from REAL cohort outputs only.

Reads  data/processed/cohort_results/<cohort>_analysis_summary.json  and  <cohort>_robustness_summary.tsv.
Refuses to write anything if no real cohort output exists, and never fills in missing numbers.
The primary (v1.1) result is reported first, whether or not it is significant.
Usage:  python scripts/make_results_report.py [cohort ...]
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "data" / "processed" / "cohort_results"


def fmt(x, nd=4):
    if x is None or (isinstance(x, float) and x != x):
        return "n/a"
    return f"{x:.{nd}g}" if isinstance(x, float) else str(x)


def main() -> int:
    cohorts = sys.argv[1:] or sorted(p.name.replace("_analysis_summary.json", "") for p in RES.glob("*_analysis_summary.json"))
    if not cohorts:
        print("No real cohort output found in", RES, "- run the analysis first. Nothing written.")
        return 1
    lines = ["# Results (auto-generated from real cohort outputs)", "",
             "Primary endpoint (pre-specified v1.1): driver vs not-listed difference in the decoy-calibrated PTM-specific score,",
             "tested with a mutation-burden x PTM-anchor stratified permutation test. Reported first, including when null.", ""]
    summaries = {}
    for c in cohorts:
        f = RES / f"{c}_analysis_summary.json"
        if not f.exists():
            print("Missing", f)
            return 1
        summaries[c] = json.loads(f.read_text(encoding="utf-8"))
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from cohort_guard import duplicate_groups
    dups = duplicate_groups({c: (s.get("input_provenance") or {}).get("maf_sha256") for c, s in summaries.items()
                             if (s.get("input_provenance") or {}).get("maf_sha256")})
    if dups:
        print("REFUSING to write a report: these cohort labels used a byte-identical MAF (same SHA-256), so they are not "
              "independent cohorts:", dups, "- re-run each cohort with its own MAF.")
        return 2
    lines += ["## 1. Primary result per cohort", "",
              "| cohort | genes tested | drivers | not-listed | driver median decoy-z | not-listed median decoy-z | stratified one-sided p |",
              "|---|---|---|---|---|---|---|"]
    for c, s in summaries.items():
        g = s["primary_group_comparison"]
        v = g.get("primary_v1_1_decoy_calibrated_stratified_test", {})
        lines.append(f"| {c} | {g.get('eligible_gene_count')} | {v.get('driver_genes_with_decoy_score', 'n/a')} | "
                     f"{v.get('not_listed_genes_with_decoy_score', 'n/a')} | {fmt(v.get('driver_median_decoy_z'))} | "
                     f"{fmt(v.get('not_listed_median_decoy_z'))} | {fmt(v.get('one_sided_p_value'))} ({v.get('status')}) |")
    lines += ["", "## 2. Secondary tests", "",
              "| cohort | burden-stratified raw-score p | v1.0 Mann-Whitney p | Cliff's delta | driver hits q<0.05 | not-listed hits q<0.05 |",
              "|---|---|---|---|---|---|"]
    for c, s in summaries.items():
        g = s["primary_group_comparison"]
        lines.append(f"| {c} | {fmt(g.get('burden_stratified_raw_score_test', {}).get('one_sided_p_value'))} | "
                     f"{fmt(g.get('one_sided_p_value'))} | {fmt(g.get('cliffs_delta'))} | "
                     f"{g.get('secondary_driver_q_lt_005_hits')} | {g.get('secondary_not_listed_q_lt_005_hits')} |")
    p_vals = [s["primary_group_comparison"].get("primary_v1_1_decoy_calibrated_stratified_test", {}).get("one_sided_p_value") for s in summaries.values()]
    if len(summaries) >= 2:
        lines += ["", "## 3. Replication", "",
                  "Cohorts must be independent (sample overlap checked outside this tool). Replication is claimed only if the primary",
                  f"result is significant in the discovery cohort AND the same-direction one-sided p is < 0.05 in the second: p-values = {[fmt(p) for p in p_vals]}."]
    lines += ["", "## 4. Threshold sensitivity (descriptive; family-corrected hits shown)", ""]
    frames = []
    for c in summaries:
        rf = RES / f"{c}_robustness_summary.tsv"
        if rf.exists():
            frames.append(pd.read_csv(rf, sep="\t"))
    if frames:
        rob = pd.concat(frames, ignore_index=True)
        lines.append(rob.to_markdown(index=False) if hasattr(rob, "to_markdown") and _has_tabulate() else rob.to_string(index=False))
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            fig, axes = plt.subplots(1, len(summaries), figsize=(5 * len(summaries), 3.6), squeeze=False)
            for ax, c in zip(axes[0], summaries):
                sub = rob[rob["cohort"] == c].pivot(index="sequence_separation_min", columns="distance_max_A", values="genes_q_lt_005_family_corrected")
                im = ax.imshow(sub.values, cmap="Blues", vmin=0)
                ax.set_xticks(range(len(sub.columns)), [f"{x:g}" for x in sub.columns]); ax.set_yticks(range(len(sub.index)), sub.index)
                ax.set_xlabel("CA-CA distance (A)"); ax.set_ylabel("min sequence separation"); ax.set_title(f"{c}: genes q<0.05 (family-corrected)")
                for (i, j), v in __import__("numpy").ndenumerate(sub.values):
                    ax.text(j, i, int(v), ha="center", va="center")
            fig.tight_layout(); out = ROOT / "threshold_sensitivity.png"; fig.savefig(out, dpi=160)
            lines += ["", f"![threshold sensitivity]({out.name})"]
        except Exception as exc:  # matplotlib optional
            lines += ["", f"(Figure not produced: {exc})"]
    else:
        lines.append("Robustness summaries not found; rerun the analysis with --robustness.")
    lines += ["", "## 5. Limits that apply to every number above", ""]
    seen = set()
    for s in summaries.values():
        for l in s.get("limitations", []):
            if l not in seen:
                seen.add(l); lines.append(f"- {l}")
    (ROOT / "RESULTS.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("Wrote", ROOT / "RESULTS.md")
    return 0


def _has_tabulate() -> bool:
    try:
        import tabulate  # noqa: F401
        return True
    except Exception:
        return False


if __name__ == "__main__":
    raise SystemExit(main())
