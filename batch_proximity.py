"""Batch 3D-proximal mutation-PTM screen over every protein in the local database (real data only).

    python batch_proximity.py brca_tcga luad_tcga --sep 8 12 20 --dist 6 8 10 12 --min-units 3 --n-perm 5000

Per study x min-separation it writes:
  results/proximity_<study>_sep<N>_dist<D>.csv  protein-level exact test, BH across proteins
  results/sites_<study>_sep<N>_dist<D>.csv     PTM-centric 3D-neighbourhood test per site, BH across ALL sites
With 2+ studies it also writes results/pooled_sep<N>_dist<D>.csv (Fisher's method per protein, BH across proteins).
Permutation columns (p_perm/q_perm) use an amino-acid-composition-matched null and are reported alongside
rather than replacing the exact hypergeometric baseline.
Pool only INDEPENDENT cohorts (different cancer types; never overlapping TCGA releases of the same cancer).
Failures and untestable cases appear in 'status', never skipped silently.
"""
import argparse, csv, sys, time
from pathlib import Path

BASE = Path(__file__).resolve().parent
sys.path.insert(0, str(BASE))
from services import analysis as A, external as X, proximity as P
from services.dbutil import db


def proteins(min_units):
    conn = db()
    rows = conn.execute("SELECT * FROM ptm_sites ORDER BY uniprot").fetchall()
    conn.close()
    by = {}
    for r in rows:
        by.setdefault(r["uniprot"], []).append(r)
    for acc, rs in by.items():
        units = A.units_from_rows(rs)
        if len(units) >= min_units:
            yield acc, rs[0]["gene"], units


def run(study, sep, dist, min_units, pause, n_perm, perm_seed, dist_mode="ca", plddt_min=P.DEFAULT_PLDDT_MIN):
    out, site_rows = [], []
    for acc, gene, units in proteins(min_units):
        row = {"gene": gene, "uniprot": acc, "n_ptm_units": len(units), "status": "ok"}
        try:
            cb = X.cbio_gene_mutations(study, gene)
            r = P.proximity_analysis(acc, units, cb, dist, sep, n_perm=n_perm, perm_seed=perm_seed,
                                     dist_mode=dist_mode, plddt_min=plddt_min)
            if not r.get("available"):
                row["status"] = r.get("reason", "unavailable")
            else:
                ck = r["cbio_check"]
                rt = r.get("retrieval") or {}
                row.update(retrieval_status=rt.get("status"), profile_id=rt.get("profile_id"),
                           sample_list_id=rt.get("sample_list_id"), dist_mode=dist_mode, plddt_min=plddt_min,
                           n_records_raw=ck.get("n_records_raw"), n_samples_mutated=ck.get("n_distinct_samples"),
                           n_sequenced=ck.get("n_sequenced"), n_missense=r["n_missense"], n_positions_all=r["n_positions_all"],
                           n_ptm_confident=r["n_ptm_confident"], n_missense_mapped=r.get("n_missense_mapped"),
                           n_eligible_residues=r["n_eligible_residues"], n_qualifying_residues=r["n_qualifying_residues"],
                           n_mutated_sequence_proximal=r.get("n_mutated_sequence_proximal"),
                           fold_enrichment=r.get("fold_enrichment"), fold_ci95=r.get("fold_ci95"),
                           p_value_reason=r.get("p_value_reason"),
                           n_synonymous=r.get("n_synonymous", 0), n_synonymous_positions=r.get("n_synonymous_positions", 0),
                           synonymous_observed=(r.get("negative_control") or {}).get("n_observed"),
                           synonymous_expected=(r.get("negative_control") or {}).get("expected"),
                           synonymous_p_value=(r.get("negative_control") or {}).get("p_value"),
                           n_positions=r["n_mutated_positions"], observed=r["n_observed"],
                           expected=r["expected"], p_value=r["p_value"], p_perm=r.get("p_perm"),
                           p_perm_note=r.get("p_perm_note"), n_permutations=r.get("n_permutations", 0),
                           context_null="unavailable: cBioPortal protein-level records do not provide sufficient genomic context for a trinucleotide-matched null",
                           analysis_note=r.get("note"),
                           hit_positions=";".join(str(h["pos"]) for h in r["hits"]))
                if r["p_value"] is None:
                    row["status"] = r["note"]
                for x in r["sites"]:
                    site_rows.append({"gene": gene, "uniprot": acc, **x,
                                      "hit_positions": ";".join(map(str, x["hit_positions"]))})
        except X.MutationRetrievalError as e:
            # unknown, NOT "zero mutations": the p-value columns stay empty and the status says which step failed
            row["status"] = f"not retrieved ({e.status}): {e}"
            row["retrieval_status"] = e.status
        except X.ExternalError as e:
            row["status"] = f"external error: {e}"
        out.append(row)
        print(f"  {gene:8s} {row['status']}", flush=True)
        time.sleep(pause)
    tested = [r for r in out if r.get("p_value") is not None]
    for r, q in zip(tested, P.bh_fdr([r["p_value"] for r in tested])):
        r["q_bh"] = q
    pt = [r for r in out if r.get("p_perm") is not None]
    for r, q in zip(pt, P.bh_fdr([r["p_perm"] for r in pt])):
        r["q_perm"] = q
    st = [r for r in site_rows if r["p_value"] is not None]
    for r, q in zip(st, P.bh_fdr([r["p_value"] for r in st])):
        r["q_bh_all_sites"] = q
    stp = [r for r in site_rows if r.get("p_perm") is not None]
    for r, q in zip(stp, P.bh_fdr([r["p_perm"] for r in stp])):
        r["q_perm_all_sites"] = q
    return out, site_rows


def write(path, rows, cols, key):
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(sorted(rows, key=lambda r: 2 if r.get(key) is None else r[key]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("studies", nargs="+")
    ap.add_argument("--sep", type=int, nargs="+", default=[8, 12, 20])
    ap.add_argument("--dist", type=float, nargs="+", default=[10.0],
                    help="One or more distance cutoffs in Angstrom (definition set by --dist-mode); several values give a sensitivity sweep.")
    ap.add_argument("--dist-mode", choices=list(P.DIST_MODES), default="ca",
                    help="ca = C-alpha to C-alpha (default); min_atom = closest atom pair (use ~4-8 A cutoffs).")
    ap.add_argument("--plddt-min", type=float, default=P.DEFAULT_PLDDT_MIN, help="Minimum pLDDT for a residue (default 70).")
    ap.add_argument("--min-units", type=int, default=3)
    ap.add_argument("--pause", type=float, default=1.0)
    ap.add_argument("--n-perm", type=int, default=5000,
                    help="Number of amino-acid-composition-matched permutations (default: 5000).")
    ap.add_argument("--perm-seed", type=int, default=1,
                    help="Random seed for reproducible permutation tests (default: 1).")
    a = ap.parse_args()
    Path("results").mkdir(exist_ok=True)
    pc = ["gene", "uniprot", "n_ptm_units", "n_records_raw", "n_samples_mutated", "n_sequenced", "n_missense", "n_positions_all", "n_ptm_confident", "n_missense_mapped", "n_eligible_residues", "n_qualifying_residues", "n_mutated_sequence_proximal", "n_synonymous", "n_synonymous_positions", "synonymous_observed", "synonymous_expected", "synonymous_p_value",
          "n_positions", "observed", "expected", "fold_enrichment", "fold_ci95", "p_value", "p_value_reason", "q_bh", "p_perm", "q_perm", "n_permutations",
          "p_perm_note", "retrieval_status", "profile_id", "sample_list_id", "dist_mode", "plddt_min", "context_null", "analysis_note", "hit_positions", "status"]
    sc = ["gene", "uniprot", "site", "ptm_type", "pos", "n_contacts", "n_universe", "n_neighbourhood", "n_mutated_in_universe", "n_mutated_in_neighbourhood",
          "expected", "fold_enrichment", "fold_ci95", "p_value", "p_value_reason", "q_bh_protein", "q_bh_all_sites", "p_perm", "q_perm_protein", "q_perm_all_sites",
          "n_permutations", "p_perm_note", "hit_positions"]
    pooled = {}
    for study in a.studies:
        for dist in a.dist:
            for sep in a.sep:
                print(f"== {study}  min_sep={sep}  dist<={dist:g} A")
                rows, srows = run(study, sep, dist, a.min_units, a.pause, a.n_perm, a.perm_seed, a.dist_mode, a.plddt_min)
                write(Path("results") / f"proximity_{study}_sep{sep}_dist{dist:g}.csv", rows, pc, "p_value")
                write(Path("results") / f"sites_{study}_sep{sep}_dist{dist:g}.csv", srows, sc, "p_value")
                for r in rows:
                    if r.get("p_value") is not None:
                        pooled.setdefault((sep, dist), {}).setdefault((r["gene"], r["uniprot"]), []).append((study, r["p_value"]))
                sig = [r["gene"] for r in rows if r.get("q_bh", 1) < 0.05]
                psig_prot = [r["gene"] for r in rows if r.get("q_perm", 1) < 0.05]
                ssig = [f"{r['gene']} {r['site']}" for r in srows if r.get("q_bh_all_sites", 1) < 0.05]
                psig = [f"{r['gene']} {r['site']}" for r in srows if r.get("q_perm_all_sites", 1) < 0.05]
                print(f"   exact protein q<0.05: {', '.join(sig) or 'none'} | permutation protein q<0.05: {', '.join(psig_prot) or 'none'} | exact site q<0.05: {', '.join(ssig) or 'none'} | permutation site q<0.05: {', '.join(psig) or 'none'}")
    if len(a.studies) > 1:
        for (sep, dist), d in pooled.items():
            rows = [{"gene": g, "uniprot": u, "n_studies": len(v), "studies": ";".join(s for s, _ in v),
                     "fisher_p": P.fisher_combine([p for _, p in v])} for (g, u), v in d.items() if len(v) >= 2]
            for r, q in zip(rows, P.bh_fdr([r["fisher_p"] for r in rows])):
                r["q_bh"] = q
            write(Path("results") / f"pooled_sep{sep}_dist{dist:g}.csv", rows, ["gene", "uniprot", "n_studies", "studies", "fisher_p", "q_bh"], "fisher_p")
            print(f"== pooled sep={sep} dist={dist:g}: {len(rows)} proteins tested in >=2 studies; q<0.05: "
                  f"{', '.join(r['gene'] for r in rows if r['q_bh'] < 0.05) or 'none'}")


if __name__ == "__main__":
    main()
