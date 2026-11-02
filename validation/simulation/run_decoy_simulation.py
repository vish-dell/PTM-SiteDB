"""Synthetic test of the PTM-specificity (decoy) control. ALL DATA HERE ARE SIMULATED.

Three tests on the same simulated genes:
  naive : hypergeometric test on distinct mutated residues (ignores mutation contexts)
  sbs96 : production p_context_null (SBS96-matched opportunity permutation)
  decoy : production p_decoy (same-gene, same-amino-acid non-PTM decoy anchors)

Scenarios (relative_risk = RR):
  null             RR = 1 everywhere
  generic_surface  surface residues mutate RR x more; PTM anchors are random S/T/Y (no PTM effect)
  ptm_specific     residues near the true PTM anchors mutate RR x more

The production function setting_gene_result() is called unmodified.
Run from the project root:   python validation/simulation/run_decoy_simulation.py
"""
import argparse, csv, math, sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import hypergeom

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "research" / "scripts"))
import annotate_and_analyze_maf as A      # production code, unmodified

OUT = Path(__file__).resolve().parent
ap = argparse.ArgumentParser()
ap.add_argument("--quick", action="store_true", help="small run to check that everything works")
ap.add_argument("--seed", type=int, default=20261101)
ap.add_argument("--genes", type=int, default=None)
ap.add_argument("--reps", type=int, default=None)
ap.add_argument("--only-null", action="store_true")
args = ap.parse_args()

CFG = dict(N=300, radius=30.0, n_anchors=6, sep=8, dist=10.0, n_perm=500, n_decoy=200,
           genes=25, reps_per_gene=16, n_mut=[20, 40], rr=[3.0, 6.0], alpha=0.05, surface_cut=0.75)
if args.quick:
    CFG.update(genes=3, reps_per_gene=4, n_perm=200, n_decoy=60, n_mut=[30], rr=[4.0])

if args.genes: CFG["genes"] = args.genes
if args.reps: CFG["reps_per_gene"] = args.reps
OTHER_AA = list("ACDEFGHIKLMNPQRVW")
CONTEXTS = [f"ctx{i:02d}" for i in range(96)]


def wilson(k, n, z=1.959963984540054):
    if n == 0:
        return (None, None)
    p = k / n; d = 1 + z * z / n; c = p + z * z / (2 * n)
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return (c - h) / d, (c + h) / d


def make_gene(rng, name):
    N = CFG["N"]; R = CFG["radius"]
    v = rng.normal(size=(N, 3)); v /= np.linalg.norm(v, axis=1, keepdims=True)
    r = R * rng.random(N) ** (1 / 3); X = v * r[:, None]; rn = r / R
    pos = np.arange(1, N + 1)
    p_sty = 0.02 + 0.30 * rn ** 2
    is_sty = rng.random(N) < p_sty
    aa = np.array([rng.choice(OTHER_AA) for _ in range(N)])
    sty_choice = rng.choice(list("STY"), size=N, p=[0.5, 0.35, 0.15])
    aa[is_sty] = sty_choice[is_sty]
    sty_idx = np.where(is_sty)[0]
    anchors = np.sort(rng.choice(sty_idx, size=CFG["n_anchors"], replace=False))
    coords = {(name, int(p)): tuple(x) for p, x in zip(pos, X)}
    ctx_w = np.exp(rng.normal(0, 0.7, size=len(CONTEXTS)))
    rows = []
    for i in range(N):
        for c in rng.choice(len(CONTEXTS), size=3, replace=False):
            rows.append((name, int(pos[i]), CONTEXTS[c], str(aa[i]), float(ctx_w[c])))
    opp = pd.DataFrame(rows, columns=["gene", "protein_pos", "context96", "ref_aa", "w"])
    near_anchor = A.near_mask(pos, pos[anchors], coords, name, CFG["sep"], CFG["dist"])
    return dict(name=name, pos=pos, coords=coords, opp=opp, anchors=pos[anchors], rn=rn, near_anchor=near_anchor)


def draw_mutations(rng, g, scenario, rr, n_mut):
    opp = g["opp"]; idx = opp["protein_pos"].to_numpy() - 1
    mult = np.ones(len(g["pos"]))
    if scenario == "generic_surface":
        mult[g["rn"] >= CFG["surface_cut"]] = rr
    elif scenario == "ptm_specific":
        mult[g["near_anchor"]] = rr
    rate = opp["w"].to_numpy() * mult[idx]
    pick = rng.choice(len(opp), size=n_mut, replace=False, p=rate / rate.sum())
    mut = opp.iloc[pick][["gene", "protein_pos", "context96"]].copy().reset_index(drop=True)
    return mut


def run_tests(g, mut, seed):
    r = A.setting_gene_result(g["name"], list(g["anchors"]), mut, g["opp"][["gene", "protein_pos", "context96", "ref_aa"]],
                              g["coords"], CFG["sep"], CFG["dist"], CFG["n_perm"], seed, True,
                              recurrence="distinct_residue", n_decoy=CFG["n_decoy"])
    p_sbs = r.get("p_context_null"); p_dec = r.get("p_decoy") if r.get("decoy_status") == "ok" else None
    mpos = np.unique(mut["protein_pos"].to_numpy())
    n_d = len(mpos); K = int(g["near_anchor"].sum()); N = len(g["pos"])
    obs = int(g["near_anchor"][mpos - 1].sum())
    p_naive = float(hypergeom.sf(obs - 1, N, K, n_d)) if K > 0 and n_d > 0 else None
    return p_naive, p_sbs, p_dec


def main():
    ss = np.random.SeedSequence(args.seed); rng_genes = np.random.default_rng(ss.spawn(1)[0])
    genes = [make_gene(rng_genes, f"SIM{i:02d}") for i in range(CFG["genes"])]
    cells = [("null", 1.0)] + ([] if args.only_null else [(s, rr) for s in ("generic_surface", "ptm_specific") for rr in CFG["rr"]])
    rows = []
    for scenario, rr in cells:
        for n_mut in CFG["n_mut"]:
            ps = {"naive": [], "sbs96": [], "decoy": []}; rng = np.random.default_rng(ss.spawn(1)[0])
            for g in genes:
                for rep in range(CFG["reps_per_gene"]):
                    mut = draw_mutations(rng, g, scenario, rr, n_mut)
                    pn, ps_, pd_ = run_tests(g, mut, int(rng.integers(1, 2 ** 31)))
                    ps["naive"].append(pn); ps["sbs96"].append(ps_); ps["decoy"].append(pd_)
            for test, vals in ps.items():
                est = [v for v in vals if v is not None]; k = sum(v <= CFG["alpha"] for v in est)
                lo, hi = wilson(k, len(est))
                rows.append(dict(scenario=scenario, relative_risk=rr, n_mutations=n_mut, test=test, replicates=len(vals),
                                 n_estimable=len(est), rate_at_0_05=(k / len(est) if est else None), ci_low=lo, ci_high=hi))
            print(f"done {scenario} RR={rr} n={n_mut}", flush=True)
    name = "decoy_simulation_quick.csv" if args.quick else ("decoy_simulation_nullcheck.csv" if args.only_null else "decoy_simulation.csv")
    with open(OUT / name, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    df = pd.DataFrame(rows)
    print("\nRate of p <= 0.05 (95% Wilson CI). SYNTHETIC DATA.")
    for (sc, rr, n), sub in df.groupby(["scenario", "relative_risk", "n_mutations"], sort=False):
        print(f"\n{sc}  RR={rr}  n_mutations={n}")
        for _, r in sub.iterrows():
            ci = f"[{r.ci_low:.3f}, {r.ci_high:.3f}]" if r.rate_at_0_05 is not None and r.ci_low is not None else "n/a"
            rate = f"{r.rate_at_0_05:.3f}" if r.rate_at_0_05 is not None else "n/a"
            print(f"  {r.test:6s} rate={rate} {ci}  estimable={r.n_estimable}/{r.replicates}")
    print("\nSaved:", OUT / name)


if __name__ == "__main__":
    main()
