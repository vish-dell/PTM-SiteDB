"""Synthetic validation of PTM-SiteDB's proximity test. ALL DATA HERE ARE SIMULATED.
Run from the project root:  python validation/simulation/run_simulation.py
"""
import json, sys, csv, math
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[2]; sys.path.insert(0, str(ROOT))
from services import proximity as P          # the production implementation, unmodified
from scipy.stats import kstest

OUT = Path(__file__).resolve().parent; CFG = json.load(open(OUT / "config.json"))
ss = np.random.SeedSequence(CFG["master_seed"])

def wilson(k, n, z=1.959963984540054):
    p = k / n; d = 1 + z*z/n; c = p + z*z/(2*n); h = z*math.sqrt(p*(1-p)/n + z*z/(4*n*n))
    return (c-h)/d, (c+h)/d

def make_cls(N, qual_idx):
    """cls in the format produced by proximity.classify: all N residues eligible; qual_idx qualify."""
    q = set(qual_idx)
    return {i+1: {"eligible": True, "qualifies": i in q} for i in range(N)}

def draw(rng, N, n, w=None):
    return rng.choice(N, size=n, replace=False, p=None if w is None else w/w.sum())

def run_p(cls, idx):
    r = P.proximity_test(cls, [int(i)+1 for i in idx], sequence=None, n_perm=0)
    return r["p_value"]          # None => not estimable (reported, never coerced)

def summarize(ps, alphas):
    est = [p for p in ps if p is not None]; n = len(est)
    row = {"reps": len(ps), "n_estimable": n}
    for a in alphas:
        k = sum(p <= a for p in est); lo, hi = wilson(k, n) if n else (None, None)
        row[f"rate@{a}"] = k/n if n else None; row[f"lo@{a}"] = lo; row[f"hi@{a}"] = hi
    return row, est

def write(name, rows):
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(OUT / name, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys); w.writeheader(); w.writerows(rows)

A = CFG["alphas"]
# ---- Experiment A: null -------------------------------------------------------------------------------
c = CFG["null"]; rows = []; pooled = {}
for N in c["N"]:
    for qf in c["qualifying_fraction"]:
        for n in c["n_mutated"]:
            if n > N: continue
            rng = np.random.default_rng(ss.spawn(1)[0]); K = max(1, round(qf*N))
            cls = make_cls(N, range(K)); ps = []
            for _ in range(c["replicates"]):
                ps.append(run_p(cls, draw(rng, N, n)))           # uniform draw: no association planted
            s, est = summarize(ps, A); s.update(experiment="A_null", N=N, K=K, n=n); rows.append(s)
            pooled[(N, K, n)] = est
write("null_fpr.csv", rows)
# ---- Experiment B: planted signal ---------------------------------------------------------------------
c = CFG["planted"]; N = c["N"]; K = round(c["qualifying_fraction"]*N); rows = []
cls = make_cls(N, range(K))
for rr in c["relative_risk"]:
    for n in c["n_mutated"]:
        rng = np.random.default_rng(ss.spawn(1)[0]); w = np.ones(N); w[:K] = rr; ps = []
        for _ in range(c["replicates"]):
            ps.append(run_p(cls, draw(rng, N, n, w)))
        s, _ = summarize(ps, A); s.update(experiment="B_planted", N=N, K=K, n=n, relative_risk=rr); rows.append(s)
write("power.csv", rows)
# ---- Alternative setting: regional mutation-rate heterogeneity aligned with qualifying residues ---------
c = CFG["regional_confounding"]; N = c["N"]; B = c["block"]; K = round(c["qualifying_fraction"]*N); rows = []
cls = make_cls(N, range(K))                       # qualifying residues form one contiguous block (a 3D patch)
for sg in c["sigma"]:
    for n in c["n_mutated"]:
        rng = np.random.default_rng(ss.spawn(1)[0]); ps = []
        for _ in range(c["replicates"]):
            seg = np.exp(rng.normal(0, sg, size=math.ceil(N/B)))     # NO true proximity effect, only regional rate
            w = np.repeat(seg, B)[:N]; ps.append(run_p(cls, draw(rng, N, n, w)))
        s, _ = summarize(ps, A); s.update(experiment="C_regional_rate_confounding", N=N, K=K, n=n, sigma=sg); rows.append(s)
write("confounding_fpr.csv", rows)
# ---- Permutation vs exact cross-check (uniform amino-acid composition) ----------------------------------
c = CFG["permutation_crosscheck"]; N = c["N"]; K = round(c["qualifying_fraction"]*N)
rng = np.random.default_rng(ss.spawn(1)[0]); cls = make_cls(N, range(K)); ex, pm = [], []
for _ in range(c["replicates"]):
    seq = "".join(rng.choice(list("ACDEFGHIKLMNPQRSTVWY"), size=N)); idx = draw(rng, N, c["n_mutated"])
    r = P.proximity_test(cls, [int(i)+1 for i in idx], sequence=seq, n_perm=c["n_perm"], perm_seed=int(rng.integers(1, 2**31)))
    if r["p_value"] is not None and r["p_perm"] is not None: ex.append(r["p_value"]); pm.append(r["p_perm"])
ex, pm = np.array(ex), np.array(pm)
write("perm_vs_exact.csv", [{"n_pairs": len(ex), "mean_abs_diff": float(np.mean(np.abs(ex-pm))), "corr": float(np.corrcoef(ex, pm)[0,1]),
                              "fpr_exact@0.05": float((ex<=.05).mean()), "fpr_perm@0.05": float((pm<=.05).mean())}])
# ---- Plots -------------------------------------------------------------------------------------------
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
allp = np.concatenate([np.array(v) for v in pooled.values()]); ps_sorted = np.sort(allp)
fig, ax = plt.subplots(1, 2, figsize=(10, 4.3))
ax[0].plot(np.linspace(0, 1, 200), [np.mean(allp <= t) for t in np.linspace(0, 1, 200)], label="empirical P(p<=t)")
ax[0].plot([0, 1], [0, 1], "k--", label="ideal (continuous uniform)"); ax[0].set_xlabel("threshold t"); ax[0].set_ylabel("fraction of null p-values <= t")
ax[0].set_title("Null calibration (SYNTHETIC)\nbelow diagonal = conservative (discrete test)"); ax[0].legend()
nr = list(csv.DictReader(open(OUT / "null_fpr.csv")))
x = [i for i, _ in enumerate(nr)]; y = [float(r["rate@0.05"]) for r in nr]
lo = [max(0.0, float(r["rate@0.05"]) - float(r["lo@0.05"])) for r in nr]; hi = [max(0.0, float(r["hi@0.05"]) - float(r["rate@0.05"])) for r in nr]
ax[1].errorbar(x, y, yerr=[lo, hi], fmt="o", ms=3); ax[1].axhline(0.05, color="k", ls="--"); ax[1].set_ylim(0, 0.08)
ax[1].set_xlabel("null scenario index (N, K, n)"); ax[1].set_ylabel("empirical FPR at alpha=0.05 (95% Wilson CI)"); ax[1].set_title("Type I error (SYNTHETIC)")
plt.tight_layout(); plt.savefig(OUT / "null_calibration.png", dpi=140); plt.close()
pr = list(csv.DictReader(open(OUT / "power.csv"))); fig, ax = plt.subplots(figsize=(6, 4.3))
for n in sorted({int(r["n"]) for r in pr}):
    rr = [r for r in pr if int(r["n"]) == n]
    ax.errorbar([float(r["relative_risk"]) for r in rr], [float(r["rate@0.05"]) for r in rr],
                yerr=[[max(0.0, float(r["rate@0.05"]) - float(r["lo@0.05"])) for r in rr], [max(0.0, float(r["hi@0.05"]) - float(r["rate@0.05"])) for r in rr]],
                marker="o", capsize=2, label=f"n={n} mutated positions")
ax.axhline(0.05, color="k", ls=":"); ax.set_xlabel("planted relative risk (qualifying vs other residues)"); ax.set_ylabel("power at alpha=0.05")
ax.set_title("Power curves (SYNTHETIC; N=300, 10% qualifying)"); ax.legend(); plt.tight_layout(); plt.savefig(OUT / "power_curves.png", dpi=140); plt.close()
ks = kstest(allp, "uniform")
print("pooled null: n_p =", len(allp), " KS vs U(0,1) stat=%.3f (expected to deviate: p-values are discrete)" % ks.statistic)
print("Done. Outputs in", OUT)
