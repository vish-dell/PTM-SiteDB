"""3D-proximal, sequence-distant cancer mutations near PTM residues, plus statistics helpers.

Question: among missense-mutated residues that are FAR from every verified PTM site in sequence (>= min_sep),
are more of them within `threshold` Angstrom of a PTM residue in the 3D structure than expected by chance?

Definitions (all configurable, defaults unchanged: 10 A, 20 residues, pLDDT >= 70, C-alpha distance):
  confident residue  : has a C-alpha and pLDDT >= plddt_min in the AlphaFold model (labelled PREDICTED).
  verified PTM       : the PTM record's residue letter equals the model residue at that UniProt position and
                       the residue is confident. Unverified PTMs are excluded and listed with the reason.
  eligible residue   : confident, not itself a PTM site, and >= min_sep from EVERY verified PTM in sequence.
                       Only eligible residues can ever be a hit, so they are the null universe (N).
  qualifying residue : eligible and within `threshold` of its nearest verified PTM (distance_mode "ca" =
                       C-alpha to C-alpha; "min_atom" = smallest distance between any two atoms).
  mutated position   : distinct missense position whose reference residue matches the model residue and that is
                       confident. Mutations on PTM residues or within min_sep of one are sequence-proximal by
                       definition; they are counted and reported but are not part of this test.
Tests: exact hypergeometric tail on (N eligible, K qualifying, n mutated-eligible, k observed), plus an
amino-acid-composition-matched permutation null drawn ONLY from eligible residues. Neither null models
trinucleotide context, replication timing or clonality, so p-values are screening evidence, not causal inference.
A statistic that cannot be validly estimated is returned as None with a specific *_reason; it is never 0 or 1.
"""
import math
import random
from collections import Counter

import numpy as np

from . import external as ext
from . import structure as S

DEFAULT_PLDDT_MIN = 70.0
MAX_PERM = 20000
DIST_MODES = ("ca", "min_atom")
SYNONYMOUS_TYPES = ("synonymous", "silent")     # cBioPortal/MAF label silent changes "Silent"


def bh_fdr(pvals):
    """Benjamini-Hochberg q-values, returned in the input order."""
    n = len(pvals)
    order = sorted(range(n), key=lambda i: pvals[i])
    q, prev = [0.0] * n, 1.0
    for rank in range(n, 0, -1):
        i = order[rank - 1]
        prev = min(prev, pvals[i] * n / rank)
        q[i] = prev
    return q


def hyper_sf(k, N, K, n):
    """P(X >= k), X ~ Hypergeometric(population N, K successes, n draws); exact, log space."""
    if not (0 <= K <= N and 0 <= n <= N):
        raise ValueError("invalid hypergeometric parameters (need 0 <= K, n <= N)")
    lo, hi = max(0, n - (N - K)), min(K, n)        # support of X; terms outside it are exactly 0
    if k <= lo:
        return 1.0
    if k > hi:
        return 0.0
    lc = lambda a, b: math.lgamma(a + 1) - math.lgamma(b + 1) - math.lgamma(a - b + 1)
    tot = lc(N, n)
    return min(1.0, sum(math.exp(lc(K, i) + lc(N - K, n - i) - tot) for i in range(k, hi + 1)))


def wilson_interval(k, n, z=1.959963984540054):
    """Approximate 95% Wilson interval for an observed proportion."""
    if n <= 0:
        return None
    phat = k / n
    den = 1 + z * z / n
    centre = (phat + z * z / (2 * n)) / den
    half = z * math.sqrt((phat * (1 - phat) / n) + z * z / (4 * n * n)) / den
    return max(0.0, centre - half), min(1.0, centre + half)


def enrichment_summary(k, n, K, N):
    """Observed/expected enrichment plus an approximate CI derived from the Wilson interval.

    The CI treats the observed mutation-position indicator as binomial and keeps the
    structural neighbourhood fraction K/N fixed. It is an effect-size summary, not a
    replacement for the exact/permutation hypothesis test.
    """
    if not (N > 0 and n > 0):
        return {"fold_enrichment": None, "fold_ci95": None}
    expected_rate = K / N
    observed_rate = k / n
    fold = (observed_rate / expected_rate) if expected_rate > 0 else None
    ci = wilson_interval(k, n)
    # A Wilson interval from fewer than 3 mutated residues is bounded by 1/expected and uninformative
    fold_ci = [ci[0] / expected_rate, ci[1] / expected_rate] if ci and expected_rate > 0 and n >= 3 else None
    return {"fold_enrichment": round(fold, 4) if fold is not None else None,
            "fold_ci95": [round(x, 4) for x in fold_ci] if fold_ci else None}


def _valid_positions(items):
    """Keep only integer positions (cBioPortal can return None / -1 / strings for odd records)."""
    out = set()
    for p in items:
        try:
            p = int(p)
        except (TypeError, ValueError):
            continue
        if p > 0:
            out.add(p)
    return out


def _atom_arrays(residues):
    """{pos: Nx3 array of heavy-atom coordinates} for min-atom distances (hydrogens skipped)."""
    out = {}
    for pos, r in residues.items():
        pts = [(x, y, z) for (name, x, y, z) in r.get("atoms", [])
               if name.lstrip("0123456789")[:1] != "H"]
        if pts:
            out[pos] = np.asarray(pts, dtype=float)
    return out


def _dist_to_set(pos, targets, ca, atoms, mode):
    """(distance, target_pos) from residue `pos` to the nearest residue in `targets`, or (None, None)."""
    best = (None, None)
    if mode == "min_atom":
        a = atoms.get(pos) if atoms else None
        if a is None:
            return best
        for t in targets:
            b = atoms.get(t)
            if b is None:
                continue
            d = float(np.sqrt(((a[:, None, :] - b[None, :, :]) ** 2).sum(-1)).min())
            if best[0] is None or d < best[0]:
                best = (d, t)
        return best
    xyz = ca.get(pos)
    for t in targets:
        if t in ca and xyz is not None:
            d = math.dist(xyz, ca[t])
            if best[0] is None or d < best[0]:
                best = (d, t)
    return best


def classify(ca, plddt, ptm_pos, threshold=10.0, min_sep=20, atoms=None, mode="ca",
             plddt_min=DEFAULT_PLDDT_MIN):
    """Per confident residue: nearest PTM in 3D and in sequence.

    `eligible` = not a PTM residue and >= min_sep from every PTM in sequence (can be a hit at all).
    `qualifies` = eligible AND nearest-PTM distance <= threshold. PTM residues and residues inside the
    sequence window carry eligible=False so that they are never part of the null universe.
    """
    if mode not in DIST_MODES:
        raise ValueError(f"distance mode must be one of {DIST_MODES}")
    conf = {p for p in ca if plddt.get(p, 0) >= plddt_min}
    ptm = sorted({p for p in ptm_pos if p in conf})
    out = {}
    if not ptm:
        return out
    ptm_set = set(ptm)
    for pos in sorted(conf):
        d, nearest = _dist_to_set(pos, ptm, ca, atoms, mode)
        if d is None:
            continue
        sep = min(abs(pos - p) for p in ptm)
        eligible = pos not in ptm_set and sep >= min_sep
        out[pos] = {"nearest_ptm_pos": nearest, "dist": round(d, 2), "min_seq_sep": sep,
                    "is_ptm": pos in ptm_set, "eligible": eligible,
                    "qualifies": bool(eligible and d <= threshold)}
    return out


def _perm_hits(sequence, universe, mutated, target, n_perm, seed):
    """Composition-matched permutation hit counts.

    Each permutation re-places the mutated positions at random UNIVERSE residues while keeping, per
    reference amino-acid letter, the observed number of distinct mutated positions. Returns
    (hit_count_array, None) or (None, reason). `target` is the set of residues that count as a hit.
    """
    n_perm = max(0, min(int(n_perm or 0), MAX_PERM))
    if n_perm == 0:
        return None, "permutation test disabled (0 permutations requested)"
    if not sequence:
        return None, "no usable protein sequence for amino-acid-composition matching"
    by_aa = {}
    for pos in sorted(universe):
        aa = _aa_at(sequence, pos)
        if aa:
            by_aa.setdefault(aa, []).append(pos)
    need = Counter(_aa_at(sequence, p) for p in mutated if _aa_at(sequence, p))
    if not need:
        return None, "no mutated position has a known reference residue"
    short = {aa: (c, len(by_aa.get(aa, []))) for aa, c in need.items() if len(by_aa.get(aa, [])) < c}
    if short:
        det = ", ".join(f"{aa}: {c} mutated vs {m} eligible" for aa, (c, m) in sorted(short.items()))
        return None, f"composition-matched null impossible, too few eligible residues of a mutated type ({det})"
    rng = np.random.default_rng(seed)
    hits = np.zeros(n_perm, dtype=np.int64)
    for aa, count in need.items():
        cand = by_aa[aa]
        in_target = np.fromiter((c in target for c in cand), dtype=bool, count=len(cand))
        if count >= len(cand):
            hits += int(in_target.sum())
            continue
        r = rng.random((n_perm, len(cand)))
        idx = np.argpartition(r, count - 1, axis=1)[:, :count]
        hits += in_target[idx].sum(axis=1)
    return hits, None


def _empirical_p(hits, k):
    return float((1 + int((hits >= k).sum())) / (len(hits) + 1))


def proximity_test(cls, mut_positions, sequence=None, n_perm=5000, perm_seed=1):
    """Exact + permutation test on distinct mutated positions within the ELIGIBLE residues.

    Returns None (never 0 or 1) for any statistic that cannot be validly estimated, each with a *_reason.
    """
    pool = {p: v for p, v in cls.items() if v.get("eligible", True)}
    valid = _valid_positions(mut_positions)
    mset = sorted(p for p in valid if p in pool)
    n_seq_proximal = sum(1 for p in valid if p in cls and not cls[p].get("eligible", True))
    n_no_coord = sum(1 for p in valid if p not in cls)
    N, K, n = len(pool), sum(v["qualifies"] for v in pool.values()), len(mset)
    k = sum(pool[p]["qualifies"] for p in mset)

    note = None
    if not cls:
        note = "No residue has confident coordinates near a confident, verified PTM residue."
    elif not N:
        note = ("Untestable: no eligible residue. Every confident residue is a PTM site or lies within the "
                "minimum sequence separation of one. Try a smaller sequence separation.")
    elif not K:
        note = ("Untestable: no eligible residue qualifies. Every eligible residue is farther than the distance "
                "cutoff from all PTM residues. Try a larger distance cutoff.")
    elif not n:
        note = "Untestable: no missense mutation falls on an eligible (confident, sequence-distant) residue."

    p_perm, perm_note, n_permutations = None, None, 0
    if N and K and n:
        hits, why = _perm_hits(sequence, pool.keys(), mset, {p for p, v in pool.items() if v["qualifies"]},
                               n_perm, perm_seed)
        if hits is None:
            perm_note = f"Permutation test not estimable: {why}."
        else:
            p_perm, n_permutations = _empirical_p(hits, k), len(hits)
    else:
        perm_note = "Permutation test not estimable: " + (note or "insufficient eligible residues or mutations") .rstrip(".") + "."

    effect = enrichment_summary(k, n, K, N)
    fold_reason = None
    if effect["fold_enrichment"] is None:
        fold_reason = ("no mutated eligible position" if not n else
                       "no qualifying residue (expected rate is zero)" if not K else "no eligible residue")
    return {"n_confident_residues": len(cls), "n_eligible_residues": N, "n_qualifying_residues": K,
            "n_mutated_positions": n, "n_observed": k,
            "n_mutated_sequence_proximal": n_seq_proximal, "n_mutated_no_coordinates": n_no_coord,
            "expected": round(n * K / N, 3) if N and n else None,
            "p_value": hyper_sf(k, N, K, n) if (N and K and n) else None,
            "p_value_reason": None if (N and K and n) else (note or "not estimable"),
            "p_perm": p_perm, "n_permutations": n_permutations,
            **effect, "fold_reason": fold_reason,
            "ci_reason": None if effect["fold_ci95"] else (fold_reason or "not estimable"),
            "ci_method": "Wilson interval on the observed qualifying fraction, expected fraction held fixed (approximate)",
            "note": note, "p_perm_note": perm_note}


def _aa_at(seq, pos):
    """Return one-letter amino acid at 1-based UniProt position, or None."""
    if not seq or not (1 <= pos <= len(seq)):
        return None
    aa = seq[pos - 1].upper()
    return None if aa == "X" else aa


def composition_matched_positions(seq, positions, pool_positions, n_perm=5000, seed=1):
    """Generate composition-matched random residue sets.

    Each permutation preserves the observed count of distinct mutated positions for each
    reference amino-acid letter, while sampling without replacement from confident residues.
    Only residues in pool_positions are eligible. If a stratum has fewer candidates than
    observed positions, the permutation is invalid and no permutations are returned.
    """
    rng = random.Random(seed)
    pool = sorted(set(pool_positions))
    by_aa = {}
    for pos in pool:
        aa = _aa_at(seq, pos)
        if aa:
            by_aa.setdefault(aa, []).append(pos)
    need = Counter(_aa_at(seq, p) for p in set(positions) if _aa_at(seq, p))
    if not need:
        return []
    if any(len(by_aa.get(aa, [])) < count for aa, count in need.items()):
        return []
    perms = []
    for _ in range(n_perm):
        chosen = []
        for aa, count in need.items():
            chosen.extend(rng.sample(by_aa[aa], count))
        perms.append(set(chosen))
    return perms


def empirical_p_from_permutations(perms, observed, score_fn=None):
    """One-sided empirical p-value with +1 correction."""
    if not perms:
        return None
    if score_fn is None:
        ge = sum(len(p) >= observed for p in perms)
    else:
        ge = sum(score_fn(p) >= observed for p in perms)
    return (1 + ge) / (len(perms) + 1)


def site_neighbourhoods(ca, plddt, units, mut_positions, radius=10.0, min_sep=20,
                        sequence=None, n_perm=5000, seed=1, atoms=None, mode="ca",
                        plddt_min=DEFAULT_PLDDT_MIN):
    """PTM-centric 3D-neighbourhood tests.

    For one PTM site the null universe is every confident residue that is >= min_sep from THAT site in
    sequence (the only residues that could ever be in its neighbourhood). The neighbourhood is the subset
    within `radius`. Mutated positions are counted only inside that universe, so mutations the neighbourhood
    definition excludes by construction cannot dilute the expectation.
    """
    conf = {q for q in ca if plddt.get(q, 0) >= plddt_min}
    M = {q for q in _valid_positions(mut_positions) if q in conf}
    out, seen = [], set()
    for u in units:
        pos = u["pos"]
        if pos not in conf or (u["ptm_type"], pos) in seen:
            continue
        seen.add((u["ptm_type"], pos))
        dists = {}
        for q in conf:
            if q == pos:
                continue
            d, _ = _dist_to_set(pos, [q], ca, atoms, mode)
            if d is not None:
                dists[q] = d
        near = [q for q, d in dists.items() if d <= radius]
        universe = {q for q in dists if abs(q - pos) >= min_sep}
        nb = {q for q in near if q in universe}
        mU = M & universe
        N, K, n = len(universe), len(nb), len(mU)
        hits = sorted(mU & nb)
        k = len(hits)
        exact_p = hyper_sf(k, N, K, n) if (N and K and n) else None
        p_reason = None
        if exact_p is None:
            p_reason = ("no eligible residue" if not N else "empty neighbourhood (no sequence-distant residue within the radius)"
                        if not K else "no mutated position in the eligible residues")
        perm_p, perm_note = None, None
        if N and K and n:
            ph, why = _perm_hits(sequence, universe, mU, nb, n_perm, seed)
            if ph is None:
                perm_note = f"Permutation test not estimable: {why}."
            else:
                perm_p = _empirical_p(ph, k)
        else:
            perm_note = f"Permutation test not estimable: {p_reason}."
        effect = enrichment_summary(k, n, K, N)
        out.append({"site": u["site"], "ptm_type": u["ptm_type"], "pos": pos, "n_contacts": len(near),
                    "n_universe": N, "n_neighbourhood": K, "n_mutated_in_universe": n,
                    "n_mutated_in_neighbourhood": k,
                    "expected": round(n * K / N, 3) if N and n else None,
                    "p_value": exact_p, "p_value_reason": p_reason, "p_perm": perm_p,
                    "n_permutations": (max(0, min(int(n_perm or 0), MAX_PERM)) if perm_p is not None else 0),
                    **effect, "hit_positions": hits, "p_perm_note": perm_note})
    tested = [o for o in out if o["p_value"] is not None]
    for o in out:
        o["q_bh_protein"] = None
        o["q_perm_protein"] = None
    for o, q in zip(tested, bh_fdr([o["p_value"] for o in tested])):
        o["q_bh_protein"] = q
    perm_tested = [o for o in out if o["p_perm"] is not None]
    for o, q in zip(perm_tested, bh_fdr([o["p_perm"] for o in perm_tested])):
        o["q_perm_protein"] = q
    return out


def fisher_combine(pvals):
    """Fisher's method for INDEPENDENT tests (e.g. different cancer types); exact chi-square tail, 2k df."""
    ps = [min(max(p, 1e-300), 1.0) for p in pvals]
    x = -sum(math.log(p) for p in ps)
    return min(1.0, math.exp(-x) * sum(x ** i / math.factorial(i) for i in range(len(ps))))


def _sequence_from_residues(rs):
    """Model sequence in UniProt numbering ('X' where the model has no residue)."""
    if not rs:
        return ""
    seq = ["X"] * max(rs)
    for pos, r in rs.items():
        if r.get("aa"):
            seq[pos - 1] = r["aa"]
    return "".join(seq)


def _ref_letter(m):
    """Reference residue of a cBioPortal mutation record (explicit `wt`, else parsed from proteinChange)."""
    wt = m.get("wt")
    if wt:
        return wt
    import re
    mm = re.match(r"^(?:p\.)?([A-Z])\d+", (m.get("protein_change") or "").strip())
    return mm.group(1) if mm else None


def _map_mutations(muts, rs, plddt, plddt_min, kinds):
    """Verify each mutation against the structure. Returns (verified Counter by pos, records, exclusions).

    A record is kept only if its protein position exists in the model, the reference residue letter equals
    the model residue (so numbering/isoform agree) and the residue is confident. Every drop has a reason.
    """
    excl, keep, by_pos = Counter(), [], Counter()
    for m in muts:
        t = (m.get("type") or "unknown")
        tl = t.lower()
        if not any(tl.startswith(k) if k == "missense" else k in tl for k in kinds):
            excl[f"not {kinds[0]} (type {t})"] += 1
            continue
        pos = m.get("pos")
        if pos is None or not isinstance(pos, int) or pos < 1:
            excl["no valid protein position"] += 1
            continue
        if m.get("pos_end") not in (None, pos):
            excl["multi-residue change (not a single-residue substitution)"] += 1
            continue
        r = rs.get(pos)
        if r is None or not r.get("ca"):
            excl["position outside the structure model or without coordinates"] += 1
            continue
        wt = _ref_letter(m)
        if wt and r.get("aa") and wt != r["aa"]:
            excl["reference residue disagrees with the model (numbering/isoform mismatch)"] += 1
            continue
        if plddt.get(pos, 0) < plddt_min:
            excl[f"low confidence (pLDDT < {plddt_min:g})"] += 1
            continue
        by_pos[pos] += 1
        keep.append(m)
    return by_pos, keep, dict(excl)


def proximity_analysis(acc, units, cbio, threshold=10.0, min_sep=20, n_perm=5000, perm_seed=1,
                       dist_mode="ca", plddt_min=DEFAULT_PLDDT_MIN):
    if len(units) < 1:
        return {"available": False, "reason": "No PTM residues recorded for this protein."}
    if dist_mode not in DIST_MODES:
        return {"available": False, "reason": f"Unknown distance mode '{dist_mode}'. Use one of {DIST_MODES}."}
    entry, text = ext.alphafold_model(acc)
    res = S.parse_pdb(text)
    chains = {}
    for (ch, num), r in res.items():
        chains.setdefault(ch, {})[num] = r
    if not chains:
        return {"available": False, "reason": "The AlphaFold file contained no parseable residues."}
    chain = "A" if "A" in chains else max(chains, key=lambda c: sum(1 for r in chains[c].values() if r["ca"]))
    rs = chains[chain]
    sequence = _sequence_from_residues(rs)
    ca = {n: r["ca"] for n, r in rs.items() if r["ca"]}
    plddt = {n: r["b"] for n, r in rs.items() if r["b"] is not None}
    atoms = _atom_arrays(rs) if dist_mode == "min_atom" else None

    # --- PTM verification: only residues whose letter matches the model and that are confident are used
    verified, ptm_excl = [], []
    for u in units:
        pos, aa, r = u["pos"], u.get("aa"), rs.get(u["pos"])
        if r is None or not r["ca"]:
            ptm_excl.append({"site": u["site"], "reason": "position outside the model or without coordinates"})
        elif aa and r.get("aa") and r["aa"] != aa:
            ptm_excl.append({"site": u["site"], "reason": f"record says {aa} but the model has {r['aa']} at {pos}"})
        elif plddt.get(pos, 0) < plddt_min:
            ptm_excl.append({"site": u["site"], "reason": f"low confidence (pLDDT {plddt.get(pos, 0):.0f} < {plddt_min:g})"})
        else:
            verified.append(u)
    ptm_pos = sorted({u["pos"] for u in verified})
    cls = classify(ca, plddt, ptm_pos, threshold, min_sep, atoms=atoms, mode=dist_mode, plddt_min=plddt_min)

    # --- mutations: retained = genuine single-residue missense that map onto the structure
    all_muts = cbio.get("mutations") or []
    miss_all = [m for m in all_muts if (m.get("type") or "").lower().startswith("missense") and m.get("pos") is not None]
    by_pos, kept, mut_excl = _map_mutations(all_muts, rs, plddt, plddt_min, ("missense",))
    stats = proximity_test(cls, list(by_pos), sequence=sequence, n_perm=n_perm, perm_seed=perm_seed)

    # Synonymous/silent calibration control (real records from the same response; not assumed neutral).
    syn_all = [m for m in all_muts if any(k in (m.get("type") or "").lower() for k in SYNONYMOUS_TYPES)
               and m.get("pos") is not None]
    syn_by_pos, _, syn_excl = _map_mutations(all_muts, rs, plddt, plddt_min, SYNONYMOUS_TYPES)
    syn_stats = proximity_test(cls, list(syn_by_pos), sequence=sequence, n_perm=n_perm, perm_seed=perm_seed)
    syn_note = None
    if not syn_all:
        syn_note = ("Control not estimable: the study returned no synonymous/silent records for this gene "
                    "(cBioPortal studies commonly omit them at import).")

    n_ptm_conf = len(ptm_pos)
    if n_ptm_conf == 0:
        stats["note"] = ("Untestable: no recorded PTM residue is both present in the AlphaFold model with a matching "
                         f"residue letter and confident (pLDDT >= {plddt_min:g}); 3D proximity cannot be assessed.")
        stats["p_value_reason"] = stats["note"]
    elif not all_muts:
        stats["note"] = ("No mutations: cBioPortal answered correctly but returned zero mutation records for this gene "
                         "in the selected study. This is a property of the cohort, not a structural problem and not "
                         "evidence of no effect; the analysis is not estimable. Try a larger study.")
        stats["p_value_reason"] = stats["note"]
    elif not miss_all:
        stats["note"] = ("No missense mutations: the study has mutation records for this gene, but none is a missense "
                         "change with a valid protein position.")
        stats["p_value_reason"] = stats["note"]
    elif not by_pos:
        stats["note"] = ("No missense mutation survived the mapping checks (see exclusion reasons below); "
                         "nothing could be tested.")
        stats["p_value_reason"] = stats["note"]

    n_mis_mismatch = sum(v for k, v in mut_excl.items() if "disagrees" in k)
    mismatch_warning = None
    if miss_all and n_mis_mismatch / len(miss_all) > 0.10:
        mismatch_warning = (f"{n_mis_mismatch} of {len(miss_all)} missense records ({n_mis_mismatch / len(miss_all):.0%}) "
                            "disagree with the model residue; the cBioPortal transcript probably differs from the "
                            "UniProt canonical sequence used by AlphaFold. Treat results with caution.")

    sites = site_neighbourhoods(ca, plddt, verified, list(by_pos), threshold, min_sep, sequence=sequence,
                                n_perm=n_perm, seed=perm_seed, atoms=atoms, mode=dist_mode, plddt_min=plddt_min)
    hits = sorted(({"pos": p, "residue": rs[p].get("aa"), "n_mutations": by_pos[p],
                    "protein_change": sorted({m["protein_change"] for m in kept
                                              if m["pos"] == p and m.get("protein_change")}),
                    "plddt": round(plddt[p], 1), "mapping_verified": True, **cls[p]}
                   for p in by_pos if p in cls and cls[p]["qualifies"]),
                  key=lambda h: (h["dist"], -h["n_mutations"]))
    mean_plddt = round(sum(plddt.values()) / len(plddt), 1) if plddt else None
    dist_def = ("C-alpha to C-alpha distance to the nearest verified PTM residue" if dist_mode == "ca" else
                "smallest distance between any heavy-atom pair of the residue and the nearest verified PTM residue")
    return {"available": True, "structure": f"AlphaFold {entry} (PREDICTED; residues with pLDDT >= {plddt_min:g} only)",
            "structure_source": "alphafold", "structure_label": "AlphaFold PREDICTED model (not experimental)",
            "structure_id": entry, "chain": chain, "mean_plddt": mean_plddt, "plddt_min": plddt_min,
            "dist_mode": dist_mode, "distance_definition": dist_def, "threshold": threshold, "min_sep": min_sep,
            "seed": perm_seed, "n_perm_requested": max(0, min(int(n_perm or 0), MAX_PERM)),
            "n_missense": len(miss_all), "n_missense_mapped": len(kept), "n_positions_all": len(by_pos),
            "n_synonymous": len(syn_all), "n_synonymous_positions": len(syn_by_pos),
            "negative_control": syn_stats, "negative_control_note": syn_note,
            "n_ptm_confident": n_ptm_conf, "n_ptm_records": len({(u["ptm_type"], u["pos"]) for u in units}),
            "ptm_exclusions": ptm_excl,
            "exclusions": {"missense_mutations": mut_excl, "synonymous_mutations": syn_excl},
            "mapping_warning": mismatch_warning,
            "sites": sites,
            "highlight": {"chain": chain,
                          "ptm": [{"site": u["site"], "pos": u["pos"], "resi": u["pos"], "aa": rs[u["pos"]].get("aa"),
                                   "plddt": round(plddt[u["pos"]], 1)} for u in verified],
                          "mutations": [{"pos": h["pos"], "resi": h["pos"], "aa": h["residue"],
                                         "change": ", ".join(h["protein_change"]), "dist": h["dist"]} for h in hits]},
            "cbio_check": {**cbio.get("check", {}), "n_sequenced": cbio.get("n_sequenced")},
            "retrieval": cbio.get("retrieval"), **stats, "hits": hits,
            "method": (f"AlphaFold C-alpha model in UniProt numbering (chain {chain}), residues with pLDDT >= {plddt_min:g}. "
                       f"Distance = {dist_def}. A residue qualifies if it is within {threshold} A of a verified PTM residue "
                       f"and >= {min_sep} residues from every verified PTM residue in sequence. Missense mutations from the "
                       "selected cBioPortal study are kept only if their reference residue matches the model, then collapsed "
                       "to distinct positions. The null universe is the eligible residues only (confident, not PTM, >= "
                       "min_sep from all PTMs). Exact hypergeometric tail plus an amino-acid-composition-matched permutation "
                       f"null (requested {max(0, min(int(n_perm or 0), MAX_PERM))}, seed {perm_seed}) drawn from eligible "
                       "residues. Neither null models trinucleotide context, clonality or mutational signatures; distance "
                       "enrichment alone does not show functional effect.")}


def composition_null(seq, positions, window, n_perm=1000, seed=1):
    """Family-wise null for hotspots: re-place modified sites at random residues of the SAME amino-acid
    letter (keeps S/T/Y/K composition, removing the bias of a uniform null), record the max window count."""
    rng = random.Random(seed)
    by_aa = {}
    for i, ch in enumerate(seq, 1):
        by_aa.setdefault(ch, []).append(i)
    need = Counter(seq[p - 1] for p in positions if 1 <= p <= len(seq))
    maxes = []
    for _ in range(n_perm):
        pts = sorted(q for aa, c in need.items() for q in rng.sample(by_aa[aa], min(c, len(by_aa[aa]))))
        best, j = 0, 0
        for i in range(len(pts)):
            while pts[i] - pts[j] >= window:
                j += 1
            best = max(best, i - j + 1)
        maxes.append(best)
    return maxes


def perm_p(maxes, k):
    return (1 + sum(m >= k for m in maxes)) / (len(maxes) + 1)