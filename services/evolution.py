"""Conservation, record-integrity audit, cancer-cohort mutation overlap and site prioritisation."""
import re
from collections import Counter, defaultdict

from . import external as ext
from .analysis import binom_sf


# ------------------------------------------------------------ alignment
def nw_map(a, b, match=2, mismatch=-1, gap=-2):
    """Needleman–Wunsch; returns list mapping index in a -> index in b (or None). Used only when the
    Ensembl source protein differs from the UniProt sequence (isoform difference)."""
    n, m = len(a), len(b)
    prev = [j * gap for j in range(m + 1)]
    tb = [bytearray(m + 1) for _ in range(n + 1)]
    for j in range(1, m + 1):
        tb[0][j] = 2
    for i in range(1, n + 1):
        cur = [i * gap] + [0] * m
        tb[i][0] = 1
        ai = a[i - 1]
        for j in range(1, m + 1):
            d = prev[j - 1] + (match if ai == b[j - 1] else mismatch)
            u, l = prev[j] + gap, cur[j - 1] + gap
            if d >= u and d >= l:
                cur[j], tb[i][j] = d, 0
            elif u >= l:
                cur[j], tb[i][j] = u, 1
            else:
                cur[j], tb[i][j] = l, 2
        prev = cur
    mp = [None] * n
    i, j = n, m
    while i > 0 or j > 0:
        t = tb[i][j]
        if i > 0 and j > 0 and t == 0:
            mp[i - 1] = j - 1
            i, j = i - 1, j - 1
        elif i > 0 and (j == 0 or t == 1):
            i -= 1
        else:
            j -= 1
    return mp


def _uniprot_to_column(uni_seq, src_aligned, memo=None):
    """Map each UniProt residue index -> alignment column of the Ensembl source sequence."""
    src = src_aligned.replace("-", "")
    cols = [c for c, ch in enumerate(src_aligned) if ch != "-"]
    if src == uni_seq:
        return cols, "exact"
    if memo is not None and src in memo:      # same Ensembl source protein for every orthologue: align once
        mp = memo[src]
    else:
        mp = nw_map(uni_seq, src)
        if memo is not None:
            memo[src] = mp
    return [cols[j] if j is not None else None for j in mp], "realigned"


def conservation(acc, gene, units, flank=5):
    uni = ext.uniprot_sequence(acc)
    orths = ext.ensembl_orthologs(gene)
    if not orths:
        raise ext.ExternalError("No orthologues with aligned sequences were returned")
    prepared, memo = [], {}
    for o in orths:
        cmap, how = _uniprot_to_column(uni, o["source_aligned"], memo)
        prepared.append((o, cmap, how))
    sites = []
    for u in units:
        i = u["pos"] - 1
        if not (0 <= i < len(uni)):
            sites.append({"site": u["site"], "ptm_type": u["ptm_type"], "pos": u["pos"], "human_residue": None,
                          "status": "position outside UniProt sequence", "orthologs": [], "n_identical": 0,
                          "n_aligned": 0})
            continue
        human = uni[i]
        status = "ok" if (u["aa"] is None or u["aa"] == human) else f"record says {u['aa']} but UniProt has {human}"
        row = {"site": u["site"], "ptm_type": u["ptm_type"], "pos": u["pos"], "human_residue": human,
               "status": status, "orthologs": [], "n_identical": 0, "n_aligned": 0}
        if status != "ok":
            sites.append(row)
            continue
        for o, cmap, how in prepared:
            col = cmap[i] if i < len(cmap) else None
            tgt = o["target_aligned"]
            if col is None or col >= len(tgt):
                res, ident, wid = None, False, None
            else:
                res = tgt[col]
                ident = res == human
                # flank identity over ±flank human residues, using each residue's own column
                same = tot = 0
                for k in range(max(0, i - flank), min(len(uni), i + flank + 1)):
                    if k == i or cmap[k] is None or cmap[k] >= len(tgt):
                        continue
                    tot += 1
                    same += tgt[cmap[k]] == uni[k]
                wid = round(same / tot, 2) if tot else None
            row["orthologs"].append({"species": o["label"], "orthology": o["orthology_type"], "residue": res,
                                     "identical": ident, "gap": res == "-", "flank_identity": wid,
                                     "protein_id": o["protein_id"], "perc_id": o["perc_id"], "mapping": how})
            if res not in (None, "-"):
                row["n_aligned"] += 1
                row["n_identical"] += ident
        sites.append(row)
    return {"n_orthologs": len(orths), "species": [o["label"] for o in orths], "sites": sites,
            "method": "Ensembl Compara orthologue alignments (real). UniProt residues are placed onto the alignment "
                      "columns of the Ensembl source protein (exact match or Needleman–Wunsch when isoforms differ). "
                      "A site is 'conserved' in a species only if the aligned residue is identical."}


# ---------------------------------------------------------- integrity
def norm_flank(fl):
    """Letters only, upper-case: drops spaces, dashes, underscores, digits and markup remnants."""
    return re.sub(r"[^A-Z]", "", (fl or "").upper())


def flank_match(seq, fl, pos=None):
    """Best ungapped placement of a stored flank on the UniProt sequence ('X' = wildcard).
    found = identity >= 0.9 (exact when 1.0). Prefers a window that covers `pos`."""
    fl = norm_flank(fl)
    n = len(fl)
    if not fl or n > len(seq):
        return None
    best = None
    for st in range(len(seq) - n + 1):
        same = sum(1 for a, b in zip(fl, seq[st:st + n]) if a == "X" or a == b)
        cov = pos is not None and st + 1 <= pos <= st + n
        key = (same, cov)
        if best is None or key > best[0]:
            best = (key, st)
            if same == n and cov:
                break
    (same, cov), st = best
    return {"start": st + 1, "identity": round(same / n, 3), "found": same / n >= 0.9,
            "exact": same == n, "covers_site": cov}


def integrity(acc, rows, catalog_len=None):
    """Compare every stored record against the UniProt sequence."""
    seq = ext.uniprot_sequence(acc)
    out, c = [], Counter()
    for r in rows:
        aa, pos = (r["site"][:1].upper(), r["position"]) if r["site"] and r["site"][0].isalpha() else (None, r["position"])
        res_ok = None
        if pos and 1 <= pos <= len(seq):
            res_ok = (seq[pos - 1] == aa) if aa else None
        elif pos:
            res_ok = False
        fm = flank_match(seq, r["flanking"], pos)
        out.append({"id": r["id"], "site": r["site"], "ptm_type": r["ptm_type"], "pmid": r["pmid"],
                    "residue_matches_uniprot": res_ok,
                    "uniprot_residue": seq[pos - 1] if pos and 1 <= pos <= len(seq) else None,
                    "stored_flank": norm_flank(r["flanking"]) or None,
                    "uniprot_window": seq[max(0, pos - 8):pos + 7] if pos and 1 <= pos <= len(seq) else None,
                    "flank_found_in_uniprot": fm["found"] if fm else None,
                    "flank_exact": fm["exact"] if fm else None,
                    "flank_identity": fm["identity"] if fm else None,
                    "flank_covers_site": fm["covers_site"] if fm else None})
        c["residue_ok" if res_ok else ("residue_bad" if res_ok is False else "residue_unchecked")] += 1
        if fm:
            c["flank_found" if fm["found"] else "flank_not_found"] += 1
            c["flank_exact"] += bool(fm["exact"])
    return {"length_uniprot": len(seq), "length_catalog": catalog_len,
            "length_matches_catalog": (len(seq) == catalog_len) if catalog_len else None,
            "n_records": len(rows), "summary": dict(c),
            "bad": [o for o in out if o["residue_matches_uniprot"] is False or o["flank_found_in_uniprot"] is False],
            "records": out}


# ------------------------------------------------ cancer-cohort overlap
def cohort_overlap(units, length, study, cb, window=0):
    """Per-PTM-site mutation counts in a user-selected cBioPortal study (no default cohort)."""
    n_seq = cb["n_sequenced"]
    muts = cb["mutations"]
    by_pos = defaultdict(list)
    for m in muts:
        by_pos[m["pos"]].append(m)
    all_samples = {m["sample"] for m in muts}
    rows, ptm_positions = [], set()
    for u in units:
        ptm_positions.add(u["pos"])
        ms = [m for off in range(-window, window + 1) for m in by_pos.get(u["pos"] + off, [])]
        at = by_pos.get(u["pos"], [])
        samples = {m["sample"] for m in ms}
        top = Counter(m["protein_change"] for m in ms).most_common(5)
        rows.append({"site": u["site"], "ptm_type": u["ptm_type"], "pos": u["pos"],
                     "n_mut_at_residue": len(at), "n_samples_at_residue": len({m["sample"] for m in at}),
                     "n_mut_window": len(ms), "n_samples_window": len(samples),
                     "freq_window": round(len(samples) / n_seq, 5) if n_seq else None,
                     "top_changes": [{"change": k, "n": v} for k, v in top],
                     "types": dict(Counter(m["type"] for m in ms)), "diseases": u["diseases"],
                     "kinases": u["kinases"], "confidence": u["best_confidence"]})
    rows.sort(key=lambda r: (-r["n_mut_at_residue"], -r["n_mut_window"]))
    total = len(muts)
    on_ptm = sum(len(by_pos.get(p, [])) for p in ptm_positions)
    p_uni = len(ptm_positions) / length if length else None
    pval = binom_sf(on_ptm, total, p_uni) if (total and p_uni) else None
    return {"study": study, "n_sequenced": n_seq, "n_mutations": total, "n_mutated_samples": len(all_samples),
            "n_mut_on_ptm_residues": on_ptm,
            "expected_on_ptm_uniform": round(total * p_uni, 2) if p_uni is not None else None,
            "p_enrichment": pval, "window": window, "sites": rows,
            "method": "Counts of mutations (all types reported by cBioPortal for this gene in the selected study) whose "
                      "protein position equals (±window) a recorded PTM residue. Enrichment null: mutations fall on "
                      "residues uniformly at random, so P(on a PTM residue)=n_PTM_residues/L (one-sided binomial). "
                      "Known driver hotspots violate this null; interpret as a screen."}


# ---------------------------------------------------------- priority
DEFAULT_W = {"confidence": 1, "in_vivo": 1, "kinase": 1, "disease": 1, "replicated": 1,
             "hotspot": 1, "crosstalk": 1}


def prioritise(units, hot, cross, weights=None):
    """Transparent additive ranking. Each criterion is 0/1 taken from stored evidence; weights default to
    equal and are user-adjustable. This is a triage aid, not a probability."""
    w = dict(DEFAULT_W)
    w.update({k: float(v) for k, v in (weights or {}).items() if k in w})
    in_hot = set()
    for r in hot.get("regions", []):
        if r["adj_p"] < 0.05:
            in_hot.update(s["pos"] for s in r["sites"])
    xt = Counter()
    for p in cross.get("pairs", []):
        if p["class"] != "same-type proximity":
            xt[(p["a"]["ptm_type"], p["a"]["pos"])] += 1
            xt[(p["b"]["ptm_type"], p["b"]["pos"])] += 1
    out = []
    for u in units:
        c = {"confidence": u["best_confidence"] == "High",
             "in_vivo": any("in vivo" in e.lower() for e in u["evidence"]),
             "kinase": bool(u["kinases"]), "disease": bool(u["diseases"]),
             "replicated": len(u["pmids"]) >= 2,
             "hotspot": u["pos"] in in_hot,
             "crosstalk": xt[(u["ptm_type"], u["pos"])] > 0}
        tot = sum(w.values()) or 1
        score = sum(w[k] for k, v in c.items() if v) / tot
        out.append({"gene": u["gene"], "uniprot": u["uniprot"], "site": u["site"], "ptm_type": u["ptm_type"],
                    "pos": u["pos"], "score": round(score, 3), "criteria": c, "kinases": u["kinases"],
                    "diseases": u["diseases"], "pmids": u["pmids"]})
    out.sort(key=lambda r: (-r["score"], r["gene"], r["pos"]))
    return out
