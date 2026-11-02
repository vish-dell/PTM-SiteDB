"""Evidence-derived analytics for PTM-SiteDB.

Design rule (project requirement): NO synthetic, simulated or assumed data.
Every function here is a deterministic computation over
  (a) records already stored in the local ptm_sites table (imported from the
      original PTM-SiteDB by sync_data.py), or
  (b) real external annotations passed in explicitly (UniProt sequences,
      EBI Proteins-API variants, AlphaFold coordinates).
Nothing is estimated from "typical" values. Statistical assumptions are stated
in the returned `method` fields so they can be quoted in a report.
"""
import math
import re
from collections import Counter, defaultdict

BLANK = {"", "—", "-", "–", "na", "n/a", "none"}
CONF_RANK = {"high": 3, "medium": 2, "low": 1}
CANCER_WORDS = re.compile(
    r"cancer|carcinoma|tumou?r|neoplas|leuk(a)?emia|lymphoma|melanoma|sarcoma|"
    r"glioma|glioblastoma|blastoma|myeloma|adenoma|malignan|polyposis", re.I)


# ----------------------------------------------------------------- helpers
def clean_kinase(k):
    k = (k or "").strip()
    return None if k.lower() in BLANK else k


def parse_site(site):
    m = re.match(r"^\s*([A-Za-z]?)\s*(\d+)", site or "")
    if not m:
        return None, None
    return (m.group(1).upper() or None), int(m.group(2))


def split_diseases(s):
    return [d.strip() for d in (s or "").split(";") if d.strip()]


def best_conf(confs):
    best = max((CONF_RANK.get(c.lower(), 0) for c in confs), default=0)
    return {3: "High", 2: "Medium", 1: "Low"}.get(best, "")


def units_from_rows(rows):
    """Collapse per-publication rows into one unit per (PTM type, residue)."""
    tmp = {}
    for r in rows:
        r = dict(r)
        aa, pos = parse_site(r.get("site"))
        if r.get("position") is not None:
            pos = r["position"]
        if pos is None:
            continue
        key = (r["uniprot"], r["ptm_type"], pos, aa)  # uniprot in key: ALL-protein calls must not merge proteins
        u = tmp.setdefault(key, {
            "gene": r["gene"], "uniprot": r["uniprot"], "ptm_type": r["ptm_type"],
            "site": r["site"], "aa": aa, "pos": pos, "kinases": set(), "diseases": set(),
            "pmids": set(), "evidence": set(), "confidences": set(), "record_ids": []})
        k = clean_kinase(r.get("kinase"))
        if k:
            u["kinases"].add(k)
        u["diseases"].update(split_diseases(r.get("disease_context")))
        if (r.get("pmid") or "").strip():
            u["pmids"].add(r["pmid"].strip())
        if (r.get("evidence") or "").strip():
            u["evidence"].add(r["evidence"].strip())
        if (r.get("confidence") or "").strip():
            u["confidences"].add(r["confidence"].strip())
        u["record_ids"].append(r["id"])
    out = []
    for u in tmp.values():
        u = dict(u)
        for k in ("kinases", "diseases", "pmids", "evidence", "confidences"):
            u[k] = sorted(u[k])
        u["best_confidence"] = best_conf(u["confidences"])
        out.append(u)
    out.sort(key=lambda x: (x["pos"], x["ptm_type"]))
    return out


def _log_binom_pmf(i, n, p):
    return (math.lgamma(n + 1) - math.lgamma(i + 1) - math.lgamma(n - i + 1)
            + i * math.log(p) + (n - i) * math.log1p(-p))


def binom_sf(k, n, p):
    """P(X >= k) for X ~ Binomial(n, p), computed in log space (no overflow)."""
    if k <= 0:
        return 1.0
    if k > n:
        return 0.0
    if p <= 0:
        return 0.0
    if p >= 1:
        return 1.0
    return min(1.0, sum(math.exp(_log_binom_pmf(i, n, p)) for i in range(k, n + 1)))


def _brief(u):
    return {"site": u["site"], "pos": u["pos"], "ptm_type": u["ptm_type"],
            "kinases": u["kinases"], "diseases": u["diseases"], "pmids": u["pmids"],
            "evidence": u["evidence"], "best_confidence": u["best_confidence"]}


# ---------------------------------------------------------- 1. hotspots
def hotspots(units, length, window=15, min_sites=3):
    """Find windows holding more modified residues than expected by chance.

    Unit of analysis = distinct modified residue positions (a residue carrying two
    PTM types counts once). Null model: the n observed positions are scattered
    uniformly over the L residues, so a window of w residues holds
    Binomial(n, w/L) of them. Windows are anchored at observed sites; overlapping
    passing windows are merged into regions. adj_p is a Bonferroni correction over
    ceil(L/w) non-overlapping windows (approximate; documented, not hidden).
    Caveat: PTM databases are curation-biased toward well-studied regions, so the
    null is a screening tool, not proof of biological enrichment.
    """
    by_pos = defaultdict(list)
    for u in units:
        by_pos[u["pos"]].append(u)
    positions = sorted(by_pos)
    n = len(positions)
    L = length or (positions[-1] if positions else 0)
    method = (f"Sliding window w={window} residues anchored at observed sites; null = uniform "
              f"placement of n={n} modified residues over L={L}; one-sided binomial tail; "
              f"Bonferroni over ceil(L/w) windows.")
    res = {"params": {"window": window, "min_sites": min_sites}, "length": L,
           "n_positions": n, "n_units": len(units), "regions": [], "profile": [], "method": method}
    if n == 0 or L == 0:
        return res
    p_win = min(1.0, window / L)
    n_tests = max(1, math.ceil(L / window))
    # density profile: modified residues in [s, s+w-1] for every start s
    prof, j0, j1 = [], 0, 0
    for s in range(1, max(1, L - window + 2)):
        while j0 < n and positions[j0] < s:
            j0 += 1
        j1 = max(j1, j0)
        while j1 < n and positions[j1] <= s + window - 1:
            j1 += 1
        prof.append([s, j1 - j0])
    res["profile"] = prof
    passing = []
    for i, s in enumerate(positions):
        inside = [q for q in positions[i:] if q <= s + window - 1]
        if len(inside) >= min_sites:
            passing.append((s, s + window - 1, len(inside), inside))
    merged = []
    for s, e, k, inside in passing:
        if merged and s <= merged[-1]["win_end"]:
            m = merged[-1]
            m["win_end"] = max(m["win_end"], e)
            m["positions"].update(inside)
            if k > m["best_k"]:
                m["best_k"], m["best_start"] = k, s
        else:
            merged.append({"win_end": e, "positions": set(inside), "best_k": k, "best_start": s,
                           "win_start": s})
    for m in merged:
        pos_list = sorted(m["positions"])
        start, end = pos_list[0], pos_list[-1]
        members = [u for q in pos_list for u in by_pos[q]]
        p_raw = binom_sf(m["best_k"], n, p_win)
        expected = n * p_win
        res["regions"].append({
            "start": start, "end": end, "span": end - start + 1,
            "n_residues": len(pos_list), "n_modifications": len(members),
            "best_window": [m["best_start"], m["best_start"] + window - 1],
            "best_window_count": m["best_k"],
            "expected_in_window": round(expected, 3),
            "fold_enrichment": round(m["best_k"] / expected, 2) if expected > 0 else None,
            "p_value": p_raw, "adj_p": min(1.0, p_raw * n_tests),
            "ptm_types": dict(Counter(u["ptm_type"] for u in members)),
            "kinases": sorted({k for u in members for k in u["kinases"]}),
            "diseases": sorted({d for u in members for d in u["diseases"]}),
            "pmids": sorted({p for u in members for p in u["pmids"]}),
            "sites": [_brief(u) for u in members],
        })
    from .proximity import bh_fdr
    for r, q in zip(res["regions"], bh_fdr([r["p_value"] for r in res["regions"]])):
        r["fdr_q"] = q
    res["regions"].sort(key=lambda r: (r["adj_p"], -r["n_residues"]))
    return res


# --------------------------------------------------------- 2. crosstalk
def crosstalk(units, max_dist=10, include_same_type=False):
    """Pairs of modified residues of DIFFERENT PTM types within max_dist residues
    (distance 0 = two PTM types on the same residue, i.e. direct competition)."""
    us = sorted(units, key=lambda u: (u["pos"], u["ptm_type"]))
    pairs = []
    for i, a in enumerate(us):
        for b in us[i + 1:]:
            d = b["pos"] - a["pos"]
            if d > max_dist:
                break
            same_type = a["ptm_type"] == b["ptm_type"]
            if same_type and not include_same_type:
                continue
            if same_type and d == 0:
                continue
            shared_k = sorted(set(a["kinases"]) & set(b["kinases"]))
            shared_d = sorted(set(a["diseases"]) & set(b["diseases"]))
            pairs.append({
                "a": _brief(a), "b": _brief(b), "distance": d,
                "class": ("same-residue competition" if d == 0 else
                          ("cross-type proximity" if not same_type else "same-type proximity")),
                "shared_kinases": shared_k, "shared_diseases": shared_d,
                "both_kinase_linked": bool(a["kinases"] and b["kinases"]),
                "pmids": sorted(set(a["pmids"]) | set(b["pmids"])),
            })
    pairs.sort(key=lambda p: (p["distance"], p["a"]["pos"]))
    type_pairs = Counter("  +  ".join(sorted([p["a"]["ptm_type"], p["b"]["ptm_type"]])) for p in pairs)
    return {"params": {"max_dist": max_dist, "include_same_type": include_same_type},
            "n_pairs": len(pairs), "type_pair_counts": dict(type_pairs), "pairs": pairs,
            "method": "All pairs of distinct (PTM type, residue) units on this protein whose "
                      "sequence separation is <= max_dist; no scoring beyond what the records state."}


# ------------------------------------------------ 3. kinase–PTM network
def kinase_network(rows, db_genes):
    """kinase -> protein -> PTM site -> disease context, strictly from records that name a kinase.
    A kinase symbol that is exactly a gene symbol in the database is flagged (cascade candidate)."""
    units = units_from_rows(rows)
    db_genes = set(db_genes)
    paths, unlinked = [], 0
    k_stats = defaultdict(lambda: {"sites": set(), "proteins": set(), "diseases": set(), "pmids": set()})
    for u in units:
        if not u["kinases"]:
            unlinked += 1
            continue
        sid = f"{u['gene']}|{u['ptm_type']}|{u['site']}"
        for k in u["kinases"]:
            paths.append({"kinase": k, "protein": u["gene"], "site_id": sid,
                          "site": f"{u['gene']} {u['site']}", "ptm_type": u["ptm_type"],
                          "diseases": u["diseases"], "pmids": u["pmids"],
                          "confidence": u["best_confidence"]})
            s = k_stats[k]
            s["sites"].add(sid)
            s["proteins"].add(u["gene"])
            s["diseases"].update(u["diseases"])
            s["pmids"].update(u["pmids"])
    kinases = sorted(({"kinase": k, "n_sites": len(v["sites"]), "n_proteins": len(v["proteins"]),
                       "diseases": sorted(v["diseases"]), "n_pmids": len(v["pmids"]),
                       "also_in_database": k in db_genes} for k, v in k_stats.items()),
                     key=lambda x: (-x["n_sites"], x["kinase"]))
    # cascades: X -> (protein K) and K -> other sites, using exact symbol matches only
    by_kinase = defaultdict(list)
    for p in paths:
        by_kinase[p["kinase"]].append(p)
    cascades = []
    for p in paths:
        g = p["protein"]
        if g in by_kinase and g != p["kinase"]:
            cascades.append({"upstream_kinase": p["kinase"], "via_site": p["site"],
                             "kinase_protein": g,
                             "downstream_sites": sorted({q["site"] for q in by_kinase[g]}),
                             "pmids": p["pmids"]})
    return {"paths": paths, "kinases": kinases, "cascades": cascades,
            "n_sites_with_kinase": len({p["site_id"] for p in paths}),
            "n_sites_without_kinase": unlinked,
            "method": "Paths are built only from records whose kinase field is filled; sites with "
                      "no kinase are counted but not drawn. Cascades require an exact symbol match "
                      "between a kinase and a protein in the database."}


# -------------------------------------------- 4/5. variants x PTM sites
def variant_ptm_overlap(units, variants, window=0):
    """Join real UniProt/EBI missense variants to PTM residues (|variant pos - site pos| <= window)."""
    if not variants:
        return {"params": {"window": window}, "hits": [], "summary": {}}
    by_pos = defaultdict(list)
    for v in variants:
        by_pos[v["pos"]].append(v)
    hits = []
    for u in units:
        for off in range(-window, window + 1):
            for v in by_pos.get(u["pos"] + off, []):
                at_res = off == 0
                hits.append({
                    "site": u["site"], "ptm_type": u["ptm_type"], "site_pos": u["pos"],
                    "site_kinases": u["kinases"], "site_diseases": u["diseases"],
                    "site_confidence": u["best_confidence"], "site_pmids": u["pmids"],
                    "variant": v["label"], "variant_pos": v["pos"], "offset": off,
                    "at_modified_residue": at_res,
                    "wildtype_matches_site_residue": (v["wt"] == u["aa"]) if (at_res and u["aa"]) else None,
                    "consequence": v.get("consequence"), "clinical": v.get("clinical", []),
                    "diseases": v.get("diseases", []), "cancer_associated": v.get("cancer", False),
                    "cosmic": v.get("cosmic", []), "dbsnp": v.get("dbsnp", []),
                    "clinvar": v.get("clinvar", []), "source": v.get("source"),
                })
    hits.sort(key=lambda h: (abs(h["offset"]), not h["cancer_associated"], h["site_pos"]))
    summ = {
        "n_variants_total": len(variants),
        "n_hits": len(hits),
        "n_at_modified_residue": sum(h["at_modified_residue"] for h in hits),
        "n_cancer_associated_hits": sum(h["cancer_associated"] for h in hits),
        "n_cosmic_hits": sum(bool(h["cosmic"]) for h in hits),
        "sites_with_hits": len({(h["site"], h["ptm_type"]) for h in hits}),
        "sites_total": len(units),
    }
    return {"params": {"window": window}, "hits": hits, "summary": summ,
            "method": "Exact positional join of PTM residues with single-residue variants from the EBI "
                      "Proteins API (UniProt curated + large-scale studies). 'Cancer-associated' means "
                      "the variant's own disease annotation matches a cancer keyword list or carries a "
                      "COSMIC cross-reference; nothing is inferred beyond that."}


# ----------------------------------------------------- 6. hypotheses
GENERIC_ASSAYS = {
    "module": ["Site-directed mutagenesis of the cluster (phospho-dead / phospho-mimetic, K->R)",
               "Site-specific antibodies or targeted MS to test co-occurrence"],
    "crosstalk": ["Mutate site A and quantify modification at site B (targeted MS / site-specific antibody)",
                  "Inhibit the annotated upstream enzyme and monitor both sites"],
    "kinase": ["In vitro kinase assay on wild-type vs. site-mutant substrate",
               "Kinase inhibitor / knockdown time-course with phospho-specific readout"],
    "validation": ["Orthogonal validation: targeted MS (PRM/SRM) or site-specific antibody"],
    "variant": ["Compare wild-type vs. variant protein for modification level and downstream function"],
    "conserved": ["Mutate the conserved residue and test function in a cell model"],
}


def _crit(name, met, applicable=True):
    return {"name": name, "met": bool(met), "applicable": applicable}


def _score(criteria):
    app = [c for c in criteria if c["applicable"]]
    met = sum(c["met"] for c in app)
    return met, len(app)


def generate_hypotheses(gene, units, hot, cross, variant_res=None, conservation=None, spatial=None):
    """Turn computed evidence into ranked, testable research QUESTIONS.

    These are prompts for experiments, not findings. `score` is simply the fraction of
    applicable evidence criteria that are met (criteria are listed on every card)."""
    out = []
    ext_var = variant_res is not None and bool(variant_res.get("hits") is not None)
    cons_map = {}
    if conservation and conservation.get("sites"):
        cons_map = {(s["site"], s["ptm_type"]): s for s in conservation["sites"]}
    n_orth = (conservation or {}).get("n_orthologs", 0)

    # H1: hotspot modules
    for r in (hot or {}).get("regions", []):
        types, kin, dis = r["ptm_types"], r["kinases"], r["diseases"]
        crit = [
            _crit("Enriched vs. uniform expectation (adj. p < 0.05)", r["adj_p"] < 0.05),
            _crit("Two or more PTM types in the region", len(types) >= 2),
            _crit("Upstream kinase annotated", bool(kin)),
            _crit("Disease/cancer context annotated", bool(dis)),
            _crit("At least one High-confidence site", any(s["best_confidence"] == "High" for s in r["sites"])),
        ]
        if ext_var:
            rng = range(r["start"], r["end"] + 1)
            crit.append(_crit("Cancer-associated variant inside the region",
                              any(h["cancer_associated"] and h["site_pos"] in rng for h in variant_res["hits"])))
        if cons_map and n_orth:
            flags = [cons_map[(s["site"], s["ptm_type"])]["n_identical"] == cons_map[(s["site"], s["ptm_type"])]["n_aligned"] >= 3
                     for s in r["sites"] if (s["site"], s["ptm_type"]) in cons_map]
            crit.append(_crit("Every mapped site identical in all orthologs that align there (>=3, gaps excluded)",
                              bool(flags) and all(flags), bool(flags)))
        tstr = ", ".join(f"{v}× {k}" for k, v in types.items())
        stmt = (f"Do the {r['n_residues']} modified residues in {gene} {r['start']}–{r['end']} ({tstr}) "
                f"act as a coordinated regulatory module"
                + (f" controlled by {', '.join(kin)}" if kin else "") + "?")
        basis = [f"{r['best_window_count']} modified residues fall in {r['best_window'][0]}–{r['best_window'][1]} "
                 f"(expected {r['expected_in_window']} under a uniform null; adj. p = {r['adj_p']:.2g}).",
                 "Sites: " + ", ".join(s["site"] for s in r["sites"]) + "."]
        out.append(_card("Hotspot module", stmt, basis, crit, "module", r["pmids"], kin, dis,
                         [s["site"] for s in r["sites"]]))

    # H2: crosstalk
    for p in (cross or {}).get("pairs", []):
        if p["class"] == "same-type proximity":
            continue
        a, b = p["a"], p["b"]
        crit = [
            _crit("Both sites High confidence", a["best_confidence"] == "High" and b["best_confidence"] == "High"),
            _crit("Both sites have an annotated kinase", p["both_kinase_linked"]),
            _crit("Shared disease context", bool(p["shared_diseases"])),
        ]
        if spatial and spatial.get("available"):
            sp = next((x for x in spatial["pairs"] if {x["pos_a"], x["pos_b"]} == {a["pos"], b["pos"]}), None)
            crit.append(_crit(f"CA–CA ≤ {spatial['threshold']} Å in the AlphaFold model (pLDDT ≥ 70 both)",
                              bool(sp and sp["ca_dist"] <= spatial["threshold"] and sp["both_confident"]),
                              sp is not None))
        if p["distance"] == 0:
            stmt = f"Do {a['ptm_type']} and {b['ptm_type']} compete for {gene} residue {a['site']}?"
        else:
            stmt = (f"Does {a['ptm_type'].lower()} at {gene} {a['site']} influence {b['ptm_type'].lower()} at "
                    f"{b['site']} ({p['distance']} residues apart), or vice versa?")
        basis = [f"{a['site']} ({a['ptm_type']}) and {b['site']} ({b['ptm_type']}) are recorded {p['distance']} residue(s) apart."]
        out.append(_card("PTM crosstalk", stmt, basis, crit, "crosstalk", p["pmids"],
                         sorted(set(a["kinases"]) | set(b["kinases"])),
                         sorted(set(a["diseases"]) | set(b["diseases"])), [a["site"], b["site"]]))

    # H3: one kinase, several sites
    by_k = defaultdict(list)
    for u in units:
        for k in u["kinases"]:
            by_k[k].append(u)
    region_pos = [range(r["start"], r["end"] + 1) for r in (hot or {}).get("regions", [])]
    for k, us in by_k.items():
        if len(us) < 2:
            continue
        pos = [u["pos"] for u in us]
        crit = [
            _crit("Sites lie in one hotspot region", any(all(p in rg for p in pos) for rg in region_pos)),
            _crit("Sites share a disease context", bool(set.intersection(*[set(u["diseases"]) for u in us]))),
            _crit("All sites High confidence", all(u["best_confidence"] == "High" for u in us)),
        ]
        stmt = (f"Does {k} act on {', '.join(u['site'] for u in us)} of {gene} in a coordinated manner "
                f"(e.g. sequentially or as a priming cascade)?")
        basis = [f"{k} is annotated as the upstream kinase for {len(us)} sites on {gene}."]
        out.append(_card("Kinase multi-site", stmt, basis, crit, "kinase",
                         sorted({p for u in us for p in u["pmids"]}), [k],
                         sorted({d for u in us for d in u["diseases"]}), [u["site"] for u in us]))

    # H4: validation priority inside significant hotspots
    for r in (hot or {}).get("regions", []):
        weak = [s for s in r["sites"] if s["best_confidence"] != "High" or "in vivo" not in [e.lower() for e in s["evidence"]]]
        if weak and r["adj_p"] < 0.05:
            crit = [_crit("Region is enriched (adj. p < 0.05)", True),
                    _crit("Region has disease context", bool(r["diseases"])),
                    _crit("Region has an annotated kinase", bool(r["kinases"]))]
            stmt = (f"Are {', '.join(s['site'] for s in weak)} in the {gene} {r['start']}–{r['end']} hotspot "
                    f"genuine? They rest on limited evidence ({', '.join(sorted({e for s in weak for e in s['evidence']} or {'n/a'}))}; "
                    f"{', '.join(sorted({s['best_confidence'] or 'n/a' for s in weak}))} confidence).")
            out.append(_card("Validation priority", stmt,
                             ["Neighbouring sites in an enriched cluster are better supported than these."],
                             crit, "validation", sorted({p for s in weak for p in s["pmids"]}), r["kinases"], r["diseases"],
                             [s["site"] for s in weak]))

    # H5: variants at / near PTM residues
    if ext_var:
        for h in variant_res["hits"]:
            if not (h["at_modified_residue"] or h["cancer_associated"]):
                continue
            crit = [
                _crit("Variant replaces the modified residue", h["at_modified_residue"]),
                _crit("Wild-type residue matches the PTM residue letter", bool(h["wildtype_matches_site_residue"]),
                      h["wildtype_matches_site_residue"] is not None),
                _crit("Cancer-associated annotation on the variant", h["cancer_associated"]),
                _crit("Somatic (COSMIC) cross-reference", bool(h["cosmic"])),
                _crit("PTM site is High confidence", h["site_confidence"] == "High"),
                _crit("Upstream kinase annotated", bool(h["site_kinases"])),
            ]
            where = "replaces" if h["at_modified_residue"] else f"lies {abs(h['offset'])} residues from"
            stmt = (f"Does {gene} {h['variant']} ({', '.join(h['diseases'][:2]) or 'variant annotation'}) alter "
                    f"{h['ptm_type'].lower()} at {h['site']}? The variant {where} the modified residue.")
            out.append(_card("Mutation–PTM–cancer", stmt,
                             [f"{h['variant']} at offset {h['offset']:+d} from {h['site']}."],
                             crit, "variant", h["site_pmids"], h["site_kinases"],
                             sorted(set(h["site_diseases"]) | set(h["diseases"])), [h["site"]]))

    # H6: conserved + disease-linked
    if cons_map and n_orth >= 3:
        for u in units:
            c = cons_map.get((u["site"], u["ptm_type"]))
            if c and c["n_aligned"] >= 3 and c["n_identical"] == c["n_aligned"] and u["diseases"]:
                crit = [_crit(f"Residue identical in all {c['n_aligned']} orthologs aligned at this site", True),
                        _crit("Disease context annotated", True),
                        _crit("High confidence", u["best_confidence"] == "High"),
                        _crit("Upstream kinase annotated", bool(u["kinases"]))]
                stmt = (f"Is the evolutionarily conserved {u['ptm_type'].lower()} at {gene} {u['site']} "
                        f"required for its role in {', '.join(u['diseases'][:2])}?")
                out.append(_card("Conserved site", stmt,
                                 [f"{u['site']} is identical in {c['n_aligned']}/{c['n_aligned']} orthologs aligned at this site."],
                                 crit, "conserved", u["pmids"], u["kinases"], u["diseases"], [u["site"]]))

    out.sort(key=lambda h: (-h["score"], -h["criteria_total"], h["statement"]))
    for i, h in enumerate(out, 1):
        h["id"] = f"{gene}-H{i:03d}"
    return out


def _card(kind, statement, basis, criteria, assay_key, pmids, kinases, diseases, sites):
    met, tot = _score(criteria)
    return {"kind": kind, "statement": statement, "basis": basis, "criteria": criteria,
            "criteria_met": met, "criteria_total": tot, "score": round(met / tot, 3) if tot else 0.0,
            "suggested_assays": GENERIC_ASSAYS[assay_key],
            "evidence": {"sites": sites, "pmids": pmids, "kinases": kinases, "diseases": diseases},
            "disclaimer": "Research question generated from database evidence; not a finding."}
