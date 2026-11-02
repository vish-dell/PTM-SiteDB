"""Residue-verified 3D distances between PTM residues (real coordinates only).

Mapping rule: a UniProt position is used only if (1) SIFTS maps it into the structure and
(2) the residue actually present at that place has the same amino-acid letter as the PTM
site (e.g. site S473 must land on a Ser/SEP). Anything else stays 'unmapped'.
"""
import math

from . import external as ext

PARENT = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
          "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
          "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
          # standard parent residues of common modified residues (PDB HETATM)
          "SEP": "S", "TPO": "T", "PTR": "Y", "MSE": "M", "ALY": "K", "M3L": "K", "MLY": "K", "MLZ": "K",
          "KCX": "K", "CSO": "C", "HYP": "P", "SNN": "N"}


def parse_pdb(text):
    """{(chain, resseq): {'aa','atoms':[(name,x,y,z)],'ca':(x,y,z)|None,'b':float|None}}"""
    res = {}
    for line in text.splitlines():
        rec = line[:6]
        if rec not in ("ATOM  ", "HETATM"):
            continue
        name3 = line[17:20].strip()
        if rec == "HETATM" and name3 not in PARENT:
            continue
        alt = line[16]
        if alt not in (" ", "A", "1"):
            continue
        try:
            an = line[12:16].strip()
            key = (line[21].strip(), int(line[22:26]))
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            b = float(line[60:66])
        except (ValueError, IndexError):
            continue
        if line[26].strip():          # insertion code: skip, keeps numbering unambiguous
            continue
        r = res.setdefault(key, {"aa": PARENT.get(name3), "atoms": [], "ca": None, "b": None})
        r["atoms"].append((an, *xyz))
        if an == "CA":
            r["ca"], r["b"] = xyz, b
    return res


def _d(p, q):
    return math.sqrt(sum((p[i] - q[i]) ** 2 for i in range(3)))


def _map_experimental(blocks, pos):
    """Yield (chain, pdb_resnum) for every SIFTS block covering UniProt position pos."""
    for b in blocks:
        try:
            us, ue = int(b["unp_start"]), int(b["unp_end"])
            if us <= pos <= ue:
                st = b["start"]
                yield str(b.get("chain_id")), int(st.get("author_residue_number", st.get("residue_number"))) + (pos - us)
        except (KeyError, TypeError, ValueError):
            continue


def _best_experimental(acc, units, max_structures=6):
    ids = ext.rcsb_entries_for(acc)
    cand = []

    def _blocks(pdb):
        try:
            return pdb, ext.sifts_mappings(pdb, acc)
        except ext.ExternalError:
            return pdb, None
    from concurrent.futures import ThreadPoolExecutor
    with ThreadPoolExecutor(max_workers=8) as ex:       # 25 sequential lookups were the main slowdown
        looked = list(ex.map(_blocks, ids))
    for pdb, blocks in looked:
        if blocks is None:
            continue
        n = sum(1 for u in units if any(True for _ in _map_experimental(blocks, u["pos"])))
        if n:
            cand.append((n, pdb, blocks))
    cand.sort(key=lambda x: -x[0])
    best = None
    for n_cov, pdb, blocks in cand[:max_structures]:
        try:
            res = parse_pdb(ext.pdb_text(pdb))
        except ext.ExternalError:
            continue
        per_chain = {}
        for u in units:
            for chain, num in _map_experimental(blocks, u["pos"]):
                r = res.get((chain, num))
                if r and r["ca"] and (u["aa"] is None or r["aa"] == u["aa"]):
                    per_chain.setdefault(chain, {})[(u["ptm_type"], u["pos"])] = (u, num)
        for chain, m in per_chain.items():
            if best is None or len(m) > len(best["mapped"]):
                best = {"source": "experimental", "id": pdb.upper(), "chain": chain, "mapped": m, "res": res}
        if best and len(best["mapped"]) == len({(u["ptm_type"], u["pos"]) for u in units}):
            break
    return best, len(ids)


def _alphafold(acc, units):
    """Find the best real AlphaFold chain for PTM coverage.

    AlphaFold normally uses chain A, but the parser deliberately does not assume that:
    if a returned model contains another chain, choose the chain that actually verifies
    the largest number of PTM residues. No residue is inferred when the amino-acid identity
    disagrees.
    """
    entry, text = ext.alphafold_model(acc)
    res = parse_pdb(text)
    chains = sorted({chain for chain, _ in res})
    candidates = []
    for chain in chains:
        m, why = {}, []
        for u in units:
            r = res.get((chain, u["pos"]))
            if r and r["ca"] and (u["aa"] is None or r["aa"] == u["aa"]):
                m[(u["ptm_type"], u["pos"])] = (u, u["pos"])
            elif not r:
                why.append(f"{u['site']}: position {u['pos']} is outside AlphaFold chain {chain}")
            else:
                why.append(f"{u['site']}: record says {u['aa']} but model has {r['aa']} at {u['pos']} on chain {chain}")
        candidates.append((len(m), chain, m, why))
    if not candidates:
        raise ext.ExternalError("AlphaFold model contains no parseable coordinate chain")
    _, chain, m, why = max(candidates, key=lambda x: (x[0], x[1] == "A"))
    return {"source": "alphafold", "id": entry, "chain": chain, "mapped": m, "res": res, "why": why}


def _coverage_summary(units, experimental, alphafold):
    keys = {(u["ptm_type"], u["pos"]): u["site"] for u in units}
    def summary(candidate):
        mapped = set(candidate["mapped"]) if candidate else set()
        return {
            "mapped": len(mapped),
            "total": len(keys),
            "sites": [keys[k] for k in sorted(mapped, key=lambda x: (x[1], x[0]))],
            "unmapped": [keys[k] for k in sorted(set(keys) - mapped, key=lambda x: (x[1], x[0]))],
            "source": candidate["source"] if candidate else None,
            "structure_id": candidate["id"] if candidate else None,
            "chain": candidate["chain"] if candidate else None,
        }
    return {"experimental": summary(experimental), "alphafold": summary(alphafold)}


def spatial_analysis(acc, units, threshold=10.0, allow_alphafold=True):
    """Return distances between PTM residues from one real coordinate model.

    Mapping policy:
      1. Prefer the best single experimental PDB chain verified by SIFTS.
      2. If AlphaFold is enabled, compare its best chain and use it when it gives
         strictly better PTM coverage (or when experimental coverage is <2).
      3. Never mix coordinates from different structures/chains for a distance.
    """
    if len(units) < 2:
        return {"available": False, "reason": "Fewer than two PTM residues are recorded for this protein."}

    notes, best_exp, best_af, n_entries = [], None, None, 0
    try:
        best_exp, n_entries = _best_experimental(acc, units)
    except ext.ExternalError as e:
        notes.append(f"Experimental structure lookup failed: {e}")

    if allow_alphafold:
        try:
            best_af = _alphafold(acc, units)
        except ext.ExternalError as e:
            notes.append(f"AlphaFold lookup failed: {e}")

    coverage = _coverage_summary(units, best_exp, best_af)
    exp_n = coverage["experimental"]["mapped"]
    af_n = coverage["alphafold"]["mapped"]

    # Prefer experimental when it is at least as complete. Otherwise use the real
    # AlphaFold model so sites absent from the experimental chain are not silently lost.
    if best_exp and exp_n >= 2 and exp_n >= af_n:
        best = best_exp
    elif best_af and af_n >= 2:
        best = best_af
        if best_exp and af_n > exp_n:
            notes.append(f"AlphaFold fallback selected because it verifies {af_n}/{len(units)} PTM sites "
                         f"versus {exp_n}/{len(units)} in the best experimental chain.")
    elif best_exp and exp_n >= 2:
        best = best_exp
    else:
        return {"available": False, "notes": notes, "n_pdb_entries": n_entries,
                "coverage": coverage,
                "reason": "Fewer than two PTM residues could be verified in any single real structure; "
                          "no distance was inferred."}

    res, mapped = best["res"], best["mapped"]
    ck = best["chain"]
    keys = sorted(mapped, key=lambda k: (k[1], k[0]))
    pairs = []
    for i in range(len(keys)):
        for j in range(i + 1, len(keys)):
            (ua, na), (ub, nb) = mapped[keys[i]], mapped[keys[j]]
            ra, rb = res[(ck, na)], res[(ck, nb)]
            ca = _d(ra["ca"], rb["ca"])
            mind = min(_d(a[1:], b[1:]) for a in ra["atoms"] for b in rb["atoms"])
            conf = True
            plddt = None
            if best["source"] == "alphafold":
                plddt = (ra["b"], rb["b"])
                conf = ra["b"] >= 70 and rb["b"] >= 70
            sep = abs(ub["pos"] - ua["pos"])
            pairs.append({"site_a": ua["site"], "type_a": ua["ptm_type"], "pos_a": ua["pos"], "resi_a": na,
                          "site_b": ub["site"], "type_b": ub["ptm_type"], "pos_b": ub["pos"], "resi_b": nb,
                          "seq_sep": sep, "ca_dist": round(ca, 2), "min_atom_dist": round(mind, 2),
                          "both_confident": conf, "plddt": plddt,
                          "cross_type": ua["ptm_type"] != ub["ptm_type"],
                          "long_range_close": sep >= 20 and ca <= threshold})
    pairs.sort(key=lambda p: p["ca_dist"])
    mapped_keys = set(mapped)
    mapping_records = []
    for u in units:
        key = (u["ptm_type"], u["pos"])
        if key in mapped:
            _, resnum = mapped[key]
            mapping_records.append({
                "site": u["site"], "ptm_type": u["ptm_type"], "uniprot_position": u["pos"],
                "residue": u.get("aa"), "structure_source": best["source"],
                "structure_id": best["id"], "chain": ck, "structure_residue_number": resnum,
                "mapping_status": "mapped",
                "experimental": best["source"] == "experimental",
                "predicted": best["source"] == "alphafold"
            })
        else:
            mapping_records.append({
                "site": u["site"], "ptm_type": u["ptm_type"], "uniprot_position": u["pos"],
                "residue": u.get("aa"), "structure_source": None, "structure_id": None, "chain": None,
                "structure_residue_number": None, "mapping_status": "unmapped",
                "experimental": False, "predicted": False
            })
    unmapped = [u["site"] for u in units if (u["ptm_type"], u["pos"]) not in mapped_keys]
    source_label = "Experimental" if best["source"] == "experimental" else "AlphaFold predicted"
    return {"available": True, "source": best["source"], "structure_id": best["id"], "chain": ck,
            "threshold": threshold, "n_mapped": len(mapped), "n_units": len(units),
            "unmapped": sorted(set(unmapped)), "mapping": mapping_records, "pairs": pairs, "notes": notes,
            "coverage": coverage,
            "n_close": sum(p["ca_dist"] <= threshold for p in pairs),
            "n_long_range_close": sum(p["long_range_close"] for p in pairs),
            "basis": ("Experimental coordinates: residue mapping via PDBe SIFTS, each residue identity verified."
                      if best["source"] == "experimental" else
                      "AlphaFold PREDICTED model, selected because it provides better real PTM coverage than the best experimental chain. Pairs with pLDDT < 70 are flagged low-confidence."),
            "note": f"Geometric proximity only — {source_label} coordinates; this does not establish functional interaction or crosstalk.",
            "provenance": {"accession": acc, "source": best["source"], "structure_id": best["id"], "chain": ck,
                           "sifts_verified": best["source"] == "experimental",
                           "predicted_model": best["source"] == "alphafold",
                           "pLDDT_threshold": 70 if best["source"] == "alphafold" else None,
                           "experimental_candidates_checked": n_entries}
    }
