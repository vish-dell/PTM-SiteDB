"""Server-side Ensembl Compara conservation for PTM sites (drop into your Flask app).

Why: calling Ensembl from the browser once per site/species trips rate limits (HTTP 429 comes back
without CORS headers, which the browser reports as just "Failed to fetch") and the spinner never ends.
Here: ONE request per protein, retries with back-off, a hard timeout, and a 24 h cache.
"""
import threading
import time

import requests
from flask import Blueprint, jsonify, request

bp = Blueprint("conservation_api", __name__)

ENSEMBL = "https://rest.ensembl.org"
SPECIES = ["pan_troglodytes", "macaca_mulatta", "mus_musculus", "rattus_norvegicus",
           "canis_lupus_familiaris", "bos_taurus", "gallus_gallus", "xenopus_tropicalis",
           "danio_rerio", "drosophila_melanogaster", "caenorhabditis_elegans"]
TTL = 24 * 3600
_cache, _lock = {}, threading.Lock()
_http = requests.Session()
_http.headers.update({"Accept": "application/json", "User-Agent": "ptm-atlas/4.0"})


def ensembl_get(path, params=None, retries=4, timeout=45):
    last = None
    for attempt in range(retries):
        try:
            r = _http.get(ENSEMBL + path, params=params, timeout=timeout)
            if r.status_code == 429:                                   # rate limited: honour Retry-After
                time.sleep(float(r.headers.get("Retry-After", 1 + attempt)))
                last = "HTTP 429 (rate limited)"
                continue
            if r.status_code in (500, 502, 503, 504):
                last = f"HTTP {r.status_code}"
                time.sleep(1 + 2 * attempt)
                continue
            r.raise_for_status()
            return r.json()
        except requests.RequestException as e:
            last = str(e)
            time.sleep(1 + attempt)
    raise RuntimeError(f"Ensembl request failed after {retries} tries: {last}")


def fetch_homologies(symbol):
    key = symbol.upper()
    with _lock:
        hit = _cache.get(key)
        if hit and time.time() - hit[0] < TTL:
            return hit[1]
    data = ensembl_get(f"/homology/symbol/homo_sapiens/{symbol}",
                       params={"type": "orthologues", "sequence": "protein", "aligned": 1,
                               "target_species": SPECIES, "content-type": "application/json"})
    homs = (data.get("data") or [{}])[0].get("homologies", [])
    with _lock:
        _cache[key] = (time.time(), homs)
    return homs


def residue_at(align_seq, pos):
    """Residue and alignment column for the pos-th (1-based) non-gap residue of an aligned sequence."""
    n = 0
    for col, ch in enumerate(align_seq):
        if ch != "-":
            n += 1
            if n == pos:
                return ch.upper(), col
    return None, None


def conservation_for_site(homs, pos, aa):
    rows, same, total = [], 0, 0
    for h in homs:
        src, tgt = h["source"], h["target"]
        s_res, col = residue_at(src["align_seq"], pos)
        if s_res is None:
            continue
        if s_res != aa:                       # Ensembl canonical numbering differs from UniProt
            return {"status": "numbering_mismatch",
                    "detail": f"Ensembl protein has {s_res} at {pos}, record says {aa}"}
        t_res = tgt["align_seq"][col].upper()
        total += 1
        ident = t_res == aa
        same += ident
        rows.append({"species": tgt["species"], "residue": t_res, "identical": ident,
                     "perc_id": tgt.get("perc_id")})
    return {"status": "ok", "n_orthologues": total, "n_identical": same,
            "fraction": round(same / total, 3) if total else None, "orthologues": rows}


@bp.route("/api/conservation/<symbol>")
def conservation(symbol):
    """GET /api/conservation/AKT1?sites=S473,T308,Y474"""
    try:
        homs = fetch_homologies(symbol)
    except Exception as e:
        return jsonify({"error": str(e)}), 502
    out = []
    for raw in filter(None, (request.args.get("sites") or "").split(",")):
        raw = raw.strip()
        try:
            aa, pos = raw[0].upper(), int("".join(c for c in raw if c.isdigit()))
        except ValueError:
            out.append({"site": raw, "status": "unreadable_label"})
            continue
        out.append({"site": raw, **conservation_for_site(homs, pos, aa)})
    return jsonify({"symbol": symbol, "n_orthologues_returned": len(homs), "sites": out})