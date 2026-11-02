"""Feature blueprint — every analysis module has its own page and its own JSON API."""
from flask import Blueprint, jsonify, render_template, request

from services import analysis as A, evolution as E, external as X, proximity as PX, structure as S
from services.catalog import PROTEIN_CATALOG
from services.dbutil import db

bp = Blueprint("features", __name__)

FEATURES = [  # (endpoint, nav label, icon, one-line purpose)
    ("hotspots", "PTM Hotspots", "🔥", "Statistically enriched PTM-dense regions"),
    ("spatial", "3D Spatial Distance", "📏", "Distances between PTM residues in real structures"),
    ("network", "Kinase–PTM Network", "🕸️", "Kinase → protein → PTM → cancer context"),
    ("conservation", "Conservation Explorer", "🧬", "Site conservation across orthologues"),
    ("mutations", "Mutation–PTM–Cancer", "🎯", "Cancer variants at or near PTM sites"),
    ("proximity", "3D-Proximal Mutations", "🧩", "Cancer mutations far in sequence but close in 3D to PTM sites"),
    ("crosstalk", "PTM Crosstalk", "🔄", "Different PTM types at nearby residues"),
    ("hypotheses", "Hypothesis Generator", "💡", "Evidence-linked research questions"),
    ("priority", "PhosphoRank", "⭐", "Transparent site prioritisation"),
    ("integrity", "Data Integrity Audit", "✅", "Every record verified against UniProt"),
]


@bp.app_context_processor
def inject_nav():
    return {"FEATURES": FEATURES}


def _proteins():
    """Catalogue proteins plus any protein present in the imported records (never hides real data)."""
    conn = db()
    counts = {r["uniprot"]: r["n"] for r in conn.execute("SELECT uniprot,COUNT(*) n FROM ptm_sites GROUP BY uniprot")}
    extra = conn.execute("SELECT uniprot, MIN(gene) g, MAX(position) L FROM ptm_sites GROUP BY uniprot").fetchall()
    conn.close()
    out = [{"gene": g, "uniprot": u, "name": n, "length": l, "sites": counts.get(u, 0)} for g, u, n, l in PROTEIN_CATALOG]
    known = {p["uniprot"] for p in out}
    out += [{"gene": r["g"], "uniprot": r["uniprot"], "name": None, "length": r["L"], "sites": counts[r["uniprot"]]}
            for r in extra if r["uniprot"] not in known]
    with_data = [p for p in out if p["sites"] > 0]  # hide proteins with no records yet
    return with_data or out


def _load(uniprot):
    conn = db()
    rows = conn.execute("SELECT * FROM ptm_sites WHERE uniprot=? ORDER BY position", (uniprot,)).fetchall()
    conn.close()
    cat = next((c for c in PROTEIN_CATALOG if c[1] == uniprot), None)
    if not cat and not rows:
        return None
    prot = ({"gene": cat[0], "uniprot": cat[1], "name": cat[2], "length": cat[3]} if cat else
            {"gene": rows[0]["gene"], "uniprot": uniprot, "name": None,
             "length": max((r["position"] or 0) for r in rows) or None})
    return prot, rows, A.units_from_rows(rows)


def _int(name, default, lo, hi):
    try:
        return max(lo, min(hi, int(request.args.get(name, default))))
    except ValueError:
        return default


def _float(name, default, lo, hi):
    try:
        return max(lo, min(hi, float(request.args.get(name, default))))
    except ValueError:
        return default


def _page(name, title, blurb):
    protos = _proteins()
    sel = request.args.get("uniprot") or next((p["uniprot"] for p in protos if p["sites"]), protos[0]["uniprot"])
    total = sum(p["sites"] for p in protos)
    info = next(f for f in FEATURES if f[0] == name)
    return render_template(f"feature_{name}.html", page=name, title=title, blurb=blurb, proteins=protos,
                           selected=sel, db_empty=(total == 0), icon=info[2])


# ------------------------------------------------------------------ pages
@bp.route("/features")
def features_index():
    conn = db()
    n = conn.execute("SELECT COUNT(*) FROM ptm_sites").fetchone()[0]
    conn.close()
    return render_template("features.html", db_empty=(n == 0))


@bp.route("/hotspots")
def hotspots_page():
    return _page("hotspots", "PTM Hotspot Analysis",
                 "Finds regions where modified residues cluster more than a uniform placement would predict.")


@bp.route("/spatial")
def spatial_page():
    return _page("spatial", "PTM Spatial-Distance Analysis",
                 "Measures 3D distances between PTM residues using real experimental (or clearly labelled predicted) coordinates.")


@bp.route("/network")
def network_page():
    return _page("network", "Kinase–PTM Network",
                 "Kinase → protein → PTM site → cancer context, built only from records that name a kinase.")


@bp.route("/conservation")
def conservation_page():
    return _page("conservation", "Conservation Explorer",
                 "Checks whether each PTM residue is identical in Ensembl Compara orthologues.")


@bp.route("/mutations")
def mutations_page():
    return _page("mutations", "Mutation–PTM–Cancer Analysis",
                 "Joins real cancer-associated variants and cBioPortal cohort mutations to PTM residues.")


@bp.route("/proximity")
def proximity_page():
    return _page("proximity", "3D-Proximal Mutation–PTM Analysis",
                 "Cancer missense mutations that sit far from every PTM site in sequence but close to one in 3D — invisible to sequence-only tools.")


@bp.route("/crosstalk")
def crosstalk_page():
    return _page("crosstalk", "PTM Crosstalk Analysis",
                 "Different PTM types on the same or neighbouring residues — candidate crosstalk, not proof.")


@bp.route("/hypotheses")
def hypotheses_page():
    return _page("hypotheses", "Research-Hypothesis Generator",
                 "Turns computed evidence into ranked, testable research questions with the evidence that supports each.")


@bp.route("/priority")
def priority_page():
    return _page("priority", "PhosphoRank — Site Prioritisation",
                 "A transparent, adjustable checklist score for choosing which sites to follow up first.")


@bp.route("/integrity")
def integrity_page():
    return _page("integrity", "Data Integrity Audit",
                 "Every stored record is checked against the UniProt sequence: residue letter and flanking window.")


# -------------------------------------------------------------------- APIs
def _guard(uniprot):
    got = _load(uniprot)
    if not got:
        return None, (jsonify({"error": "Protein not found"}), 404)
    if not got[1]:
        return None, (jsonify({"empty": True, "message": "No PTM records for this protein in the local database. "
                               "Run sync_data.py (or import_tsv.py) first."}), 200)
    return got, None


def _ext_error(e):
    return jsonify({"available": False, "reason": str(e)}), 200


@bp.route("/api/f/hotspots/<uniprot>")
def api_hotspots(uniprot):
    got, err = _guard(uniprot)
    if err:
        return err
    prot, rows, units = got
    w = _int("window", 15, 5, 60)
    res = A.hotspots(units, prot["length"], w, _int("min_sites", 3, 2, 12))
    try:
        seq = X.uniprot_sequence(uniprot)
        maxes = PX.composition_null(seq, sorted({u["pos"] for u in units}), w)
        for r in res["regions"]:
            r["perm_p_fwer"] = PX.perm_p(maxes, r["best_window_count"])
        res["method"] += (" Also: family-wise permutation p (1000 re-placements keeping each residue letter's "
                          "composition; statistic = max sites in any window).")
    except X.ExternalError:
        res["method"] += " Composition-preserving permutation p not available (UniProt sequence unreachable)."
    return jsonify({"protein": prot, **res})


@bp.route("/api/f/crosstalk/<uniprot>")
def api_crosstalk(uniprot):
    got, err = _guard(uniprot)
    if err:
        return err
    prot, rows, units = got
    return jsonify({"protein": prot, **A.crosstalk(units, _int("max_dist", 10, 0, 50),
                                                   request.args.get("same_type") == "1")})


@bp.route("/api/f/network/<uniprot>")
def api_network(uniprot):
    """uniprot='ALL' returns the database-wide network."""
    conn = db()
    if uniprot == "ALL":
        rows = conn.execute("SELECT * FROM ptm_sites").fetchall()
    else:
        rows = conn.execute("SELECT * FROM ptm_sites WHERE uniprot=?", (uniprot,)).fetchall()
    genes = [r["gene"] for r in conn.execute("SELECT DISTINCT gene FROM ptm_sites")]
    conn.close()
    if not rows:
        return jsonify({"empty": True, "message": "No PTM records available. Run sync_data.py first."})
    return jsonify(A.kinase_network(rows, genes))


@bp.route("/api/f/spatial/<uniprot>")
def api_spatial(uniprot):
    got, err = _guard(uniprot)
    if err:
        return err
    prot, rows, units = got
    try:
        res = S.spatial_analysis(uniprot, units, _float("threshold", 10.0, 4, 30),
                                 request.args.get("alphafold", "1") == "1")
    except X.ExternalError as e:
        return _ext_error(e)
    return jsonify({"protein": prot, **res})


@bp.route("/api/f/spatial/<uniprot>/model")
def api_spatial_model(uniprot):
    """The exact coordinate file (cached real PDB/AlphaFold text) behind a spatial result, for the 3D view."""
    src, sid = request.args.get("source"), request.args.get("id", "")
    try:
        if src == "experimental" and sid.isalnum() and len(sid) == 4:
            return X.pdb_text(sid), 200, {"Content-Type": "text/plain"}
        if src == "alphafold":
            return X.alphafold_model(uniprot)[1], 200, {"Content-Type": "text/plain"}
    except X.ExternalError as e:
        return str(e), 502
    return "bad request", 400


@bp.route("/api/f/conservation/<uniprot>")
def api_conservation(uniprot):
    got, err = _guard(uniprot)
    if err:
        return err
    prot, rows, units = got
    try:
        return jsonify({"protein": prot, "available": True, **E.conservation(uniprot, prot["gene"], units)})
    except X.ExternalError as e:
        return _ext_error(e)


@bp.route("/api/f/variants/<uniprot>")
def api_variants(uniprot):
    got, err = _guard(uniprot)
    if err:
        return err
    prot, rows, units = got
    try:
        variants = X.ebi_variants(uniprot)
    except X.ExternalError as e:
        return _ext_error(e)
    return jsonify({"protein": prot, "available": True, **A.variant_ptm_overlap(units, variants, _int("window", 0, 0, 10))})


@bp.route("/api/cbio/studies")
def api_cbio_studies():
    try:
        return jsonify({"available": True, "studies": X.cbio_studies()})
    except X.ExternalError as e:
        return _ext_error(e)


@bp.route("/api/cbio/suggest/<uniprot>")
def api_cbio_suggest(uniprot):
    got = _load(uniprot)
    if not got:
        return jsonify({"error": "Protein not found"}), 404
    try:
        return jsonify({"available": True, "gene": got[0]["gene"], **X.cbio_suggest_studies(got[0]["gene"])})
    except X.ExternalError as e:
        return _ext_error(e)


@bp.route("/api/f/cohort/<uniprot>")
def api_cohort(uniprot):
    got, err = _guard(uniprot)
    if err:
        return err
    prot, rows, units = got
    study = request.args.get("study", "").strip()
    if not study:
        return jsonify({"available": False, "reason": "Select a cBioPortal study first; no cohort is assumed."})
    try:
        cb = X.cbio_gene_mutations(study, prot["gene"])
    except X.ExternalError as e:
        return _ext_error(e)
    return jsonify({"protein": prot, "available": True,
                    **E.cohort_overlap(units, prot["length"], study, cb, _int("window", 0, 0, 10))})


@bp.route("/api/f/proximity/<uniprot>")
def api_proximity(uniprot):
    got, err = _guard(uniprot)
    if err:
        return err
    prot, rows, units = got
    study = request.args.get("study", "").strip()
    if not study:
        return jsonify({"available": False, "reason": "Select a cBioPortal study first; no cohort is assumed."})
    mode = request.args.get("dist_mode", "ca")
    if mode not in PX.DIST_MODES:
        return jsonify({"available": False, "reason": f"Unknown distance mode '{mode}'. Use one of {list(PX.DIST_MODES)}."}), 200
    try:
        cb = X.cbio_gene_mutations(study, prot["gene"])
    except X.MutationRetrievalError as e:
        # The answer is UNKNOWN here (not "zero mutations"): say exactly which retrieval step failed.
        return jsonify({"available": False, "reason": str(e), "retrieval_status": e.status,
                        "retrieval": e.retrieval, "study": study,
                        "notes": ["Mutation retrieval failed; no mutation data from any other study was substituted."]}), 200
    except X.ExternalError as e:
        return _ext_error(e)
    try:
        res = PX.proximity_analysis(uniprot, units, cb, _float("threshold", 10.0, 1, 30), _int("min_sep", 20, 5, 200),
                                    n_perm=_int("n_perm", 5000, 0, PX.MAX_PERM), perm_seed=_int("seed", 1, 0, 2**31 - 1),
                                    dist_mode=mode, plddt_min=_float("plddt_min", PX.DEFAULT_PLDDT_MIN, 0, 100))
    except X.ExternalError as e:
        return jsonify({"available": False, "reason": f"Structure retrieval failed: {e}", "study": study,
                        "retrieval": cb.get("retrieval"),
                        "notes": ["Mutations were retrieved, but the AlphaFold model could not be; nothing was inferred."]}), 200
    return jsonify({"protein": prot, "study": study, **res})


@bp.route("/api/f/integrity/<uniprot>")
def api_integrity(uniprot):
    got, err = _guard(uniprot)
    if err:
        return err
    prot, rows, units = got
    try:
        return jsonify({"protein": prot, "available": True, **E.integrity(uniprot, rows, prot["length"])})
    except X.ExternalError as e:
        return _ext_error(e)


@bp.route("/api/f/priority/<uniprot>")
def api_priority(uniprot):
    weights = {k: request.args.get(k) for k in E.DEFAULT_W if request.args.get(k) not in (None, "")}
    if uniprot == "ALL":
        conn = db()
        ups = [r[0] for r in conn.execute("SELECT DISTINCT uniprot FROM ptm_sites")]
        conn.close()
    else:
        ups = [uniprot]
    out = []
    for u in ups:
        got = _load(u)
        if not got or not got[1]:
            continue
        prot, rows, units = got
        hot = A.hotspots(units, prot["length"])
        cross = A.crosstalk(units, 10)
        out += E.prioritise(units, hot, cross, weights)
    if not out:
        return jsonify({"empty": True, "message": "No PTM records available. Run sync_data.py first."})
    out.sort(key=lambda r: (-r["score"], r["gene"], r["pos"]))
    return jsonify({"weights": {**E.DEFAULT_W, **{k: float(v) for k, v in weights.items()}}, "sites": out[:500],
                    "n_sites": len(out)})


@bp.route("/api/f/hypotheses/<uniprot>")
def api_hypotheses(uniprot):
    got, err = _guard(uniprot)
    if err:
        return err
    prot, rows, units = got
    hot = A.hotspots(units, prot["length"])
    cross = A.crosstalk(units, _int("max_dist", 10, 0, 50))
    layers = {"records": {"used": True}}
    spatial = variants = cons = None
    if request.args.get("spatial") == "1":
        try:
            spatial = S.spatial_analysis(uniprot, units)
            layers["spatial"] = {"used": bool(spatial.get("available")), "reason": spatial.get("reason")}
        except X.ExternalError as e:
            layers["spatial"] = {"used": False, "reason": str(e)}
    if request.args.get("variants") == "1":
        try:
            variants = A.variant_ptm_overlap(units, X.ebi_variants(uniprot), 2)
            layers["variants"] = {"used": True}
        except X.ExternalError as e:
            layers["variants"] = {"used": False, "reason": str(e)}
    if request.args.get("conservation") == "1":
        try:
            cons = E.conservation(uniprot, prot["gene"], units)
            layers["conservation"] = {"used": True}
        except X.ExternalError as e:
            layers["conservation"] = {"used": False, "reason": str(e)}
    hyps = A.generate_hypotheses(prot["gene"], units, hot, cross, variants, cons, spatial)
    return jsonify({"protein": prot, "layers": layers, "n": len(hyps), "hypotheses": hyps,
                    "policy": "Hypotheses are research questions derived from the listed evidence. They are not results."})

