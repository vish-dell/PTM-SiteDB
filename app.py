
from flask import Flask, render_template, request, jsonify, Response, redirect, url_for
import sqlite3, csv, io, json, re, math
from collections import Counter
from pathlib import Path


from services.catalog import PROTEIN_CATALOG
from services.dbutil import db, ensure_schema, DB
from features import bp as features_bp

app = Flask(__name__)
from conservation_api import bp as conservation_bp

app.register_blueprint(conservation_bp)
app.register_blueprint(features_bp)

def init_db():
    # Schema only. The database is filled exclusively by sync_data.py / import_tsv.py from the real source.
    ensure_schema()

def rows_from_db(params=None):
    conn=db()
    sql="SELECT * FROM ptm_sites WHERE 1=1"
    vals=[]
    p=params or {}
    for key,col in [("q",None),("gene","gene"),("ptm_type","ptm_type"),("confidence","confidence"),
                    ("evidence","evidence"),("kinase","kinase"),("disease","disease_context")]:
        v=(p.get(key) or "").strip()
        if not v: continue
        if key=="q":
            sql += " AND (gene LIKE ? OR uniprot LIKE ? OR site LIKE ? OR ptm_type LIKE ? OR kinase LIKE ? OR disease_context LIKE ? OR pmid LIKE ?)"
            vals += [f"%{v}%"]*7
        else:
            sql += f" AND {col} LIKE ?"
            vals.append(f"%{v}%")
    sql += " ORDER BY gene, position, id"
    rows=conn.execute(sql,vals).fetchall()
    conn.close()
    return rows

@app.route("/")
def home():
    conn=db()
    stats={
        "sites": conn.execute("SELECT COUNT(*) FROM ptm_sites").fetchone()[0],
        "proteins": conn.execute("SELECT COUNT(DISTINCT uniprot) FROM ptm_sites").fetchone()[0],
        "catalogue_proteins": len(PROTEIN_CATALOG),
        "kinases": conn.execute("SELECT COUNT(DISTINCT kinase) FROM ptm_sites WHERE kinase IS NOT NULL AND trim(kinase)!='' AND kinase!='—'").fetchone()[0],
        "publications": conn.execute("SELECT COUNT(DISTINCT pmid) FROM ptm_sites WHERE pmid IS NOT NULL AND trim(pmid)!=''").fetchone()[0],
    }
    conn.close()
    return render_template("home.html", stats=stats, db_empty=(stats["sites"]==0))

@app.route("/search")
def search():
    rows=rows_from_db(request.args)
    return render_template("search.html", rows=rows, args=request.args)

@app.route("/browse")
def browse():
    view = request.args.get("view", "proteins")
    if view == "sites":
        rows = rows_from_db({})
        return render_template("browse.html", view=view, rows=rows,
                                total_proteins=len(PROTEIN_CATALOG))
    conn = db()
    counts = {r["uniprot"]: r["n"] for r in
              conn.execute("SELECT uniprot, COUNT(*) n FROM ptm_sites GROUP BY uniprot")}
    conn.close()
    proteins = [{"gene": g, "uniprot": u, "name": n, "length": l, "sites": counts.get(u, 0)}
                for g, u, n, l in PROTEIN_CATALOG]
    return render_template("browse.html", view=view, proteins=proteins,
                            total_proteins=len(PROTEIN_CATALOG))

@app.route("/statistics")
def statistics():
    conn=db()
    stats={
        "sites": conn.execute("SELECT COUNT(*) FROM ptm_sites").fetchone()[0],
        "proteins": conn.execute("SELECT COUNT(DISTINCT uniprot) FROM ptm_sites").fetchone()[0],
        "catalogue_proteins": len(PROTEIN_CATALOG),
        "kinases": conn.execute("SELECT COUNT(DISTINCT kinase) FROM ptm_sites WHERE kinase!='' AND kinase!='—'").fetchone()[0],
        "publications": conn.execute("SELECT COUNT(DISTINCT pmid) FROM ptm_sites WHERE pmid!=''").fetchone()[0],
    }
    types=conn.execute("SELECT ptm_type,COUNT(*) n FROM ptm_sites GROUP BY ptm_type ORDER BY n DESC").fetchall()
    genes=conn.execute("SELECT gene,COUNT(*) n FROM ptm_sites GROUP BY gene ORDER BY n DESC").fetchall()
    evidence=conn.execute("""SELECT evidence,COUNT(*) n FROM ptm_sites
                              WHERE evidence IS NOT NULL AND trim(evidence)!=''
                              GROUP BY evidence ORDER BY n DESC""").fetchall()
    confidence=conn.execute("""SELECT confidence,COUNT(*) n FROM ptm_sites
                                WHERE confidence IS NOT NULL AND trim(confidence)!=''
                                GROUP BY confidence ORDER BY n DESC""").fetchall()
    kinase_rows=conn.execute("""SELECT kinase FROM ptm_sites
                                 WHERE kinase IS NOT NULL AND trim(kinase)!='' AND kinase!='—'""").fetchall()
    disease_rows=conn.execute("""SELECT disease_context FROM ptm_sites
                                  WHERE disease_context IS NOT NULL AND trim(disease_context)!=''""").fetchall()
    conn.close()

    kinase_counts = Counter(r["kinase"].strip() for r in kinase_rows)
    kinases = [{"kinase": k, "n": n} for k, n in kinase_counts.most_common(10)]

    disease_counts = Counter()
    for r in disease_rows:
        for d in r["disease_context"].split(";"):
            d = d.strip()
            if d:
                disease_counts[d] += 1
    diseases = [{"disease": d, "n": n} for d, n in disease_counts.most_common(6)]

    return render_template("statistics.html", stats=stats, types=types, genes=genes,
                            evidence=evidence, confidence=confidence,
                            kinases=kinases, diseases=diseases)

@app.route("/protein/<uniprot>")
def protein(uniprot):
    conn=db()
    rows=conn.execute("SELECT * FROM ptm_sites WHERE uniprot=? ORDER BY position", (uniprot,)).fetchall()
    conn.close()
    cat = next((c for c in PROTEIN_CATALOG if c[1] == uniprot), None)
    if rows:
        prot = dict(rows[0])
        if cat:
            prot["name"] = cat[2]
            prot["length"] = cat[3]
        else:
            prot["name"] = None
            prot["length"] = max((r["position"] or 0) for r in rows) or None
    else:
        if not cat:
            return "Protein not found", 404
        prot = {"gene": cat[0], "uniprot": cat[1], "name": cat[2], "length": cat[3]}

    pmdata = {
        "gene": prot.get("gene"), "uniprot": prot.get("uniprot"), "length": prot.get("length"),
        "ptms": [{"id": r["id"], "site": r["site"], "ptm_type": r["ptm_type"], "position": r["position"],
                  "kinase": r["kinase"], "disease_context": r["disease_context"], "evidence": r["evidence"],
                  "confidence": r["confidence"], "pmid": r["pmid"]} for r in rows],
    }
    return render_template("protein.html", rows=rows, protein=prot, pmdata=pmdata)

@app.route("/viewer/<uniprot>")
def viewer(uniprot):
    conn = db()
    rows = conn.execute("SELECT * FROM ptm_sites WHERE uniprot=? ORDER BY position", (uniprot,)).fetchall()
    conn.close()
    cat = next((c for c in PROTEIN_CATALOG if c[1] == uniprot), None)
    if not rows and not cat:
        return "Protein not found", 404
    if rows:
        protein_info = dict(rows[0])
        if cat:
            protein_info.update({"gene": cat[0], "uniprot": cat[1], "name": cat[2], "length": cat[3]})
        else:
            protein_info.setdefault("name", None)
            protein_info["length"] = max((r["position"] or 0 for r in rows), default=0) or None
    else:
        protein_info = {"gene": cat[0], "uniprot": cat[1], "name": cat[2], "length": cat[3]}
    return render_template("viewer.html", protein=protein_info, ptms=[dict(r) for r in rows])

@app.route("/api/search")
def api_search():
    rows=rows_from_db(request.args)
    return jsonify([dict(r) for r in rows])

@app.route("/api/protein/<uniprot>")
def api_protein(uniprot):
    conn=db()
    rows=conn.execute("SELECT * FROM ptm_sites WHERE uniprot=? ORDER BY position",(uniprot,)).fetchall()
    conn.close()
    if not rows: return jsonify({"error":"Protein not found"}),404
    return jsonify({"protein":dict(rows[0]),"ptms":[dict(r) for r in rows]})

@app.route("/api/protein/<uniprot>/structure")
def api_structure(uniprot):
    return jsonify({"uniprot":uniprot,"policy":{
        "experimental_first":True,
        "alphafold_fallback":True,
        "never_fake_mapping":True
    }})

@app.route("/api/statistics")
def api_statistics():
    conn=db()
    out={
        "sites":conn.execute("SELECT COUNT(*) FROM ptm_sites").fetchone()[0],
        "proteins":conn.execute("SELECT COUNT(DISTINCT uniprot) FROM ptm_sites").fetchone()[0],
        "catalogue_proteins":len(PROTEIN_CATALOG),
        "kinases":conn.execute("SELECT COUNT(DISTINCT kinase) FROM ptm_sites WHERE kinase!='' AND kinase!='—'").fetchone()[0],
        "publications":conn.execute("SELECT COUNT(DISTINCT pmid) FROM ptm_sites WHERE pmid!=''").fetchone()[0],
        "ptm_types":[dict(r) for r in conn.execute("SELECT ptm_type,COUNT(*) n FROM ptm_sites GROUP BY ptm_type ORDER BY n DESC")]
    }
    conn.close()
    return jsonify(out)

@app.route("/api/download")
def api_download():
    rows=rows_from_db({})
    output=io.StringIO()
    fields=["gene","uniprot","ptm_type","site","position","flanking","kinase","disease_context","evidence","confidence","pmid"]
    w=csv.DictWriter(output,fieldnames=fields,delimiter="\t")
    w.writeheader()
    for r in rows: w.writerow({f:r[f] for f in fields})
    return Response(output.getvalue(),mimetype="text/tab-separated-values",
                    headers={"Content-Disposition":"attachment; filename=ptm_sitedb.tsv"})

@app.route("/api/download.json")
def api_download_json():
    rows=rows_from_db({})
    fields=["gene","uniprot","ptm_type","site","position","flanking","kinase","disease_context","evidence","confidence","pmid"]
    data=[{f:r[f] for f in fields} for r in rows]
    return Response(json.dumps(data, indent=2), mimetype="application/json",
                    headers={"Content-Disposition":"attachment; filename=ptm_sitedb.json"})

@app.route("/api/download.sql")
def api_download_sql():
    conn=db()
    rows=conn.execute("SELECT * FROM ptm_sites ORDER BY id").fetchall()
    conn.close()
    def esc(v):
        if v is None: return "NULL"
        if isinstance(v, int): return str(v)
        return "'" + str(v).replace("'", "''") + "'"
    fields=["gene","uniprot","ptm_type","site","position","flanking","kinase","disease_context","evidence","confidence","pmid","source_url","protein_id"]
    lines=["-- PTM-SiteDB SQL dump", "-- Generated from the local ptm_sites table. No records invented.", "",
    """CREATE TABLE IF NOT EXISTS ptm_sites (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    gene TEXT NOT NULL, uniprot TEXT NOT NULL, ptm_type TEXT NOT NULL,
    site TEXT NOT NULL, position INTEGER, flanking TEXT, kinase TEXT,
    disease_context TEXT, evidence TEXT, confidence TEXT, pmid TEXT,
    source_url TEXT, protein_id TEXT,
    UNIQUE(uniprot, ptm_type, site, pmid)
);""", ""]
    cols=", ".join(fields)
    for r in rows:
        vals=", ".join(esc(r[f]) for f in fields)
        lines.append(f"INSERT INTO ptm_sites ({cols}) VALUES ({vals});")
    return Response("\n".join(lines), mimetype="application/sql",
                    headers={"Content-Disposition":"attachment; filename=ptm_sitedb.sql"})

@app.route("/api/download/protein/<gene>.<fmt>")
def api_download_protein(gene, fmt):
    conn=db()
    rows=conn.execute("SELECT * FROM ptm_sites WHERE gene=? ORDER BY position",(gene,)).fetchall()
    conn.close()
    fields=["gene","uniprot","ptm_type","site","position","flanking","kinase","disease_context","evidence","confidence","pmid"]
    if fmt == "json":
        data=[{f:r[f] for f in fields} for r in rows]
        return Response(json.dumps(data, indent=2), mimetype="application/json",
                        headers={"Content-Disposition":f"attachment; filename={gene}_ptm_sites.json"})
    output=io.StringIO()
    w=csv.DictWriter(output,fieldnames=fields,delimiter="\t")
    w.writeheader()
    for r in rows: w.writerow({f:r[f] for f in fields})
    return Response(output.getvalue(),mimetype="text/tab-separated-values",
                    headers={"Content-Disposition":f"attachment; filename={gene}_ptm_sites.tsv"})

@app.route("/download")
def download():
    conn=db()
    total=conn.execute("SELECT COUNT(*) FROM ptm_sites").fetchone()[0]
    conn.close()
    genes=sorted({g for g,_,_,_ in PROTEIN_CATALOG})
    return render_template("download.html", total=total, genes=genes)

@app.route("/help")
def help_page():
    return render_template("help.html")

@app.route("/connection-check")
def connection_check():
    """Local diagnostics for database and research-source readiness."""
    conn = db()
    site_count = int(conn.execute("SELECT COUNT(*) FROM ptm_sites").fetchone()[0])
    gene_count = int(conn.execute("SELECT COUNT(DISTINCT gene) FROM ptm_sites").fetchone()[0])
    last_sync = conn.execute("SELECT synced_at, source, pages, records_parsed, records_new FROM sync_log ORDER BY id DESC LIMIT 1").fetchone()
    conn.close()
    research_root = Path(__file__).resolve().parent / "research" / "data" / "processed"
    research_files = {
        "dbPTM source table": research_root / "dbptm_phosphorylation_source.txt",
        "IntOGen driver table": research_root / "intogen_drivers_source.tsv",
        "PTM panel": research_root / "benchmark_panel.csv",
        "source manifest": research_root / "source_manifest.json",
    }
    research_ready = {k: {"exists": v.exists(), "bytes": v.stat().st_size if v.exists() else 0}
                      for k, v in research_files.items()}
    return render_template("connection.html", db_path=str(DB), site_count=site_count, gene_count=gene_count,
                           last_sync=dict(last_sync) if last_sync else None, research_ready=research_ready,
                           research_root=str(research_root))

@app.route("/api/health")
def api_health():
    conn = db()
    try:
        sites = int(conn.execute("SELECT COUNT(*) FROM ptm_sites").fetchone()[0])
        genes = int(conn.execute("SELECT COUNT(DISTINCT gene) FROM ptm_sites").fetchone()[0])
        last = conn.execute("SELECT synced_at, source, pages, records_parsed, records_new FROM sync_log ORDER BY id DESC LIMIT 1").fetchone()
    finally:
        conn.close()
    return jsonify({"ok": True, "database_path": str(DB), "ptm_records": sites, "distinct_genes": genes,
                    "database_populated": sites > 0, "last_sync": dict(last) if last else None})

init_db()

# The original PTM-SiteDB Flask app stays at the root with its original templates/CSS.
# The expanded cohort workflow is a separate Flask application mounted below /research,
# preventing route collisions and keeping the original database/search/Analysis Suite intact.
from werkzeug.middleware.dispatcher import DispatcherMiddleware
from research.app import app as research_app

# Keep the public Flask object and its test client; wrap its WSGI handler to dispatch
# /research/* requests to the independent research app.
app.wsgi_app = DispatcherMiddleware(app.wsgi_app, {"/research": research_app})

if __name__ == "__main__":
    app.run(debug=False, host="127.0.0.1", port=5000, use_reloader=False)
