import sqlite3, re, csv, sys, time, collections, urllib.request, urllib.parse
import xml.etree.ElementTree as ET
from pathlib import Path

DB = "instance/ptm.db"
EXPECT = {"P04637":"TP53","P31749":"AKT1","P35222":"CTNNB1","P00533":"EGFR","P28482":"MAPK1",
          "P15056":"BRAF","P38398":"BRCA1","O14757":"CHEK1","P42336":"PIK3CA"}
GENE_WORDS = {"P04637":["p53","tp53"],"P31749":["akt","pkb"],"P35222":["catenin"],"P00533":["egfr","epidermal growth factor"],
              "P28482":["erk","mapk","raf"],"P15056":["raf"],"P38398":["brca1"],"O14757":["chk1","chek1","checkpoint","atm"],
              "P42336":["pi3k","p110","phosphoinositide","pik3ca"]}
GENERIC = ["phospho","proteom","kinase","signal","mass spectrom","cancer","tumor","cell cycle","atlas","enrich"]
VALID_TYPES = {"Phosphorylation","Acetylation","Ubiquitination","Methylation","SUMOylation","Glycosylation"}
results = []   # (section, status, message)
def log(sec, ok, msg): results.append((sec, "PASS" if ok else "FAIL", msg)); print(("PASS " if ok else "FAIL ") + f"[{sec}] {msg}")
def note(sec, msg): results.append((sec, "NOTE", msg)); print(f"NOTE [{sec}] {msg}")

c = sqlite3.connect(DB); c.row_factory = sqlite3.Row
rows = [dict(r) for r in c.execute("select * from ptm_sites")]
print(f"Database: {DB} | records {len(rows)}\n")

# ---------- 1. COUNTS ----------
uniq = {(r["uniprot"], r["site"]) for r in rows}
prots = {r["uniprot"] for r in rows}
kin = {r["kinase"] for r in rows if r["kinase"] and len(r["kinase"].strip()) > 1}
pm_all = set()
for r in rows: pm_all |= set(re.findall(r"\d+", r["pmid"] or ""))
note("COUNTS", f"records={len(rows)} unique_sites={len(uniq)} proteins={len(prots)} kinases={len(kin)} publications={len(pm_all)}")

# ---------- 2. SITES (offline part) ----------
bad = []
for r in rows:
    p = []
    if EXPECT.get(r["uniprot"]) != r["gene"]: p.append("accession/gene mismatch")
    m = re.fullmatch(r"([A-Z])(\d+)", r["site"] or "")
    if not m: p.append("bad site format")
    else:
        if r["position"] is not None and int(r["position"]) != int(m.group(2)): p.append("position != site number")
        if r["ptm_type"] == "Phosphorylation" and m.group(1) not in "STY": p.append("phospho not on S/T/Y")
    if r["ptm_type"] not in VALID_TYPES: p.append("invalid ptm_type")
    if p: bad.append((r["gene"], r["site"], p))
log("SITES", not bad, f"format/position/type check: {len(bad)} problem(s)")
for g, s, p in bad: print("     ", g, s, p)
dup = collections.Counter((r["uniprot"], r["site"], r["ptm_type"], r["pmid"]) for r in rows)
d = [k for k, v in dup.items() if v > 1]
log("SITES", not d, f"exact duplicate rows (same protein+site+type+PMID): {len(d)}")
note("SITES", f"{len(rows)-len(uniq)} rows are extra PMIDs for an already-listed site (allowed)")

# ---------- 3. PROTEINS / UNIPROT SEQUENCE ----------
def fetch(url, t=60): return urllib.request.urlopen(url, timeout=t).read().decode()
seqs = {}
try:
    for a in prots:
        seqs[a] = "".join(fetch(f"https://rest.uniprot.org/uniprotkb/{a}.fasta").splitlines()[1:])
    log("PROTEINS", all(re.fullmatch(r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}", a) for a in prots), "accession format valid")
    note("PROTEINS", "lengths: " + ", ".join(f"{EXPECT.get(a,a)}={len(s)}" for a, s in sorted(seqs.items())))
    sb = []
    for r in rows:
        m = re.fullmatch(r"([A-Z])(\d+)", r["site"] or ""); s = seqs.get(r["uniprot"], "")
        if not m: continue
        aa, n = m.group(1), int(m.group(2))
        if n > len(s): sb.append((r["gene"], r["site"], "beyond protein length")); continue
        if s[n-1] != aa: sb.append((r["gene"], r["site"], f"UniProt has {s[n-1]}")); continue
        f = (r["flanking"] or "").replace("-", "").replace("_", "").upper()
        if f and f not in s: sb.append((r["gene"], r["site"], "flank not in sequence"))
    log("SITES", not sb, f"residue + flank vs UniProt sequence: {len(sb)} problem(s)")
    for x in sb: print("     ", x)
except Exception as e:
    note("PROTEINS", f"UniProt check skipped (offline?): {e}")

# ---------- 4. KINASES ----------
placeholders = collections.Counter(r["kinase"] for r in rows if r["kinase"] is not None and (len(r["kinase"].strip()) <= 1 or r["kinase"].strip().lower() in ("none","nan","null","n/a")))
note("KINASES", f"placeholder values (blank/dash/None/nan): {dict(placeholders)}")
real = [r for r in rows if r["kinase"] and len(r["kinase"].strip()) > 1]
case = collections.defaultdict(set)
for r in real: case[r["kinase"].strip().lower()].add(r["kinase"])
mixed = {k: v for k, v in case.items() if len(v) > 1}
spaces = [r["kinase"] for r in real if r["kinase"] != r["kinase"].strip() or "  " in r["kinase"]]
numeric = [r["kinase"] for r in real if r["kinase"].strip().isdigit()]
log("KINASES", not mixed and not spaces and not numeric, f"mixed-case dupes={len(mixed)}, extra spaces={len(spaces)}, numeric={len(numeric)}")
TYR_ONLY = {"FGFR1","SRC","ABL1","EGFR","JAK2","FYN","LCK"}
mis = [(r["gene"], r["site"], r["kinase"]) for r in real if r["kinase"] in TYR_ONLY and r["site"][0] in "ST"]
log("KINASES", not mis, f"tyrosine-kinase on S/T site: {len(mis)}")
for x in mis: print("     ", x)
note("KINASES", "records with kinase: " + str(len(real)) + " -> " + ", ".join(sorted({r['kinase'] for r in real})))

# ---------- 5. PMIDs ----------
pm_bad = []
for r in rows:
    raw = r["pmid"] or ""
    ids = re.findall(r"\d+", raw)
    if not ids: pm_bad.append((r["gene"], r["site"], "no PMID"))
    if len(ids) != len(set(ids)): pm_bad.append((r["gene"], r["site"], "duplicate PMID in cell"))
    if re.search(r"\s", raw.strip()) and ";" not in raw and "," not in raw: pm_bad.append((r["gene"], r["site"], "space in PMID"))
    for p in ids:
        if not (1 <= int(p) <= 45_000_000): pm_bad.append((r["gene"], r["site"], f"implausible PMID {p}"))
log("PMIDS", not pm_bad, f"format check: {len(pm_bad)} problem(s)")
for x in pm_bad: print("     ", x)
try:
    pm = sorted(pm_all); info = {}
    for i in range(0, len(pm), 100):
        u = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + urllib.parse.urlencode({"db":"pubmed","id":",".join(pm[i:i+100]),"retmode":"xml"})
        root = ET.fromstring(fetch(u).encode())
        for art in root.findall(".//PubmedArticle"):
            pid = art.findtext(".//MedlineCitation/PMID"); t = art.find(".//ArticleTitle")
            info[pid] = {"title": "".join(t.itertext()) if t is not None else "",
                         "retracted": any("Retract" in (x.text or "") for x in art.findall(".//PublicationType"))}
        time.sleep(0.4)
    missing = [p for p in pm if p not in info]
    log("PMIDS", not missing, f"PMIDs not found on PubMed: {missing}")
    retr = [p for p, v in info.items() if v["retracted"]]
    log("PMIDS", not retr, f"retracted PMIDs: {retr}")
    off = []
    for r in rows:
        for p in re.findall(r"\d+", r["pmid"] or ""):
            t = info.get(p, {}).get("title", "").lower()
            if t and not any(k in t for k in GENERIC + GENE_WORDS.get(r["uniprot"], [])): off.append((r["gene"], r["site"], p, t[:70]))
    note("PMIDS", f"title looks unrelated to protein or phosphorylation (review by hand): {len(off)}")
    for x in off: print("     ", x)
except Exception as e:
    note("PMIDS", f"PubMed check skipped (offline?): {e}")

# ---------- 6. VERIFICATION REPORT ----------
vr = Path("verification_report.csv")
if vr.exists():
    v = list(csv.DictReader(open(vr, encoding="utf-8-sig")))
    verd = collections.Counter(x["verdict"].split(":")[0] for x in v)
    note("VERIFY", f"verification_report.csv entries={len(v)} verdicts={dict(verd)}")
    log("VERIFY", len(v) == len(uniq) and set(verd) <= {"VERIFIED_EXPERIMENTAL_STRUCTURE"}, "report matches database and every site is structure-covered")
else:
    note("VERIFY", "verification_report.csv not found - run batch_verify.py instance/ptm.db")

# ---------- 7. APP PAGES + COUNTS ----------
try:
    sys.path.insert(0, ".")
    import app as A
    cl = A.app.test_client()
    pages = ["/", "/browse", "/browse?view=sites", "/statistics", "/download", "/connection-check", "/help",
             "/api/statistics", "/api/health", "/api/search?q=S473", "/search?q=TP53"]
    for a in sorted(prots): pages += [f"/protein/{a}", f"/viewer/{a}", f"/api/protein/{a}"]
    try:
        for rule in A.app.url_map.iter_rules():
            if rule.rule.startswith("/features") or rule.rule.startswith("/api/f/"):
                if not rule.arguments: pages.append(rule.rule)
    except Exception: pass
    errs = []
    for p in dict.fromkeys(pages):
        try:
            resp = cl.get(p)
            if resp.status_code >= 400: errs.append((p, resp.status_code))
        except Exception as e: errs.append((p, f"EXC {type(e).__name__}: {str(e)[:60]}"))
    log("PAGES", not errs, f"{len(set(pages))} routes requested, {len(errs)} error(s) (pages needing internet may fail offline)")
    for x in errs: print("     ", x)
    st = cl.get("/api/statistics").get_json()
    log("PAGES", st.get("sites") == len(uniq) and st.get("proteins") == len(prots) and st.get("kinases") == len(kin) and st.get("publications") == len(pm_all),
        f"/api/statistics {({k: st.get(k) for k in ('sites','proteins','kinases','publications')})} vs database {len(uniq)}/{len(prots)}/{len(kin)}/{len(pm_all)}")
    home = cl.get("/").get_data(as_text=True)
    log("PAGES", all(str(n) in home for n in (len(uniq), len(prots), len(kin), len(pm_all))), "home page shows the real counts")
    leak = [k for k in ("pmid_old","pmids_dbptm","pmid_score","pmid_count") if k in cl.get("/api/search?q=S473").get_data(as_text=True)]
    log("PAGES", not leak, f"internal columns leaking through /api/search: {leak}")
except Exception as e:
    note("PAGES", f"app import/test failed: {type(e).__name__}: {e}")

# ---------- 8. CODE: templates/static references ----------
tdir = Path("templates")
refs = set()
for f in list(Path(".").glob("*.py")) + list(Path("services").glob("*.py")):
    refs |= set(re.findall(r'render_template\(\s*["\']([^"\']+)', f.read_text(encoding="utf-8", errors="ignore")))
miss_t = sorted(t for t in refs if not (tdir / t).exists() and not Path("research/templates", t).exists())
log("CODE", not miss_t, f"templates referenced but missing: {miss_t}")
static_refs = set()
for f in tdir.rglob("*.html"):
    static_refs |= set(re.findall(r"static/([\w./-]+\.(?:js|css|png|svg))", f.read_text(encoding="utf-8", errors="ignore")))
miss_s = sorted(s for s in static_refs if not (Path("static") / s).exists())
log("CODE", not miss_s, f"static files referenced but missing: {miss_s}")
app_src = "".join(p.read_text(encoding="utf-8", errors="ignore") for p in [Path("app.py"), Path("features.py")] if p.exists())
routes = set(re.findall(r'route\(\s*["\']([^"\'<]+)', app_src))
js_calls = set()
for f in Path("static/js").glob("*.js"):
    js_calls |= set(re.findall(r'fetch\(\s*[`"\'](/(?:api|features)[^`"\'?${]*)', f.read_text(encoding="utf-8", errors="ignore")))
miss_api = sorted(u for u in js_calls if not any(r.startswith(u.rstrip("/")) or u.startswith(r.rstrip("/")) for r in routes))
note("CODE", f"JS endpoints with no obvious matching route (review by hand): {miss_api}")
mojibake = []
for f in list(Path(".").glob("*.py")) + list(Path("static").rglob("*.js")) + list(tdir.rglob("*.html")):
    if "\u00e2\u20ac" in f.read_text(encoding="utf-8", errors="ignore"): mojibake.append(str(f))
log("CODE", not mojibake, f"files containing mis-decoded characters (e.g. 'â€”'): {mojibake}")
css = Path("static/css/main.css").read_text(encoding="utf-8", errors="ignore") if Path("static/css/main.css").exists() else ""
log("CODE", not re.search(r"(?<![.\w>-])table\s*\{[^}]*min-width\s*:\s*1000px", css), "no global table{min-width:1000px} rule")

# ---------- SUMMARY ----------
fails = [r for r in results if r[1] == "FAIL"]
print("\n" + "="*60)
print(f"SUMMARY: {sum(r[1]=='PASS' for r in results)} passed, {len(fails)} failed, {sum(r[1]=='NOTE' for r in results)} notes")
for s, _, m in fails: print("  FAIL:", s, "-", m)
with open("full_check_report.csv", "w", newline="", encoding="utf-8-sig") as fh:
    w = csv.writer(fh); w.writerow(["section","status","message"]); w.writerows(results)
print("Report saved to full_check_report.csv")
