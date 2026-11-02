"""Live fetchers for REAL public annotations, cached on disk (instance/cache).

Nothing here fabricates data: if a source is unreachable or returns nothing usable, the
functions raise ExternalError (or return an explicit 'unavailable' structure) and the
caller shows that state to the user instead of guessing.

Sources: UniProt REST, Ensembl REST (Compara orthologues), EBI Proteins API (variants),
RCSB Search + files, PDBe SIFTS, AlphaFold DB, cBioPortal public API.
"""
import hashlib
import json
import os
import re
import time

import requests

from .dbutil import BASE

CACHE = BASE / "instance" / "cache"
UA = {"User-Agent": "PTM-Phosphoproteome-Cancer/2.0 (academic project)"}
TIMEOUT = 40

MODEL_SPECIES = ["mus_musculus", "rattus_norvegicus", "macaca_mulatta", "canis_lupus_familiaris",
                 "bos_taurus", "sus_scrofa", "gallus_gallus", "xenopus_tropicalis", "danio_rerio",
                 "drosophila_melanogaster", "caenorhabditis_elegans"]
SPECIES_LABEL = {"mus_musculus": "Mouse", "rattus_norvegicus": "Rat", "macaca_mulatta": "Rhesus macaque",
                 "canis_lupus_familiaris": "Dog", "bos_taurus": "Cow", "sus_scrofa": "Pig",
                 "gallus_gallus": "Chicken", "xenopus_tropicalis": "Xenopus", "danio_rerio": "Zebrafish",
                 "drosophila_melanogaster": "Fruit fly", "caenorhabditis_elegans": "C. elegans"}


class ExternalError(RuntimeError):
    pass


def _key(*parts):
    return hashlib.sha1(json.dumps(parts, sort_keys=True, default=str).encode()).hexdigest()


def fetch(url, params=None, headers=None, method="GET", body=None, ttl=7 * 86400, text=False, timeout=None):
    """GET/POST with a JSON-on-disk cache. Failures are never cached."""
    CACHE.mkdir(parents=True, exist_ok=True)
    f = CACHE / (_key(method, url, params, body, text) + ".json")
    if f.exists():
        try:
            c = json.loads(f.read_text(encoding="utf-8"))
            if time.time() - c["t"] < ttl:
                return c["d"]
        except Exception:
            pass
    h = dict(UA)
    h.update(headers or {})
    r = None
    for attempt in range(3):      # transient server errors (e.g. cBioPortal 503) are retried with backoff
        try:
            if method == "POST":
                r = requests.post(url, params=params, json=body, headers=h, timeout=timeout or TIMEOUT)
            else:
                r = requests.get(url, params=params, headers=h, timeout=timeout or TIMEOUT)
        except requests.RequestException as e:
            if attempt == 2 or isinstance(e, requests.Timeout) and attempt >= 1:
                raise ExternalError(f"{url.split('/')[2]} unreachable: {e.__class__.__name__}"
                                    + (" (server too slow; a successful result is cached, so retry or try a smaller study)"
                                       if isinstance(e, requests.Timeout) else ""))
        else:
            if r.status_code not in (429, 500, 502, 503, 504) or attempt == 2:
                break
        time.sleep(1.5 * (attempt + 1))
    if r.status_code == 404:
        raise ExternalError(f"{url.split('/')[2]}: not found (404)")
    if not r.ok:
        raise ExternalError(f"{url.split('/')[2]}: HTTP {r.status_code}")
    try:
        data = r.text if text else r.json()
    except ValueError:
        raise ExternalError(f"{url.split('/')[2]}: unreadable response")
    f.write_text(json.dumps({"t": time.time(), "d": data}), encoding="utf-8")
    return data


# ------------------------------------------------------------------ UniProt
def uniprot_sequence(acc):
    txt = fetch(f"https://rest.uniprot.org/uniprotkb/{acc}.fasta", text=True, ttl=30 * 86400)
    seq = "".join(l.strip() for l in txt.splitlines() if l and not l.startswith(">"))
    if not seq:
        raise ExternalError("UniProt returned an empty sequence")
    return seq


# ------------------------------------------------------------------ Ensembl
def ensembl_orthologs(symbol, species=None):
    """Compara orthologues of a human gene with Compara-aligned protein sequences.

    One request for 11 species with aligned sequences is heavy and Ensembl often answers HTTP 500,
    so species are requested in small batches. A failed batch is skipped (and reported only if
    NOTHING came back), so one slow species no longer blanks the whole page.
    """
    species = species or MODEL_SPECIES
    out, errors, seen = [], [], set()
    for i in range(0, len(species), 3):
        batch = species[i:i + 3]
        params = [("type", "orthologues"), ("sequence", "protein"), ("aligned", "1"),
                  ("content-type", "application/json")] + [("target_species", s) for s in batch]
        try:
            d = fetch(f"https://rest.ensembl.org/homology/symbol/homo_sapiens/{symbol}", params=params,
                      ttl=30 * 86400)
            homs = d["data"][0]["homologies"]
        except ExternalError as e:
            errors.append(str(e))
            continue
        except (KeyError, IndexError, TypeError):
            continue
        for h in homs:
            s, t = h.get("source", {}), h.get("target", {})
            sa, ta = s.get("align_seq") or s.get("seq"), t.get("align_seq") or t.get("seq")
            key = (t.get("species"), t.get("protein_id"))
            if not sa or not ta or key in seen:
                continue
            seen.add(key)
            out.append({"species": t.get("species"), "label": SPECIES_LABEL.get(t.get("species"), t.get("species")),
                        "orthology_type": h.get("type"), "protein_id": t.get("protein_id"),
                        "gene_id": t.get("id"), "perc_id": t.get("perc_id"), "source_protein_id": s.get("protein_id"),
                        "source_aligned": sa, "target_aligned": ta})
    if not out:
        raise ExternalError("Ensembl returned no orthologue data" + (f" ({errors[0]})" if errors else " for this gene"))
    return out


# -------------------------------------------------------------- EBI variants
def ebi_variants(acc):
    d = fetch(f"https://www.ebi.ac.uk/proteins/api/variation/{acc}", headers={"Accept": "application/json"},
              ttl=14 * 86400)
    from .analysis import CANCER_WORDS
    out = []
    for f in d.get("features", []):
        b, e = str(f.get("begin", "")), str(f.get("end", ""))
        if not b.isdigit() or b != e:
            continue
        wt, alt = f.get("wildType"), f.get("alternativeSequence")
        if not wt or not alt or len(wt) != 1:
            continue
        xrefs = f.get("xrefs") or []
        cosmic = [x.get("id") for x in xrefs if "cosmic" in (x.get("name") or "").lower()]
        dbsnp = [x.get("id") for x in xrefs if (x.get("name") or "").lower() == "dbsnp"]
        clinvar = [x.get("id") for x in xrefs if "clinvar" in (x.get("name") or "").lower()]
        diseases = sorted({a.get("name") for a in (f.get("association") or []) if a.get("disease") and a.get("name")})
        clin = sorted({c.get("type") for c in (f.get("clinicalSignificances") or []) if c.get("type")})
        cancer = bool(cosmic) or any(CANCER_WORDS.search(x or "") for x in diseases)
        out.append({"pos": int(b), "wt": wt, "alt": alt, "label": f"{wt}{b}{alt}",
                    "consequence": f.get("consequenceType"), "clinical": clin, "diseases": diseases,
                    "cancer": cancer, "cosmic": cosmic, "dbsnp": dbsnp, "clinvar": clinvar,
                    "source": f.get("sourceType")})
    return out


# ------------------------------------------------------------ structures
def rcsb_entries_for(acc, rows=25):
    q = {"query": {"type": "terminal", "service": "text", "parameters": {
        "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
        "operator": "exact_match", "value": acc}},
        "return_type": "entry", "request_options": {"paginate": {"start": 0, "rows": rows}}}
    d = fetch("https://search.rcsb.org/rcsbsearch/v2/query", method="POST", body=q, ttl=14 * 86400)
    return [x["identifier"] for x in d.get("result_set", []) if x.get("identifier")]


def sifts_mappings(pdb, acc):
    d = fetch(f"https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb.lower()}", ttl=60 * 86400)
    blocks = (((d.get(pdb.lower()) or {}).get("UniProt") or {}).get(acc) or {}).get("mappings") or []
    return blocks


def pdb_text(pdb):
    return fetch(f"https://files.rcsb.org/download/{pdb.upper()}.pdb", text=True, ttl=180 * 86400)


def alphafold_entry(acc):
    """The AlphaFold DB API record for exactly this accession (never just "the first row").

    The API can return several rows (isoforms / fragments of very long proteins). Prefer the row whose
    uniprotAccession equals the requested accession and that has a coordinate file; among those, the lowest
    fragment number. Anything else would silently give coordinates of a different protein or numbering.
    """
    base = (acc or "").split("-")[0]
    d = fetch(f"https://alphafold.ebi.ac.uk/api/prediction/{base}", ttl=30 * 86400)
    if not d or not isinstance(d, list):
        raise ExternalError("AlphaFold DB has no model for this accession")
    rows = [r for r in d if isinstance(r, dict) and r.get("pdbUrl")]
    if not rows:
        raise ExternalError("AlphaFold DB entry has no PDB file")
    exact = [r for r in rows if (r.get("uniprotAccession") or base) == base]
    if not exact:
        raise ExternalError(f"AlphaFold DB returned models only for other accessions than {base}")

    def frag(r):
        m = re.search(r"-F(\d+)$", str(r.get("entryId") or ""))
        return int(m.group(1)) if m else 1
    return sorted(exact, key=frag)[0]


def alphafold_model(acc):
    e = alphafold_entry(acc)
    text = fetch(e["pdbUrl"], text=True, ttl=180 * 86400)
    if not isinstance(text, str) or "ATOM" not in text:
        raise ExternalError("AlphaFold DB returned a file without coordinate records")
    return e.get("entryId") or acc, text


# --------------------------------------------------------------- cBioPortal
CBIO = os.environ.get("CBIO_API", "https://www.cbioportal.org/api")


def cbio_studies():
    d = fetch(f"{CBIO}/studies", params={"projection": "SUMMARY"}, headers={"Accept": "application/json"},
              ttl=86400)
    return sorted(({"id": s["studyId"], "name": s.get("name"), "cancer_type": s.get("cancerTypeId"),
                    "samples": s.get("allSampleCount")} for s in d), key=lambda s: s["name"] or "")


class MutationRetrievalError(ExternalError):
    """A cBioPortal retrieval problem with a machine-readable `status`.

    A study that simply has no mutation in the gene is NOT an error (it returns n_records_raw == 0).
    These are the cases where the answer is unknown, so downstream code must report "not retrieved"
    and never treat it as "zero mutations":
      study_not_found | gene_not_found | no_mutation_profile | no_sample_list |
      request_failed | malformed_response | incomplete
    """

    def __init__(self, message, status="request_failed", retrieval=None):
        super().__init__(message)
        self.status = status
        self.retrieval = retrieval or {}


_PROT_CHANGE = re.compile(r"^(?:p\.)?([A-Z\*])(\d+)([A-Z\*]|fs.*|del|_.*)?$")


def cbio_gene_mutations(study_id, symbol, page_size=50000, max_pages=20):
    """Real mutation records for one gene in one cBioPortal study, with a full retrieval audit trail.

    Verified step by step (each failure has its own status, see MutationRetrievalError):
      1. the study exists;
      2. the study has a MUTATION_EXTENDED molecular profile (its real id is discovered, not assumed to be
         "<study>_mutations");
      3. the study has a sequenced sample list (discovered; "<study>_sequenced" or category
         all_cases_with_mutation_data) and its sample count is read;
      4. the gene symbol resolves to an Entrez id;
      5. mutation pages are requested until a short page arrives. A full final page, a repeated page
         (server ignoring the page number) or a response that is not a list is reported as incomplete/malformed.
    A valid response with zero records returns an empty list and retrieval["status"] == "ok_no_mutations".
    """
    if not re.fullmatch(r"[A-Za-z0-9_\-]+", study_id or ""):
        raise MutationRetrievalError("Invalid study id", "study_not_found")
    if not re.fullmatch(r"[A-Za-z0-9_\-\.]+", symbol or ""):
        raise MutationRetrievalError("Invalid gene symbol", "gene_not_found")
    H = {"Accept": "application/json"}
    ret = {"study_id": study_id, "gene": symbol, "profile_id": None, "profile_candidates": [],
           "sample_list_id": None, "n_sequenced": None, "entrez_gene_id": None, "page_size": page_size,
           "pages": 0, "complete": False, "n_malformed_records": 0, "status": "pending", "warnings": []}

    def fail(msg, status):
        ret["status"] = status
        raise MutationRetrievalError(msg, status, ret)

    # 1. study
    try:
        fetch(f"{CBIO}/studies/{study_id}", headers=H, ttl=86400)
    except ExternalError as e:
        fail(f"cBioPortal study '{study_id}' does not exist" if "404" in str(e) else str(e),
             "study_not_found" if "404" in str(e) else "request_failed")

    # 2. mutation profile
    try:
        profiles = fetch(f"{CBIO}/studies/{study_id}/molecular-profiles", headers=H, ttl=86400)
    except ExternalError as e:
        fail(f"Could not list molecular profiles for '{study_id}': {e}", "request_failed")
    if not isinstance(profiles, list):
        fail("cBioPortal returned a malformed molecular-profile list", "malformed_response")
    cands = sorted(p.get("molecularProfileId") for p in profiles
                   if isinstance(p, dict) and p.get("molecularAlterationType") == "MUTATION_EXTENDED"
                   and p.get("molecularProfileId"))
    ret["profile_candidates"] = cands
    if not cands:
        fail(f"Study '{study_id}' has no mutation (MUTATION_EXTENDED) profile, so mutations cannot be retrieved "
             "from it. This is different from 'zero mutations'.", "no_mutation_profile")
    profile = f"{study_id}_mutations" if f"{study_id}_mutations" in cands else cands[0]
    if len(cands) > 1:
        ret["warnings"].append(f"{len(cands)} mutation profiles found ({', '.join(cands)}); used {profile}.")
    ret["profile_id"] = profile

    # 3. sequenced sample list
    try:
        lists = fetch(f"{CBIO}/studies/{study_id}/sample-lists", headers=H, ttl=86400)
    except ExternalError as e:
        fail(f"Could not list sample lists for '{study_id}': {e}", "request_failed")
    if not isinstance(lists, list):
        fail("cBioPortal returned a malformed sample-list response", "malformed_response")
    chosen = next((l for l in lists if isinstance(l, dict) and l.get("sampleListId") == f"{study_id}_sequenced"), None) \
        or next((l for l in lists if isinstance(l, dict) and l.get("category") == "all_cases_with_mutation_data"), None)
    if not chosen:
        fail(f"Study '{study_id}' has no sequenced sample list (cannot define the sequenced cohort).", "no_sample_list")
    ret["sample_list_id"] = chosen.get("sampleListId")
    n_seq = chosen.get("sampleCount")
    if n_seq is None and isinstance(chosen.get("sampleIds"), list):
        n_seq = len(chosen["sampleIds"])
    ret["n_sequenced"] = n_seq

    # 4. gene
    try:
        g = fetch(f"{CBIO}/genes/{symbol}", headers=H, ttl=90 * 86400)
    except ExternalError as e:
        fail(f"cBioPortal does not know gene symbol '{symbol}'" if "404" in str(e) else str(e),
             "gene_not_found" if "404" in str(e) else "request_failed")
    entrez = (g or {}).get("entrezGeneId") if isinstance(g, dict) else None
    if not entrez:
        fail("cBioPortal has no Entrez id for this gene", "gene_not_found")
    ret["entrez_gene_id"] = int(entrez)

    # 5. mutations, paginated
    records, first_keys, seen_pages = [], None, set()
    for page in range(max_pages):
        try:
            batch = fetch(f"{CBIO}/molecular-profiles/{profile}/mutations",
                          params={"sampleListId": ret["sample_list_id"], "entrezGeneId": entrez,
                                  "projection": "SUMMARY", "pageSize": page_size, "pageNumber": page},
                          headers=H, ttl=14 * 86400, timeout=150)
        except ExternalError as e:
            fail(f"Mutation request failed on page {page}: {e}", "request_failed")
        if not isinstance(batch, list):
            fail("cBioPortal returned a mutation response that is not a list (API schema changed?)", "malformed_response")
        ret["pages"] = page + 1
        sig = tuple((m.get("sampleId"), m.get("startPosition"), m.get("proteinChange"), m.get("mutationType"))
                    for m in batch[:20] if isinstance(m, dict))
        if page > 0 and batch and sig in seen_pages:
            ret["warnings"].append("Server returned the same page twice (pageNumber ignored); data may be truncated.")
            ret["status"] = "incomplete"
            raise MutationRetrievalError("cBioPortal pagination did not advance; the mutation list may be incomplete",
                                         "incomplete", ret)
        seen_pages.add(sig)
        records.extend(batch)
        if len(batch) < page_size:
            ret["complete"] = True
            break
    if not ret["complete"]:
        fail(f"Stopped after {max_pages} pages of {page_size}; the mutation list is incomplete", "incomplete")

    from collections import Counter
    out, skipped, by_type = [], Counter(), Counter()
    for m in records:
        if not isinstance(m, dict) or not m.get("sampleId"):
            ret["n_malformed_records"] += 1
            continue
        by_type[m.get("mutationType") or "unknown"] += 1
        pc = (m.get("proteinChange") or "").strip()
        mm = _PROT_CHANGE.match(pc)
        start = m.get("proteinPosStart")
        end = m.get("proteinPosEnd")
        try:
            start = int(start) if start not in (None, "", -1) else None
            end = int(end) if end not in (None, "", -1) else start
        except (TypeError, ValueError):
            start = end = None
        if start is None and mm:
            start = end = int(mm.group(2))
        if not start or start < 1:
            skipped["no_protein_position"] += 1
            continue
        wt = mm.group(1) if mm else None
        alt = mm.group(3) if mm and mm.group(3) and len(mm.group(3)) == 1 else None
        out.append({"pos": int(start), "pos_end": int(end or start), "wt": wt, "alt": alt, "protein_change": pc,
                    "type": m.get("mutationType"), "sample": m.get("sampleId"), "patient": m.get("patientId")})
    ret["status"] = "ok" if records else "ok_no_mutations"
    check = {"n_records_raw": len(records), "n_skipped_no_position": sum(skipped.values()),
             "skipped_reasons": dict(skipped),
             "n_distinct_samples": len({m.get("sampleId") for m in records if isinstance(m, dict)}),
             "by_type": dict(by_type),
             "profile_id": profile, "sample_list_id": ret["sample_list_id"],
             "api_url": f"{CBIO}/molecular-profiles/{profile}/mutations?sampleListId={ret['sample_list_id']}"
                        f"&entrezGeneId={int(entrez)}&projection=SUMMARY&pageSize={page_size}&pageNumber=0"}
    return {"n_sequenced": n_seq, "mutations": out, "check": check, "retrieval": ret}


def cbio_suggest_studies(symbol, limit=10, min_samples=300):
    """Large public cBioPortal studies in which `symbol` has real missense mutations.

    A study with no mutation in the gene (small cohorts often have none) makes every mutation-based
    analysis untestable. This only reports counts returned by cBioPortal; nothing is estimated.
    """
    from concurrent.futures import ThreadPoolExecutor
    keys = ("pan_can_atlas", "msk_", "impact", "genie", "mixed_pan", "metastatic")
    big = [s for s in cbio_studies()
           if (s.get("samples") or 0) >= min_samples and any(k in s["id"] for k in keys)]
    big.sort(key=lambda s: -(s.get("samples") or 0))
    big = big[:limit]

    def one(st):
        try:
            cb = cbio_gene_mutations(st["id"], symbol)
        except ExternalError:
            return None
        n = sum(1 for m in cb["mutations"] if (m.get("type") or "").lower().startswith("missense"))
        return {"id": st["id"], "name": st["name"], "samples": st["samples"], "n_missense": n} if n else None

    with ThreadPoolExecutor(max_workers=4) as pool:
        found = [r for r in pool.map(one, big) if r]
    found.sort(key=lambda r: -r["n_missense"])
    return {"checked": len(big), "studies": found}

