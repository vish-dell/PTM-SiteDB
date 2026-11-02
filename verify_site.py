import sys, json, re, urllib.request

acc, site = sys.argv[1], sys.argv[2]          # e.g. P04637 S15
aa, pos = site[0].upper(), int(site[1:])

def get(url, data=None):
    req = urllib.request.Request(url, data=data, headers={"Content-Type": "application/json"})
    return urllib.request.urlopen(req, timeout=60).read().decode()

seq = "".join(get(f"https://rest.uniprot.org/uniprotkb/{acc}.fasta").splitlines()[1:])
print(f"UniProt residue {pos} = {seq[pos-1]}  (record says {aa})  flank: {seq[max(0,pos-8):pos+7]}")

term = lambda a, v: {"type": "terminal", "service": "text",
                     "parameters": {"attribute": a, "operator": "exact_match", "value": v}}
q = {"query": {"type": "group", "logical_operator": "and", "nodes": [
        term("rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_name", "UniProt"),
        term("rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession", acc)]},
     "return_type": "entry",
     "request_options": {"paginate": {"start": 0, "rows": 300}, "results_content_type": ["experimental"]}}
ids = [r["identifier"] for r in json.loads(get("https://search.rcsb.org/rcsbsearch/v2/query", json.dumps(q).encode())).get("result_set", [])]
print(len(ids), "experimental entries found")

for pdb in ids:
    try:
        m = json.loads(get(f"https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb.lower()}"))[pdb.lower()]["UniProt"][acc]["mappings"]
    except Exception:
        continue
    for b in m:
        if not (b["unp_start"] <= pos <= b["unp_end"]):
            continue
        chain = b["chain_id"]
        resi = b["start"]["author_residue_number"] + (pos - b["unp_start"])
        txt = get(f"https://files.rcsb.org/download/{pdb}.pdb")
        resn = next((l[17:20] for l in txt.splitlines()
                     if l[:6] in ("ATOM  ", "HETATM") and l[21] == chain
                     and l[22:26].strip() == str(resi) and l[26] == " "), None)
        missing = re.search(rf"REMARK 465\s+\w{{3}}\s+{chain}\s+{resi}\s*$", txt, re.M)
        print(f"{pdb} chain {chain} resi {resi}: {resn or ('MISSING (REMARK 465)' if missing else 'no coordinates')}")