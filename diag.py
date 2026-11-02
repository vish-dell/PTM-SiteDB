import json, os, sqlite3
acc = "P27361"
p = f"verify_cache/search_{acc}.json"
ids = [r["identifier"] for r in json.load(open(p))["result_set"]] if os.path.exists(p) else []
print(acc, "RCSB entries:", len(ids), ids[:20])
for i in ids:
    f = f"verify_cache/sifts_{i}.json"
    if not os.path.exists(f):
        print(i, "no SIFTS file cached"); continue
    try:
        m = json.load(open(f))[i.lower()]["UniProt"][acc]["mappings"]
        print(i, [(b["unp_start"], b["unp_end"], b["chain_id"]) for b in m])
    except Exception as e:
        print(i, "SIFTS has no mapping for", acc, e)
con = sqlite3.connect("instance/ptm.db")
print()
for r in con.execute("select gene, site, pmid, evidence, confidence from ptm_sites where (uniprot='P31749' and site='S124') or (uniprot='O14757' and site='S365') or (uniprot='P27361' and site in ('T202','Y204')) or (uniprot='P04637' and site='S315')"):
    print(r)