"""Why don't PTM sites map onto the 3D structure?  READ-ONLY diagnostic (never edits the database).

For every protein in instance/ptm.db it fetches the UniProt sequence + signal/propeptide/initiator-Met features and tests
ONLY principled numbering candidates (offset 0 = UniProt numbering, +signal-peptide length, +1 for a removed initiator Met,
+propeptide/transit end), not arbitrary shifts. A record "supports" an offset when the residue letter at position+offset
equals the recorded letter. Letter agreement alone is weak evidence (~5-10% by chance), so each protein gets a binomial
p-value; treat a verdict as a hypothesis to confirm against the source paper, not as a correction.
Usage (internet needed):  python validation\\diagnose_site_mapping.py [--db instance\\ptm.db]
"""
import argparse, collections, json, sqlite3, sys, urllib.request
from pathlib import Path
from scipy.stats import binomtest

def candidate_offsets(features):
    off = {0: "UniProt numbering (as recorded)"}
    for f in features:
        t, end = f.get("type", ""), (f.get("location", {}).get("end", {}) or {}).get("value")
        if end is None: continue
        if t == "Signal": off[end] = f"signal peptide removed (+{end})"
        elif t == "Initiator methionine": off[end] = "initiator Met removed (+1)"
        elif t in ("Propeptide", "Transit peptide"): off[end] = f"{t.lower()} removed (+{end})"
    return off

def diagnose(records, seq, offsets):
    """records: [(site, position)]. Returns {offset: (n_match, n_total, p_value)}."""
    freq = collections.Counter(seq); n = len(seq); out = {}
    for off in offsets:
        hits, p0s = 0, []
        for site, pos in records:
            q = pos + off
            if 1 <= q <= n and seq[q - 1] == site[0]: hits += 1
            p0s.append(max(freq[site[0]] / n, 1e-6))
        p0 = sum(p0s) / len(p0s)
        out[off] = (hits, len(records), binomtest(hits, len(records), p0, alternative="greater").pvalue)
    return out

def per_record(records, seq, offsets):
    rows = []
    for site, pos in records:
        ok = [o for o in offsets if 1 <= pos + o <= len(seq) and seq[pos + o - 1] == site[0]]
        rows.append((site, "BEYOND PROTEIN END" if pos > len(seq) and 0 not in ok else ("ok as recorded" if 0 in ok else (f"letter fits offset {ok}" if ok else "fits no candidate offset"))))
    return rows

def fetch(acc):
    with urllib.request.urlopen(f"https://rest.uniprot.org/uniprotkb/{acc}.json", timeout=60) as r:
        j = json.load(r)
    return j["sequence"]["value"], j.get("features", []), j.get("entryType", "")

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--db", default=str(Path(__file__).resolve().parents[1] / "instance" / "ptm.db")); a = ap.parse_args()
    c = sqlite3.connect(a.db); by = collections.defaultdict(list)
    for acc, gene, site, pos in c.execute("select uniprot,gene,site,position from ptm_sites"): by[(acc, gene)].append((site.strip(), int(pos)))
    for (acc, gene), recs in sorted(by.items()):
        try: seq, feats, _ = fetch(acc)
        except Exception as e: print(f"{gene} {acc}: could not fetch UniProt ({e})"); continue
        offs = candidate_offsets(feats); res = diagnose(recs, seq, offs)
        print(f"\n{gene} {acc}  length={len(seq)}  records={len(recs)}")
        for o, (h, t, p) in sorted(res.items(), key=lambda x: -x[1][0]): print(f"   offset {o:+4d}  {offs[o]:45s} {h}/{t} letters agree  (binomial p={p:.2g})")
        for s, v in per_record(recs, seq, offs): 
            if v != "ok as recorded": print("     ", s, "->", v)