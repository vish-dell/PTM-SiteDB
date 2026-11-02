import sqlite3, re, urllib.request, collections
EXPECT = {"P04637":"TP53","P31749":"AKT1","P35222":"CTNNB1","P00533":"EGFR","P28482":"MAPK1",
          "P15056":"BRAF","P38398":"BRCA1","O14757":"CHEK1","P42336":"PIK3CA"}
def seq(acc):
    url = f"https://rest.uniprot.org/uniprotkb/{acc}.fasta"
    t = urllib.request.urlopen(url, timeout=60).read().decode()
    return "".join(t.splitlines()[1:])
c = sqlite3.connect("instance/ptm.db")
rows = c.execute("select id,gene,uniprot,ptm_type,site,position,flanking,pmid from ptm_sites").fetchall()
seqs = {a: seq(a) for a in {r[2] for r in rows}}
print("sequence lengths:", {a: len(s) for a, s in seqs.items()})
stat = collections.defaultdict(lambda: [0, set(), 0])
problems = []
for id_, gene, up, typ, site, pos, flank, pmid in rows:
    p = []
    if EXPECT.get(up) != gene: p.append(f"accession/gene mismatch ({up} vs {gene})")
    m = re.fullmatch(r"([A-Z])(\d+)", site or "")
    if not m: p.append("bad site format")
    else:
        aa, n = m.group(1), int(m.group(2))
        if pos is not None and int(pos) != n: p.append(f"position {pos} != {n}")
        s = seqs.get(up, "")
        if n > len(s): p.append("beyond protein length")
        elif s[n-1] != aa: p.append(f"UniProt has {s[n-1]} at {n}")
        else:
            if flank:
                f = flank.replace("-", "").replace("_", "").upper()
                if f and f not in s: p.append("flank not found in sequence")
        if aa not in "STY": p.append("not S/T/Y")
    if typ != "Phosphorylation": p.append("PTM type " + str(typ))
    ids = re.findall(r"\d+", pmid or "")
    if not ids: p.append("no PMID")
    if len(ids) != len(set(ids)): p.append("duplicate PMID in cell")
    st = stat[gene]; st[0] += 1; st[1].add(site)
    if p: st[2] += 1; problems.append((gene, site, p))
print(f"{'gene':8}{'records':>8}{'unique sites':>14}{'with problems':>15}")
for g, (n, u, bad) in sorted(stat.items()): print(f"{g:8}{n:>8}{len(u):>14}{bad:>15}")
print("TOTAL records:", len(rows), "| unique sites:", len({(r[2], r[4]) for r in rows}))
print("PROBLEMS:", len(problems))
for g, s, p in problems: print(" ", g, s, "|", "; ".join(p))
