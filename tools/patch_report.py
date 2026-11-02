from pathlib import Path
p = Path("validation/FINAL_REPORT.md")
s = p.read_text(encoding="utf-8")
a = "has only **67 PTM records** (51 phospho, 11 acetyl, 5 ubiquitin)"
b = "now has **30 verified phosphorylation sites** (cleaned 25 Oct 2026; originally 67 records: 51 phospho, 11 acetyl, 5 ubiquitin; 37 rows kept in ptm_sites_quarantine)"
c = "67 original-site records (DB)"
d = "30 verified phosphosites (DB, after cleanup)"
assert s.count(a) == 1 and s.count(c) == 1, "text not found exactly once - nothing changed"
p.write_text(s.replace(a, b).replace(c, d), encoding="utf-8")
print("patched")
