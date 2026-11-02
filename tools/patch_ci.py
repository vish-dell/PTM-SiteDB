from pathlib import Path
p = Path("services/proximity.py")
s = p.read_text(encoding="utf-8")
old = "    fold_ci = [ci[0] / expected_rate, ci[1] / expected_rate] if ci and expected_rate > 0 else None"
new = ("    # A Wilson interval from fewer than 3 mutated residues is bounded by 1/expected and uninformative\n"
       "    fold_ci = [ci[0] / expected_rate, ci[1] / expected_rate] if ci and expected_rate > 0 and n >= 3 else None")
assert s.count(old) == 1, "line not found exactly once - nothing changed"
p.write_text(s.replace(old, new), encoding="utf-8")
print("patched")
