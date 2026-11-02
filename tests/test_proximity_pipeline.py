"""End-to-end check of proximity_analysis using a tiny hand-built structure (test fixture only, never app data)."""
import os, sys, unittest
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from services import proximity as P


def pdb_line(i, x, b):
    return f"ATOM  {i:5d}  CA  ALA A{i:4d}    {x:8.3f}{0:8.3f}{0:8.3f}  1.00{b:6.2f}           C"


class Pipeline(unittest.TestCase):
    def test_hit_found_and_reported(self):
        # residue 1 = PTM; residue 60 sits 6 A away (3D-close, 59 apart); residue 61 is 80 A away
        xs = {1: 0.0, 2: 40.0, 60: 6.0, 61: 80.0}
        text = "\n".join(pdb_line(i, x, 90.0) for i, x in xs.items())
        P.ext.alphafold_model = lambda acc: ("AF-TEST-F1", text)
        cb = {"mutations": [{"pos": 60, "protein_change": "A60T", "type": "Missense_Mutation"},
                            {"pos": 60, "protein_change": "A60V", "type": "Missense_Mutation"},
                            {"pos": 61, "protein_change": "A61T", "type": "Missense_Mutation"},
                            {"pos": 2, "protein_change": "A2*", "type": "Nonsense_Mutation"},
                            {"pos": 60, "protein_change": "A60A", "type": "Synonymous_Mutation"}]}
        res = P.proximity_analysis("TEST", [{"pos": 1, "site": "S1", "ptm_type": "Phosphorylation"}], cb, threshold=10, min_sep=20)
        self.assertTrue(res["available"])
        self.assertEqual([h["pos"] for h in res["hits"]], [60])
        self.assertEqual(res["hits"][0]["n_mutations"], 2)
        self.assertEqual(res["n_missense"], 3)          # nonsense excluded
        self.assertEqual(res["n_observed"], 1)
        self.assertEqual(res["n_qualifying_residues"], 1)
        self.assertEqual(res["n_synonymous"], 1)
        self.assertEqual(res["negative_control"]["n_observed"], 1)
        # Eligible residues (confident, not the PTM, >= min_sep from it) = {60, 61}: N=2, K=1 (60), and both are
        # mutated (n=2), so every eligible residue is drawn and P(X>=1) = 1 exactly: the data cannot discriminate.
        # (The former 0.5 used N=4, wrongly counting the PTM residue and sequence-near residue 2.)
        self.assertAlmostEqual(res["p_value"], 1.0)


if __name__ == "__main__":
    unittest.main()
