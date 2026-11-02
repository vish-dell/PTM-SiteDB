import os, sys, unittest
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from services import structure as S
from services import external as X


def pdb_line(serial, chain, resi, aa3, x, b=90.0):
    return f"ATOM  {serial:5d}  CA  {aa3:>3s} {chain}{resi:4d}    {x:8.3f}{0:8.3f}{0:8.3f}  1.00{b:6.2f}           C"


class StructureFallbackTests(unittest.TestCase):
    def test_alphafold_selects_chain_with_best_real_ptm_coverage(self):
        text = "\n".join([
            pdb_line(1, "A", 1, "SER", 0),
            pdb_line(2, "B", 1, "SER", 0),
            pdb_line(3, "B", 2, "THR", 5),
        ])
        old = X.alphafold_model
        try:
            X.alphafold_model = lambda acc: ("AF-TEST-F1", text)
            units = [
                {"site": "S1", "pos": 1, "aa": "S", "ptm_type": "Phosphorylation"},
                {"site": "T2", "pos": 2, "aa": "T", "ptm_type": "Phosphorylation"},
            ]
            r = S._alphafold("PTEST", units)
            self.assertEqual(r["chain"], "B")
            self.assertEqual(len(r["mapped"]), 2)
        finally:
            X.alphafold_model = old

    def test_spatial_reports_experimental_and_alphafold_coverage(self):
        text = "\n".join([
            pdb_line(1, "A", 1, "SER", 0),
            pdb_line(2, "A", 2, "THR", 5),
            pdb_line(3, "A", 3, "LYS", 10),
        ])
        old_af = X.alphafold_model
        old_rcsb = X.rcsb_entries_for
        old_sifts = X.sifts_mappings
        old_pdb = X.pdb_text
        try:
            X.rcsb_entries_for = lambda acc: ["TEST"]
            X.sifts_mappings = lambda pdb, acc: [
                {"unp_start": 1, "unp_end": 1, "chain_id": "A", "start": {"author_residue_number": 1}}
            ]
            X.pdb_text = lambda pdb: text
            X.alphafold_model = lambda acc: ("AF-TEST-F1", text)
            units = [
                {"site": "S1", "pos": 1, "aa": "S", "ptm_type": "Phosphorylation"},
                {"site": "T2", "pos": 2, "aa": "T", "ptm_type": "Phosphorylation"},
                {"site": "K3", "pos": 3, "aa": "K", "ptm_type": "Acetylation"},
            ]
            r = S.spatial_analysis("PTEST", units, threshold=10, allow_alphafold=True)
            self.assertTrue(r["available"])
            self.assertEqual(r["source"], "alphafold")
            self.assertEqual(r["n_mapped"], 3)
            self.assertEqual(r["coverage"]["experimental"]["mapped"], 1)
            self.assertEqual(r["coverage"]["alphafold"]["mapped"], 3)
            self.assertTrue(any("AlphaFold fallback selected" in n for n in r["notes"]))
        finally:
            X.alphafold_model = old_af
            X.rcsb_entries_for = old_rcsb
            X.sifts_mappings = old_sifts
            X.pdb_text = old_pdb


if __name__ == "__main__":
    unittest.main()
