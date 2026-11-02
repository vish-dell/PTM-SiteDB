"""Tests for benchmark_report and the batch filename convention (tiny hand-built fixtures, never app data)."""
import os, re, sys, unittest
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import benchmark_report as B


class Benchmark(unittest.TestCase):
    def setUp(self):
        self.project = [
            {"gene": "A", "uniprot": "P1", "site": "S1", "q_perm_all_sites": "0.01"},
            {"gene": "B", "uniprot": "P2", "site": "S2", "q_perm_all_sites": "0.40"},
            {"gene": "C", "uniprot": "P3", "site": "S3", "q_perm_all_sites": "0.02"},
            {"gene": "D", "uniprot": "P4", "site": "S4", "q_perm_all_sites": "0.60"},
            {"gene": "E", "uniprot": "P5", "site": "S5", "q_perm_all_sites": ""},
        ]
        self.bench = [
            {"gene": "A", "reference_result": "positive"},   # TP
            {"gene": "B", "reference_result": "positive"},   # FN
            {"gene": "C", "reference_result": "negative"},   # FP
            {"gene": "D", "reference_result": "negative"},   # TN
            {"gene": "E", "reference_result": "positive"},   # untestable: no q
            {"gene": "Z", "reference_result": "positive"},   # absent from project
            {"gene": "A", "reference_result": "maybe"},      # invalid label
        ]

    def test_counts_and_metrics(self):
        rows = B.compare(self.project, self.bench)
        s = B.summarize(rows)
        self.assertEqual((s["TP"], s["FN"], s["FP"], s["TN"]), (1, 1, 1, 1))
        self.assertEqual(s["n_not_comparable"], 2)       # untested rows are never counted as negatives
        self.assertEqual(s["n_invalid_reference"], 1)
        self.assertEqual((s["precision"], s["recall"], s["specificity"]), (0.5, 0.5, 0.5))

    def test_no_division_by_zero(self):
        s = B.summarize(B.compare(self.project, [{"gene": "Z", "reference_result": "positive"}]))
        self.assertIsNone(s["precision"]); self.assertIsNone(s["recall"])


class Filenames(unittest.TestCase):
    def test_batch_and_sensitivity_agree_on_names(self):
        src = open(os.path.join(os.path.dirname(__file__), "..", "batch_proximity.py"), encoding="utf-8").read()
        self.assertIn("sep{sep}_dist{dist:g}.csv", src)
        sens = open(os.path.join(os.path.dirname(__file__), "..", "sensitivity_report.py"), encoding="utf-8").read()
        self.assertTrue(re.search(r"proximity_\{args\.study\}_sep\{sep\}_dist\{dist:g\}\.csv", sens))


if __name__ == "__main__":
    unittest.main()
