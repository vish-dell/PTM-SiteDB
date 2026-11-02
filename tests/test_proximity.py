import math, os, sys, unittest
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from services import proximity as P


class T(unittest.TestCase):
    def test_bh_matches_hand_calc(self):
        q = P.bh_fdr([0.01, 0.04, 0.03, 0.005])
        self.assertTrue(all(abs(a - b) < 1e-12 for a, b in zip(q, [0.02, 0.04, 0.04, 0.02])))

    def test_hypergeometric_exact(self):
        # N=10, K=4, n=3: P(X>=2) = (C(4,2)C(6,1)+C(4,3))/C(10,3) = 40/120
        self.assertAlmostEqual(P.hyper_sf(2, 10, 4, 3), 40 / 120)
        self.assertEqual(P.hyper_sf(0, 10, 4, 3), 1.0)
        self.assertEqual(P.hyper_sf(5, 10, 4, 3), 0.0)

    def test_classify_3d_close_seq_far(self):
        ca = {1: (0, 0, 0), 2: (50, 0, 0), 60: (5, 0, 0), 61: (80, 0, 0)}
        c = P.classify(ca, {k: 90 for k in ca}, [1], threshold=10, min_sep=20)
        self.assertTrue(c[60]["qualifies"])       # 5 A away, 59 residues apart
        self.assertFalse(c[2]["qualifies"])       # far in 3D
        self.assertFalse(c[1]["qualifies"])       # the PTM itself (sep 0)

    def test_untestable_gives_none_not_p1(self):
        ca = {1: (0, 0, 0), 2: (50, 0, 0), 60: (80, 0, 0)}
        c = P.classify(ca, {k: 90 for k in ca}, [1], threshold=10, min_sep=20)
        r = P.proximity_test(c, [2, 60])
        self.assertIsNone(r["p_value"])
        self.assertIn("no eligible residue", r["note"])

    def test_low_confidence_excluded(self):
        ca = {1: (0, 0, 0), 60: (5, 0, 0)}
        c = P.classify(ca, {1: 90, 60: 50}, [1])
        self.assertNotIn(60, c)

    def test_composition_null_keeps_letters_and_is_seeded(self):
        seq = "S" * 30 + "A" * 70
        a = P.composition_null(seq, [5, 10, 15], 15, n_perm=200, seed=3)
        self.assertEqual(a, P.composition_null(seq, [5, 10, 15], 15, n_perm=200, seed=3))
        self.assertTrue(all(1 <= m <= 3 for m in a))   # 3 sites in total; never more than 3 in a window
        self.assertAlmostEqual(P.perm_p([3, 2, 1, 1], 3), 2 / 5)


if __name__ == "__main__":
    unittest.main()


class Extra(unittest.TestCase):
    def test_fisher_matches_closed_form(self):
        import math
        self.assertAlmostEqual(P.fisher_combine([0.05, 0.05]), 0.0025 * (1 - math.log(0.0025)), places=9)
        self.assertAlmostEqual(P.fisher_combine([0.3]), 0.3, places=9)

    def test_site_neighbourhood(self):
        ca = {1: (0, 0, 0), 2: (3, 0, 0), 60: (6, 0, 0), 70: (50, 0, 0), 80: (90, 0, 0)}
        pl = {k: 90 for k in ca}
        r = P.site_neighbourhoods(ca, pl, [{"site": "S1", "ptm_type": "Phos", "pos": 1}], [60, 70], 10, 20)[0]
        self.assertEqual(r["n_contacts"], 2)               # residues 2 and 60 within 10 A
        self.assertEqual(r["n_neighbourhood"], 1)          # residue 2 is sequence-near -> excluded
        self.assertEqual(r["hit_positions"], [60])
        # Null universe = residues that could ever qualify: not the PTM itself, >= min_sep away in sequence
        # -> {60, 70, 80}: N=3, K=1 (60), n=2 mutated-eligible (60, 70), k=1.
        # P(X>=1) = 1 - C(2,2)/C(3,2) = 2/3. (The former 0.4 used N=5, counting the PTM residue and the
        # sequence-near residue 2, which can never qualify; verified independently with scipy below.)
        self.assertAlmostEqual(r["p_value"], 2 / 3)
        from scipy.stats import hypergeom
        self.assertAlmostEqual(r["p_value"], hypergeom.sf(0, 3, 1, 2))


class HypergeometricBoundaries(unittest.TestCase):
    """hyper_sf must equal scipy's inclusive upper tail P(X >= k) = hypergeom.sf(k - 1, N, K, n)."""

    def test_matches_scipy_exhaustively_on_small_grid(self):
        from scipy.stats import hypergeom
        for N in range(1, 13):
            for K in range(0, N + 1):
                for n in range(0, N + 1):
                    for k in range(0, n + 2):
                        self.assertAlmostEqual(P.hyper_sf(k, N, K, n), hypergeom.sf(k - 1, N, K, n), places=10,
                                               msg=(k, N, K, n))

    def test_inclusive_tail_and_edges(self):
        self.assertEqual(P.hyper_sf(0, 10, 3, 4), 1.0)            # P(X>=0) is always 1
        self.assertEqual(P.hyper_sf(4, 10, 3, 4), 0.0)            # k above support (K=3)
        self.assertAlmostEqual(P.hyper_sf(1, 4, 1, 2), 0.5)       # 1 - C(3,2)/C(4,2)
        self.assertAlmostEqual(P.hyper_sf(1, 2, 1, 2), 1.0)       # all eligible residues drawn: uninformative
        self.assertAlmostEqual(P.hyper_sf(2, 6, 3, 3), 0.5, places=12)  # (C(3,2)C(3,1)+C(3,3))/C(6,3) = 10/20

    def test_invalid_parameters_raise(self):
        with self.assertRaises(ValueError):
            P.hyper_sf(1, 5, 6, 2)
