import os, sys, unittest
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from services import evolution as E

SEQ = "MSDVAIVKEGWLHKRGEYIKTWRPRYFLLKNDGTFIGYKERPQDVDQREAPLNNFSVAQCQLMKTERPRPNTFIIRCLQWTTVIERTFHVETPEEREEWTTAIQTVADGLKKQEEEEMDFRSGSPSDNSGAEEMEVSLAKPKHRVTMNEFEYLKLLGKGTFGKVILVKEKATGRYYAMKILKKEVIVAKDEVAHTLTENRVLQNSRHPFLTALKYSFQTHDRLCFVMEYANGGELFFHLSRERVFSEDRARFYGAEIVSALDYLHSEKNVVYRDLKLENLMLDKDGHIKITDFGLCKEGIKDGATMKTFCGTPEYLAPEVLEDNDYGRAVDWWGLGVVMYEMMCGRLPFYNQDHEKLFELILMEEIRFPRTLGPEAKSLLSGLLKKDPKQRLGGGSEDAKEIMQHRFFAGIVWQHVYEKKLSPPFKPQVTSETDTRYFDEEFTAQMITITPPDQDDSMECVDSERRPHFPQFSYSASGTA"


class T(unittest.TestCase):
    def test_spaces_and_separators_are_ignored(self):
        w = SEQ[465:480]
        for fl in (w, " ".join(w), w[:7] + "-" + w[7] + "_" + w[8:]):
            m = E.flank_match(SEQ, fl, 473)
            self.assertTrue(m["found"] and m["exact"], fl)

    def test_wildcard_and_one_mismatch(self):
        w = SEQ[465:480]
        self.assertTrue(E.flank_match(SEQ, "X" + w[1:], 473)["exact"])
        bad = w[:3] + ("A" if w[3] != "A" else "G") + w[4:]
        m = E.flank_match(SEQ, bad, 473)
        self.assertTrue(m["found"] and not m["exact"])

    def test_unrelated_flank_not_found(self):
        self.assertFalse(E.flank_match(SEQ, "WWWWWWWWWWWWWWW", 473)["found"])
        self.assertIsNone(E.flank_match(SEQ, "", 473))


if __name__ == "__main__":
    unittest.main()
