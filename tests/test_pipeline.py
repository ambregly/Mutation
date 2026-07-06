"""Tests de bout en bout sur les donnees d'exemple (bibliotheque standard seule)."""

import os
import sys
import tempfile
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from mutmatch import ods                       # noqa: E402
from mutmatch.vcf import read_vcf              # noqa: E402
from mutmatch.references import read_reference, _parse_coord  # noqa: E402
from mutmatch.match import FilterConfig, merge_all, match_known  # noqa: E402

EX = os.path.join(ROOT, "examples")


class TestOdsRoundTrip(unittest.TestCase):
    def test_write_then_read(self):
        rows = [["a", "b"], ["1", "x:2"], ["3", "4"]]
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "t.ods")
            ods.write_table(path, rows)
            back = ods.read_table(path)
        self.assertEqual(back, rows)


class TestParseCoord(unittest.TestCase):
    def test_chr_pos_ref_alt(self):
        # format reel de la colonne Position de Results_patientJB
        self.assertEqual(
            _parse_coord("19-33301387-C-CGGAAGATGCCCCG"),
            ("19", 33301387, "C", "CGGAAGATGCCCCG"),
        )

    def test_chr_pos(self):
        self.assertEqual(_parse_coord("13:28034314"), ("13", 28034314, "", ""))
        self.assertEqual(_parse_coord("chr13:28034314"), ("13", 28034314, "", ""))

    def test_position_seule(self):
        self.assertEqual(_parse_coord("28035000"), ("", 28035000, "", ""))

    def test_non_coordonnee(self):
        self.assertEqual(_parse_coord("CEBPA:c.1015_1027dup"), ("", None, "", ""))


class TestVcf(unittest.TestCase):
    def test_sample_and_info(self):
        vf = read_vcf(os.path.join(EX, "vcf", "JB_01_1_fast.gz.results.vcf"))
        self.assertEqual(vf.sample, "JB_01")
        self.assertEqual(vf.pair, "1")
        self.assertEqual(len(vf.variants), 3)
        v = vf.variants[0]
        self.assertEqual(v.pos, 28034317)
        self.assertEqual(v.m, 17)
        self.assertEqual(v.wt, 55123)
        self.assertTrue(v.is_dup)
        self.assertAlmostEqual(v.vaf, 0.000308)


class TestPipeline(unittest.TestCase):
    def setUp(self):
        self.vcfs = [
            read_vcf(os.path.join(EX, "vcf", "JB_01_1_fast.gz.results.vcf")),
            read_vcf(os.path.join(EX, "vcf", "JB_01_2_fast.gz.results.vcf")),
        ]
        self.known = (
            read_reference(os.path.join(EX, "Fichier_CHU.ods"))
            + read_reference(os.path.join(EX, "Results_patientJB.ods"))
        )

    def test_merge_pairs(self):
        per = merge_all(self.vcfs, FilterConfig())
        self.assertIn("JB_01", per)
        variants = {mv.pos: mv for mv in per["JB_01"]}
        # l'ITD est present dans les deux paires
        itd = variants[28034317]
        self.assertEqual(itd.found_in_pairs, "1,2")
        self.assertEqual(itd.total_m, 32)
        self.assertTrue(itd.is_dup)

    def test_known_matches(self):
        per = merge_all(self.vcfs, FilterConfig())
        matches = match_known(per, self.known)
        by_query = {m.known.query: m for m in matches}
        # correspondance exacte via les colonnes ref/alt du Fichier_CHU
        self.assertEqual(by_query["FLT3-ITD-76"].match_level, "exact")
        # correspondance exacte via la coordonnee chr-pos-ref-alt de Results
        self.assertEqual(by_query["FLT3-ITD-76;p.Xaa"].match_level, "exact")
        # correspondance par position seule
        self.assertEqual(by_query["FLT3-ins36"].match_level, "position")
        # mutation connue non detectee
        self.assertEqual(by_query["FLT3-D835"].match_level, "none")

    def test_filter_min_m(self):
        per = merge_all(self.vcfs, FilterConfig(min_m=5))
        # seul l'ITD (M=15/17) passe M>=5
        positions = {mv.pos for mv in per["JB_01"]}
        self.assertEqual(positions, {28034317})

    def test_require_both_pairs(self):
        per = merge_all(self.vcfs, FilterConfig(require_both_pairs=True))
        positions = {mv.pos for mv in per["JB_01"]}
        self.assertEqual(positions, {28034317})


if __name__ == "__main__":
    unittest.main()
