"""Tests de bout en bout sur les donnees d'exemple (bibliotheque standard seule)."""

import os
import sys
import tempfile
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from mutmatch import ods                       # noqa: E402
from mutmatch.vcf import read_vcf              # noqa: E402
from mutmatch.references import (               # noqa: E402
    read_reference, read_any, read_clinical_wide, _parse_coord, norm_hgvs,
    mutation_type,
)
from mutmatch.match import FilterConfig, merge_all, match_known  # noqa: E402
from mutmatch.report import synthesis_table     # noqa: E402

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
        self.known = read_reference(os.path.join(EX, "Fichier_CHU.ods"))

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
        # ITD FLT3 : correspondance exacte (chr+pos+ref+alt)
        self.assertEqual(by_query["FLT3:c.1747_1794dup"].match_level, "exact")
        # seule la position coincide -> 'position'
        self.assertEqual(by_query["FLT3:c.pos"].match_level, "position")
        # CEBPA (chr19) hors de portee du VCF FLT3 -> non detecte
        self.assertEqual(by_query["CEBPA:c.1015_1027dup"].match_level, "none")
        # SNV FLT3 absent du VCF -> non detecte
        self.assertEqual(by_query["FLT3:c.2503G>T"].match_level, "none")

    def test_filter_min_m(self):
        per = merge_all(self.vcfs, FilterConfig(min_m=5))
        # seul l'ITD (M=15/17) passe M>=5
        positions = {mv.pos for mv in per["JB_01"]}
        self.assertEqual(positions, {28034317})

    def test_require_both_pairs(self):
        per = merge_all(self.vcfs, FilterConfig(require_both_pairs=True))
        positions = {mv.pos for mv in per["JB_01"]}
        self.assertEqual(positions, {28034317})


class TestClinicalWide(unittest.TestCase):
    def test_read_any_detects_clinical(self):
        kind, items = read_any(os.path.join(EX, "Results_patientJB.ods"))
        self.assertEqual(kind, "clinical")
        genes = {m.gene for m in items}
        self.assertEqual(genes, {"FLT3", "CEBPA"})
        flt3 = [m for m in items if m.hgvs_c == "c.1747_1794dup"][0]
        self.assertEqual(flt3.sample_id, "JB_01")
        self.assertEqual(flt3.hgvs, "FLT3:c.1747_1794dup")
        self.assertEqual(flt3.hgvs_p, "p.G583_E598dup")

    def test_read_any_detects_coord(self):
        kind, items = read_any(os.path.join(EX, "Fichier_CHU.ods"))
        self.assertEqual(kind, "coord")

    def test_norm_hgvs_join(self):
        self.assertEqual(
            norm_hgvs("FLT3 : c.1747_1794dup; p.G583_E598dup"),
            norm_hgvs("FLT3:c.1747_1794dup"),
        )

    def test_mutation_type(self):
        self.assertEqual(mutation_type("c.1747_1794dup"), "duplication")
        self.assertEqual(mutation_type("c.2503G>T"), "SNV")
        self.assertEqual(mutation_type("c.863_864insCCTG"), "insertion")
        self.assertEqual(mutation_type("c.282_297del"), "deletion")
        self.assertEqual(mutation_type("c.284_294delinsGT"), "delins")

    def test_synthesis(self):
        per = merge_all(self.vcfs, FilterConfig())
        matches = match_known(per, self.known)
        clinical = read_clinical_wide(os.path.join(EX, "Results_patientJB.ods"))
        vcf_chroms = {mv.chrom for v in per.values() for mv in v}
        rows = synthesis_table(clinical, matches, vcf_chroms)
        header = rows[0]
        idx = {h: i for i, h in enumerate(header)}
        by_gene_c = {(r[idx["gene"]], r[idx["hgvs_c"]]): r for r in rows[1:]}

        # duplication FLT3 : cible de FiLT3r et detectee
        itd = by_gene_c[("FLT3", "c.1747_1794dup")]
        self.assertEqual(itd[idx["type_mutation"]], "duplication")
        self.assertEqual(itd[idx["cible_filt3r"]], "oui")
        self.assertEqual(itd[idx["detecte"]], "oui")
        self.assertEqual(itd[idx["caryotype"]], "46,XX[40]")
        self.assertIn("confirmee", itd[idx["commentaire"]])

        # SNV FLT3 : hors perimetre FiLT3r (pas une duplication)
        snv = by_gene_c[("FLT3", "c.2503G>T")]
        self.assertEqual(snv[idx["type_mutation"]], "SNV")
        self.assertEqual(snv[idx["cible_filt3r"]], "non")

        # CEBPA : autre gene, hors perimetre
        cebpa = by_gene_c[("CEBPA", "c.1015_1027dup")]
        self.assertEqual(cebpa[idx["cible_filt3r"]], "non")
        self.assertEqual(cebpa[idx["chrom"]], "19")
        self.assertIn("hors perimetre", cebpa[idx["commentaire"]])

    def setUp(self):
        self.vcfs = [
            read_vcf(os.path.join(EX, "vcf", "JB_01_1_fast.gz.results.vcf")),
            read_vcf(os.path.join(EX, "vcf", "JB_01_2_fast.gz.results.vcf")),
        ]
        self.known = read_reference(os.path.join(EX, "Fichier_CHU.ods"))


if __name__ == "__main__":
    unittest.main()
