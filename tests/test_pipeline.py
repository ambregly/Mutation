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
from mutmatch.report import synthesis_table, novel_table  # noqa: E402
from mutmatch.match import MergedVariant, drop_recurrent  # noqa: E402

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


class TestPlot(unittest.TestCase):
    def test_read_exons_bed(self):
        import plot_itd
        exons = plot_itd.read_exons(os.path.join(EX, "flt3_exons.example.bed"))
        self.assertEqual(len(exons), 2)
        self.assertEqual(exons[0]["chrom"], "13")
        self.assertEqual(exons[0]["start"], 28034100)
        self.assertEqual(exons[0]["color"], "#F6C177")

    def test_read_exons_fasta(self):
        import plot_itd
        with tempfile.TemporaryDirectory() as d:
            fa = os.path.join(d, "ref.fa")
            with open(fa, "w") as fh:
                fh.write(">13:28033760-28034429:-1 FLT3\nACGT\n")
            exons = plot_itd.read_exons(fa)
        self.assertEqual(exons[0]["start"], 28033760)
        self.assertEqual(exons[0]["end"], 28034429)

    def test_radius_modes(self):
        import plot_itd
        # M : echelle log -> 1000 reads donne un plus gros rayon que 10
        self.assertGreater(plot_itd._radius(1000, "m"), plot_itd._radius(10, "m"))
        # VAF : plus la VAF est haute, plus le point est gros
        self.assertGreater(plot_itd._radius(0.02, "vaf"),
                           plot_itd._radius(0.001, "vaf"))

    def test_build_svg_runs(self):
        import plot_itd
        variants = plot_itd._read(_mini_tsv(), dup_only=False)
        svg = plot_itd.build_svg(variants, size_by="m")
        self.assertIn("<svg", svg)

    def test_read_missed(self):
        import plot_itd
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "synthese.tsv")
            with open(path, "w") as fh:
                fh.write("sample\tpos\thgvs_c\tcible_filt3r\tdetecte\n")
                fh.write("JB_15\t28034135\tc.1745_1783dup\toui\tnon\n")   # a marquer
                fh.write("JB_05\t28034124\tc.1747_1794dup\toui\toui\n")   # detectee: ignore
                fh.write("JB_01\t\tc.x\tnon\t?\n")                        # sans pos: ignore
            missed = plot_itd.read_missed(path)
        self.assertEqual(len(missed), 1)
        self.assertEqual(missed[0]["sample"], "JB_15")
        self.assertEqual(missed[0]["pos"], 28034135)


def _mini_tsv():
    """Ecrit un mini variants_par_echantillon.tsv et renvoie son chemin."""
    d = tempfile.mkdtemp()
    path = os.path.join(d, "v.tsv")
    with open(path, "w") as fh:
        fh.write("sample\tchrom\tpos\tref\talt\tsvlen\tDUP\ttrouve_dans_paires\t"
                 "nb_paires\tM_R1\tVAF_R1\tWT_R1\tM_R2\tVAF_R2\tWT_R2\tM_total\t"
                 "VAF_max\tconnu\tniveau_correspondance\tquery_connue\tpatient_id\n")
        fh.write("JB_01\t13\t28034124\tA\tATTC\t48\toui\t1,2\t2\t320\t0.013\t100\t"
                 "427\t0.017\t100\t747\t0.017\toui\texact\tFLT3:x\tJB_01\n")
    return path


class TestTriage(unittest.TestCase):
    def _mv(self, sample, pos, alt, svlen, dup, pairs, vaf, m):
        mv = MergedVariant(sample=sample, chrom="13", pos=pos, ref="C", alt=alt,
                           svlen=svlen, is_dup=dup)
        for p in pairs:
            mv.per_pair[p] = {"m": m, "wt": 10000, "vaf": vaf, "filter": "PASS"}
        return mv

    def test_triage_labels(self):
        # variant recurrent (present dans 6 echantillons) a basse VAF -> artefact
        per = {}
        for i in range(6):
            s = "JB_%02d" % i
            per[s] = [self._mv(s, 100, "CGG", 3, True, ["1", "2"], 0.0002, 2)]
        # un candidat ITD net chez un seul patient
        per["JB_99"] = [self._mv("JB_99", 28034100, "C" + "GAT" * 6, 18,
                                 True, ["1", "2"], 0.02, 800)]
        rows = novel_table(per, candidate_min_vaf=0.01, artifact_min_samples=5)
        idx = {h: i for i, h in enumerate(rows[0])}
        by_pos = {r[idx["pos"]]: r for r in rows[1:]}
        self.assertTrue(by_pos["100"][idx["triage"]].startswith("artefact"))
        self.assertEqual(by_pos["28034100"][idx["triage"]], "candidat ITD")
        # tri : le candidat (VAF 0.02) doit etre en premier
        self.assertEqual(rows[1][idx["pos"]], "28034100")

    def test_drop_recurrent(self):
        # meme variant (pos 100) chez 6 echantillons -> artefact ; pos 500 unique
        per = {}
        for i in range(6):
            s = "JB_%02d" % i
            per[s] = [self._mv(s, 100, "CGG", 3, True, ["1", "2"], 0.0002, 2)]
        per["JB_00"].append(
            self._mv("JB_00", 500, "CAAA", 3, True, ["1", "2"], 0.02, 800))
        out, removed = drop_recurrent(per, min_samples=5)
        self.assertEqual(removed, 6)  # le variant recurrent chez 6 patients
        # le variant unique (pos 500) est conserve
        positions = {mv.pos for vs in out.values() for mv in vs}
        self.assertEqual(positions, {500})

    def test_drop_recurrent_disabled(self):
        per = {"JB_01": [self._mv("JB_01", 100, "CGG", 3, True, ["1"], 0.5, 9)]}
        out, removed = drop_recurrent(per, min_samples=0)
        self.assertEqual(removed, 0)

    def test_low_vaf_is_noise(self):
        per = {"JB_01": [self._mv("JB_01", 200, "CA", 1, False, ["1"], 0.0001, 1)]}
        rows = novel_table(per, candidate_min_vaf=0.01, artifact_min_samples=5)
        idx = {h: i for i, h in enumerate(rows[0])}
        self.assertIn("bruit", rows[1][idx["triage"]])
        self.assertEqual(rows[1][idx["cadre_lecture"]], "hors-cadre")


if __name__ == "__main__":
    unittest.main()
