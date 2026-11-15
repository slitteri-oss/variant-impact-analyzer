import unittest
import sys
import os

# Agregamos 'src' al path para importar las funciones de analyzer.py
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../src")))

from analyzer import translate_dna, classify_variant


class TestVariantAnalyzer(unittest.TestCase):

    def test_translate_dna_standard(self):
        # ATG (M) + GCT (A) + TAA (*)
        self.assertEqual(translate_dna("ATGGCTTAA"), "MA*")

    def test_translate_dna_incomplete_codon(self):
        # Ignora nucleótidos sobrantes que no completan un triplete
        self.assertEqual(translate_dna("ATGGCTTA"), "MA")

    def test_classify_synonymous(self):
        ref_aa, sample_aa, v_type, _ = classify_variant("TTC", "TTT")
        self.assertEqual((ref_aa, sample_aa, v_type), ("F", "F", "SYNONYMOUS"))

    def test_classify_missense_non_conservative(self):
        ref_aa, sample_aa, v_type, impact = classify_variant("GGC", "CGG")
        self.assertEqual((ref_aa, sample_aa, v_type), ("G", "R", "MISSENSE"))
        self.assertIn("Non-conservative", impact)

    def test_classify_nonsense(self):
        ref_aa, sample_aa, v_type, _ = classify_variant("TGG", "TAG")
        self.assertEqual((ref_aa, sample_aa, v_type), ("W", "*", "NONSENSE"))


if __name__ == "__main__":
    unittest.main()
