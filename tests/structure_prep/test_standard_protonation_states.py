"""
Regression test for standard ("textbook") protonation-state assignment.

ProtonationStateAnalyzer.get_default_recommendations(standard_states=True)
ignores the PROPKA per-residue analysis entirely and assigns canonical
physiological-pH states: acids (ASP/GLU) deprotonated; bases/thiol/phenol
(LYS/ARG/CYS/TYR) protonated; HIS as the epsilon tautomer (HIE); and — the
case PROPKA cannot reach on its own because it does not model non-standard
groups — heme propionate PRN as the deprotonated acid PRD.

The test deliberately sets each residue's analysis fields (protonated /
probability) to the OPPOSITE of the standard answer, so a passing test proves
the standard path overrides the analysis rather than coincidentally agreeing
with it.
"""
import unittest

from proprep.structure_prep.protonation_worker import ProtonationStateAnalyzer


def _res(restype, chain, number, protonated, probability):
    return dict(type=restype, chain=chain, number=number,
                protonated=protonated, probability=probability, pKa=7.0)


class StandardProtonationStatesTest(unittest.TestCase):
    def setUp(self):
        # Instance without running __init__ — get_default_recommendations is a
        # pure function of self.results.
        self.analyzer = ProtonationStateAnalyzer.__new__(ProtonationStateAnalyzer)
        # Each residue's analysis is set to the OPPOSITE of the standard state.
        self.analyzer.results = {
            "A_10_ASP": _res("ASP", "A", 10, protonated=True,  probability=0.95),
            "A_20_GLU": _res("GLU", "A", 20, protonated=True,  probability=0.95),
            "A_30_HIS": _res("HIS", "A", 30, protonated=True,  probability=0.95),
            "A_40_CYS": _res("CYS", "A", 40, protonated=False, probability=0.05),
            "A_50_LYS": _res("LYS", "A", 50, protonated=False, probability=0.05),
            "A_60_TYR": _res("TYR", "A", 60, protonated=False, probability=0.05),
            "A_70_ARG": _res("ARG", "A", 70, protonated=False, probability=0.05),
            "A_80_PRN": _res("PRN", "A", 80, protonated=True,  probability=0.95),
        }

    def test_standard_states_override_analysis(self):
        recs = self.analyzer.get_default_recommendations(standard_states=True)

        self.assertEqual(recs["A_10_ASP"], "ASP")  # acid deprotonated
        self.assertEqual(recs["A_20_GLU"], "GLU")  # acid deprotonated
        self.assertEqual(recs["A_30_HIS"], "HIE")  # epsilon tautomer
        self.assertEqual(recs["A_40_CYS"], "CYS")  # neutral thiol
        self.assertEqual(recs["A_50_LYS"], "LYS")  # protonated base
        self.assertEqual(recs["A_60_TYR"], "TYR")  # protonated phenol
        self.assertEqual(recs["A_70_ARG"], "ARG")  # protonated base (passthrough)
        self.assertEqual(recs["A_80_PRN"], "PRD")  # heme propionate, deprotonated

    def test_prn_to_prd_only_in_standard_mode(self):
        # The PROPKA-based default leaves PRN as PRN (it has no PRD/PRP notion);
        # only standard mode renames it to the deprotonated acid.
        default = self.analyzer.get_default_recommendations()
        self.assertEqual(default["A_80_PRN"], "PRN")

        standard = self.analyzer.get_default_recommendations(standard_states=True)
        self.assertEqual(standard["A_80_PRN"], "PRD")

    def test_standard_differs_from_propka_for_titratables(self):
        # Sanity: with the opposite-of-standard analysis above, the PROPKA path
        # produces the protonated/charged forms, confirming the two paths really
        # diverge (so the override in test 1 is meaningful).
        default = self.analyzer.get_default_recommendations()
        self.assertEqual(default["A_10_ASP"], "ASH")  # protonated acid
        self.assertEqual(default["A_50_LYS"], "LYN")  # neutral lysine
        self.assertEqual(default["A_30_HIS"], "HIP")  # doubly protonated


if __name__ == "__main__":
    unittest.main()
