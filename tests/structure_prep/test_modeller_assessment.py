"""
Regression tests for the rigorous MODELLER model-assessment reporting in
ModellerInterface (structure_completeness.py).

The old reporting normalized MOLPDF by residue count and bucketed the result
into invented "excellent/good/acceptable/poor" tiers with prose verdicts
("near X-ray quality", "ready for MD"). That heuristic was replaced with:

  * Global scores read straight from the MODELLER API (no stdout scraping):
    normalized DOPE z-score (Shen & Sali 2006) and GA341 (John & Sali 2003),
    with raw DOPE / molpdf shown for provenance only.
  * A per-residue assessment of the rebuilt region: each MODELLER-built residue
    is compared to the same chain's experimentally-resolved residues via a
    within-structure normalized-DOPE z-score (self-calibrating, no universal
    cutoff), flagged at z > BUILT_RESIDUE_Z_THRESHOLD.

The pure-logic tests run anywhere; the per-residue tests require a working
MODELLER install + license and are skipped otherwise.
"""
import io
import os
import unittest

from rich.console import Console

from proprep.structure_prep import structure_completeness as sc
from proprep.structure_prep.structure_completeness import ModellerInterface

# A small real protein shipped in the repo (BPTI, 58 residues, blank chain id).
_BPTI = os.path.join(
    os.path.dirname(sc.__file__),
    "..", "pb_titrate", "tests", "bpti_asp3", "5pti_cleaned.pdb",
)


def _modeller_ready():
    """True if MODELLER imports AND a license lets us build an Environ."""
    if not sc.HAS_MODELLER:
        return False
    key_file = os.path.expanduser("~/.proprep/modeller_key")
    if "KEY_MODELLER" not in os.environ and os.path.exists(key_file):
        with open(key_file) as fh:
            os.environ["KEY_MODELLER"] = fh.read().strip()
    try:
        import modeller
        modeller.Environ()
        return True
    except Exception:
        return False


MODELLER_READY = _modeller_ready()


class ExtractAssessmentTest(unittest.TestCase):
    """_extract_assessment reads mdl.outputs without stdout parsing."""

    def setUp(self):
        self.mi = ModellerInterface(Console(file=io.StringIO()))

    def test_reads_all_scores(self):
        class FakeMdl:
            # GA341 is returned by MODELLER as a list whose [0] is the score.
            outputs = [{
                "Normalized DOPE score": -1.42,
                "GA341 score": [0.95, 0.0, 0.0],
                "DOPE score": -2543.7,
                "molpdf": 1234.5,
            }]

        a = self.mi._extract_assessment(FakeMdl())
        self.assertAlmostEqual(a["normalized_dope"], -1.42)
        self.assertAlmostEqual(a["ga341"], 0.95)
        self.assertAlmostEqual(a["dope"], -2543.7)
        self.assertAlmostEqual(a["molpdf"], 1234.5)

    def test_missing_outputs_degrade_to_none(self):
        # No outputs attribute, empty list, and missing keys must not raise.
        for mdl in (object(), type("M", (), {"outputs": []})()):
            a = self.mi._extract_assessment(mdl)
            self.assertIsNone(a["normalized_dope"])
            self.assertIsNone(a["ga341"])


class InterpretationBandTest(unittest.TestCase):
    """Calibrated bands match the published thresholds."""

    def setUp(self):
        self.mi = ModellerInterface(Console(file=io.StringIO()))

    def test_normalized_dope_bands(self):
        self.assertIn("Native-like", self.mi._interpret_normalized_dope(-1.0)[0])
        self.assertIn("Native-like", self.mi._interpret_normalized_dope(-2.3)[0])
        self.assertIn("Borderline", self.mi._interpret_normalized_dope(-0.5)[0])
        self.assertIn("Likely incorrect", self.mi._interpret_normalized_dope(0.4)[0])

    def test_ga341_bands(self):
        self.assertIn("Reliable", self.mi._interpret_ga341(0.7)[0])
        self.assertIn("Reliable", self.mi._interpret_ga341(1.0)[0])
        self.assertIn("Low", self.mi._interpret_ga341(0.69)[0])

    def test_display_never_crashes_without_modeller(self):
        # Renders the global table even when no scores are present.
        from Bio.PDB import PDBParser
        struct = PDBParser(QUIET=True).get_structure("s", _BPTI)
        assessment = {"normalized_dope": None, "ga341": None,
                      "dope": None, "molpdf": None}
        self.mi._display_quality_assessment(assessment, struct, None)


@unittest.skipUnless(MODELLER_READY, "MODELLER unavailable or unlicensed")
class BuiltRegionAssessmentTest(unittest.TestCase):
    """Per-residue rebuilt-region assessment against the live MODELLER API."""

    def setUp(self):
        import modeller
        self.mi = ModellerInterface(Console(file=io.StringIO()))
        self.env = modeller.Environ()
        self.env.io.hetatm = True
        self.env.io.atom_files_directory = [os.path.dirname(os.path.abspath(_BPTI))]
        self.env.libs.topology.read(file="${LIB}/top_heav.lib")
        self.env.libs.parameters.read(file="${LIB}/par.lib")

    def test_built_residues_scored_against_chain_baseline(self):
        built = {("", 20), ("", 21), ("", 22)}
        res = self.mi._assess_built_region(self.env, _BPTI, built)
        self.assertNotIn("error", res)
        self.assertEqual(res["n_built"], 3)
        # Every built residue gets a z-score (58-residue baseline available).
        for r in res["residues"]:
            self.assertIsNotNone(r["z"])

    def test_chain_mismatch_falls_back_to_number_matching(self):
        # Real chain id is "" — wrong chain "A" must still match by number.
        res = self.mi._assess_built_region(self.env, _BPTI, {("A", 20), ("A", 21)})
        self.assertNotIn("error", res)
        self.assertEqual(res["n_built"], 2)

    def test_failure_is_swallowed_not_raised(self):
        # A nonexistent file must yield {'error': ...}, never an exception.
        res = self.mi._assess_built_region(self.env, "/no/such/model.pdb",
                                           {("A", 1)})
        self.assertIn("error", res)


class HomologyModelerAssessmentTest(unittest.TestCase):
    """The homology-modeling path (full de novo model) reports the same
    calibrated GLOBAL scores. There is no rebuilt-region baseline here — the
    whole model is built — so only Pass 1 applies."""

    def test_extract_and_bands(self):
        from proprep.structure_prep.homology_modeler import HomologyModeler as H

        class FakeMdl:
            outputs = [{
                "Normalized DOPE score": -1.21,
                "GA341 score": [0.92, 0.0],
                "DOPE score": -3201.4,
                "molpdf": 987.6,
            }]

        a = H._extract_assessment(FakeMdl())
        self.assertAlmostEqual(a["normalized_dope"], -1.21)
        self.assertAlmostEqual(a["ga341"], 0.92)
        self.assertIsNone(H._extract_assessment(object())["ga341"])
        self.assertIn("Native-like", H._interpret_normalized_dope(-1.2)[0])
        self.assertIn("Likely incorrect", H._interpret_normalized_dope(0.3)[0])
        self.assertIn("Reliable", H._interpret_ga341(0.7)[0])

    def test_display_handles_missing_scores(self):
        from proprep.structure_prep.homology_modeler import HomologyModeler as H

        none_scores = {"normalized_dope": None, "ga341": None,
                       "dope": None, "molpdf": None}
        # No structure and no scores must still render without raising.
        H.display_quality_assessment(none_scores, None, Console(file=io.StringIO()))


if __name__ == "__main__":
    unittest.main()
