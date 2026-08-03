"""
Regression test for PROPKA determinant parsing.

ProtonationStateAnalyzer._parse_propka_determinants extracts the per-residue
breakdown that explains each predicted pKa (desolvation + side-chain/backbone/
Coulombic determinants), and — critically for cofactor-heavy systems — flags
which contributions come from ligands (partner shown as RESNAME ATOM CHAIN
rather than RESNAME RESNUM CHAIN).

The fixture below is copied verbatim from real propka3 .pka output (fixed-width
columns matter — the parser locates columns from the dashed separator line).
"""
import os
import tempfile
import unittest

from proprep.structure_prep.protonation_worker import ProtonationStateAnalyzer

# Verbatim slice of a propka3 .pka detailed-determinant table. Column alignment
# is significant; do not reflow.
FIXTURE = (
    "---------  -----   ------   ---------------------    --------------    --------------    --------------\n"
    "                            DESOLVATION  EFFECTS       SIDECHAIN          BACKBONE        COULOMBIC\n"
    " RESIDUE    pKa    BURIED     REGULAR      RE        HYDROGEN BOND     HYDROGEN BOND      INTERACTION\n"
    "---------  -----   ------   ---------   ---------    --------------    --------------    --------------\n"
    "\n"
    "ASP  86 A   3.25    29 %    0.83  363   0.22    0   -0.35 LYS  95 A   -0.70 CYS  87 A   -0.01 LYS 610 A\n"
    "ASP  86 A                                            0.00 XXX   0 X    0.00 XXX   0 X   -0.01 FMN  N5 A\n"
    "ASP  86 A                                            0.00 XXX   0 X    0.00 XXX   0 X   -0.53 LYS  95 A\n"
    "\n"
    "GLU 150 A   4.49     0 %    0.29  279   0.00    0   -0.28 FMN O3' A    0.00 XXX   0 X   -0.01 LYS  38 A\n"
    "\n"
    "--------------------------------------------------------------------------------------------------------\n"
    "SUMMARY OF THIS PREDICTION\n"
)


class PropkaDeterminantParseTest(unittest.TestCase):
    def setUp(self):
        # Instance without running __init__ — we only exercise the pure parser.
        self.analyzer = ProtonationStateAnalyzer.__new__(ProtonationStateAnalyzer)
        fd, self.path = tempfile.mkstemp(suffix=".pka")
        with os.fdopen(fd, "w") as fh:
            fh.write(FIXTURE)

    def tearDown(self):
        os.remove(self.path)

    def test_parses_determinants_and_flags_ligands(self):
        d = self.analyzer._parse_propka_determinants(self.path)

        # Keyed as chain_resnum_type (matches custom_pkas / results keys).
        self.assertIn("A_86_ASP", d)
        self.assertIn("A_150_GLU", d)

        asp = d["A_86_ASP"]
        self.assertAlmostEqual(asp["pKa"], 3.25)
        self.assertAlmostEqual(asp["desolvation_regular"], 0.83)
        self.assertAlmostEqual(asp["desolvation_re"], 0.22)
        # Null "XXX" cells dropped; 5 real determinants for Asp86
        # (sidechain LYS95, backbone CYS87, Coulombic LYS610/FMN-N5/LYS95).
        self.assertEqual(len(asp["determinants"]), 5)

        # The FMN N5 contribution must be present and flagged as a ligand.
        fmn = [x for x in asp["determinants"]
               if x["partner_resname"] == "FMN" and x["partner_id"] == "N5"]
        self.assertEqual(len(fmn), 1)
        self.assertTrue(fmn[0]["is_ligand"])
        self.assertEqual(fmn[0]["kind"], "coulombic")

        # Protein partners are NOT flagged as ligands.
        lys = [x for x in asp["determinants"] if x["partner_resname"] == "LYS"]
        self.assertTrue(lys and all(not x["is_ligand"] for x in lys))

        # Glu150's side-chain determinant is the FMN O3' cofactor oxygen.
        glu = d["A_150_GLU"]
        sc = [x for x in glu["determinants"] if x["kind"] == "sidechain"]
        self.assertEqual(len(sc), 1)
        self.assertEqual((sc[0]["partner_resname"], sc[0]["partner_id"]), ("FMN", "O3'"))
        self.assertTrue(sc[0]["is_ligand"])
        self.assertAlmostEqual(sc[0]["value"], -0.28)

    def test_missing_table_returns_empty(self):
        fd, empty = tempfile.mkstemp(suffix=".pka")
        with os.fdopen(fd, "w") as fh:
            fh.write("no determinant table here\n")
        try:
            self.assertEqual(self.analyzer._parse_propka_determinants(empty), {})
        finally:
            os.remove(empty)


if __name__ == "__main__":
    unittest.main()
