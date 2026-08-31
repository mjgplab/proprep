"""Regression: parse_propka_summary must not drop residues numbered >= 1000.

ProPKA writes the SUMMARY "Group" column fixed-width, so the separator between
the residue name and its number disappears once the number reaches four
digits::

    ASP   4 A     2.34       3.80      <- 3-digit, space present
    LYS1150 A    11.83      10.50      <- 4-digit, space absent

_SUMMARY_RE originally required ``\\s+`` there, so every titratable residue
numbered >= 1000 failed to match and was silently omitted from the PB-vs-ProPKA
comparison. On a 9YUQ cytochrome filament (6 protomers, residues 1-1188) that
dropped the entire sixth protomer -- 15 protein sites -- and the loss was
invisible because propionates legitimately have no ProPKA value either, so the
tally just read "53 with no ProPKA value".

The lines below are verbatim from a real run
(pb_titrate_work/propka_compare/min_for_propka.pka).
"""
from proprep.pb_titrate.propka_compare import parse_propka_summary


SUMMARY = """\
propka3 output

SUMMARY OF THIS PREDICTION
       Group      pKa  model-pKa   ligand atom-type
   ASP   4 A     2.34       3.80                      
   ASP 122 A     8.47       3.80                      
   LYS 996 A    11.53      10.50                      
   ASP1012 A     3.41       3.80                      
   TYR1045 A    11.02      10.00                      
   LYS1150 A    11.83      10.50                      
   LYS1165 A    12.34      10.50                      
   LYS1188 A    10.53      10.50                      
   N+    1 A     7.87       8.00                      
   C-  198 A     3.20       3.20                      
   HCO  NB A    -5.43       5.00                NAR   
   HCOND12 A     1.35       5.00                NAR   
   HCONE21 A    -1.66       5.00                NAR   
--------------------------------------------------------
Free energy of folding
"""


def _parse(tmp_path):
    p = tmp_path / "min_for_propka.pka"
    p.write_text(SUMMARY)
    return parse_propka_summary(p)


def test_four_digit_residue_numbers_are_parsed(tmp_path):
    """The bug: these five were dropped entirely."""
    out = _parse(tmp_path)
    assert out[("ASP", 1012)] == 3.41
    assert out[("TYR", 1045)] == 11.02
    assert out[("LYS", 1150)] == 11.83
    assert out[("LYS", 1165)] == 12.34
    assert out[("LYS", 1188)] == 10.53


def test_three_digit_residue_numbers_still_parsed(tmp_path):
    """No regression on the format that already worked."""
    out = _parse(tmp_path)
    assert out[("ASP", 4)] == 2.34
    assert out[("ASP", 122)] == 8.47
    assert out[("LYS", 996)] == 11.53


def test_termini_and_ligand_rows_still_skipped(tmp_path):
    """\\s* must not let the non-residue rows through.

    Termini ("N+", "C-") fail on [A-Z]{2,3}; heme ligand rows ("HCO  NB",
    "HCOND12", "HCONE21") fail because (\\d+) cannot match a group column
    starting with a letter, with or without intervening whitespace.
    """
    out = _parse(tmp_path)
    assert all(rn not in ("N+", "C-", "HCO", "HC") for rn, _ in out)
    # HCOND12 must not be misread as HCO + some number.
    assert ("HCO", 12) not in out
    assert ("HCO", 21) not in out
    assert len(out) == 8   # exactly the 8 standard residues above
