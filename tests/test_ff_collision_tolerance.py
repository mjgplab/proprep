"""Regression: the FF-collision detector must treat rounding-level numeric
differences as harmless duplication, not conflicts.

When a ligand also coordinates a metal, the MCPB site frcmod re-lists that
ligand's GAFF bonds rounded to one decimal (e.g. `ca-ca 354.2`), while the
standalone ligand frcmod keeps antechamber's `ca-ca 354.25`. These are the SAME
parameter; exact `!=` comparison flagged all of them and tripped the Topology
Generator's collision safety net (`_verify_no_ff_collisions_or_abort`), blocking
the build. The comparison now uses a numeric tolerance while still catching
genuinely different values (force constants, phases, periodicities).
"""

from proprep.ff_compat.matrix import _values_equivalent, _dict_diffs


def test_rounding_is_equivalent():
    # One-decimal rounding of a bond (K, r0) is not a conflict.
    assert _values_equivalent((354.25, 1.399), (354.2, 1.3990))
    assert _values_equivalent((332.45, 1.383), (332.4, 1.3830))
    # Scalars (MASS) too.
    assert _values_equivalent(12.01, 12.01)


def test_real_differences_still_flagged():
    assert not _values_equivalent((300.0, 1.5), (350.0, 1.5))      # different K
    assert not _values_equivalent((10.0, 0.0, 2), (10.0, 180.0, 2))  # phase 0 vs 180
    assert not _values_equivalent((1.1, 180.0, 2), (1.1, 180.0, 3))  # periodicity


def test_non_numeric_requires_exact_match():
    # Library unit-type names must match exactly (no tolerance concept).
    assert _values_equivalent("HD1", "HD1")
    assert not _values_equivalent("HD1", "HE1")


def test_dict_diffs_ignores_rounding_but_keeps_conflicts():
    a = {("ca", "ca"): (354.25, 1.399), ("c", "n"): (317.5, 1.397)}
    b = {("ca", "ca"): (354.2, 1.3990), ("c", "n"): (490.0, 1.335)}
    diffs = _dict_diffs(a, b)
    keys = {d[0] for d in diffs}
    assert ("ca", "ca") not in keys      # rounding -> not a conflict
    assert ("c", "n") in keys            # genuinely different -> flagged


def test_length_mismatch_is_a_difference():
    assert not _values_equivalent((1.0, 2.0), (1.0, 2.0, 3.0))


# --------------------------------------------------------------------------- #
# parse_set: a metal site's frcmod field is a LIST (bonded + ligand GAFF)
# --------------------------------------------------------------------------- #
#
# discover_forcefield_files returns `frcmod` as a list for a metal site (the
# MCPB bonded frcmod plus each organic ligand's own GAFF frcmod). parse_set once
# took a single path and Path(list) crashed the collision safety net; it now
# parses and merges every frcmod in the list into one signature.

from proprep.ff_compat.parser import parse_set


def _frcmod(path, mass_line, nonb_line):
    path.write_text(
        f"remark\nMASS\n{mass_line}\n\nBOND\n\nANGLE\n\nDIHE\n\n"
        f"IMPROPER\n\nNONBON\n{nonb_line}\n"
    )


def test_parse_set_merges_list_of_frcmods(tmp_path):
    bonded = tmp_path / "site_1_bonded.frcmod"
    _frcmod(bonded, "M1 54.94", "  M1  1.20  0.010")   # metal type
    ligand = tmp_path / "e4z.frcmod"
    _frcmod(ligand, "ca 12.01", "  ca  1.90  0.086")   # ligand GAFF type

    sig = parse_set("dimn", [str(bonded), str(ligand)], [])

    # Both frcmods' types are present — the list was merged, not dropped.
    assert "M1" in sig.nonb and "ca" in sig.nonb
    # Both source files recorded.
    assert any(s.endswith("site_1_bonded.frcmod") for s in sig.source_files)
    assert any(s.endswith("e4z.frcmod") for s in sig.source_files)


def test_parse_set_still_accepts_single_frcmod_string(tmp_path):
    solo = tmp_path / "solo.frcmod"
    _frcmod(solo, "ca 12.01", "  ca  1.90  0.086")
    sig = parse_set("solo", str(solo), [])  # bare string, backward-compat
    assert "ca" in sig.nonb
