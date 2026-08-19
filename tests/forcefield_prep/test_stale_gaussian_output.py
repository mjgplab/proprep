"""
A Gaussian log describing a different model than its input must be refused.

Reported after step 14 fitted RESP at -1 for a site whose large_resp.gjf says
-2. Step 12 had been re-run, rebuilding the large model, but Gaussian had not
been re-run, so large_resp.log was left over from the previous model:

    site 1   large_resp.gjf  89 atoms   large_resp.log  ran on 104 atoms

The charge disagreement is the visible symptom, not the defect. The ESP grid in
that log belongs to the old model, so a RESP fit against it is meaningless under
any charge constraint -- and it fails silently, because the output still looks
like charges.

The test is on CONTENT, not timestamps. Modification times are not evidence
about what a file contains: copying a log back from a cluster rewrites them,
and re-running step 12 rewrites an input that may be byte-identical to the one
that ran. Site 2 of the same run is exactly that case -- regenerated input,
identical geometry, valid log -- and an mtime test called it stale.

Compare against the log's ``Symbolic Z-matrix`` echo, which is the input
verbatim. The ``Input orientation`` table is NOT interchangeable: Gaussian
re-centers and reorients the molecule there, so comparing it against the .gjf
reports rigid-body motion as a difference.
"""

import pytest

from proprep.forcefield_prep.structure_preprocessor import StructurePreprocessor


ATOMS = [
    ("N", -36.960, -14.480, -45.924),
    ("C", -37.673, -15.521, -46.439),
    ("FE", -38.838, -15.236, -47.072),
]


def _gjf(atoms, charge=-2, mult=11, radii=True):
    body = "\n".join(f"{e}   -1  {x:.8f}  {y:.8f}  {z:.8f}" for e, x, y, z in atoms)
    text = (f"%chk=large_resp.chk\n"
            f"#P B3LYP/6-31G* Pop(MK,ReadRadii)\n\n"
            f"large model resp\n\n"
            f"{charge} {mult}\n{body}\n")
    if radii:
        # Looks exactly like atom lines; must not be read as geometry.
        text += "\nFe 1.383\nFe 1.383\n\n"
    return text


def _log(atoms, charge=-2, mult=11, symbolic=True):
    body = "\n".join(f" {e}                    -1   {x}    {y}    {z} "
                     for e, x, y, z in atoms)
    head = " Entering Gaussian System\n"
    if symbolic:
        head += (f" Symbolic Z-matrix:\n"
                 f" Charge = {charge} Multiplicity = {mult}\n{body}\n\n")
    return head + (f" Charge = {charge} Multiplicity = {mult}\n"
                   f" NAtoms=   {len(atoms)}\n"
                   f" Normal termination of Gaussian 16\n")


def _models(tmp_path, gjf_atoms, log_atoms, **kw):
    d = tmp_path / "models"
    d.mkdir()
    (d / "large_resp.gjf").write_text(_gjf(gjf_atoms, **kw))
    (d / "large_resp.log").write_text(_log(log_atoms))
    (d / "large.pdb").write_text("REMARK Total charge: -2.0\nEND\n")
    return d


# --------------------------------------------------------------------------- #
# the guard
# --------------------------------------------------------------------------- #

def test_a_matching_model_passes(tmp_path):
    models = _models(tmp_path, ATOMS, ATOMS)

    assert StructurePreprocessor._stale_gaussian_output(models) is None


def test_a_different_atom_count_is_stale(tmp_path):
    """Site 1 of the reported run: 89 atoms in, 104 in the log."""
    models = _models(tmp_path, ATOMS[:2], ATOMS)

    reason = StructurePreprocessor._stale_gaussian_output(models)

    assert reason is not None
    assert "2 atoms" in reason and "3" in reason


def test_a_substituted_residue_is_stale(tmp_path):
    """Same atom count, different elements -- e.g. a gap residue swapped."""
    swapped = [("O", *ATOMS[0][1:]), *ATOMS[1:]]
    models = _models(tmp_path, swapped, ATOMS)

    reason = StructurePreprocessor._stale_gaussian_output(models)

    assert reason is not None and "element" in reason


def test_a_moved_model_is_stale(tmp_path):
    """Same atoms, different coordinates: the model was rebuilt."""
    moved = [(e, x + 1.7, y, z) for e, x, y, z in ATOMS]
    models = _models(tmp_path, moved, ATOMS)

    reason = StructurePreprocessor._stale_gaussian_output(models)

    assert reason is not None and "coordinates" in reason


def test_trailing_zeros_are_not_a_difference(tmp_path):
    """
    The .gjf writes -36.96000000 and the log echoes -36.96. Gaussian strips
    trailing zeros rather than rounding, so real files agree exactly; this
    pins that the two spellings parse to the same number.
    """
    models = _models(tmp_path, ATOMS, ATOMS)

    assert "-36.96000000" in (models / "large_resp.gjf").read_text()
    assert "-36.96 " in (models / "large_resp.log").read_text()
    assert StructurePreprocessor._stale_gaussian_output(models) is None


def test_a_difference_below_tolerance_passes(tmp_path):
    """Headroom for a format that does round, without admitting a real move."""
    nudged = [(e, x + 5e-4, y, z) for e, x, y, z in ATOMS]
    models = _models(tmp_path, ATOMS, nudged)

    assert StructurePreprocessor._stale_gaussian_output(models) is None


# --------------------------------------------------------------------------- #
# what must NOT be reported as stale
# --------------------------------------------------------------------------- #

def test_a_regenerated_but_identical_input_passes(tmp_path):
    """
    Site 2 of the reported run. Step 12 rewrote the input with byte-identical
    geometry and Gaussian was never re-run; the log is still valid. An mtime
    test failed this.
    """
    models = _models(tmp_path, ATOMS, ATOMS)
    # Rewriting the file is exactly what step 12 does.
    (models / "large_resp.gjf").write_text(_gjf(ATOMS))

    assert StructurePreprocessor._stale_gaussian_output(models) is None


def test_timestamps_are_not_consulted(tmp_path):
    """A log copied back from a cluster carries a fresh mtime; content rules."""
    import os

    models = _models(tmp_path, ATOMS[:2], ATOMS)
    os.utime(models / "large_resp.log", (2_000_000_000,) * 2)   # far future

    assert StructurePreprocessor._stale_gaussian_output(models) is not None


def test_readradii_entries_are_not_read_as_geometry(tmp_path):
    """
    The MK radii block sits after the geometry and looks just like it
    ("Fe 1.383"), so reading to end-of-file would inflate the input and report
    every cluster site as stale.
    """
    models = _models(tmp_path, ATOMS, ATOMS, radii=True)

    assert len(StructurePreprocessor._gjf_geometry(models / "large_resp.gjf")) == 3
    assert StructurePreprocessor._stale_gaussian_output(models) is None


# --------------------------------------------------------------------------- #
# unreadable inputs must not be guessed at
# --------------------------------------------------------------------------- #

def test_a_missing_log_is_not_staleness(tmp_path):
    """That is 'run Gaussian', which the caller already reports."""
    models = _models(tmp_path, ATOMS, ATOMS)
    (models / "large_resp.log").unlink()

    assert StructurePreprocessor._stale_gaussian_output(models) is None


def test_a_missing_gjf_is_not_staleness(tmp_path):
    models = _models(tmp_path, ATOMS, ATOMS)
    (models / "large_resp.gjf").unlink()

    assert StructurePreprocessor._stale_gaussian_output(models) is None


def test_a_log_without_the_input_echo_is_not_staleness(tmp_path):
    """No Symbolic Z-matrix means nothing to compare, not a mismatch."""
    d = tmp_path / "models"
    d.mkdir()
    (d / "large_resp.gjf").write_text(_gjf(ATOMS))
    (d / "large_resp.log").write_text(_log(ATOMS, symbolic=False))

    assert StructurePreprocessor._stale_gaussian_output(d) is None


def test_a_missing_directory_is_not_an_error(tmp_path):
    assert StructurePreprocessor._stale_gaussian_output(tmp_path / "nope") is None


# --------------------------------------------------------------------------- #
# the geometry parsers
# --------------------------------------------------------------------------- #

def test_both_files_use_the_same_layout(tmp_path):
    """
    element, frozen-atom flag, x, y, z -- in the .gjf AND the log echo. Reading
    the flag as a coordinate is what made an identical model look 1.7 A apart.
    """
    d = _models(tmp_path, ATOMS, ATOMS)

    assert (StructurePreprocessor._gjf_geometry(d / "large_resp.gjf")
            == StructurePreprocessor._log_input_geometry(d / "large_resp.log")
            == [("N", -36.960, -14.480, -45.924),
                ("C", -37.673, -15.521, -46.439),
                ("FE", -38.838, -15.236, -47.072)])


def test_geometry_without_a_freeze_flag_parses():
    parsed = StructurePreprocessor._parse_gaussian_geometry(
        ["N  -36.96  -14.48  -45.924\n"])

    assert parsed == [("N", -36.96, -14.48, -45.924)]


def test_element_case_does_not_matter():
    parsed = StructurePreprocessor._parse_gaussian_geometry(
        ["fe  -1  1.0  2.0  3.0\n"])

    assert parsed[0][0] == "FE"


def test_unparseable_files_return_none(tmp_path):
    assert StructurePreprocessor._gjf_geometry(tmp_path / "nope.gjf") is None
    assert StructurePreprocessor._log_input_geometry(tmp_path / "nope.log") is None


# --------------------------------------------------------------------------- #
# the charge the stale log would have supplied
# --------------------------------------------------------------------------- #

def test_the_stale_log_disagrees_with_the_input(tmp_path):
    """Without the guard the log's charge wins -- the reported -1 vs -2."""
    d = tmp_path / "models"
    d.mkdir()
    (d / "large_resp.gjf").write_text(_gjf(ATOMS[:2], charge=-2))
    (d / "large_resp.log").write_text(_log(ATOMS, charge=-1))

    assert StructurePreprocessor._esp_charge_multiplicity(d)[0] == -1
    assert StructurePreprocessor._stale_gaussian_output(d) is not None


def test_a_current_log_is_still_preferred(tmp_path):
    """The guard must not disturb the normal path."""
    models = _models(tmp_path, ATOMS, ATOMS)

    assert StructurePreprocessor._esp_charge_multiplicity(models) == (-2, 11)
    assert StructurePreprocessor._stale_gaussian_output(models) is None
