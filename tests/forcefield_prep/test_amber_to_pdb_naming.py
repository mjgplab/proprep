"""
The intermediate PDB must keep the topology's atom names.

ProPrep converts tLEaP's prmtop/rst7 back to PDB and reloads that file against
the same libraries on the next pass. ambpdb translates to PDB v3 conventions by
default, so an externally supplied library using older names had them rewritten
in between: O1P/O2P became OP1/OP2.

The first (demetallated) build succeeded, because it loads the ORIGINAL
coordinates. The second failed on a file tLEaP had never seen:

    FATAL:  Atom .R<FAD 1311>.A<OP2 85> does not have a type.

The topology was correct throughout. Only the PDB in between was translated,
and it is the one that gets reloaded -- so it is an internal artifact that
should speak the libraries' naming, not the PDB standard's. ``-aatm`` writes
names as the topology holds them.

Measured on the reported system: 2 atoms out of 20479 change, exactly the two
that broke, with byte-identical column layout and the element column intact.
"""

import subprocess
from pathlib import Path

import pytest
from rich.console import Console

from proprep.forcefield_prep.structure_preprocessor import StructurePreprocessor


def _preprocessor():
    p = StructurePreprocessor.__new__(StructurePreprocessor)
    p.console = Console(quiet=True)
    return p


class _Result:
    def __init__(self, stdout="ATOM      1  N   MET A   1\n", stderr=""):
        self.stdout, self.stderr = stdout, stderr


# --------------------------------------------------------------------------- #
# the flag
# --------------------------------------------------------------------------- #

def test_the_naming_flag_is_passed(monkeypatch, tmp_path):
    """Without it ambpdb rewrites names the next tLEaP pass has to match."""
    captured = {}

    def fake_run(cmd, **kwargs):
        captured["cmd"] = cmd
        return _Result()

    monkeypatch.setattr(subprocess, "run", fake_run)
    _preprocessor()._convert_amber_to_pdb("p.prmtop", "c.rst7", tmp_path / "o.pdb")

    assert captured["cmd"][0] == "ambpdb"
    assert "-aatm" in captured["cmd"]


def test_the_output_is_written(monkeypatch, tmp_path):
    monkeypatch.setattr(subprocess, "run", lambda cmd, **k: _Result("ATOM  X\n"))
    out = tmp_path / "o.pdb"

    _preprocessor()._convert_amber_to_pdb("p", "c", out)

    assert out.read_text() == "ATOM  X\n"


# --------------------------------------------------------------------------- #
# failures must stop the step, not be written over silently
# --------------------------------------------------------------------------- #

def test_a_missing_ambpdb_raises(monkeypatch, tmp_path):
    """
    The caller writes this file and then reads it back. Continuing would fail
    later on a missing file instead of here on the real cause.
    """
    def boom(cmd, **kwargs):
        raise FileNotFoundError("ambpdb")

    monkeypatch.setattr(subprocess, "run", boom)

    with pytest.raises(RuntimeError, match="AmberTools"):
        _preprocessor()._convert_amber_to_pdb("p", "c", tmp_path / "o.pdb")


def test_a_failing_ambpdb_raises_with_its_message(monkeypatch, tmp_path):
    def boom(cmd, **kwargs):
        raise subprocess.CalledProcessError(
            1, cmd, stderr="Error: could not read topology\n")

    monkeypatch.setattr(subprocess, "run", boom)

    with pytest.raises(RuntimeError, match="could not read topology"):
        _preprocessor()._convert_amber_to_pdb("p", "c", tmp_path / "o.pdb")


def test_empty_output_raises(monkeypatch, tmp_path):
    monkeypatch.setattr(subprocess, "run", lambda cmd, **k: _Result(stdout="  \n"))

    with pytest.raises(RuntimeError, match="no output"):
        _preprocessor()._convert_amber_to_pdb("p", "c", tmp_path / "o.pdb")


def test_no_file_is_left_behind_on_failure(monkeypatch, tmp_path):
    monkeypatch.setattr(subprocess, "run", lambda cmd, **k: _Result(stdout=""))
    out = tmp_path / "o.pdb"

    with pytest.raises(RuntimeError):
        _preprocessor()._convert_amber_to_pdb("p", "c", out)

    assert not out.exists()


def test_there_is_no_second_converter(monkeypatch, tmp_path):
    """
    cpptraj used to be a fallback. Both ship with AmberTools, so if ambpdb is
    missing cpptraj is too -- and it writes a third naming convention.
    """
    calls = []

    def record(cmd, **kwargs):
        calls.append(cmd[0] if isinstance(cmd, list) else cmd)
        raise FileNotFoundError("ambpdb")

    monkeypatch.setattr(subprocess, "run", record)

    with pytest.raises(RuntimeError):
        _preprocessor()._convert_amber_to_pdb("p", "c", tmp_path / "o.pdb")

    assert calls == ["ambpdb"]


# --------------------------------------------------------------------------- #
# against real AmberTools, if present
# --------------------------------------------------------------------------- #

def _ambpdb_available():
    try:
        subprocess.run(["ambpdb", "--help"], capture_output=True, check=False)
        return True
    except FileNotFoundError:
        return False


needs_amber = pytest.mark.skipif(not _ambpdb_available(),
                                 reason="ambpdb not on PATH")


@needs_amber
def test_the_flag_is_what_preserves_v2_names(tmp_path):
    """
    Pins the actual behaviour rather than trusting the flag's description: the
    default translates the phosphate oxygens and -aatm does not.
    """
    help_text = subprocess.run(["ambpdb", "--help"], capture_output=True,
                               text=True).stdout + subprocess.run(
        ["ambpdb", "--help"], capture_output=True, text=True).stderr

    assert "-aatm" in help_text, "ambpdb no longer offers -aatm"
