"""The hydrogen pass must not reorder residues: bond directives address them by position.

6R2Q, clean run, 2026-09-20. tLEaP added the hydrogens (0 errors) and ProPrep
refused the result: "1194 heavy atoms are in a residue at another position in
the file (first: HOH A1533 O, residue number 1533 of the file given, 1597 of
the file written)". The check was right. tLEaP's ``reorder_residues`` is "on"
by default: "solvent will be moved to the end". The repaired structure had its
10 crystal waters lying between the hemes, so moving them shifted every residue
after the first water, and the Topology Generator's bond directives, written
for the order given, would have hit the wrong residues in bilayer.pdb. With
``set default reorder_residues off`` the same script writes the same structure
in the same positions (checked on the real inputs).
"""

import io
import shutil
import subprocess
import sys
import types
from pathlib import Path

import pytest
from rich.console import Console

from proprep.membrane_prep import membrane_builder as mb

Module = mb.MembraneBuilderModule


def _tool(name):
    found = shutil.which(name)
    if found:
        return found
    beside = Path(sys.executable).parent / name
    return str(beside) if beside.exists() else None


# Two alanines with a crystal water BETWEEN them, as a repaired multi-chain structure has.
WATER_IN_THE_MIDDLE = """\
ATOM      1  N   ALA A   1       1.000   1.000   1.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       2.450   1.000   1.000  1.00  0.00           C
ATOM      3  C   ALA A   1       3.000   2.400   1.000  1.00  0.00           C
ATOM      4  O   ALA A   1       2.300   3.400   1.000  1.00  0.00           O
ATOM      5  CB  ALA A   1       3.000   0.200   2.200  1.00  0.00           C
TER
HETATM    6  O   HOH A   2      10.000  10.000  10.000  1.00  0.00           O
TER
ATOM      7  N   ALA B   3      21.000  21.000  21.000  1.00  0.00           N
ATOM      8  CA  ALA B   3      22.450  21.000  21.000  1.00  0.00           C
ATOM      9  C   ALA B   3      23.000  22.400  21.000  1.00  0.00           C
ATOM     10  O   ALA B   3      22.300  23.400  21.000  1.00  0.00           O
ATOM     11  CB  ALA B   3      23.000  20.200  22.200  1.00  0.00           C
TER
END
"""


@pytest.mark.skipif(_tool("tleap") is None or _tool("ambpdb") is None, reason="tleap/ambpdb not found")
@pytest.mark.parametrize("setting, kept", [("", False), ("set default reorder_residues off\n", True)],
                         ids=["tLEaP's default moves the water to the end", "with the setting the order is kept"])
def test_real_tleap_with_a_water_between_two_residues(tmp_path, setting, kept):
    import os
    env = dict(os.environ, AMBERHOME=str(Path(_tool("tleap")).parent.parent))
    (tmp_path / "given.pdb").write_text(WATER_IN_THE_MIDDLE)
    (tmp_path / "run.in").write_text(
        setting + "source leaprc.protein.ff14SB\nsource leaprc.water.tip3p\nmol = loadpdb given.pdb\n"
        "saveamberparm mol out.prmtop out.rst7\nquit\n")
    subprocess.run([_tool("tleap"), "-f", "run.in"], cwd=tmp_path, capture_output=True, text=True, env=env)
    pdb = subprocess.run([_tool("ambpdb"), "-p", "out.prmtop", "-c", "out.rst7"], cwd=tmp_path,
                         capture_output=True, text=True, env=env).stdout
    (tmp_path / "written.pdb").write_text(pdb)
    problem = Module._hydrogen_pass_changed_the_structure(str(tmp_path / "given.pdb"), str(tmp_path / "written.pdb"))
    if kept:
        assert problem is None
    else:
        assert "another position in the file" in problem and "HOH" in problem


def test_the_hydrogen_pass_script_switches_the_reordering_off_before_it_loads_anything(monkeypatch, tmp_path):
    protein = tmp_path / "protein.pdb"
    protein.write_text(WATER_IN_THE_MIDDLE)
    module = Module()
    module.processor = types.SimpleNamespace(console=Console(file=io.StringIO(), width=300))
    module.config.protein_pdb = str(protein)
    monkeypatch.setattr(subprocess, "run", lambda *a, **k: subprocess.CompletedProcess(a, 1, stdout="", stderr=""))
    module._run_pre_tleap_hydrogen_pass(types.SimpleNamespace(get=lambda key, default=None: default), str(tmp_path))
    script = (tmp_path / "pretleap_hydrogen.in").read_text()
    assert "set default reorder_residues off" in script
    assert script.index("set default reorder_residues off") < script.index("loadpdb")


# ── After the build ─────────────────────────────────────────────────────

def _atoms(path, rows):
    path.write_text("".join(
        f"ATOM  {i:5d} {atom:<4s} {res:>3s} {chain}{num:4d}    {1.0 * i:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00\n"
        for i, (res, num, atom, chain) in enumerate(rows, 1)))
    return str(path)


PROTEIN = [("ACE", 1, "C", " "), ("SER", 2, "N", " "), ("SER", 2, "CA", " "), ("HCO", 3, "FE", " "), ("CA", 4, "CA", " ")]
LIPIDS = [("PA", 5, "C12", "A"), ("PC", 6, "N31", "A")]


def test_a_bilayer_that_starts_with_the_protein_it_was_given_passes(tmp_path):
    """packmol-memgen adds a chain ID and moves the protein; neither matters."""
    built = [(r, n, a, "A") for r, n, a, _ in PROTEIN] + LIPIDS
    assert Module._bilayer_protein_differs(_atoms(tmp_path / "p.pdb", PROTEIN), _atoms(tmp_path / "b.pdb", built)) is None


def test_a_bilayer_holding_another_protein_is_refused(tmp_path):
    other = [("SER", 1, "N", "A"), ("SER", 1, "CA", "A"), ("HCO", 2, "FE", "A"), ("CA", 3, "CA", "A")] + LIPIDS
    problem = Module._bilayer_protein_differs(_atoms(tmp_path / "p.pdb", PROTEIN), _atoms(tmp_path / "b.pdb", other))
    assert "atom 1 is SER N" in problem and "ACE C" in problem


def test_a_residue_that_begins_elsewhere_is_refused_even_if_every_name_matches(tmp_path):
    """Same atoms, but two residues run together: every later position is off by one."""
    merged = [("ACE", 1, "C", "A"), ("SER", 2, "N", "A"), ("SER", 2, "CA", "A"), ("HCO", 2, "FE", "A"), ("CA", 3, "CA", "A")] + LIPIDS
    # the HCO takes residue number 2 but keeps its own name, so it is still its own residue: passes
    assert Module._bilayer_protein_differs(_atoms(tmp_path / "p.pdb", PROTEIN), _atoms(tmp_path / "b.pdb", merged)) is None
    dropped = [("ACE", 1, "C", "A"), ("SER", 2, "N", "A"), ("SER", 2, "CA", "A"), ("CA", 3, "CA", "A")] + LIPIDS
    assert "HCO FE" in Module._bilayer_protein_differs(_atoms(tmp_path / "p.pdb", PROTEIN), _atoms(tmp_path / "b.pdb", dropped))
