"""
Regression tests: pure inorganic metal clusters (Fe2S2, Fe4S4, ...) are a
distinct triage category (F), NOT organometallic (C).

An organometallic cofactor (heme, MoCo pterin) has a carbon-bearing organic
scaffold that is parameterized separately with the metal removed. A pure
inorganic cluster has no organic fragment — nothing to hand to the
small-molecule parameterizer — so the WHOLE residue (metals + bridging
sulfides) is withheld from the standard-FF tLEaP pass and owned by MCPB, then
reinserted as ONE residue afterward.

The discriminator is: multi-atom, contains a metal, contains no carbon.
"""

from pathlib import Path

from rich.console import Console
from Bio.PDB import PDBParser

from proprep.forcefield_prep.structure_preprocessor import StructurePreprocessor


def _bare_instance(tmp_path):
    inst = StructurePreprocessor.__new__(StructurePreprocessor)
    inst.console = Console()

    class _WS:
        def __init__(self):
            self.d = {}

        def get(self, k, default=None):
            return self.d.get(k, default)

        def set(self, k, v):
            self.d[k] = v

    inst.workspace = _WS()
    inst._output_dir = Path(tmp_path)
    return inst


class _Atom:
    def __init__(self, name, element):
        self.name = name
        self.element = element


class _Residue:
    def __init__(self, resname, atoms):
        self.resname = resname
        self._atoms = atoms

    def get_atoms(self):
        return iter(self._atoms)


def test_pure_cluster_discriminator(tmp_path):
    inst = _bare_instance(tmp_path)

    fes = _Residue("FES", [_Atom("FE1", "FE"), _Atom("FE2", "FE"),
                           _Atom("S1", "S"), _Atom("S2", "S")])
    sf4 = _Residue("SF4", [_Atom(f"FE{i}", "FE") for i in range(1, 5)] +
                          [_Atom(f"S{i}", "S") for i in range(1, 5)])
    heme = _Residue("HEM", [_Atom("FE", "FE"), _Atom("NA", "N"),
                            _Atom("C1A", "C"), _Atom("C2A", "C")])
    mte = _Residue("MTE", [_Atom("MO", "MO"), _Atom("C1", "C"),
                           _Atom("N1", "N"), _Atom("S1", "S")])
    zn = _Residue("ZN", [_Atom("ZN", "ZN")])

    assert inst._is_pure_metal_cluster(fes) is True
    assert inst._is_pure_metal_cluster(sf4) is True
    assert inst._is_pure_metal_cluster(heme) is False   # has carbon scaffold
    assert inst._is_pure_metal_cluster(mte) is False    # has carbon scaffold
    assert inst._is_pure_metal_cluster(zn) is False     # single atom → category D


_PDB = """ATOM      1  N   ALA A   1      0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1      1.500   0.000   0.000  1.00  0.00           C
ATOM      3  C   ALA A   1      2.000   1.400   0.000  1.00  0.00           C
ATOM      4  N   GLY A   2      3.300   1.600   0.000  1.00  0.00           N
ATOM      5  CA  GLY A   2      4.000   2.900   0.000  1.00  0.00           C
HETATM    6 FE1  FES A 100      10.000  10.000  10.000  1.00  0.00          FE
HETATM    7 FE2  FES A 100      12.000  10.000  10.000  1.00  0.00          FE
HETATM    8  S1  FES A 100      11.000  11.500  10.000  1.00  0.00           S
HETATM    9  S2  FES A 100      11.000   8.500  10.000  1.00  0.00           S
END
"""


def test_cluster_round_trip_removed_whole_reinserted_as_one_residue(tmp_path):
    inst = _bare_instance(tmp_path)
    pdb = tmp_path / "input.pdb"
    pdb.write_text(_PDB)
    inst._pdb_file = str(pdb)

    # 1. Extract: every atom of the cluster, sharing one cluster_id.
    atoms = inst._extract_cluster_atoms_from_key("A:100:FES")
    assert len(atoms) == 4
    assert all(a.cluster_id == "A:100:FES" for a in atoms)
    assert sorted(a.element for a in atoms) == ["Fe", "Fe", "S", "S"]

    # 2. Remove: all four atoms gone from the metal-free structure.
    metal_free = inst._remove_metals_from_structure(atoms, str(pdb))
    mfs = PDBParser(QUIET=True).get_structure("mf", str(metal_free))
    assert not [a for a in mfs[0].get_atoms()
                if a.get_parent().get_resname().strip() == "FES"]

    # 3. Reinsert: exactly ONE FES residue holding all four atoms.
    final = inst._insert_metals(metal_free, atoms)
    s = PDBParser(QUIET=True).get_structure("f", str(final))
    fes = [r for r in s[0].get_residues() if r.get_resname().strip() == "FES"]
    assert len(fes) == 1
    assert sorted(a.name.strip() for a in fes[0].get_atoms()) == \
        ["FE1", "FE2", "S1", "S2"]

    # Protein survives; reinsertion map records the cluster's move.
    assert len([r for r in s[0].get_residues()
                if r.get_resname().strip() in ("ALA", "GLY")]) == 2
    assert ("A", 100) in inst.workspace.get("preprocessing_metal_reinsertion_map")
