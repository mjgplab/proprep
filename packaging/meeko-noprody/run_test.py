"""Run by conda-build after installing the package (see meta.yaml)."""
import importlib

from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy

# A ligand prepared end to end: atom types, torsion tree, PDBQT text.
mol = Chem.AddHs(Chem.MolFromSmiles("OCC(=O)O"))
AllChem.EmbedMolecule(mol, randomSeed=1)
setup = MoleculePreparation().prepare(mol)[0]
text, ok, error = PDBQTWriterLegacy.write_string(setup)
assert ok, error
types = [line.split()[-1] for line in text.splitlines() if line.startswith("ATOM")]
assert types.count("OA") == 3 and types.count("HD") == 2, types       # glycolic acid's O and polar H
assert "TORSDOF 3" in text, text
print("ligand prepared:", len(types), "atoms, types", sorted(set(types)))

# Without prody, the one module that needs it says so.
try:
    importlib.import_module("meeko.covalentbuilder")
except ImportError as missing:
    assert "prody" in str(missing), missing
    print("meeko.covalentbuilder needs prody, as expected")
else:
    raise AssertionError("meeko.covalentbuilder imported without prody")
