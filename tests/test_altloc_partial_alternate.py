"""An alternate location that models only part of a residue used to be offered
as an ordinary choice, and picking it dropped every atom it did not model.
In PDB 8ZST the C-terminal LEU 129 has alternates A and B with all eight atoms
and alternate C with only the five side-chain atoms. Choosing C removed the
backbone N, C and O, so tLEaP rebuilt them and closed the chain with a 3.79
Angstrom peptide bond. Regression for that, and for the labels the picker
never asked about, which used to survive into the written PDB."""

from Bio.PDB import PDBParser, PDBIO

from proprep.structure_prep import altloc_picker
from proprep.structure_prep.structure_completeness import (
    AltlocSelector,
    StructureCompletenessModule,
)

# LEU 129 as 8ZST models it: A and B complete, C side chain only. Plus a lone
# label on a neighbouring atom that has no partner, which the picker never
# offers, and a disordered water, which it never asks about either.
_PDB = """\
ATOM      1  N  AARG A 128      10.000  10.000  10.000  0.60 10.00           N
ATOM      2  CA AARG A 128      11.000  10.000  10.000  0.60 10.00           C
ATOM      3  C  AARG A 128      12.000  10.000  10.000  0.60 10.00           C
ATOM      4  O  AARG A 128      12.500  10.000  10.000  1.00 10.00           O
ATOM      5  N  ALEU A 129      13.000  10.000  10.000  0.60 10.00           N
ATOM      6  CA ALEU A 129      14.000  10.000  10.000  0.60 10.00           C
ATOM      7  C  ALEU A 129      15.000  10.000  10.000  0.60 10.00           C
ATOM      8  O  ALEU A 129      15.500  10.000  10.000  0.60 10.00           O
ATOM      9  CB ALEU A 129      14.000  11.000  10.000  0.60 10.00           C
ATOM     10  N  BLEU A 129      13.100  10.100  10.000  0.30 10.00           N
ATOM     11  CA BLEU A 129      14.100  10.100  10.000  0.30 10.00           C
ATOM     12  C  BLEU A 129      15.100  10.100  10.000  0.30 10.00           C
ATOM     13  O  BLEU A 129      15.600  10.100  10.000  0.30 10.00           O
ATOM     14  CB BLEU A 129      14.100  11.100  10.000  0.30 10.00           C
ATOM     15  CB CLEU A 129      14.200  11.200  10.000  0.10 10.00           C
HETATM   16  O  AHOH A 200      20.000  20.000  20.000  0.50 10.00           O
END
"""


def _structure(tmp_path):
    pdb = tmp_path / "in.pdb"
    pdb.write_text(_PDB)
    return PDBParser(QUIET=True).get_structure("t", str(pdb))


def _written(tmp_path, structure, selector):
    out = tmp_path / "out.pdb"
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(out), selector)
    return PDBParser(QUIET=True).get_structure("o", str(out))


def _atom_names(structure, resnum):
    return {a.get_name() for r in structure.get_residues() if r.id[1] == resnum for a in r}


def test_fill_source_picks_the_highest_occupancy_alternate_per_atom():
    atom_letters = {"N": {"A": 0.6, "B": 0.3}, "CB": {"A": 0.6, "B": 0.3, "C": 0.1}}
    assert altloc_picker.fill_source(atom_letters) == {
        "N": "A",
        "CB": "A",
    }


def test_fill_source_breaks_ties_alphabetically_so_replays_match():
    assert altloc_picker.fill_source(
        {"N": {"B": 0.5, "A": 0.5}}
    ) == {"N": "A"}


def _atom_letters(structure, resnum):
    """atom name -> {alternate: occupancy}, as the picker builds it."""
    out = {}
    for residue in structure.get_residues():
        if residue.id[1] != resnum:
            continue
        for atom in residue:
            if atom.is_disordered() and hasattr(atom, "child_dict"):
                for letter, alt in atom.child_dict.items():
                    if letter.strip():
                        out.setdefault(atom.get_name(), {})[letter.strip()] = alt.occupancy
            elif atom.altloc.strip():
                out.setdefault(atom.get_name(), {})[atom.altloc.strip()] = atom.occupancy
    return out


def test_choosing_a_partial_alternate_keeps_the_backbone(tmp_path):
    structure = _structure(tmp_path)
    plan = altloc_picker.fill_plan(
        _atom_letters(structure, 129), "C"
    )
    assert plan == {"N": "A", "CA": "A", "C": "A", "O": "A"}
    kept = _written(
        tmp_path, structure, AltlocSelector({("A", 129): "C"}, {("A", 129): plan})
    )
    assert _atom_names(kept, 129) == {"N", "CA", "C", "O", "CB"}, (
        "the backbone the chosen alternate does not model must be filled in"
    )


def test_without_a_fill_plan_a_partial_alternate_still_loses_atoms(tmp_path):
    # Documents the old behaviour, so the fill plan is not quietly dropped.
    structure = _structure(tmp_path)
    kept = _written(tmp_path, structure, AltlocSelector({("A", 129): "C"}))
    assert _atom_names(kept, 129) == {"CB"}


def test_a_complete_alternate_needs_no_fill(tmp_path):
    structure = _structure(tmp_path)
    assert altloc_picker.fill_plan(
        _atom_letters(structure, 129), "A"
    ) == {}
    kept = _written(tmp_path, structure, AltlocSelector({("A", 129): "A"}))
    assert _atom_names(kept, 129) == {"N", "CA", "C", "O", "CB"}


def test_labels_the_picker_never_asked_about_are_cleared(tmp_path):
    structure = _structure(tmp_path)
    cleared = StructureCompletenessModule._clear_altloc_labels(structure)
    assert cleared > 0
    remaining = [
        (a.get_parent().get_resname(), a.get_name())
        for a in structure.get_atoms()
        if a.altloc.strip()
    ]
    assert remaining == [], f"still labelled: {remaining}"
