"""Re-snap of redox-site bonds during structure sync.

Structure repair (MODELLER) rebuilds sidechains and can relabel symmetry-
equivalent atoms — a Glu/Asp carboxylate whose OE1/OE2 (OD1/OD2) names swap is
the canonical case seen on 7K0W. A bond defined before repair stores the
coordinating atom's *position*; carrying its *name* forward would bind the metal
to the wrong (renamed) oxygen, which then lacks the metal-ligand parameters.
RedoxSiteSync._find_bonded_atom_coords resolves each bond endpoint to the atom
physically nearest its pre-repair coordinate instead of by name.
"""

import numpy as np
import pytest
from Bio.PDB.StructureBuilder import StructureBuilder
from rich.console import Console

from proprep.structure_prep.structure_completeness import (
    RedoxSiteSync,
    ResidueIdentity,
)


def _glu_structure(near_name, far_name):
    """A single Glu with one oxygen at the coordinating (near) position and one
    far, named per the arguments so we can simulate a repair relabel."""
    sb = StructureBuilder()
    sb.init_structure("s"); sb.init_model(0); sb.init_chain("A"); sb.init_seg(" ")
    sb.init_residue("GLU", " ", 104, " ")
    near = np.array((-0.22, 0.92, 7.98), dtype=float)   # coordinating O
    far = np.array((0.30, -0.10, 6.10), dtype=float)
    sb.init_atom(near_name, near, 0.0, 1.0, " ", near_name, 1, "O")
    sb.init_atom(far_name, far, 0.0, 1.0, " ", far_name, 2, "O")
    sb.init_atom("CD", np.array((0.0, 0.4, 7.0)), 0.0, 1.0, " ", "CD", 3, "C")
    return sb.get_structure()


P = (-0.22, 0.92, 7.98)   # pre-repair coordinating-O position (what the bond stored)
FAR = (0.30, -0.10, 6.10)


def test_resnap_picks_nearest_atom_through_relabel():
    # Repair renamed the coordinating oxygen OE1 -> OE2; the bond stored OE1's
    # position P. Name-based lookup would return the far atom (now 'OE1');
    # position-based must return the near atom (now 'OE2').
    sync = RedoxSiteSync(Console())
    resid = ResidueIdentity("A", 104, "GLU", " ")
    struct = _glu_structure(near_name="OE2", far_name="OE1")

    coords, name = sync._find_bonded_atom_coords(struct, resid, P, "O", "OE1")
    assert name == "OE2"
    assert tuple(round(c, 2) for c in coords) == (-0.22, 0.92, 7.98)


def test_resnap_no_relabel_keeps_name():
    # No relabel: the coordinating atom is still OE1 at P. Resolver returns OE1.
    sync = RedoxSiteSync(Console())
    resid = ResidueIdentity("A", 104, "GLU", " ")
    struct = _glu_structure(near_name="OE1", far_name="OE2")

    _, name = sync._find_bonded_atom_coords(struct, resid, P, "O", "OE1")
    assert name == "OE1"


def test_resnap_respects_element_filter():
    # The nearest atom of a different element must be ignored: anchoring on the
    # metal-facing oxygen, a nearby carbon must not be chosen.
    sync = RedoxSiteSync(Console())
    resid = ResidueIdentity("A", 104, "GLU", " ")
    struct = _glu_structure(near_name="OE2", far_name="OE1")
    # CD (carbon) sits near the anchor too, but element 'O' must exclude it.
    _, name = sync._find_bonded_atom_coords(struct, resid, P, "O", "OE1")
    assert name in ("OE1", "OE2")  # an oxygen, never CD


def test_resnap_falls_back_to_name_when_nothing_close():
    # If no same-element atom is within tolerance (a big rebuild), keep the old
    # name rather than snapping to a distant atom.
    sync = RedoxSiteSync(Console())
    resid = ResidueIdentity("A", 104, "GLU", " ")
    struct = _glu_structure(near_name="OE2", far_name="OE1")
    far_anchor = (50.0, 50.0, 50.0)   # nothing within _RESNAP_TOL
    coords, name = sync._find_bonded_atom_coords(struct, resid, far_anchor, "O", "OE1")
    assert name == "OE1"              # fell back to the name-based lookup
    assert coords == FAR              # OE1 is the far atom in this structure
