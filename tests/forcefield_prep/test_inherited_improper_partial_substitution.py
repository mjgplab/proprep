"""Regression: inherited dihedral/improper substitution must generate PARTIAL
substitutions, not only fully-substituted terms.

When a metal ligates a carboxylate, one of its two oxygens is retyped (O2 -> Y*)
while the other stays O2. The standard improper that keeps the -COO- group planar
is 2C-O2-CO-O2 (O2 in TWO positions). The metal-free prmtop carries that term with
concrete types; MCPB must re-declare it for the retyped residue.

The bug: substitution options were built from renamed types only, so product()
over 2C-O2-CO-O2 produced only the fully-substituted 2C-Y2-CO-Y5 (which never
occurs) and never the mixed 2C-Y2-CO-O2 / 2C-O2-CO-Y5 that tleap actually needs
-> "No sp2 improper torsion term for 2C-Y2-CO-O2". Fix: include the original type
alongside its renamed variants so partial substitutions are enumerated.
"""

from proprep.forcefield_prep.metal_site_parameterizer import MetalSiteWorkflowManager
from proprep.forcefield_prep.mcpb.frcmod_builder import FrcmodBuilder


class _Console:
    def print(self, *a, **k):
        pass


class _Provider:
    def __init__(self, dihedrals, impropers):
        self.dihedral_parameters = dihedrals
        self.improper_parameters = impropers


def _inst():
    inst = MetalSiteWorkflowManager.__new__(MetalSiteWorkflowManager)
    inst.console = _Console()
    return inst


def _type_assignments():
    # Two carboxylate oxygens retyped to different Y* variants (two ligating
    # carboxylates), plus their sp3/sp2 carbons left standard.
    return {
        (0.0, 0.0, 0.0): {"renamed": True, "original_type": "O2", "renamed_type": "Y2"},
        (1.0, 0.0, 0.0): {"renamed": True, "original_type": "O2", "renamed_type": "Y5"},
        (2.0, 0.0, 0.0): {"renamed": False, "original_type": "2C", "renamed_type": "2C"},
        (3.0, 0.0, 0.0): {"renamed": False, "original_type": "CO", "renamed_type": "CO"},
    }


def test_partial_improper_substitutions_are_generated():
    inst = _inst()
    fb = FrcmodBuilder()
    prov = _Provider(
        dihedrals={},
        impropers={("2C", "O2", "CO", "O2"): [(10.5, 180.0, 2)]},
    )
    inst._add_inherited_dihedral_parameters(fb, _type_assignments(), prov)

    emitted = {tuple(d["types"]) for d in fb.impropers}
    # The mixed terms tleap needs (one O2 stays original) must be present.
    assert ("2C", "Y2", "CO", "O2") in emitted, emitted
    assert ("2C", "O2", "CO", "Y5") in emitted, emitted
    # And they inherit the parent value.
    for d in fb.impropers:
        if tuple(d["types"]) == ("2C", "Y2", "CO", "O2"):
            assert d["pk"] == 10.5 and d["pn"] == 2


def test_partial_proper_substitutions_are_generated():
    inst = _inst()
    fb = FrcmodBuilder()
    # A proper dihedral with the parent type in two positions.
    prov = _Provider(
        dihedrals={("O2", "CO", "2C", "O2"): [(1.1, 180.0, 2)]},
        impropers={},
    )
    inst._add_inherited_dihedral_parameters(fb, _type_assignments(), prov)

    emitted = set()
    for d in fb.dihedrals:
        emitted.add(tuple(d["types"]))
    # A partial substitution keeping one original O2 must appear (in some
    # canonical/reversed orientation).
    def present(q):
        return q in emitted or tuple(reversed(q)) in emitted
    assert present(("Y2", "CO", "2C", "O2")) or present(("O2", "CO", "2C", "Y2")), emitted


def test_all_original_combo_is_skipped():
    inst = _inst()
    fb = FrcmodBuilder()
    prov = _Provider(
        dihedrals={},
        impropers={("2C", "O2", "CO", "O2"): [(10.5, 180.0, 2)]},
    )
    inst._add_inherited_dihedral_parameters(fb, _type_assignments(), prov)
    emitted = {tuple(d["types"]) for d in fb.impropers}
    # The unchanged parent term is not re-emitted (tleap already has it).
    assert ("2C", "O2", "CO", "O2") not in emitted
