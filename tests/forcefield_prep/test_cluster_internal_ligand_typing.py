"""
Regression test: pure-cluster internal atoms (Fe-S bridging sulfides, a Mo-S-O
core) are typed as metal ligands and get UNIQUE Y* types, not the generic
element type.

They coordinate the metal within the residue, so they never produce a
'coordinate' bond and used to fall through as "S"/"O". Sharing a standard type
can't carry per-atom bonded parameters and is inconsistent with the Y-typed
external ligands. A pure cluster residue (multi-atom, has a metal, no carbon)
has no organic part, so every non-metal atom is coordinating core.
"""

from rich.console import Console

from proprep.forcefield_prep.metal_site_parameterizer import MetalSiteWorkflowManager


def _wf():
    wf = MetalSiteWorkflowManager.__new__(MetalSiteWorkflowManager)
    wf.console = Console()
    return wf


class _Atom:
    def __init__(self, chain, resid, resname, name, element, coords):
        self.chain = chain
        self.resid = resid
        self.resname = resname
        self.atom_name = name
        self.element = element
        self.coords = coords


class _Site:
    def __init__(self, atoms):
        self.atoms = atoms


def _c(i):
    return (float(i), 0.0, 0.0)


def test_cluster_internal_atoms_flagged_metals_and_heme_excluded():
    wf = _wf()
    atoms = [
        _Atom('A', 1311, 'FES', 'FE1', 'Fe', _c(1)),
        _Atom('A', 1311, 'FES', 'FE2', 'Fe', _c(2)),
        _Atom('A', 1311, 'FES', 'S1', 'S', _c(3)),
        _Atom('A', 1311, 'FES', 'S2', 'S', _c(4)),
        _Atom('A', 1312, 'MOS', 'MO', 'Mo', _c(5)),
        _Atom('A', 1312, 'MOS', 'S', 'S', _c(6)),
        _Atom('A', 1312, 'MOS', 'O1', 'O', _c(7)),
        _Atom('A', 1312, 'MOS', 'O2', 'O', _c(8)),
        # heme: carbon-bearing organometallic → NOT a pure cluster
        _Atom('A', 200, 'HEM', 'FE', 'Fe', _c(10)),
        _Atom('A', 200, 'HEM', 'NA', 'N', _c(11)),
        _Atom('A', 200, 'HEM', 'C1A', 'C', _c(12)),
    ]
    metal_coords = {_c(1), _c(2), _c(5), _c(10)}

    flagged = wf._cluster_internal_ligand_coords(_Site(atoms), metal_coords)
    names = sorted(f"{a.resname}:{a.atom_name}" for a in atoms if a.coords in flagged)

    assert names == ['FES:S1', 'FES:S2', 'MOS:O1', 'MOS:O2', 'MOS:S']
    # metals never flagged
    assert all(a.coords not in flagged for a in atoms if a.element in ('Fe', 'Mo'))
    # heme (has carbon) fully excluded
    assert all(a.coords not in flagged for a in atoms if a.resname == 'HEM')


def test_bridging_sulfides_get_unique_Y_types():
    wf = _wf()

    def mk(center=False, lig=False):
        return {'is_center': center, 'is_metal_ligand': lig,
                'renamed_type': None, 'renamed': False}

    ta = {
        _c(1): mk(center=True), _c(2): mk(center=True),          # 2 Fe -> M1, M2
        _c(9): mk(lig=True), _c(10): mk(lig=True),               # 2 Cys SG
        _c(11): mk(lig=True), _c(12): mk(lig=True),              # 2 Cys SG
        _c(3): mk(lig=True), _c(4): mk(lig=True),                # 2 bridging S (the fix)
    }
    wf._apply_systematic_renaming(redox_site=None, type_assignments=ta)

    metals = sorted(a['renamed_type'] for a in ta.values() if a['is_center'])
    ligands = sorted(a['renamed_type'] for a in ta.values() if a['is_metal_ligand'])
    assert metals == ['M1', 'M2']
    assert ligands == ['Y1', 'Y2', 'Y3', 'Y4', 'Y5', 'Y6']
    # no generic 'S' survives; all types unique
    all_types = metals + ligands
    assert len(all_types) == len(set(all_types))
    assert 'S' not in all_types
