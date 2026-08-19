"""
The parameterization prompt selects SITES; each goes to its own parameterizer.

The prompt used to mean "combine these residues into one unit" — machinery for
a modified amino acid covalently bound to a cofactor, which is now handled by
defining that pair as a site in the Redox Site Detector. Residues are grouped
into sites there, so the prompt has no grouping left to do.

For metal sites the old routing went further: it passed combined_residues[0] to
MCPB and let the checklist re-detect and parameterize EVERY metal site in the
structure, so the selection changed nothing. Selecting one of two equivalent
Fe2S2 clusters was impossible, even though reusing the first cluster's
parameters on the second is exactly what the emitted reuse transformer is for.

Honoring the selection removes the property that made the M*/Y* atom-type
numbering safe for free — one pass over every site could always start at zero.
So the offsets are now seeded: from the workspace within a session, and from
the fingerprints of earlier runs on disk in a fresh one.
"""

import io

import pytest
from rich.console import Console

from proprep.forcefield_prep.forcefield_parameterizer import ForcefieldParameterizer
from proprep.forcefield_prep.metal_site_parameterizer import (
    mcpb_type_index, MCPB_METAL_TYPE_NAMES, MCPB_LIGAND_TYPE_NAMES,
)
from proprep.forcefield_prep.structure_preprocessor import StructurePreprocessor


class _Res:
    def __init__(self, name, chain_id, resid):
        self.name = name
        self.chain_id = chain_id
        self.resid = resid


def _unit(category, members, site_id=None, metals=None, ligands=None):
    return {"category": category, "members": members, "site_id": site_id,
            "metals": metals if metals is not None else [],
            "ligands": ligands if ligands is not None else []}


def _metal_unit(site_id, resname, resid, ligand=None):
    metal = _Res(resname, "A", resid)
    members = [metal]
    ligands = []
    if ligand:
        lig = _Res(*ligand)
        members.append(lig)
        ligands.append(lig)
    return _unit("metal_site", members, site_id, [metal], ligands)


# The 4UHX layout: two Fe2S2 clusters and a Mo cofactor with an MTE ligand.
def _fes_a():
    return _metal_unit("site_1", "FES", 3001)


def _moco():
    return _metal_unit("site_2", "MOS", 3004, ligand=("MTE", "A", 3003))


def _fes_b():
    return _metal_unit("site_3", "FES", 3002)


def _fad():
    return _unit("small_molecule", [_Res("FAD", "A", 3006)])


def _param(monkeypatch=None):
    p = ForcefieldParameterizer.__new__(ForcefieldParameterizer)
    p.console = Console(file=io.StringIO(), width=200)
    p.processor = None
    p._workspace_writes = {}
    p.update_workspace = lambda k, v: p._workspace_writes.__setitem__(k, v)
    return p


# --------------------------------------------------------------------------- #
# routing
# --------------------------------------------------------------------------- #

def test_each_selected_site_goes_to_its_own_parameterizer():
    p = _param()
    metal_calls, other_calls = [], []
    p._parameterize_metal_sites = lambda units, all_units=None: metal_calls.append(units)
    p._route_unit = lambda unit: other_calls.append(unit)

    selected = [_fes_a(), _fad()]
    p._parameterize_selected_units(selected, selected)

    assert len(metal_calls) == 1
    assert [u["site_id"] for u in metal_calls[0]] == ["site_1"]
    assert [u["category"] for u in other_calls] == ["small_molecule"]


def test_several_metal_sites_go_in_one_mcpb_pass():
    """MCPB numbers M*/Y* across sites, so selected sites share one run."""
    p = _param()
    metal_calls = []
    p._parameterize_metal_sites = lambda units, all_units=None: metal_calls.append(units)
    p._route_unit = lambda unit: None

    selected = [_fes_a(), _moco(), _fes_b()]
    p._parameterize_selected_units(selected, selected)

    assert len(metal_calls) == 1, "metal sites must not be run one checklist each"
    assert [u["site_id"] for u in metal_calls[0]] == ["site_1", "site_2", "site_3"]


def test_only_selected_metal_sites_are_recorded():
    """The two-Fe2S2 case: parameterize one, reuse it on the other."""
    p = _param()
    p._parameterize_metal_site = lambda res, announce=True: None

    all_units = [_fes_a(), _moco(), _fes_b()]
    p._parameterize_metal_sites([all_units[1]], all_units)

    assert p._workspace_writes["mcpb_selected_site_ids"] == ["site_2"]


def test_unselected_metal_sites_are_named():
    p = _param()
    p._parameterize_metal_site = lambda res, announce=True: None

    all_units = [_fes_a(), _moco(), _fes_b()]
    p._parameterize_metal_sites([all_units[0]], all_units)

    out = p.console.file.getvalue()
    assert "1 site selected" in out
    assert "FES (A:3002)" in out, "the unselected cluster should be named"
    assert "MOS (A:3004)" in out
    assert "reuse transformer" in out


def test_one_failing_site_does_not_cancel_the_rest():
    p = _param()
    p._parameterize_metal_sites = lambda units, all_units=None: None
    routed = []

    def flaky(unit):
        routed.append(unit["members"][0].name)
        if unit["members"][0].name == "FAD":
            raise RuntimeError("antechamber exploded")

    p._route_unit = flaky
    selected = [_fad(), _unit("small_molecule", [_Res("MLI", "A", 3010)])]
    p._parameterize_selected_units(selected, selected)

    assert routed == ["FAD", "MLI"]
    assert "failed" in p.console.file.getvalue()


def test_combine_machinery_is_gone():
    """Grouping belongs to the Redox Site Detector, not this prompt."""
    for gone in ("_handle_combined_residue_selection", "_ask_combine_or_separate",
                 "_resolve_category_conflict"):
        assert not hasattr(ForcefieldParameterizer, gone), f"{gone} still present"


# --------------------------------------------------------------------------- #
# atom-type numbering
# --------------------------------------------------------------------------- #

def test_type_index_is_a_position_not_a_digit():
    """Y9 -> 8 and YA -> 9: the label stops matching the count past nine."""
    assert mcpb_type_index("Y1") == 0
    assert mcpb_type_index("Y9") == 8
    assert mcpb_type_index("YA") == 9
    assert mcpb_type_index("YB") == 10
    assert mcpb_type_index("M1") == 0
    assert mcpb_type_index("M3") == 2
    assert mcpb_type_index("SG") is None
    assert mcpb_type_index("") is None


def _preprocessor(tmp_path, workspace=None, interactive=False):
    sp = StructurePreprocessor.__new__(StructurePreprocessor)
    sp._output_dir = tmp_path / "metal_site_params_FES_A_3001"
    sp._output_dir.mkdir(parents=True, exist_ok=True)
    sp.console = Console(file=io.StringIO(), width=200)
    sp.workspace = workspace
    sp.processor = None
    sp._interactive = interactive
    return sp


def _write_fingerprint(site_dir, mapped):
    models = site_dir / "models"
    models.mkdir(parents=True, exist_ok=True)
    lines = [f"{i}-RES-AT   {i}  X  -> {name}" for i, name in enumerate(mapped, 1)]
    (models / "standard.fingerprint").write_text("\n".join(lines) + "\n")


def test_offsets_recovered_from_a_previous_run_on_disk(tmp_path):
    """The fresh-session case: nothing in the workspace, everything on disk."""
    sp = _preprocessor(tmp_path)
    _write_fingerprint(sp._output_dir / "site_1",
                       ["M1", "M2", "Y1", "Y2", "Y3", "Y4", "Y5", "Y6"])
    _write_fingerprint(sp._output_dir / "site_2",
                       ["M3", "Y7", "Y8", "Y9", "YA", "YB"])

    metal_next, ligand_next = sp._seed_mcpb_type_offsets()

    # 3 metals used (M1-M3) and 11 ligating atoms (Y1-YB) -> next M4 / YC.
    assert (metal_next, ligand_next) == (3, 11)
    assert MCPB_METAL_TYPE_NAMES[metal_next] == "M4"
    assert MCPB_LIGAND_TYPE_NAMES[ligand_next] == "YC"


def test_prior_run_in_a_sibling_directory_is_found(tmp_path):
    """A second site's run gets its own metal_site_params_* directory."""
    sibling = tmp_path / "metal_site_params_MOS_A_3004"
    _write_fingerprint(sibling / "site_2", ["M3", "Y7", "Y8"])
    sp = _preprocessor(tmp_path)

    metal_next, ligand_next = sp._seed_mcpb_type_offsets()

    assert (metal_next, ligand_next) == (3, 8)


class _WS:
    def __init__(self, d=None):
        self.d = dict(d or {})

    def get(self, k, default=None):
        return self.d.get(k, default)

    def set(self, k, v):
        self.d[k] = v


def test_workspace_and_disk_take_the_higher_mark(tmp_path):
    """Same-session state and on-disk history must not undercut each other."""
    ws = _WS({"mcpb_type_offsets": {"metal": 5, "ligand": 2}})
    sp = _preprocessor(tmp_path, workspace=ws)
    _write_fingerprint(sp._output_dir / "site_1", ["M1", "Y1", "Y2", "Y3", "Y4"])

    metal_next, ligand_next = sp._seed_mcpb_type_offsets()

    assert metal_next == 5     # workspace ahead of disk
    assert ligand_next == 4    # disk ahead of workspace


def test_no_history_starts_at_zero(tmp_path):
    sp = _preprocessor(tmp_path, workspace=_WS())

    assert sp._seed_mcpb_type_offsets() == (0, 0)


def test_rerunning_a_site_can_reuse_its_own_names(tmp_path):
    """Re-deriving a site should replace its entry, not strand its old types."""
    sp = _preprocessor(tmp_path, interactive=False)
    site_dir = sp._output_dir / "site_2"
    _write_fingerprint(site_dir, ["M3", "Y7", "Y8", "Y9", "YA", "YB"])

    # Non-interactive declines the reuse (it cannot ask), but the prior naming
    # is still reported so the run is not silent about it.
    assert sp._offer_prior_type_reuse(site_dir, "site_2") is None
    out = sp.console.file.getvalue()
    assert "parameterized before" in out
    assert "M3" in out and "YB" in out


def test_reuse_returns_the_sites_own_starting_positions(tmp_path, monkeypatch):
    sp = _preprocessor(tmp_path, interactive=True)
    site_dir = sp._output_dir / "site_2"
    _write_fingerprint(site_dir, ["M3", "Y7", "Y8", "Y9", "YA", "YB"])

    monkeypatch.setattr(
        "proprep.forcefield_prep.structure_preprocessor.confirm_with_context",
        lambda processor, prompt, **kw: True)

    assert sp._offer_prior_type_reuse(site_dir, "site_2") == (2, 6)  # M3, Y7


def test_no_prior_output_means_no_reuse_offer(tmp_path):
    sp = _preprocessor(tmp_path, interactive=True)

    assert sp._offer_prior_type_reuse(sp._output_dir / "site_9", "site_9") is None
    assert sp.console.file.getvalue() == ""
