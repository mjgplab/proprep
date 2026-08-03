"""Regression: a cached tLEaP template freezes its `mol = loadpdb <path>` line at
generation time. When the template is reused (user answers "no" to Reconfigure),
the priority structure may have advanced — most importantly, MCPB Force Field
Integration writes the renamed ``prepared_structure.pdb`` only after mcpb-4, so a
template cached earlier still points at the metal-free / un-renamed preprocessing
PDB. `_reconcile_template_loadpdb` must rewrite the stale load path (and the
saveamberparm output prefix) instead of silently building from the old structure.
"""

from proprep.tleap_prep.tleap_input_generator import TLeapInputGenerator


class _Console:
    def print(self, *a, **k):
        pass


def _gen(priority_pdb):
    g = TLeapInputGenerator.__new__(TLeapInputGenerator)

    class _P:
        console = _Console()
    g.processor = _P()
    g._select_priority_pdb_file = lambda silent=False: priority_pdb
    g._saved = None
    g._save_tleap_template = lambda t: setattr(g, "_saved", t)
    return g


_TEMPLATE = """# ProPrep-generated tLEaP Input File
source leaprc.protein.ff14SB

# === LOAD STRUCTURE ===
mol = loadpdb /work/site/remapped.pdb

check mol
saveamberparm mol remapped.prmtop remapped.rst7
quit"""


def test_rewrites_stale_loadpdb_and_prefix():
    g = _gen("/work/site/prepared_structure.pdb")
    out = g._reconcile_template_loadpdb(_TEMPLATE)

    assert "loadpdb /work/site/prepared_structure.pdb" in out
    assert "remapped.pdb" not in out
    # Output prefix tracks the new basename.
    assert "saveamberparm mol prepared_structure.prmtop prepared_structure.rst7" in out
    # Change is persisted back to the workspace.
    assert g._saved == out


def test_noop_when_already_current():
    g = _gen("/work/site/remapped.pdb")
    out = g._reconcile_template_loadpdb(_TEMPLATE)
    assert out == _TEMPLATE
    assert g._saved is None  # nothing rewritten, nothing persisted


def test_noop_when_no_priority_pdb():
    g = _gen(None)
    out = g._reconcile_template_loadpdb(_TEMPLATE)
    assert out == _TEMPLATE
    assert g._saved is None


def test_relative_vs_absolute_same_file_is_noop(tmp_path):
    # Same file expressed as a relative path must not trigger a rewrite.
    pdb = tmp_path / "prepared_structure.pdb"
    pdb.write_text("ATOM\n")
    template = _TEMPLATE.replace(
        "/work/site/remapped.pdb", str(pdb)
    ).replace("remapped.prmtop", "prepared_structure.prmtop").replace(
        "remapped.rst7", "prepared_structure.rst7"
    )
    g = _gen(str(pdb))
    out = g._reconcile_template_loadpdb(template)
    assert out == template
    assert g._saved is None
