"""
Tests for the auto-launched 3D viewer at the redox-site bond-definition prompt.

The prepared structure tLEaP emits has blank chain IDs and globally renumbered
residues, so residue selections must be by bare residue number (chain-qualified
only when a chain survived). Each numbered residue is coloured with its table
row's palette index so "row [N]" maps to a colour in 3D. A viewer failure must
never propagate out of the helper.
"""

import proprep.structure_prep.viewer_coordinator as vc
from proprep.forcefield_prep.structure_preprocessor import StructurePreprocessor


class _RecordingViewer:
    def __init__(self):
        self.shown = None
        self.cleared = False
        self.highlights = []      # (selection, color, label)
        self.focused = None

    def show_structure(self, f, force=False):
        self.shown = f

    def clear_annotations(self):
        self.cleared = True

    def highlight(self, selection, *, style="ball+stick", color="element",
                  label=None, focused=False, force=False, **kw):
        self.highlights.append((selection, color, label))

    def focus_on(self, selection):
        self.focused = selection


class _Console:
    def print(self, *a, **k):
        pass


def _bare_preprocessor():
    p = object.__new__(StructurePreprocessor)
    p.console = _Console()
    p._final_pdb = None
    return p


def test_blank_chain_uses_bare_resid_nonblank_uses_chain(tmp_path, monkeypatch):
    rec = _RecordingViewer()
    monkeypatch.setattr(vc, "viewer", rec)

    pdb = tmp_path / "prepared.pdb"
    pdb.write_text("END\n")

    # (chain, resname, resid, icode) — mix of blank (post-tLEaP) and 'A'
    residue_list = [
        ("", "HID", 45, ""),
        ("", "E4Z", 182, ""),
        ("A", "MN", 185, ""),
    ]
    p = _bare_preprocessor()
    ok = p._launch_bond_definition_viewer(residue_list, str(pdb))
    assert ok is True

    assert rec.shown == str(pdb)
    assert rec.cleared is True

    sels = [h[0] for h in rec.highlights]
    assert sels == ["45", "182", "(:A and 185)"]

    # Palette colour matches 1-based table row.
    colors = [h[1] for h in rec.highlights]
    assert colors == ["palette:1", "palette:2", "palette:3"]

    # Labels carry the row number + residue for legibility.
    labels = [h[2] for h in rec.highlights]
    assert labels == ["[1] HID45", "[2] E4Z182", "[3] MN185"]

    # Focus spans the whole site.
    assert rec.focused == "45 or 182 or (:A and 185)"


def test_missing_structure_returns_false_without_launch(tmp_path, monkeypatch):
    rec = _RecordingViewer()
    monkeypatch.setattr(vc, "viewer", rec)

    p = _bare_preprocessor()
    ok = p._launch_bond_definition_viewer(
        [("", "MN", 185, "")], str(tmp_path / "does_not_exist.pdb")
    )
    assert ok is False
    assert rec.shown is None  # never tried to launch


def test_viewer_exception_is_swallowed(tmp_path, monkeypatch):
    class _Boom:
        def show_structure(self, *a, **k):
            raise RuntimeError("viewer crashed")
    monkeypatch.setattr(vc, "viewer", _Boom())

    pdb = tmp_path / "prepared.pdb"
    pdb.write_text("END\n")

    p = _bare_preprocessor()
    # Must not raise — the bond prompt has to survive a dead viewer.
    ok = p._launch_bond_definition_viewer([("", "MN", 185, "")], str(pdb))
    assert ok is False
