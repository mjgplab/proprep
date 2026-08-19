"""
The viewer can re-read a structure file that was edited in place.

``show_structure(path)`` is a no-op when the path has not changed — right for
idempotent callers, wrong once the FILE has been rewritten. Adding a hydrogen to
a cluster residue edits the structure in place, so the viewer kept rendering the
state before the edit.

``refresh_structure()`` re-serves the same path. It must never start a viewer or
open a tab: it is called after an edit regardless of whether anyone is looking.
"""

import pytest

from proprep.structure_prep.viewer_coordinator import ViewerCoordinator


class _FakeViewer:
    def __init__(self):
        self.launches = []
        self.selected_structures = []
        self.annotation_config = {}
        self.viewer_config = {}
        self.shape_config = {}

    def _launch_viewer(self, open_browser=False):
        self.launches.append(open_browser)


def _coordinator(running, viewer=None):
    c = ViewerCoordinator()
    c._viewer = viewer or _FakeViewer()
    c.is_running = lambda: running
    return c


def test_refresh_reserves_the_file_without_opening_a_tab():
    v = _FakeViewer()
    c = _coordinator(True, v)

    c.refresh_structure()

    assert v.launches == [False], "a refresh must not pop a browser tab"


def test_refresh_does_nothing_when_no_viewer_is_running():
    v = _FakeViewer()
    c = _coordinator(False, v)

    c.refresh_structure()

    assert v.launches == [], "a refresh must never start a viewer"


def test_refresh_keeps_the_structure_list_and_annotations():
    """The path is unchanged, so overlays stay valid."""
    v = _FakeViewer()
    v.selected_structures = ["/tmp/structure.pdb"]
    v.annotation_config = {"cluster_h_focus": {"selection": ":A and 3004"}}
    c = _coordinator(True, v)

    c.refresh_structure()

    assert v.selected_structures == ["/tmp/structure.pdb"]
    assert "cluster_h_focus" in v.annotation_config


def test_refresh_survives_a_broken_viewer():
    """Viewer plumbing must never bubble into the caller's flow."""
    class _Broken(_FakeViewer):
        def _launch_viewer(self, open_browser=False):
            raise RuntimeError("no display")

    c = _coordinator(True, _Broken())

    c.refresh_structure()   # must not raise


# --------------------------------------------------------------------------- #
# the caller
# --------------------------------------------------------------------------- #

from proprep.forcefield_prep.structure_preprocessor import StructurePreprocessor  # noqa: E402


class _RecordingViewer:
    def __init__(self, boom=False):
        self.calls = []
        self.boom = boom

    def _record(self, name, *a, **k):
        if self.boom:
            raise RuntimeError("viewer gone")
        self.calls.append((name, a, k))

    def show_structure(self, *a, **k):
        self._record("show_structure", *a, **k)

    def refresh_structure(self, *a, **k):
        self._record("refresh_structure", *a, **k)

    def highlight(self, *a, **k):
        self._record("highlight", *a, **k)

    def unhighlight(self, *a, **k):
        self._record("unhighlight", *a, **k)

    def focus_on(self, *a, **k):
        self._record("focus_on", *a, **k)


def _preprocessor(monkeypatch, viewer):
    sp = StructurePreprocessor.__new__(StructurePreprocessor)
    sp._pdb_file = "/tmp/structure.pdb"
    monkeypatch.setattr(
        "proprep.structure_prep.viewer_coordinator.viewer", viewer, raising=False)
    return sp


def test_prompt_focuses_the_cluster_without_refreshing(monkeypatch):
    v = _RecordingViewer()
    sp = _preprocessor(monkeypatch, v)

    sp._focus_viewer_on_cluster(":A and 3004")

    names = [c[0] for c in v.calls]
    assert "refresh_structure" not in names, "nothing has changed on disk yet"
    assert names.index("highlight") < names.index("focus_on")
    assert ("focus_on", (":A and 3004",), {}) in v.calls


def test_after_adding_a_hydrogen_the_file_is_re_read(monkeypatch):
    v = _RecordingViewer()
    sp = _preprocessor(monkeypatch, v)

    sp._focus_viewer_on_cluster(":A and 3004", refresh=True)

    names = [c[0] for c in v.calls]
    assert "refresh_structure" in names
    # Re-read before re-highlighting, so the overlay lands on the new contents.
    assert names.index("refresh_structure") < names.index("highlight")


def test_a_dead_viewer_never_interrupts_the_prompt(monkeypatch):
    sp = _preprocessor(monkeypatch, _RecordingViewer(boom=True))

    sp._focus_viewer_on_cluster(":A and 3004", refresh=True)   # must not raise
