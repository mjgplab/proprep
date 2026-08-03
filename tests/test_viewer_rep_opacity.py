"""Translucency support for the NGL viewer's representations.

The browser exposes opacity as a per-representation slider; these tests
cover the Python half — that ``viewer_config['rep_opacity']`` reaches the
per-structure representation dicts the browser consumes, and that
``ViewerCoordinator.set_opacity`` writes that key.
"""

import pytest

from proprep.structure_prep.interactive_structure_viewer import (
    InteractiveStructureViewer,
)


def _viewer(viewer_config=None):
    v = InteractiveStructureViewer()
    v.selected_structures = ["/tmp/fake_structure.pdb"]
    v.viewer_config = viewer_config or {}
    return v


def _reps_by_id(config):
    return {r["id"]: r for r in config["structures"][0]["representations"]}


def test_default_reps_carry_no_opacity_key():
    """Untouched reps stay opaque — no key means NGL's default of 1.0."""
    reps = _reps_by_id(_viewer()._build_viewer_config())
    assert "default_protein" in reps
    for rep in reps.values():
        assert "opacity" not in rep


def test_rep_opacity_applies_to_named_rep_only():
    config = _viewer({"rep_opacity": {"default_protein": 0.35}})._build_viewer_config()
    reps = _reps_by_id(config)
    assert reps["default_protein"]["opacity"] == pytest.approx(0.35)
    assert "opacity" not in reps["default_ligands"]


def test_rep_opacity_accepts_multiple_reps():
    config = _viewer({
        "rep_opacity": {"default_protein": 0.4, "default_waters": 0.2},
    })._build_viewer_config()
    reps = _reps_by_id(config)
    assert reps["default_protein"]["opacity"] == pytest.approx(0.4)
    assert reps["default_waters"]["opacity"] == pytest.approx(0.2)


@pytest.mark.parametrize(
    "supplied,expected",
    [(-0.5, 0.0), (0.0, 0.0), (1.0, 1.0), (4.2, 1.0), ("0.6", 0.6)],
)
def test_rep_opacity_is_clamped_and_coerced(supplied, expected):
    """Out-of-range values clamp rather than reaching NGL and misrendering."""
    config = _viewer(
        {"rep_opacity": {"default_protein": supplied}}
    )._build_viewer_config()
    reps = _reps_by_id(config)
    assert reps["default_protein"]["opacity"] == pytest.approx(expected)


def test_unknown_rep_id_is_ignored():
    config = _viewer({"rep_opacity": {"no_such_rep": 0.5}})._build_viewer_config()
    for rep in _reps_by_id(config).values():
        assert "opacity" not in rep


def test_rep_opacity_survives_hidden_reps():
    """Opacity is independent of visibility — a hidden rep keeps its value."""
    config = _viewer({
        "show_waters": False,
        "rep_opacity": {"default_waters": 0.25},
    })._build_viewer_config()
    waters = _reps_by_id(config)["default_waters"]
    assert waters["visible"] is False
    assert waters["opacity"] == pytest.approx(0.25)


class _StubViewer:
    """Minimal stand-in for the viewer the coordinator drives."""

    def __init__(self):
        self.selected_structures = ["/tmp/fake_structure.pdb"]
        self.viewer_config = {}
        self.annotation_config = {}
        self.updates = 0

    def update_annotations(self, *args, **kwargs):
        self.updates += 1
        return True


def _coordinator_with_stub(stub):
    from proprep.structure_prep import viewer_coordinator as vc_mod

    coord = vc_mod.ViewerCoordinator()
    coord._viewer = stub
    coord._ensure_viewer = lambda: stub
    coord.is_running = lambda: True
    return coord


def test_set_opacity_writes_rep_opacity_and_pushes():
    stub = _StubViewer()
    _coordinator_with_stub(stub).set_opacity("protein", 0.3)
    assert stub.viewer_config["rep_opacity"] == {"default_protein": 0.3}
    assert stub.updates == 1


def test_set_opacity_accepts_full_rep_id():
    stub = _StubViewer()
    _coordinator_with_stub(stub).set_opacity("default_ligands", 0.5)
    assert stub.viewer_config["rep_opacity"] == {"default_ligands": 0.5}


def test_set_opacity_accumulates_across_calls():
    stub = _StubViewer()
    coord = _coordinator_with_stub(stub)
    coord.set_opacity("protein", 0.3)
    coord.set_opacity("waters", 0.1)
    assert stub.viewer_config["rep_opacity"] == {
        "default_protein": 0.3,
        "default_waters": 0.1,
    }


def test_set_opacity_clamps():
    stub = _StubViewer()
    coord = _coordinator_with_stub(stub)
    coord.set_opacity("protein", 2.0)
    assert stub.viewer_config["rep_opacity"]["default_protein"] == 1.0
    coord.set_opacity("protein", -1.0)
    assert stub.viewer_config["rep_opacity"]["default_protein"] == 0.0


def test_set_opacity_without_structure_is_a_noop():
    stub = _StubViewer()
    stub.selected_structures = []
    _coordinator_with_stub(stub).set_opacity("protein", 0.3)
    assert "rep_opacity" not in stub.viewer_config
    assert stub.updates == 0


def test_coordinator_output_round_trips_into_the_viewer_config():
    """What set_opacity writes is what _build_viewer_config consumes."""
    stub = _StubViewer()
    _coordinator_with_stub(stub).set_opacity("protein", 0.3)

    reps = _reps_by_id(_viewer(stub.viewer_config)._build_viewer_config())
    assert reps["default_protein"]["opacity"] == pytest.approx(0.3)
