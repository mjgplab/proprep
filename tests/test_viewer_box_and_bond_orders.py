"""Structure viewer: a box shape (the docking search box) and drawn bond orders.

The viewer drew only spheres, so a docking search box had no way to be
shown; it is drawn as twelve solid edges, never a translucent solid, in a
colour the page picks by background by computed contrast (the user is
visually impaired). NGL draws every bond single unless multipleBond is set,
so a ligand's double and aromatic bonds were invisible even from an SDF;
highlight(..., multiple_bond=True) asks for them, and structure_index lets
an overlay target a ligand shown as a second structure.
"""

import re
from pathlib import Path

from proprep.structure_prep import viewer_coordinator as vc_mod
from proprep.structure_prep.interactive_structure_viewer import InteractiveStructureViewer

ROOT = Path(__file__).resolve().parents[1]
TEMPLATE = (ROOT / "src/proprep/structure_prep/templates/ngl_viewer.html").read_text()


def _contrast(a, b):
    def lum(h):
        rgb = [int(h.lstrip("#")[i:i + 2], 16) / 255 for i in (0, 2, 4)]
        lin = [c / 12.92 if c <= 0.03928 else ((c + 0.055) / 1.055) ** 2.4 for c in rgb]
        return 0.2126 * lin[0] + 0.7152 * lin[1] + 0.0722 * lin[2]
    hi, lo = sorted((lum(a), lum(b)), reverse=True)
    return (hi + 0.05) / (lo + 0.05)


class _Viewer:
    def __init__(self, structures):
        self.selected_structures = structures
        self.shape_config = {}
        self.annotation_config = {}
        self.viewer_config = {}
        self.pushed = 0

    def update_annotations(self, *a, **k):
        self.pushed += 1
        return True

    def _launch_viewer(self, open_browser=True):
        pass


def _coordinator(monkeypatch, structures=("/x/prot.pdb",)):
    c = vc_mod.ViewerCoordinator()
    fake = _Viewer(list(structures))
    monkeypatch.setattr(c, "_ensure_viewer", lambda: fake)
    monkeypatch.setattr(c, "is_running", lambda: True)
    monkeypatch.setattr(c, "_owns_live_server", lambda v: True)
    monkeypatch.setattr(vc_mod, "_is_web_shell_mode", lambda: False)
    return c, fake


def _built(fake, tmp_path):
    v = InteractiveStructureViewer.__new__(InteractiveStructureViewer)
    pdb = tmp_path / "prot.pdb"; pdb.write_text("END\n")
    lig = tmp_path / "lig.sdf"; lig.write_text("\n")
    v.selected_structures = [str(pdb), str(lig)]
    v.available_annotations = {}
    v.annotation_config = fake.annotation_config
    v.viewer_config = {}
    v.shape_config = fake.shape_config
    v.trajectory_files = {}
    v.density_by_file = {}
    v.scene_override = None
    v.processor = None
    return v._build_viewer_config()


def test_box_is_drawn_as_twelve_unlit_opaque_edges():
    """Measured on rendered pixels in headless Chrome: shaded cylinders at opacity 0.5 came out at
    3.1-3.4:1 against the background; opaque wide lines keep the chosen colour."""
    assert "shape.addWideline(p1, p2, rgb, BOX_LINE_WIDTH, shapeData.label || 'box')" in TEMPLATE
    assert "} else if (shapeData.type === 'box') {" in TEMPLATE


def test_box_colours_clear_seven_to_one_on_every_background():
    light, dark = re.search(r"BOX_COLOURS = \{ light: '(#[0-9a-f]{6})', dark: '(#[0-9a-f]{6})' \}", TEMPLATE).groups()
    for bg in ("#1e1e1e", "#000000"):
        assert _contrast(light, bg) >= 7.0
    assert _contrast(dark, "#ffffff") >= 7.0
    assert "if (lastShapes.some(sh => sh.type === 'box')) applyShapesFromConfig({ shapes: lastShapes });" in TEMPLATE


def test_show_box_reaches_the_page_config(monkeypatch, tmp_path):
    c, fake = _coordinator(monkeypatch)
    c.show_box((1.0, 2.0, 3.0), (20.0, 22.0, 24.0), label="docking_box")
    assert fake.pushed == 1
    shapes = _built(fake, tmp_path)["shapes"]
    assert shapes == [{"label": "docking_box", "type": "box", "coords": [0.0, 0.0, 0.0], "radius": None,
                       "color": "auto", "opacity": 1.0, "center": [1.0, 2.0, 3.0], "size": [20.0, 22.0, 24.0]}]
    c.show_box(None, label="docking_box")
    assert "docking_box" not in fake.shape_config


def test_bond_orders_and_the_ligand_structure_can_be_asked_for(monkeypatch, tmp_path):
    assert "params.multipleBond = 'symmetric';" in TEMPLATE
    c, fake = _coordinator(monkeypatch, ("/x/prot.pdb", "/x/lig.sdf"))
    c.highlight("all", label="ligand", multiple_bond=True, structure_index=1)
    cfg = _built(fake, tmp_path)
    ligand_reps = [r for r in cfg["structures"][1]["representations"] if r["id"] == "ann_ligand"]
    assert len(ligand_reps) == 1 and ligand_reps[0]["multiple_bond"] is True
    assert not [r for r in cfg["structures"][0]["representations"] if r["id"] == "ann_ligand"]


def test_highlight_defaults_are_unchanged(monkeypatch, tmp_path):
    c, fake = _coordinator(monkeypatch)
    c.highlight(":A and 56", label="site")
    assert fake.annotation_config["site"]["structure_indices"] == [0]
    assert "multiple_bond" not in fake.annotation_config["site"]


def test_default_ligand_representations_draw_bond_orders(tmp_path):
    """A second overlay on the same atoms hid the double bonds behind the default single sticks
    (seen in headless Chrome); the defaults draw bond orders themselves. A PDB file has none, so
    only SDF/mol2 structures change."""
    v = InteractiveStructureViewer.__new__(InteractiveStructureViewer)
    pdb = tmp_path / "prot.pdb"; pdb.write_text("END\n")
    v.selected_structures = [str(pdb)]
    v.available_annotations = {}; v.annotation_config = {}; v.viewer_config = {}; v.shape_config = {}
    v.trajectory_files = {}; v.density_by_file = {}; v.scene_override = None; v.processor = None
    reps = {r["id"]: r for r in v._build_viewer_config()["structures"][0]["representations"]}
    assert reps["default_ligands"].get("multiple_bond") is True
    assert reps["default_nonstandard"].get("multiple_bond") is True
    assert "multiple_bond" not in reps["default_waters"]
