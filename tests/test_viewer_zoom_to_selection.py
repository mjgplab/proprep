"""Zoom button on each representation row of the NGL viewer.

Every row in the viewer's sidebar, annotations from the PDB Filter or the
Protonation State Analyzer included, carries a Zoom button that frames the
row's selection. The feature lives in the browser template, so, as in
``test_viewer_dihedral_measurement.py``, the behavioural tests extract the
shipped ``zoomToRep`` and run it under Node against stand-ins for NGL and the
page. They skip when Node is unavailable.

Checked against the real page in headless Chrome when written (6R2Q, HEC A
901): the view centre landed on the selection's atom centre and the camera
distance went from 251 to 34 A; a selection matching nothing left both as
they were.
"""

import json
import re
import shutil
import subprocess
from pathlib import Path

import pytest

TEMPLATE = (
    Path(__file__).parent.parent
    / "src" / "proprep" / "structure_prep" / "templates" / "ngl_viewer.html"
)

needs_node = pytest.mark.skipif(shutil.which("node") is None, reason="Node not installed")


@pytest.fixture(scope="module")
def template_text():
    return TEMPLATE.read_text()


# ── Wiring ──────────────────────────────────────────────────────────────

def test_every_rep_row_gets_a_zoom_button(template_text):
    row = template_text[template_text.index("function buildRepRow("):
                        template_text.index("function buildStyleOptions(")]
    assert 'class="zoom-btn"' in row
    assert "zoomToRep(' + structIdx + '," in row
    # Not only on annotation rows: the button is outside any isAnnotation branch.
    assert "isAnnotation ?" not in row[row.index('class="zoom-btn"') - 200:
                                       row.index('class="zoom-btn"')]


def test_zoom_button_is_a_word_with_its_own_high_contrast_style(template_text):
    """The icon buttons are #666 on #353535; Zoom must not inherit that."""
    assert ">Zoom</button>" in template_text
    css = re.search(r"\.zoom-btn \{(.*?)\}", template_text, re.S).group(1)
    assert "color: #ffffff" in css and "background: #2d2d2d" in css


# ── Behaviour, run under Node ───────────────────────────────────────────

def _run_zoom(template_text, selection, atom_count, rep_id="ann_site"):
    fn = re.search(r"( {8}function zoomToRep\(.*?\n {8}\}\n)", template_text, re.S).group(1)
    duration = re.search(r"const ZOOM_DURATION_MS = (\d+);", template_text).group(1)
    script = f"""
        const calls = [], flagged = [];
        const ZOOM_DURATION_MS = {duration};
        const NGL = {{ Selection: function (s) {{ this.string = s; }} }};
        const repState = {{ 0: {{ ann_site: {{ config: {{ selection: {json.dumps(selection)}, visible: false }} }} }} }};
        const components = {{ 0: {{
            structure: {{ getView: (sel) => ({{ atomCount: {atom_count} }}) }},
            autoView: (sel, ms) => calls.push([sel, ms]),
        }} }};
        const document = {{ getElementById: (id) => ({{ classList: {{
            add: (c) => flagged.push(id + ':' + c), remove: () => {{}} }} }}) }};
        const setTimeout = () => {{}};
        {fn}
        zoomToRep(0, {json.dumps(rep_id)});
        console.log(JSON.stringify({{ calls, flagged }}));
    """
    out = subprocess.run(["node", "-e", script], capture_output=True, text=True, check=True)
    return json.loads(out.stdout)


@needs_node
def test_zoom_frames_the_rows_current_selection(template_text):
    """The selection as typed in the row now, hidden row or not, with the glide."""
    result = _run_zoom(template_text, "[HEC] and :A and 901", atom_count=43)
    assert result["calls"] == [["[HEC] and :A and 901", 600]]
    assert result["flagged"] == []


@needs_node
def test_selection_matching_no_atom_leaves_the_view_and_says_so(template_text):
    result = _run_zoom(template_text, ":Z and 99999", atom_count=0)
    assert result["calls"] == []
    assert sorted(result["flagged"]) == ["sel-0-ann_site:error", "zoom-0-ann_site:error"]


@needs_node
def test_unknown_row_is_ignored(template_text):
    result = _run_zoom(template_text, "protein", atom_count=10, rep_id="ann_removed")
    assert result == {"calls": [], "flagged": []}
