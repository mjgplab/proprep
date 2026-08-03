"""Dihedral measurement in the NGL viewer.

The feature lives entirely in the browser template, so these tests work on
``ngl_viewer.html`` itself rather than on a Python object:

* structural checks that the mode is wired up (button, mode table, NGL
  representation, list styling), which need nothing but the file;
* numerical checks that run the shipped ``dihedralAngle`` under Node against
  geometries whose torsion is known by construction. Extracting the function
  from the template rather than restating the maths in Python is the point:
  a Python reimplementation could agree with itself while the code the user
  actually runs is wrong.

The numerical tests skip when Node is unavailable. Node is a convenience for
checking browser code from the test suite, not a ProPrep dependency.
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


@pytest.fixture(scope="module")
def template_text():
    return TEMPLATE.read_text()


# ── Wiring ──────────────────────────────────────────────────────────────

def test_dihedral_mode_button_exists(template_text):
    assert 'id="dihedral-mode-btn"' in template_text
    assert "toggleMeasurementMode('dihedral')" in template_text


def test_dihedral_declared_in_mode_table(template_text):
    """Four atoms, not three — the count drives when the pick completes."""
    table = re.search(r"dihedral:\s*\{[^}]*\}", template_text)
    assert table, "dihedral missing from MEASUREMENT_MODES"
    assert "atoms: 4" in table.group(0)


def test_mode_table_covers_every_button(template_text):
    """Every mode button must have a table entry.

    toggleMeasurementMode and clearModeButtons drive the buttons by iterating
    the table's keys, so a button without an entry is dead and an entry
    without a button throws on a null element.
    """
    buttons = set(re.findall(r'id="(\w+)-mode-btn"', template_text))
    table = re.search(
        r"const MEASUREMENT_MODES = \{(.*?)\n        \};", template_text, re.S
    )
    assert table, "MEASUREMENT_MODES table not found"
    declared = set(re.findall(r"^\s*(\w+):\s*\{", table.group(1), re.M))
    assert buttons == declared


def test_atom_count_not_hardcoded_by_ternary(template_text):
    """The old 'distance ? 2 : 3' form yields 3 for a dihedral."""
    assert "measurementMode === 'distance' ? 2 : 3" not in template_text


def test_uses_ngl_dihedral_representation(template_text):
    assert "addRepresentation('dihedral'" in template_text
    assert "atomQuad" in template_text


def test_measurement_list_styles_dihedral(template_text):
    assert ".measurement-item.dihedral" in template_text


# ── Pick marker ─────────────────────────────────────────────────────────

def test_pick_marker_uses_a_fixed_radius(template_text):
    """Not scaled off the van der Waals radius.

    ``scale`` on a spacefill multiplies the vdW radius, which made the marker
    both too large (0.68 A on a carbon, against a 1.54 A C-C bond, so markers
    on bonded atoms buried their neighbours) and element-dependent, so an H
    pick and an Fe pick looked nothing alike.
    """
    marker = re.search(
        r"addRepresentation\('spacefill',\s*\{(.*?)\}\);", template_text, re.S
    )
    assert marker, "pick marker representation not found"
    body = marker.group(1)
    assert "radiusType: 'size'" in body
    assert "radiusSize: PICK_MARKER_RADIUS" in body
    assert "scale:" not in body, "scale would reintroduce vdW-proportional sizing"


def test_pick_marker_radius_clears_a_bonded_neighbour(template_text):
    """Two markers on bonded atoms must not meet.

    The shortest bond worth measuring is around 1.0 A (X-H), so the marker
    diameter has to stay under that for the pair to read as two spheres.
    """
    m = re.search(r"const PICK_MARKER_RADIUS = ([\d.]+);", template_text)
    assert m, "PICK_MARKER_RADIUS not declared"
    radius = float(m.group(1))
    assert 0 < radius * 2 < 1.0


# ── Geometry ────────────────────────────────────────────────────────────

def _extract_dihedral_fn(text):
    start = text.index("function dihedralAngle(")
    depth, i = 0, text.index("{", start)
    for j in range(i, len(text)):
        if text[j] == "{":
            depth += 1
        elif text[j] == "}":
            depth -= 1
            if depth == 0:
                return text[start:j + 1]
    raise AssertionError("could not extract dihedralAngle from the template")


def _run_dihedral(template_text, quads):
    """Evaluate the shipped dihedralAngle on each quad of xyz triples."""
    node = shutil.which("node")
    if node is None:
        pytest.skip("node not available; skipping browser-side geometry checks")

    script = (
        _extract_dihedral_fn(template_text)
        + "\nconst pt = a => ({x: a[0], y: a[1], z: a[2]});"
        + "\nconst quads = " + json.dumps(quads) + ";"
        + "\nconsole.log(JSON.stringify(quads.map("
        + "q => dihedralAngle(pt(q[0]), pt(q[1]), pt(q[2]), pt(q[3])))));"
    )
    out = subprocess.run(
        [node, "-e", script], capture_output=True, text=True, check=True
    )
    return json.loads(out.stdout)


# Torsion about the z axis: atom1 sits at angle 0 in the xy plane and atom4 at
# angle theta, so the dihedral is theta by construction.
def _quad(theta_deg):
    import math
    t = math.radians(theta_deg)
    return [[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 1.0],
            [math.cos(t), math.sin(t), 1.0]]


@pytest.mark.parametrize("theta", [0, 60, 90, 120, 180, -60, -90, -120])
def test_known_torsions(template_text, theta):
    got = _run_dihedral(template_text, [_quad(theta)])[0]
    assert got == pytest.approx(theta, abs=1e-6)


def test_sign_distinguishes_gauche_conformers(template_text):
    """The reason for atan2: acos would collapse these onto one value."""
    plus, minus = _run_dihedral(template_text, [_quad(60), _quad(-60)])
    assert plus == pytest.approx(60.0, abs=1e-6)
    assert minus == pytest.approx(-60.0, abs=1e-6)
    assert plus != minus


def test_result_stays_within_180(template_text):
    values = _run_dihedral(template_text, [_quad(t) for t in range(-179, 180, 7)])
    assert all(-180.0 <= v <= 180.0 for v in values)


def test_degenerate_axis_does_not_produce_nan(template_text):
    """Atoms 2 and 3 coincident leave no axis to measure about."""
    quad = [[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
    got = _run_dihedral(template_text, [quad])[0]
    assert got == 0
