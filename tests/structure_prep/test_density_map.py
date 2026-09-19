"""A structure against its electron density.

The geometry is tested on a synthetic crystal: a skewed (triclinic) cell, 80 atoms lying
OUTSIDE the origin cell as deposited models do, and a periodic map made of a Gaussian on each
atom, written as PDBe writes its maps (one whole unit cell from the origin).
`PROPREP_NETWORK_TESTS=1` adds 1LTZ from PDBe, whose strongest peak must be its iron."""

import io
import json
import os
import struct
import types
import urllib.error
import urllib.request

import numpy as np
import pytest
import requests
from rich.console import Console

from proprep.structure_prep import density_map as dm

CELL = (40.0, 36.0, 30.0, 80.0, 105.0, 70.0)
SAMPLING = (64, 56, 48)
# 80 atoms in a blob centred in a NEIGHBOURING cell, so none lies inside the map as fetched.
# Enough atoms that the mean density over them means something: over a handful, one atom
# landing on a peak decides the average.
ATOMS_FRACTIONAL = np.array([-0.45, 0.35, 1.40]) + np.random.default_rng(7).uniform(-0.22, 0.22, size=(80, 3))


def _orthogonalization():
    return np.linalg.inv(dm.fractionalization_matrix(CELL))


def _atoms_cartesian():
    return ATOMS_FRACTIONAL @ _orthogonalization().T


def _whole_cell_map():
    """Periodic density: a Gaussian (width 0.7 A) on every atom and on all its images."""
    axes = [np.arange(n) / n for n in SAMPLING]
    fx, fy, fz = np.meshgrid(*axes, indexing="ij")
    grid_fractional = np.stack([fx, fy, fz], axis=-1)
    orthogonalization = _orthogonalization()
    grid = np.zeros(SAMPLING, dtype=np.float32)
    for atom in ATOMS_FRACTIONAL:
        delta = grid_fractional - atom
        delta -= np.rint(delta)                                  # nearest periodic image
        distance_squared = ((delta @ orthogonalization.T) ** 2).sum(axis=-1)
        grid += np.exp(-distance_squared / (2 * 0.7 ** 2)).astype(np.float32)
    grid -= grid.mean()
    return dm.DensityMap(grid=grid, cell=CELL, sampling=SAMPLING, start=(0, 0, 0),
                         mean=float(grid.mean()), sigma=float(grid.std()))


def _moved(coordinates):
    c, s = np.cos(np.pi / 6), np.sin(np.pi / 6)
    return coordinates @ np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]]).T + np.array([7.0, -3.0, 11.0])


def _write_pdb(path, coordinates, header=True, rename=False):
    scale = dm.fractionalization_matrix(CELL)
    with open(path, "w") as f:
        if header:
            f.write(f"HEADER    OXIDOREDUCTASE                          01-JAN-20   9XYZ              \n")
            f.write("CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1\n" % CELL)
            for i in range(3):
                f.write("SCALE%d    %10.6f%10.6f%10.6f     %10.5f\n" % (i + 1, *scale[i], 0.0))
        for i, (x, y, z) in enumerate(coordinates, 1):
            resname, atom = ("HIE", " CB ") if rename else ("HIS", " CA ")
            f.write(f"ATOM  {i:5d} {atom} {resname} A{i:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           C  \n")
        f.write("END\n")
    return str(path)


# ---------------------------------------------------------------------------
# CCP4 format
# ---------------------------------------------------------------------------

def test_a_map_written_and_read_back_is_the_same(tmp_path):
    density = _whole_cell_map()
    density.start = (-7, 3, 11)
    path = str(tmp_path / "m.ccp4")
    dm.write_ccp4(density, path)
    back = dm.read_ccp4(path)
    assert np.array_equal(back.grid, density.grid)
    assert back.start == (-7, 3, 11) and back.sampling == SAMPLING
    assert np.allclose(back.cell, CELL)


def test_a_map_stored_with_another_axis_order_reads_to_the_same_grid(tmp_path):
    """mapc/mapr/maps = 2,1,3: the fast file axis follows b, the next follows a."""
    density = _whole_cell_map()
    path = str(tmp_path / "m.ccp4")
    dm.write_ccp4(density, path)
    raw = bytearray(open(path, "rb").read())
    nx, ny, nz = SAMPLING
    struct.pack_into("<3i", raw, 0, ny, nx, nz)                   # columns follow b, rows follow a
    struct.pack_into("<3i", raw, 16, 5, 4, 6)                     # file-axis starts: b, a, c
    struct.pack_into("<3i", raw, 64, 2, 1, 3)
    raw[1024:] = np.ascontiguousarray(density.grid.transpose(2, 0, 1), dtype="<f4").tobytes()   # [c, a, b]
    open(path, "wb").write(bytes(raw))

    back = dm.read_ccp4(path)
    assert np.array_equal(back.grid, density.grid)
    assert back.start == (4, 5, 6)                                # as a, b, c


def test_what_is_not_a_map_is_refused(tmp_path):
    path = tmp_path / "x.ccp4"
    path.write_bytes(b"<html>not found</html>" * 100)
    with pytest.raises(dm.DensityMapError):
        dm.read_ccp4(str(path))


# ---------------------------------------------------------------------------
# Re-cutting: the map is one cell at the origin, the model is elsewhere
# ---------------------------------------------------------------------------

def test_every_atom_lies_outside_the_map_as_fetched():
    assert ((ATOMS_FRACTIONAL < 0) | (ATOMS_FRACTIONAL >= 1)).any(axis=1).all()


def test_the_recut_box_holds_the_atoms_and_copies_the_periodic_density_exactly():
    density = _whole_cell_map()
    atoms = _atoms_cartesian()
    region = dm.recut_around(density, atoms, margin=3.0)

    low = np.floor(ATOMS_FRACTIONAL.min(axis=0) * SAMPLING)
    high = np.ceil(ATOMS_FRACTIONAL.max(axis=0) * SAMPLING)
    assert (np.array(region.start) < low).all()
    assert (np.array(region.start) + region.grid.shape - 1 > high).all()
    assert min(region.start) < 0                                   # it reaches into neighbouring cells

    rng = np.random.default_rng(1)
    for _ in range(200):
        i, j, k = (rng.integers(0, n) for n in region.grid.shape)
        source = [(index + start) % n for index, start, n in zip((i, j, k), region.start, SAMPLING)]
        assert region.grid[i, j, k] == density.grid[tuple(source)]
    assert (region.mean, region.sigma) == (density.mean, density.sigma)   # whole-cell statistics travel


def test_in_a_skewed_cell_the_density_is_on_the_atoms(tmp_path):
    density = _whole_cell_map()
    atoms = _atoms_cartesian()
    region = dm.recut_around(density, atoms, margin=3.0)
    path = str(tmp_path / "cut.ccp4")
    dm.write_ccp4(region, path)
    back = dm.read_ccp4(path)
    back.mean, back.sigma = density.mean, density.sigma

    on_atoms = dm.density_at(back, atoms)
    assert on_atoms.min() > 5.0                                    # sigma: each atom sits on its peak
    assert np.allclose(on_atoms, dm.density_at(density, atoms))    # the cut agrees with the periodic wrap
    elsewhere = _moved(atoms)
    assert abs(dm.density_at(density, elsewhere).mean()) < 0.5     # sampled through the periodic cell
    with pytest.raises(dm.DensityMapError, match="outside the map"):
        dm.density_at(back, elsewhere)                             # the cut box has edges

    # where NGL will put the strongest voxel (its getMatrix: cell basis / sampling, applied to
    # index + NXSTART..) is an atom, to within the grid spacing
    peak = np.array(np.unravel_index(back.grid.argmax(), back.grid.shape)) + np.array(back.start)
    position = (peak / np.array(SAMPLING)) @ _orthogonalization().T
    assert np.linalg.norm(atoms - position, axis=1).min() < density.spacing


def test_a_partial_map_cannot_be_extended():
    density = _whole_cell_map()
    density.grid = density.grid[:10]
    with pytest.raises(dm.DensityMapError, match="whole unit cell"):
        dm.recut_around(density, _atoms_cartesian(), margin=2.0)


def test_an_entry_in_a_non_standard_orientation_is_refused():
    density = _whole_cell_map()
    atoms = _atoms_cartesian()
    standard = dm.fractionalization_matrix(CELL)
    dm.check_standard_frame(density, atoms, standard, np.zeros(3))            # standard: fine
    dm.check_standard_frame(density, atoms, None, None)
    turn = np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    with pytest.raises(dm.DensityMapError, match="standard crystal orientation"):
        dm.check_standard_frame(density, atoms, standard @ turn, np.zeros(3))


# ---------------------------------------------------------------------------
# Is the structure still where it was deposited?
# ---------------------------------------------------------------------------

def test_frame_check_is_by_position_not_by_name(tmp_path):
    atoms = _atoms_cartesian()
    deposited = _write_pdb(tmp_path / "dep.pdb", atoms)
    renamed = _write_pdb(tmp_path / "renamed.pdb", atoms[:60], header=False, rename=True)
    moved = _write_pdb(tmp_path / "moved.pdb", _moved(atoms))

    same = dm.frame_check(renamed, deposited, grid_spacing=0.6)
    assert same.same_frame and same.median_distance == 0.0 and same.coincident_fraction == 1.0
    assert not dm.frame_check(moved, deposited, grid_spacing=0.6).same_frame


def test_the_map_itself_catches_a_moved_structure_that_kept_its_header(tmp_path):
    """Compared with itself a file always agrees; its atoms in the density do not."""
    density = _whole_cell_map()
    map_path = str(tmp_path / "9xyz_2fofc.ccp4")
    dm.write_ccp4(density, map_path)
    atoms = _atoms_cartesian()
    deposited = _write_pdb(tmp_path / "dep.pdb", atoms)
    moved = _write_pdb(tmp_path / "moved.pdb", _moved(atoms))       # CRYST1, SCALE and HEADER intact

    assert dm.frame_check(moved, moved, grid_spacing=0.6).same_frame               # the hole
    assert dm.mean_density_at_atoms(map_path, deposited, deposited) > dm.IN_DENSITY_SIGMA
    assert dm.mean_density_at_atoms(map_path, moved, moved) < dm.IN_DENSITY_SIGMA  # closed


def test_pdb_frame_records_are_read_and_an_nmr_cell_is_recognisable(tmp_path):
    crystal = _write_pdb(tmp_path / "x.pdb", _atoms_cartesian())
    cell, scale, shift, idcode = dm.read_pdb_frame(crystal)
    assert idcode == "9XYZ" and np.allclose(cell, CELL)
    assert np.allclose(scale, dm.fractionalization_matrix(CELL), atol=1e-6) and np.allclose(shift, 0)

    nmr = tmp_path / "nmr.pdb"
    nmr.write_text("HEADER    DE NOVO PROTEIN                         01-JAN-20   2RVD              \n"
                   "CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1           1\n"
                   "SCALE1      1.000000  0.000000  0.000000        0.00000\n"
                   "SCALE2      0.000000  1.000000  0.000000        0.00000\n"
                   "SCALE3      0.000000  0.000000  1.000000        0.00000\n"
                   "ATOM      1  CA  GLY A   1       0.000   0.000   0.000  1.00  0.00           C  \n")
    cell, scale, _, idcode = dm.read_pdb_frame(str(nmr))
    assert idcode == "2RVD" and min(cell[:3]) <= 1.0 and scale is None


# ---------------------------------------------------------------------------
# Fetching
# ---------------------------------------------------------------------------

class _Response:
    def __init__(self, status, content=b""):
        self.status_code, self.content = status, content


def test_fetch_no_maps_unreachable_and_cached(tmp_path, monkeypatch):
    monkeypatch.setattr(dm.requests, "get", lambda url, timeout=None: _Response(404))
    assert dm.fetch_density_maps("2RVD", str(tmp_path)) is None

    def down(url, timeout=None):
        raise requests.exceptions.ConnectionError("no route to host")
    monkeypatch.setattr(dm.requests, "get", down)
    with pytest.raises(dm.DensityMapError, match="could not be reached"):
        dm.fetch_density_maps("1LTZ", str(tmp_path))

    asked = []
    monkeypatch.setattr(dm.requests, "get", lambda url, timeout=None: asked.append(url) or _Response(200, b"x" * 4096))
    paths = dm.fetch_density_maps("1LTZ", str(tmp_path))
    assert asked == ["https://www.ebi.ac.uk/pdbe/entry-files/1ltz.ccp4", "https://www.ebi.ac.uk/pdbe/entry-files/1ltz_diff.ccp4"]
    assert dm.fetch_density_maps("1LTZ", str(tmp_path)) == paths and len(asked) == 2      # second time from disk


# ---------------------------------------------------------------------------
# Viewer: which structures get density, and what the page is told
# ---------------------------------------------------------------------------

@pytest.fixture
def project(tmp_path, monkeypatch):
    """A project with the synthetic entry 9XYZ 'downloaded' and its maps 'on PDBe'."""
    density = _whole_cell_map()
    atoms = _atoms_cartesian()
    files = {
        "deposited": _write_pdb(tmp_path / "9XYZ.pdb", atoms),
        "filtered": _write_pdb(tmp_path / "filtered.pdb", atoms[:60], header=False, rename=True),
        "moved": _write_pdb(tmp_path / "aligned.pdb", _moved(atoms)),
    }

    def fake_fetch(pdb_id, out_dir, timeout=120):
        if pdb_id != "9XYZ":
            return None
        paths = {}
        for kind in dm.MAP_KINDS:
            paths[kind] = os.path.join(out_dir, f"9xyz_{kind}.ccp4")
            dm.write_ccp4(density, paths[kind])
        return paths
    monkeypatch.setattr(dm, "fetch_density_maps", fake_fetch)

    from proprep.structure_prep.interactive_structure_viewer import InteractiveStructureViewer
    console = Console(width=300, record=True, file=io.StringIO(), force_terminal=False)
    workspace = {"rcsb_pdb_file": files["deposited"], "rcsb_pdb_files": [files["deposited"]]}
    viewer = InteractiveStructureViewer()
    viewer.processor = types.SimpleNamespace(console=console, _get_workspace=lambda: workspace)
    return viewer, files, console, density


def test_density_goes_to_structures_still_in_the_deposited_frame(project):
    viewer, files, console, density = project
    assert viewer._attach_density(files["deposited"], 3.0)
    assert viewer._attach_density(files["filtered"], 3.0)          # renamed, no header: found by position
    assert not viewer._attach_density(files["moved"], 3.0)         # kept its header: refused by the map
    out = " ".join(console.export_text().split())                  # the console wraps long lines
    assert "whatever its header says" in out and "The level required is 1.0 sigma" in out


def test_an_nmr_structure_is_told_why(project, tmp_path):
    viewer, _, console, _ = project
    nmr = tmp_path / "2RVD.pdb"
    nmr.write_text("HEADER    DE NOVO PROTEIN                         01-JAN-20   2RVD              \n"
                   "CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1           1\n"
                   "ATOM      1  CA  GLY A   1       0.000   0.000   0.000  1.00  0.00           C  \n")
    assert not viewer._attach_density(str(nmr), 3.0)
    out = " ".join(console.export_text().split())
    assert "2RVD has no unit cell" in out and "9XYZ" not in out     # not blamed on an unrelated entry


def test_the_page_is_given_whole_cell_statistics_and_a_rewritten_file_loses_its_maps(project):
    viewer, files, _, density = project
    viewer._attach_density(files["deposited"], 3.0)
    viewer._attach_density(files["filtered"], 3.0)
    viewer.selected_structures = [files["moved"], files["filtered"], files["deposited"]]

    by_index = viewer._density_for_index()
    assert sorted(by_index) == [1, 2]                              # by file, whatever the viewer's order
    config = viewer._build_viewer_config()
    entry = next(s for s in config["structures"] if s["index"] == 2)["density"]
    assert entry["pdb_id"] == "9XYZ"
    assert [m["url"] for m in entry["maps"]] == ["/density/2/2fofc", "/density/2/fofc"]
    assert entry["maps"][0]["sigma"] == pytest.approx(density.sigma)      # the cell's, not the box's
    recut = dm.read_ccp4(by_index[2]["maps"][0]["path"])
    assert recut.sigma != pytest.approx(density.sigma)
    assert "density" not in next(s for s in config["structures"] if s["index"] == 0)

    with open(files["filtered"], "a") as f:                        # rewritten in place after the maps were fitted
        f.write("REMARK moved\n")
    assert sorted(viewer._density_for_index()) == [2]


def test_the_server_serves_the_maps_by_structure_index(project):
    from proprep.structure_prep.viewer_server import ViewerServer
    viewer, files, _, _ = project
    viewer._attach_density(files["deposited"], 3.0)
    viewer.selected_structures = [files["deposited"]]
    record = viewer._density_for_index()[0]
    server = ViewerServer(config={}, structure_files=[files["deposited"]], port=8790,
                          density_files={0: {m["kind"]: m["path"] for m in record["maps"]}})
    assert server.start(open_browser=False)
    try:
        base = server.get_url().rsplit("/", 1)[0]
        body = urllib.request.urlopen(f"{base}/density/0/2fofc").read()
        assert body[208:212] == b"MAP " and len(body) == os.path.getsize(record["maps"][0]["path"])
        for missing in ("/density/0/nope", "/density/5/2fofc"):
            with pytest.raises(urllib.error.HTTPError) as error:
                urllib.request.urlopen(base + missing)
            assert error.value.code == 404
    finally:
        server.stop()


def test_the_page_contours_at_absolute_levels():
    """NGL's own 'sigma' mode would take the statistics of the re-cut box, not of the cell."""
    from proprep.utils.paths import get_package_dir
    page = (get_package_dir() / "structure_prep" / "templates" / "ngl_viewer.html").read_text()
    block = page[page.index("async function attachDensity"):page.index("function setDensityVisible")]
    assert "isolevelType: 'value'" in block and "isolevelType: 'sigma'" not in block
    assert block.count("d.mean + d.sign * n * d.sigma") == 1 and "map.mean + sign *" in block


def test_the_page_shows_the_whole_map_first_and_a_box_belongs_to_the_model():
    """Seen in the browser by mgp: a 24 A cube mid-protein that stayed mid-screen on panning. NGL pins
    a volume's box to the centre of the VIEW (setBox on every frame); ProPrep anchors it on the model."""
    from proprep.utils.paths import get_package_dir
    page = (get_package_dir() / "structure_prep" / "templates" / "ngl_viewer.html").read_text()
    assert 'id="density-box" min="0" max="60" step="1" value="0"' in page          # whole map to begin with
    assert '<span id="density-box-label" style="min-width:64px">whole map</span>' in page
    assert "r.viewer.signals.ticked.remove(r.setBox, r)" in page                    # the box leaves the view centre
    assert "setParameters({ boxCenter: new NGL.Vector3(densityCentre.x" in page     # and sits on a point of the model
    assert "r.viewer.signals.ticked.add(r.setBox, r)" in page                       # 'follows the view' puts it back
    assert "typeof r.setBox !== 'function'" in page                                 # an NGL without it: said, not assumed
    # a selection NGL cannot read matches every atom; that is refused, not centred on
    assert "n === structure.atomCount" in page


# ---------------------------------------------------------------------------
# Highlighting the mesh near chosen atoms (mgp is visually impaired: contrast is a requirement)
# ---------------------------------------------------------------------------

def _page():
    from proprep.utils.paths import get_package_dir
    return (get_package_dir() / "structure_prep" / "templates" / "ngl_viewer.html").read_text()


def test_the_pages_highlight_code_agrees_with_brute_force(tmp_path):
    """Runs the functions the page ships, cut out of it, under Node. They are plain JavaScript:
    no NGL, no network."""
    import shutil
    import subprocess
    node = shutil.which("node")
    if not node:
        pytest.skip("node is not installed")
    page = _page()
    core = page[page.index("// <density-highlight-core>"):page.index("// </density-highlight-core>")]
    script = tmp_path / "highlight.cjs"
    script.write_text(core + r"""
let seed = 12345; function rnd() { seed = (seed * 1103515245 + 12345) % 2147483648; return seed / 2147483648; }
const atoms = [], vertices = [];
for (let i = 0; i < 40; i++) atoms.push({ x: rnd() * 30 - 15, y: rnd() * 30 - 15, z: rnd() * 30 - 15 });   // negatives too
for (let i = 0; i < 20000; i++) vertices.push([rnd() * 40 - 20, rnd() * 40 - 20, rnd() * 40 - 20]);
const out = { wrong: {}, near: {} };
for (const radius of [1.0, 2.0, 3.5]) {
  const grid = buildDensityAtomGrid(atoms, radius); let wrong = 0, near = 0;
  for (const v of vertices) {
    let brute = false;
    for (const a of atoms) { const d = (a.x-v[0])**2 + (a.y-v[1])**2 + (a.z-v[2])**2; if (d <= radius * radius) { brute = true; break; } }
    if (nearDensityAtom(grid, radius, v[0], v[1], v[2]) !== brute) wrong++;
    if (brute) near++;
  }
  out.wrong[radius] = wrong; out.near[radius] = near;
}
const grid = buildDensityAtomGrid([{ x: 0, y: 0, z: 0 }], 2);
const state = { on: true, grid: grid, radius: 2, highlight: 0xffe600, rest: 0.4 };
out.colours = {
  nearMesh: densityVertexColour(state, 0x6aa9ff, true, 1, 0, 0),
  farMesh: densityVertexColour(state, 0x6aa9ff, true, 9, 0, 0),
  nearDifference: densityVertexColour(state, 0x5fd35f, false, 1, 0, 0),
  farUndimmed: densityVertexColour(Object.assign({}, state, { rest: 1 }), 0x6aa9ff, true, 9, 0, 0),
  off: densityVertexColour(Object.assign({}, state, { on: false }), 0x6aa9ff, true, 1, 0, 0),
  nothingChosen: densityVertexColour(Object.assign({}, state, { grid: new Map() }), 0x6aa9ff, true, 9, 0, 0),
};
out.contrast = [1, 0.4].map(f => +densityContrast(0xffe600, densityDimmed(0x6aa9ff, f)).toFixed(1));
out.blackOnWhite = +densityContrast(0x000000, 0xffffff).toFixed(1);
console.log(JSON.stringify(out));
""")
    result = subprocess.run([node, str(script)], capture_output=True, text=True, timeout=120)
    assert result.returncode == 0, result.stderr
    out = json.loads(result.stdout)

    assert out["wrong"] == {"1": 0, "2": 0, "3.5": 0}             # the spatial grid never disagrees
    assert 0 < out["near"]["1"] < out["near"]["2"] < out["near"]["3.5"]
    colours = out["colours"]
    assert colours["nearMesh"] == 0xffe600                         # the highlight
    assert colours["farMesh"] == 0x2a4466                          # the rest at 40%: same hue, darker
    assert colours["nearDifference"] == 0x5fd35f                   # green and red keep their meaning
    assert colours["farUndimmed"] == colours["off"] == 0x6aa9ff
    assert colours["nothingChosen"] == 0x6aa9ff                    # nothing chosen: nothing dimmed
    # WCAG ratios. Yellow on the undimmed mesh is below the 3:1 asked of graphics, which is why the
    # rest is dimmed and why 40% is the default: the brightest setting that reaches 7:1
    assert out["blackOnWhite"] == 21.0
    assert out["contrast"][0] < 3.0 and out["contrast"][1] >= 7.0


def test_the_highlight_is_one_mesh_coloured_by_position_with_the_contrast_shown():
    page = _page()
    assert "this.positionColor = function (v)" in page and "colorScheme: schemeId" in page
    assert 'id="density-rest" min="10" max="100" step="5" value="40"' in page
    assert 'id="density-colour-hl" value="#ffe600"' in page and 'id="density-colour-mesh" value="#6aa9ff"' in page
    assert "Contrast, highlight to the rest: ' + contrast.toFixed(1)" in page     # a number, not only a look
    # not by pointer precision alone: the ligand list and the selection field set the highlight too
    assert "centreDensityAt(t, t.label, 0); setDensityChosen(t.atoms);" in page
    assert "setDensityChosen(matched);" in page and 'value="residue"> whole residues' in page
    # clicks mean something else while measuring, and other structures are not in this map's frame
    assert "if (measurementMode || !pickingProxy || densityStructIdx === null) { return; }" in page
    assert "pickingProxy.component !== components[densityStructIdx]" in page


# ---------------------------------------------------------------------------
# Live
# ---------------------------------------------------------------------------

@pytest.mark.skipif(not os.environ.get("PROPREP_NETWORK_TESTS"), reason="set PROPREP_NETWORK_TESTS=1 to fetch from PDBe and RCSB")
def test_live_1ltz_the_strongest_peak_is_the_iron(tmp_path):
    pdb = str(tmp_path / "1LTZ.pdb")
    urllib.request.urlretrieve("https://files.rcsb.org/download/1LTZ.pdb", pdb)
    assert dm.fetch_density_maps("2RVD", str(tmp_path)) is None     # NMR

    paths = dm.fetch_density_maps("1LTZ", str(tmp_path))
    prepared = dm.prepare_maps_for_structure(pdb, paths, str(tmp_path), margin=5.0)
    recut = dm.read_ccp4(prepared[0]["path"])
    iron = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in open(pdb) if l.startswith("HETATM") and l[17:20] == " FE"])
    cell, scale, shift, _ = dm.read_pdb_frame(pdb)

    peak = np.array(np.unravel_index(recut.grid.argmax(), recut.grid.shape)) + np.array(recut.start)
    position = (peak / np.array(recut.sampling)) @ np.linalg.inv(dm.fractionalization_matrix(cell)).T
    assert np.linalg.norm(position - iron[0]) < recut.spacing       # 0.30 A; NGL's own parser puts it there too
    assert dm.mean_density_at_atoms(paths["2fofc"], pdb, pdb) > 3.0
