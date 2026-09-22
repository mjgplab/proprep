"""The placement ProPrep is about to build is compared with the OPM database's.

Checked live on 6R2Q when written: PPM3's placement is 1.1 degrees from OPM's
membrane normal and 0.4 A from its depth, the same way up; MEMEMBED's is 56.1
degrees and 85 A away; OPM answers "no entry" for a soluble protein (1UBQ).

The prepared protein below is renumbered from 1, has no chain IDs, renamed
residues and added hydrogens, as tLEaP and ProPrep leave it, so nothing in it
can be paired with OPM's file by name. Its atoms have not moved, and that is
what the comparison rests on.
"""

import io
import types

import numpy as np
import pytest
from rich.console import Console

from proprep.membrane_prep import membrane_builder as mb
from proprep.membrane_prep import opm_check as oc


def rot(axis, degrees):
    a, (x, y, z) = np.radians(degrees), np.array(axis, float) / np.linalg.norm(axis)
    k = np.array([[0, -z, y], [z, 0, -x], [-y, x, 0]])
    return np.eye(3) + np.sin(a) * k + (1 - np.cos(a)) * k @ k


# Deposited model: chain, residue name, residue number, atom, element, x, y, z
DEPOSITED = [
    ("A", "HIS", 71, "N", "N", 4.171, 112.259, 9.750), ("A", "HIS", 71, "CA", "C", 5.402, 111.560, 10.210),
    ("A", "HIS", 71, "NE2", "N", 8.950, 109.120, 12.480), ("A", "CYS", 112, "SG", "S", 70.176, 10.988, 5.768),
    ("B", "LEU", 300, "CA", "C", 33.420, 61.730, -14.905), ("B", "LEU", 300, "CB", "C", 34.100, 60.480, -15.522),
    ("C", "HEC", 901, "FE", "FE", 23.639, 107.015, 3.850), ("C", "HEC", 901, "NA", "N", 25.100, 106.000, 4.900),
]
PREPARED_NAMES = ["HIO", "HIO", "HIO", "CYX", "LEU", "LEU", "HCO", "HCO"]
TO_OPM = (rot([1, 2, -0.5], 218.4), np.array([12.5, -40.25, 7.0]))       # deposited frame -> OPM's frame


def _line(i, chain, resname, resid, name, element, xyz):
    return (f"ATOM  {i:5d} {name:<4s} {resname:>3s} {chain}{resid:4d}    "
            f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00  0.00          {element:>2s}\n")


def _deposited(path, move=None, header=True):
    text = "HEADER    ELECTRON TRANSPORT                      18-MAR-19   6R2Q              \n" if header else ""
    for i, (c, rn, ri, n, el, *xyz) in enumerate(DEPOSITED, 1):
        xyz = np.array(xyz)
        text += _line(i, c, rn, ri, n, el, xyz @ move[0] + move[1] if move else xyz)
    path.write_text(text + "END\n")
    return str(path)


def _prepared(path, move=None):
    """Renumbered, chainless, renamed, with a hydrogen after every atom."""
    text, i = "", 0
    for k, ((c, rn, ri, n, el, *xyz), new_name) in enumerate(zip(DEPOSITED, PREPARED_NAMES)):
        for atom, element, shift in ((n, el, 0.0), (f"H{k}", "H", 0.9)):
            i += 1
            p = np.array(xyz) + shift
            text += _line(i, " ", new_name, k + 1, atom, element, p @ move[0] + move[1] if move else p)
    path.write_text(text + "END\n")
    return str(path)


def _files(tmp_path, placement):
    opm = _deposited(tmp_path / "opm_6r2q.pdb", move=TO_OPM, header=False)
    (tmp_path / "opm_6r2q.pdb").write_text("REMARK      1/2 of bilayer thickness:   11.3\n"
                                           + (tmp_path / "opm_6r2q.pdb").read_text())
    return (_prepared(tmp_path / "placed.pdb", move=placement), _prepared(tmp_path / "protein_with_h.pdb"),
            _deposited(tmp_path / "6R2Q.pdb"), opm)


def _then(first, second):
    """Apply ``first`` then ``second`` (each a rotation and a shift, for row vectors)."""
    return first[0] @ second[0], first[1] @ second[0] + second[1]


def test_the_same_placement_reads_as_identical(tmp_path):
    result = oc.cross_check(*_files(tmp_path, TO_OPM))
    assert result["tilt_deg"] < 0.05 and result["flipped"] is False
    assert abs(result["centre_z"] - result["reference_centre_z"]) < 0.01
    assert result["atoms_in_deposited_frame"] == 8 and result["opm_half_thickness"] == 11.3


def test_spinning_about_the_normal_and_sliding_in_the_plane_change_nothing(tmp_path):
    same = _then(TO_OPM, (rot([0, 0, 1], 77.0), np.array([31.0, -18.0, 0.0])))
    result = oc.cross_check(*_files(tmp_path, same))
    assert result["tilt_deg"] < 0.05 and result["flipped"] is False
    assert abs(result["centre_z"] - result["reference_centre_z"]) < 0.01


def test_a_tilt_and_a_depth_error_are_measured(tmp_path):
    off = _then(TO_OPM, (rot([1, 0, 0], 25.0), np.zeros(3)))
    off = (off[0], off[1] + np.array([0.0, 0.0, -6.5]))
    result = oc.cross_check(*_files(tmp_path, off))
    assert result["tilt_deg"] == pytest.approx(25.0, abs=0.1) and result["flipped"] is False
    # 6.5 A lower, on top of what tilting about x does to the centre.
    tilted = _then(TO_OPM, (rot([1, 0, 0], 25.0), np.zeros(3)))
    expected = np.array([d[5:] for d in DEPOSITED]).mean(0) @ tilted[0] + tilted[1]
    assert result["centre_z"] == pytest.approx(expected[2] - 6.5 + 0.9 / 2 * tilted[0][:, 2].sum(), abs=0.05)


def test_upside_down_is_said_in_words(tmp_path):
    flipped = _then(TO_OPM, (rot([1, 0, 0], 180.0), np.zeros(3)))
    result = oc.cross_check(*_files(tmp_path, flipped))
    assert result["flipped"] is True and result["tilt_deg"] < 0.05
    assert abs(result["centre_z"] - result["reference_centre_z"]) < 0.01      # compared turned back over
    text = "\n".join(oc.report_lines("6r2q", result))
    assert "OTHER WAY UP" in text and "N-terminus orientation" in text
    assert "same way up" in "\n".join(oc.report_lines("6r2q", oc.cross_check(*_files(tmp_path, TO_OPM))))


def _repaired(path, move=None):
    """The prepared protein after a MODELLER repair with ProPrep's freeze, as seen live on a
    6R2Q fragment (4-residue gap, MODELLER 10.8): every chain renamed and every residue
    renumbered (B 400 became A 1), a built loop that the deposited model does not have,
    and the NAMES of equivalent atoms swapped (OD1/OD2, NH1/NH2, ...) while no resolved
    atom moved: 487 of 488 were still at exactly their deposited coordinates. Without the
    freeze, none were."""
    swapped = {"N": "NE2", "NE2": "N"}                     # labels exchanged, atoms where they were
    text, i = "", 0
    for k, (c, rn, ri, n, el, *xyz) in enumerate(DEPOSITED):
        i += 1
        p = np.array(xyz)
        text += _line(i, "A", rn, 500 + k, swapped.get(n, n) if rn == "HIS" else n, el,
                      p @ move[0] + move[1] if move else p)
        if k == 3:                                          # a loop MODELLER built, between two residues
            for j, loop_xyz in enumerate(([40.0, 60.0, 1.0], [41.2, 60.9, 1.4], [42.5, 60.1, 2.2])):
                i += 1
                q = np.array(loop_xyz)
                text += _line(i, "A", "GLY", 600 + j, "CA", "C", q @ move[0] + move[1] if move else q)
    path.write_text(text + "END\n")
    return str(path)


def test_a_structure_repaired_by_modeller_is_still_compared(tmp_path):
    _, _, loaded, opm = _files(tmp_path, TO_OPM)
    unplaced = _repaired(tmp_path / "repaired.pdb")
    placed = _repaired(tmp_path / "repaired_oriented.pdb", move=TO_OPM)
    assert oc.find_loaded_entry({"rcsb_pdb_file": loaded}, unplaced) == (loaded, "6R2Q")
    result = oc.cross_check(placed, unplaced, loaded, opm)
    assert result["atoms_in_deposited_frame"] == len(DEPOSITED)       # the built loop is not among them
    assert result["tilt_deg"] < 0.05 and result["flipped"] is False
    assert abs(result["centre_z"] - result["reference_centre_z"]) < 0.01


def test_a_structure_moved_since_loading_is_not_compared(tmp_path):
    placed, _, loaded, opm = _files(tmp_path, TO_OPM)
    moved = _prepared(tmp_path / "aligned.pdb", move=(rot([0, 1, 0], 40.0), np.array([3.0, 0, 0])))
    with pytest.raises(oc.OpmCheckError, match="has been moved since"):
        oc.cross_check(placed, moved, loaded, opm)


def test_an_opm_file_of_another_model_is_not_compared(tmp_path):
    placed, unplaced, loaded, _ = _files(tmp_path, TO_OPM)
    other = _deposited(tmp_path / "other.pdb", move=(np.diag([1.0, 1.0, 1.4]), np.zeros(3)), header=False)
    with pytest.raises(oc.OpmCheckError, match="not the same model"):
        oc.cross_check(placed, unplaced, loaded, other)


def test_the_entry_is_found_by_position_and_named_by_its_header(tmp_path):
    _, unplaced, loaded, _ = _files(tmp_path, TO_OPM)
    elsewhere = _deposited(tmp_path / "1abc.pdb", move=(rot([0, 0, 1], 10.0), np.array([50.0, 0, 0])))
    workspace = {"rcsb_pdb_files": [elsewhere, loaded], "local_pdb_file": str(tmp_path / "missing.pdb")}
    assert oc.find_loaded_entry(workspace, unplaced) == (loaded, "6R2Q")
    assert oc.find_loaded_entry({"rcsb_pdb_file": elsewhere}, unplaced) is None
    assert oc.find_loaded_entry(None, unplaced) is None


# ── Asking OPM ──────────────────────────────────────────────────────────

def _requests(monkeypatch, status=200, content=b"ATOM      1  N   LYS A   1\n", error=None):
    import requests
    asked = []

    def get(url, timeout):
        asked.append(url)
        if error:
            raise error
        return types.SimpleNamespace(status_code=status, content=content)

    monkeypatch.setattr(requests, "get", get)
    return asked


def test_entry_is_downloaded_once_and_reused(monkeypatch, tmp_path):
    asked = _requests(monkeypatch)
    path = oc.fetch_opm_entry("6R2Q", str(tmp_path))
    assert path.endswith("opm_6r2q.pdb") and asked == [oc.OPM_URL.format(pdb_id="6r2q")]
    assert oc.fetch_opm_entry("6r2q", str(tmp_path)) == path and len(asked) == 1


def test_an_entry_opm_does_not_hold_is_none_not_an_error(monkeypatch, tmp_path):
    _requests(monkeypatch, status=404)
    assert oc.fetch_opm_entry("1UBQ", str(tmp_path)) is None


def test_no_network_is_said_briefly(monkeypatch, tmp_path):
    import requests
    _requests(monkeypatch, error=requests.exceptions.ConnectionError("HTTPSConnectionPool(...) " * 20))
    with pytest.raises(oc.OpmCheckError) as raised:
        oc.fetch_opm_entry("6R2Q", str(tmp_path))
    assert "internet connection" in str(raised.value) and "HTTPSConnectionPool" not in str(raised.value)


# ── In the build ────────────────────────────────────────────────────────

def _module():
    module = mb.MembraneBuilderModule()
    out = io.StringIO()
    module.processor = types.SimpleNamespace(console=Console(file=out, width=300))
    return module, out


def test_the_build_reports_the_comparison(monkeypatch, tmp_path):
    placed, unplaced, loaded, opm = _files(tmp_path, TO_OPM)
    monkeypatch.setattr(oc, "fetch_opm_entry", lambda pdb_id, out_dir: opm)
    module, out = _module()
    module._cross_check_with_opm({"rcsb_pdb_file": loaded}, placed, unplaced, str(tmp_path))
    said = out.getvalue()
    assert "OPM has 6R2Q" in said and "same way up" in said and "membrane normal: 0.0°" in said


@pytest.mark.parametrize("fetch, expected", [
    (lambda pdb_id, out_dir: None, "OPM has no entry for 6R2Q"),
    (lambda pdb_id, out_dir: (_ for _ in ()).throw(oc.OpmCheckError("OPM could not be reached")), "Not compared with OPM"),
])
def test_nothing_about_opm_can_stop_a_build(monkeypatch, tmp_path, fetch, expected):
    placed, unplaced, loaded, _ = _files(tmp_path, TO_OPM)
    monkeypatch.setattr(oc, "fetch_opm_entry", fetch)
    module, out = _module()
    assert module._cross_check_with_opm({"rcsb_pdb_file": loaded}, placed, unplaced, str(tmp_path)) is None
    assert expected in out.getvalue()


def test_it_can_be_switched_off(monkeypatch, tmp_path):
    placed, unplaced, loaded, _ = _files(tmp_path, TO_OPM)
    monkeypatch.setattr(oc, "fetch_opm_entry", lambda *a: pytest.fail("OPM was asked"))
    module, out = _module()
    module.config.opm_cross_check = False
    module._cross_check_with_opm({"rcsb_pdb_file": loaded}, placed, unplaced, str(tmp_path))
    assert out.getvalue() == ""
