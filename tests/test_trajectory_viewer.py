"""Playing an MD trajectory in the structure viewer.

NGL reads Amber NetCDF in the browser; the server streams the file at
/trajectory/<index>, the viewer config names it per structure, and one
cpptraj run writes a matching first-frame PDB + NetCDF (optionally
stripped of solvent) so the atoms line up.
"""

from __future__ import annotations

import os
import stat
import urllib.request
from pathlib import Path

import pytest

from proprep.md_prep.trajectory_view import cpptraj_input, write_view_files, SOLVENT_MASK
from proprep.structure_prep.interactive_structure_viewer import InteractiveStructureViewer
from proprep.structure_prep.viewer_coordinator import ViewerCoordinator
from proprep.structure_prep.viewer_server import ViewerServer

ROOT = Path(__file__).resolve().parents[1]


def test_cpptraj_input_images_then_optionally_strips():
    s = cpptraj_input("p.prmtop", "prod.nc", "v.pdb", "v.nc", strip_solvent=True)
    lines = s.splitlines()
    assert lines[:3] == ["parm p.prmtop", "trajin prod.nc", "autoimage"]
    assert lines[3] == f"strip {SOLVENT_MASK}"
    assert "trajout v.pdb pdb onlyframes 1" in lines and "trajout v.nc netcdf" in lines
    assert "strip" not in cpptraj_input("p", "t", "a", "b", strip_solvent=False)


def test_write_view_files_names_outputs_and_fails_loudly(tmp_path):
    fake = tmp_path / "cpptraj"
    fake.write_text("#!/bin/sh\n# write whatever trajout names the script asks for\n"
                    "for f in $(grep '^trajout' \"$2\" | awk '{print $2}'); do echo x > \"$f\"; done\n")
    fake.chmod(fake.stat().st_mode | stat.S_IXUSR)
    pdb, nc = write_view_files("p.prmtop", str(tmp_path / "prod.nc"), str(tmp_path / "viewer"),
                               strip_solvent=True, cpptraj=str(fake))
    assert pdb.name == "prod_view_stripped.pdb" and nc.name == "prod_view_stripped.nc"
    assert pdb.exists() and nc.exists()
    pdb2, _ = write_view_files("p.prmtop", str(tmp_path / "prod.nc"), str(tmp_path / "viewer"),
                               strip_solvent=False, cpptraj=str(fake))
    assert pdb2.name == "prod_view.pdb"

    broken = tmp_path / "broken"; broken.write_text("#!/bin/sh\necho boom; exit 3\n")
    broken.chmod(broken.stat().st_mode | stat.S_IXUSR)
    with pytest.raises(RuntimeError, match="exit 3"):
        write_view_files("p", str(tmp_path / "prod.nc"), str(tmp_path / "v2"), cpptraj=str(broken))


def test_viewer_config_names_the_trajectory_per_structure():
    v = InteractiveStructureViewer()
    v.selected_structures = ["/tmp/a_view.pdb", "/tmp/other.pdb"]
    v.trajectory_files = {0: "/tmp/a_view.nc"}
    structures = v._build_viewer_config()["structures"]
    assert structures[0]["trajectory"] == {"url": "/trajectory/0", "name": "a_view.nc", "ext": "nc"}
    assert "trajectory" not in structures[1]


def test_coordinator_show_trajectory_binds_and_launches(monkeypatch):
    coord = ViewerCoordinator()
    v = InteractiveStructureViewer()
    coord._viewer = v
    launched = []
    monkeypatch.setattr(v, "_launch_viewer", lambda open_browser=True: launched.append(open_browser) or True)
    monkeypatch.setattr(coord, "is_running", lambda: False)
    coord.show_trajectory("/tmp/x.pdb", "/tmp/x.nc", show_waters=True, force=True)
    assert v.selected_structures == ["/tmp/x.pdb"]
    assert v.trajectory_files == {0: "/tmp/x.nc"}
    assert v.viewer_config["show_waters"] is True
    assert launched == [True]
    # a later plain structure view drops the trajectory
    monkeypatch.setattr(coord, "is_running", lambda: True)
    coord.show_structure("/tmp/y.pdb")
    assert v.trajectory_files == {}


def test_server_streams_trajectory_bytes(tmp_path):
    pdb = tmp_path / "s.pdb"; pdb.write_text("ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\nEND\n")
    nc = tmp_path / "s.nc"; payload = os.urandom(70000); nc.write_bytes(payload)
    server = ViewerServer(config={"structures": []}, structure_files=[str(pdb)],
                          port=8999, trajectory_files={0: str(nc)})
    server.port = server.find_available_port(8999)
    assert server.start(open_browser=False)
    try:
        url = f"http://127.0.0.1:{server.port}/trajectory/0?v=1"
        with urllib.request.urlopen(url, timeout=5) as r:
            assert r.headers["Content-Type"] == "application/x-netcdf"
            assert r.read() == payload
        with pytest.raises(urllib.error.HTTPError) as exc:
            urllib.request.urlopen(f"http://127.0.0.1:{server.port}/trajectory/1", timeout=5)
        assert exc.value.code == 404
    finally:
        server.stop()


def test_template_has_player_and_attach_hook():
    html = (ROOT / "src/proprep/structure_prep/templates/ngl_viewer.html").read_text()
    for needle in ('id="trajectory-section"', "NGL.autoLoad(url", "addTrajectory(frames",
                   "NGL.TrajectoryPlayer", "attachTrajectory(idx, structInfo)", "seekTrajectory",
                   "stepTrajectory(", "setTrajectoryParam('timeout'", "setTrajectoryParam('mode'",
                   "setTrajectoryParam('direction'", "setTrajectoryParam('step'",
                   "setTrajectorySuperpose(", "setTrajectoryOption('removePbc'",
                   "setTrajectoryOption('centerPbc'", "e.key === ' '", "isTypingInField"):
        assert needle in html, needle


BENCH = Path("/private/tmp/claude-501/-Users-mjgp-Desktop-Software-proprep/5f31c069-8233-4513-842d-b218f5ab85c4/scratchpad/bench/workshop")


@pytest.mark.skipif(not (BENCH / "prod.nc").exists() or not os.environ.get("AMBERHOME"),
                    reason="needs the local benchmark trajectory and AmberTools")
def test_real_cpptraj_writes_matching_pdb_and_nc(tmp_path):
    from proprep.md_prep.trajectory_view import frame_count
    pdb, nc = write_view_files(str(BENCH / "prmtop"), str(BENCH / "prod.nc"), str(tmp_path), strip_solvent=True)
    atoms = sum(1 for l in pdb.read_text().splitlines() if l.startswith(("ATOM", "HETATM")))
    assert atoms == 138 and frame_count(str(nc)) == 300
