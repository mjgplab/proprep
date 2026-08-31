"""proprep-library: snapshot, reset and restore the user library without ever
deleting anything."""
import json
from pathlib import Path

import pytest

from proprep.utils.user_library_snapshots import UserLibrary, main, LIBRARY_PARTS


def _seed(root: Path):
    ff = root / "forcefield_params" / "specialized_residues" / "heme" / "x"
    ff.mkdir(parents=True); (ff / "a.frcmod").write_text("A"); (ff / "a.lib").write_text("L")
    tr = root / "transformers"; tr.mkdir()
    (tr / "t_one.json").write_text("{}"); (tr / "__pycache__").mkdir(); (tr / "__pycache__" / "x.pyc").write_text("")
    (root / "settings.json").write_text("{}")            # not library: must survive
    (root / "feedback_private_key").write_text("secret")  # not library: must survive


def test_snapshot_copies_library_only_and_skips_caches(tmp_path):
    _seed(tmp_path); lib = UserLibrary(tmp_path)
    dest = lib.create_snapshot("fig1", note="lysozyme practice")
    assert (dest / "forcefield_params" / "specialized_residues" / "heme" / "x" / "a.frcmod").exists()
    assert (dest / "transformers" / "t_one.json").exists()
    assert not (dest / "transformers" / "__pycache__").exists()
    assert not (dest / "settings.json").exists()
    m = json.loads((dest / "manifest.json").read_text())
    assert m["name"] == "fig1" and m["note"] == "lysozyme practice"
    assert m["parts"]["forcefield_params"]["files"] == 2 and m["parts"]["transformers"]["files"] == 1
    with pytest.raises(FileExistsError):
        lib.create_snapshot("fig1")


def test_reset_snapshots_then_empties_and_leaves_settings_alone(tmp_path):
    _seed(tmp_path); lib = UserLibrary(tmp_path)
    saved = lib.reset()
    assert saved is not None and saved.name.startswith("before_reset_")
    assert lib.is_empty()
    for part in LIBRARY_PARTS:
        assert (tmp_path / part).is_dir()                      # directories kept for code paths
    assert (tmp_path / "settings.json").read_text() == "{}"
    assert (tmp_path / "feedback_private_key").read_text() == "secret"
    assert (saved / "transformers" / "t_one.json").exists()     # nothing lost
    assert lib.reset() is None                                  # empty: nothing to snapshot


def test_restore_brings_a_snapshot_back_and_saves_the_current_state(tmp_path):
    _seed(tmp_path); lib = UserLibrary(tmp_path)
    lib.create_snapshot("fig1")
    lib.reset()
    (tmp_path / "transformers" / "t_two.json").write_text("{}")  # new work after reset
    saved = lib.restore("fig1")
    assert saved.name.startswith("before_restore_")
    assert (saved / "transformers" / "t_two.json").exists()      # the new work was kept
    assert (tmp_path / "transformers" / "t_one.json").exists()   # the snapshot is back
    assert not (tmp_path / "transformers" / "t_two.json").exists()
    with pytest.raises(FileNotFoundError):
        lib.restore("nope")


def test_cli_round_trip(tmp_path, capsys):
    _seed(tmp_path)
    assert main(["--root", str(tmp_path), "status"]) == 0
    assert "forcefield_params" in capsys.readouterr().out
    assert main(["--root", str(tmp_path), "snapshot", "clean", "-m", "fresh install"]) == 0
    assert main(["--root", str(tmp_path), "reset", "-y"]) == 0
    assert UserLibrary(tmp_path).is_empty()
    assert main(["--root", str(tmp_path), "list"]) == 0
    out = capsys.readouterr().out
    assert "clean" in out and "before_reset_" in out and "fresh install" in out
    assert main(["--root", str(tmp_path), "restore", "clean", "-y"]) == 0
    assert (tmp_path / "transformers" / "t_one.json").exists()
    assert main(["--root", str(tmp_path), "restore", "missing", "-y"]) == 1
