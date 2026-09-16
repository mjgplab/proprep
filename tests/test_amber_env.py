"""bootstrap_amber_env: a bare console-script launch must still find tleap."""

from __future__ import annotations

import os
import stat
import sys
from pathlib import Path

import pytest

from proprep.utils import amber_env


def _make_amber_tree(root: Path) -> Path:
    (root / "dat" / "leap" / "parm").mkdir(parents=True)
    (root / "bin").mkdir()
    tleap = root / "bin" / "tleap"
    tleap.write_text("#!/bin/sh\n")
    tleap.chmod(tleap.stat().st_mode | stat.S_IXUSR)
    python = root / "bin" / "python3.12"
    python.write_text("")
    return python


@pytest.fixture
def bundled(tmp_path, monkeypatch):
    """A fake ProPrep prefix whose interpreter is the running one."""
    root = tmp_path / "ProPrep"
    python = _make_amber_tree(root)
    monkeypatch.setattr(sys, "executable", str(python))
    monkeypatch.setattr(sys, "prefix", str(tmp_path / "elsewhere"))
    monkeypatch.setattr(sys, "base_prefix", str(tmp_path / "elsewhere"))
    return root


def test_bare_launch_sets_amberhome_and_prepends_bin(bundled):
    env = {"PATH": "/usr/bin:/bin"}
    assert amber_env.bootstrap_amber_env(env) == str(bundled)
    assert env["AMBERHOME"] == str(bundled)
    assert env["PATH"].split(os.pathsep) == [str(bundled / "bin"), "/usr/bin", "/bin"]


def test_idempotent(bundled):
    env = {"PATH": "/usr/bin"}
    amber_env.bootstrap_amber_env(env)
    once = dict(env)
    amber_env.bootstrap_amber_env(env)
    assert env == once
    assert env["PATH"].count(str(bundled / "bin")) == 1


def test_bundled_bin_moved_to_front_if_buried(bundled):
    env = {"PATH": f"/opt/stale-amber/bin:{bundled / 'bin'}"}
    amber_env.bootstrap_amber_env(env)
    assert env["PATH"].split(os.pathsep) == [str(bundled / "bin"), "/opt/stale-amber/bin"]


def test_existing_valid_amberhome_is_respected(bundled, tmp_path):
    other = tmp_path / "amber24"
    _make_amber_tree(other)
    env = {"AMBERHOME": str(other), "PATH": "/usr/bin"}
    assert amber_env.bootstrap_amber_env(env) == str(other)
    assert env["AMBERHOME"] == str(other)
    assert env["PATH"].split(os.pathsep)[0] == str(other / "bin")
    assert str(bundled / "bin") not in env["PATH"]


def test_stale_amberhome_is_replaced(bundled, tmp_path):
    env = {"AMBERHOME": str(tmp_path / "gone"), "PATH": "/usr/bin"}
    amber_env.bootstrap_amber_env(env)
    assert env["AMBERHOME"] == str(bundled)


def test_no_amber_tree_leaves_env_alone(tmp_path, monkeypatch):
    monkeypatch.setattr(sys, "executable", str(tmp_path / "venv" / "bin" / "python"))
    monkeypatch.setattr(sys, "prefix", str(tmp_path / "venv"))
    monkeypatch.setattr(sys, "base_prefix", str(tmp_path / "venv"))
    env = {"PATH": "/usr/bin"}
    assert amber_env.bootstrap_amber_env(env) is None
    assert env == {"PATH": "/usr/bin"}


def test_sys_prefix_fallback_when_interpreter_is_wrapped(tmp_path, monkeypatch):
    root = tmp_path / "ProPrep"
    _make_amber_tree(root)
    monkeypatch.setattr(sys, "executable", str(tmp_path / "wrapper" / "bin" / "python"))
    monkeypatch.setattr(sys, "prefix", str(root))
    monkeypatch.setattr(sys, "base_prefix", str(root))
    env: dict = {}
    assert amber_env.bootstrap_amber_env(env) == str(root)
    assert env["PATH"] == str(root / "bin")


def test_symlinked_interpreter_resolves_to_real_tree(bundled, tmp_path):
    venv_bin = tmp_path / "venv" / "bin"
    venv_bin.mkdir(parents=True)
    link = venv_bin / "python"
    link.symlink_to(bundled / "bin" / "python3.12")
    sys.executable = str(link)
    env: dict = {}
    assert amber_env.bootstrap_amber_env(env) == str(bundled)


def test_default_operates_on_os_environ(bundled, monkeypatch):
    monkeypatch.delenv("AMBERHOME", raising=False)
    monkeypatch.setenv("PATH", "/usr/bin")
    amber_env.bootstrap_amber_env()
    assert os.environ["AMBERHOME"] == str(bundled)
    assert os.environ["PATH"].startswith(str(bundled / "bin"))


# ---------------------------------------------------------------------------
# wiring: both console-script entry points bootstrap before doing anything
# ---------------------------------------------------------------------------

def test_proprep_cli_bootstraps_first(monkeypatch):
    calls = []
    monkeypatch.setattr(amber_env, "bootstrap_amber_env", lambda env=None: calls.append("cli"))
    monkeypatch.setattr(sys, "argv", ["proprep", "--version"])
    from proprep import main as cli
    with pytest.raises(SystemExit) as exc:
        cli.main()
    assert exc.value.code == 0
    assert calls == ["cli"]


def test_proprep_web_bootstraps_before_uvicorn(monkeypatch):
    calls = []
    monkeypatch.setattr(amber_env, "bootstrap_amber_env", lambda env=None: calls.append("web"))
    import uvicorn
    monkeypatch.setattr(uvicorn, "run", lambda *a, **k: calls.append("uvicorn"))
    monkeypatch.setattr(sys, "argv", ["proprep-web", "--no-browser", "--strict-port", "--port", "0"])
    from proprep.web import __main__ as web
    assert web.main() == 0
    assert calls == ["web", "uvicorn"]
