"""Release tooling: the version-lockstep guard and the public-snapshot export.

These scripts encode the release procedure (docs/RELEASE_PROCEDURE.md); the
tests keep them honest about which files carry the version and what the
public snapshot contains.
"""

from __future__ import annotations

import re
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
GUARD = ROOT / "scripts" / "check_version_lockstep.sh"
EXPORT = ROOT / "scripts" / "export_public_snapshot.sh"

PINNED_FILES = [
    "pyproject.toml",
    "setup.py",
    "recipe/meta.yaml",
    "install_proprep.sh",
    "update_proprep_in_ambertools.sh",
    "CITATION.cff",
    "constructor/construct.yaml",
]


def run(cmd, **kw):
    return subprocess.run(cmd, capture_output=True, text=True, **kw)


def test_guard_reports_every_pinned_file_and_agrees():
    res = run(["bash", str(GUARD)])
    assert res.returncode == 0, res.stdout + res.stderr
    for f in PINNED_FILES:
        assert f in res.stdout, f"{f} not reported"
    # constructor carries two pins and both must be listed
    assert "construct.yaml (version:)" in res.stdout and "construct.yaml (proprep=)" in res.stdout
    m = re.search(r"OK: all version pins agree on (\d+\.\d+\.\d+)", res.stdout)
    assert m, res.stdout
    version = m.group(1)
    # asserting the agreed value works both ways
    assert run(["bash", str(GUARD), version]).returncode == 0
    assert run(["bash", str(GUARD), "0.0.0"]).returncode == 1


def test_guard_fails_when_one_pin_drifts(tmp_path):
    """Copy the repo's pinned files into a scratch tree and drift one pin."""
    scratch = tmp_path / "repo"
    (scratch / "scripts").mkdir(parents=True)
    shutil.copy(GUARD, scratch / "scripts" / GUARD.name)
    for f in PINNED_FILES:
        dst = scratch / f
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(ROOT / f, dst)
    upd = scratch / "update_proprep_in_ambertools.sh"
    upd.write_text(re.sub(r'^PROPREP_VERSION="[^"]+"', 'PROPREP_VERSION="1.14.0"', upd.read_text(), count=1, flags=re.M))
    res = run(["bash", str(scratch / "scripts" / GUARD.name)])
    assert res.returncode == 1
    assert "NOT in lockstep" in res.stdout


def test_documented_pins_match_the_guard():
    """docs/RELEASE_PROCEDURE.md must list exactly the files the guard checks."""
    doc = (ROOT / "docs" / "RELEASE_PROCEDURE.md").read_text()
    for f in PINNED_FILES:
        assert f"`{f}`" in doc, f"{f} missing from the release procedure's pin table"


@pytest.mark.skipif(
    run(["git", "-C", str(ROOT), "rev-parse", "-q", "--verify", "refs/tags/v1.16.0"]).returncode != 0
    or shutil.which("zip") is None,
    reason="needs the v1.16.0 tag and the zip tool",
)
def test_export_public_snapshot_applies_the_rules(tmp_path):
    res = run(["bash", str(EXPORT), "v1.16.0", str(tmp_path)])
    assert res.returncode == 0, res.stdout + res.stderr
    snap = tmp_path / "proprep-v1.16.0-clean"
    assert snap.is_dir() and (tmp_path / "proprep-v1.16.0-clean.zip").is_file()
    for excluded in ("docs", "examples", "prototypes", "tools", "README.public.md"):
        assert not (snap / excluded).exists(), f"{excluded} must not be in the public snapshot"
    for kept in ("src/proprep/main.py", "tests", "recipe/meta.yaml", "install_proprep.sh",
                 "update_proprep_in_ambertools.sh", "CITATION.cff", ".zenodo.json", "LICENSE", "README.md"):
        assert (snap / kept).exists(), f"{kept} missing from the public snapshot"
    public_readme = run(["git", "-C", str(ROOT), "show", "v1.16.0:README.public.md"]).stdout
    assert (snap / "README.md").read_text() == public_readme
    # exactly the tracked files minus the exclusions (723 at 1.16.0, matched
    # file-for-file against github.com/mjgplab/proprep on 2026-08-29)
    tracked = run(["git", "-C", str(ROOT), "ls-tree", "-r", "--name-only", "v1.16.0"]).stdout.split()
    expected = [
        p for p in tracked
        if not p.startswith(("docs/", "examples/", "prototypes/", "tools/")) and p != "README.public.md"
    ]
    exported = sorted(str(p.relative_to(snap)) for p in snap.rglob("*") if p.is_file())
    assert exported == sorted(expected)


def test_export_refuses_unknown_tag(tmp_path):
    res = run(["bash", str(EXPORT), "v0.0.0", str(tmp_path)])
    assert res.returncode == 2 and "does not exist" in res.stderr
