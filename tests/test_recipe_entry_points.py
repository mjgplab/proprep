"""Every console script in pyproject.toml must be a conda entry point too.

The conda package is ``noarch: python``. Conda writes a launcher at install
time, with the user's interpreter in its shebang, only for the commands listed
under ``build: entry_points`` in recipe/meta.yaml. A command that is missing
there still ships, as the launcher pip wrote on the build machine, pinned to
that machine's interpreter (e.g. ``#!.../bin/python3.13`` in an environment
that has Python 3.12), and fails with "bad interpreter".
"""

import re
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def _pyproject_scripts() -> dict[str, str]:
    with open(ROOT / "pyproject.toml", "rb") as fh:
        return tomllib.load(fh)["project"]["scripts"]


def _recipe_entry_points() -> dict[str, str]:
    text = (ROOT / "recipe" / "meta.yaml").read_text()
    block = re.search(r"^  entry_points:\n((?:    - .*\n)+)", text, re.MULTILINE)
    assert block, "build: entry_points block not found in recipe/meta.yaml"
    entries = {}
    for line in block.group(1).splitlines():
        name, _, target = line.strip()[2:].partition("=")
        entries[name.strip()] = target.strip()
    return entries


def test_every_console_script_is_a_conda_entry_point():
    scripts = _pyproject_scripts()
    entry_points = _recipe_entry_points()
    missing = sorted(set(scripts) - set(entry_points))
    assert not missing, f"add to recipe/meta.yaml entry_points: {missing}"


def test_entry_points_match_pyproject_targets():
    scripts = _pyproject_scripts()
    for name, target in _recipe_entry_points().items():
        assert scripts.get(name) == target, f"{name}: recipe {target!r} vs pyproject {scripts.get(name)!r}"
