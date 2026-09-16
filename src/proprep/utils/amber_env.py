"""Make the bundled AmberTools reachable from a bare console-script launch.

ProPrep ships AmberTools inside its own prefix: ``bin/`` holds the Python
interpreter next to ``tleap``, ``antechamber``, ``cpptraj`` and the rest,
and ``dat/leap`` holds their data. That prefix only becomes an Amber
environment through activation (``conda activate`` or ``source amber.sh``),
which exports ``AMBERHOME`` and prepends ``bin/`` to ``PATH``. A setuptools
console script such as ``~/ProPrep/bin/proprep-web`` skips activation: its
shebang pins the interpreter and nothing else runs, so every
``subprocess.run(["tleap", ...])`` fails with "not found in PATH" and every
``os.environ["AMBERHOME"]`` lookup raises.

:func:`bootstrap_amber_env` repairs the process environment from the running
interpreter. Each entry point calls it once, before anything that shells out.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path
from typing import MutableMapping, Optional


def _is_amber_root(root: Path) -> bool:
    """True when ``root`` holds both the leap data and a runnable tleap."""
    return (root / "dat" / "leap" / "parm").is_dir() and os.access(
        root / "bin" / "tleap", os.X_OK
    )


def _candidate_roots() -> list[Path]:
    """Prefixes that might be the bundled Amber tree, most specific first.

    The resolved interpreter comes first so a venv layered on the conda env
    still finds the conda tree; ``sys.prefix``/``sys.base_prefix`` cover an
    interpreter launched through a wrapper. Duplicates are dropped.
    """
    seen: list[Path] = []
    for candidate in (
        Path(sys.executable).resolve().parent.parent,
        Path(sys.prefix),
        Path(sys.base_prefix),
    ):
        if candidate not in seen:
            seen.append(candidate)
    return seen


def find_amber_home(env: Optional[MutableMapping[str, str]] = None) -> Optional[Path]:
    """The Amber tree this process should use, or None.

    An already-set, valid ``AMBERHOME`` wins: it is a deliberate choice
    (a cluster's full Amber with pmemd.cuda, say) and ``PATH`` must agree
    with it, not with the bundled tree. Otherwise the bundled tree is
    located from the interpreter.
    """
    env = os.environ if env is None else env
    existing = env.get("AMBERHOME")
    if existing and _is_amber_root(Path(existing)):
        return Path(existing)
    for root in _candidate_roots():
        if _is_amber_root(root):
            return root
    return None


def bootstrap_amber_env(env: Optional[MutableMapping[str, str]] = None) -> Optional[str]:
    """Set ``AMBERHOME`` and put its ``bin/`` first on ``PATH``.

    Operates on ``os.environ`` by default, or on any mapping (a child
    process's env dict). Prepends rather than appends so the tree
    ``AMBERHOME`` names shadows any stale Amber elsewhere on ``PATH``:
    ``dat/leap`` must match the binaries. Idempotent. Returns the
    ``AMBERHOME`` in effect, or None when no Amber tree can be found, in
    which case the environment is left untouched.
    """
    env = os.environ if env is None else env
    root = find_amber_home(env)
    if root is None:
        return None
    env["AMBERHOME"] = str(root)
    bindir = str(root / "bin")
    parts = [p for p in env.get("PATH", "").split(os.pathsep) if p]
    if not parts or parts[0] != bindir:
        parts = [bindir] + [p for p in parts if p != bindir]
        env["PATH"] = os.pathsep.join(parts)
    return env["AMBERHOME"]
