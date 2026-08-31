"""Supply the MODELLER license key at run time instead of install time.

MODELLER reads its key from ``modeller.config``, a two-line module the conda
package writes at install time: its post-link script bakes in ``$KEY_MODELLER``
if that was exported, otherwise leaves the placeholder ``'XXXX'``.
``modeller/__init__.py`` then does ``from modeller import config`` and pushes
``config.license`` and ``config.install_dir`` into the C library once, on
first import.

That makes the key a property of the *installation*. For distributables
(PyInstaller bundle, ``constructor`` installer) we want it to be a property of
the *user*, so that no shipped artifact ever contains a key. This module
answers the ``modeller.config`` import first, from one of:

  1. the ``KEY_MODELLER`` environment variable
  2. ``~/.proprep/modeller_key`` (first non-blank, non-comment line)

Precedence and safety:

* Outside a frozen bundle, if neither source exists **nothing is changed**:
  whatever the install baked in is used, so users of ``install_proprep.sh``
  (which asks for the key at install time) see no difference.
* Inside a PyInstaller bundle the real ``modeller.config`` is deliberately not
  packaged (it would carry the build machine's key), so the stand-in is always
  installed there, with ``'XXXX'`` when no runtime key is found.
* A runtime key always wins over a baked one, so a wrong key can be fixed by
  editing the file rather than reinstalling.

It must run before the first ``import modeller`` in the process. ``main.py``
calls :func:`configure_modeller_license` at import time and the PyInstaller
runtime hook (``hook-modeller-config.py``) calls it before the entry script.
"""

from __future__ import annotations

import glob
import importlib.util
import os
import sys
import types
from pathlib import Path
from typing import Optional, Tuple

KEY_ENV_VAR = "KEY_MODELLER"
KEY_FILE_RELATIVE = Path(".proprep") / "modeller_key"
PLACEHOLDER_KEY = "XXXX"

# Return values of configure_modeller_license(), for logging/tests.
STATUS_CONFIGURED = "configured"
STATUS_NO_RUNTIME_KEY = "no-runtime-key"
STATUS_NO_MODELLER = "no-modeller"
STATUS_ALREADY_IMPORTED = "already-imported"


def key_file_path() -> Path:
    """Location of the per-user key file (computed per call so HOME changes apply)."""
    return Path.home() / KEY_FILE_RELATIVE


def read_runtime_key() -> Optional[str]:
    """Return the runtime key, or None if neither source provides one."""
    key = os.environ.get(KEY_ENV_VAR, "").strip()
    if key:
        return key
    try:
        with open(key_file_path(), encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if line and not line.startswith("#"):
                    return line
    except OSError:
        pass
    return None


def runtime_key_source() -> Optional[str]:
    """Human-readable description of where the runtime key came from, or None."""
    if os.environ.get(KEY_ENV_VAR, "").strip():
        return f"${KEY_ENV_VAR}"
    if read_runtime_key() is not None:
        return str(key_file_path())
    return None


def _is_frozen() -> bool:
    return bool(getattr(sys, "frozen", False)) and hasattr(sys, "_MEIPASS")


def _frozen_install_dir() -> Optional[str]:
    """MODELLER data directory bundled by proprep.spec (``modeller-<ver>/``)."""
    matches = glob.glob(os.path.join(getattr(sys, "_MEIPASS", ""), "modeller-*"))
    return matches[0] if matches else None


def _installed_config() -> Optional[Tuple[Optional[str], Optional[str]]]:
    """Read ``(install_dir, license)`` from the installed ``modeller/config.py``
    WITHOUT importing the modeller package (find_spec does not execute it)."""
    try:
        spec = importlib.util.find_spec("modeller")
    except (ImportError, ValueError):
        return None
    if spec is None or not spec.submodule_search_locations:
        return None
    for location in spec.submodule_search_locations:
        cfg = Path(location) / "config.py"
        if cfg.is_file():
            namespace: dict = {}
            try:
                exec(compile(cfg.read_text(encoding="utf-8"), str(cfg), "exec"), namespace)
            except Exception:
                return None
            return namespace.get("install_dir"), namespace.get("license")
    return None


def _is_our_stand_in(module: object) -> bool:
    return bool(getattr(module, "__proprep_runtime_key__", False))


def configure_modeller_license() -> str:
    """Pre-seed ``sys.modules['modeller.config']`` from the runtime key.

    Returns one of the ``STATUS_*`` constants. Safe to call repeatedly.
    """
    already = sys.modules.get("modeller")
    if already is not None and not _is_our_stand_in(sys.modules.get("modeller.config")):
        # modeller was imported before we ran; its C library has already read
        # whatever config it found. Nothing we set now would take effect.
        return STATUS_ALREADY_IMPORTED

    key = read_runtime_key()
    frozen = _is_frozen()
    if key is None and not frozen:
        return STATUS_NO_RUNTIME_KEY

    if frozen:
        install_dir = _frozen_install_dir()
    else:
        found = _installed_config()
        if found is None:
            return STATUS_NO_MODELLER
        install_dir = found[0]
    if not install_dir:
        return STATUS_NO_MODELLER

    config = types.ModuleType("modeller.config")
    config.license = key or PLACEHOLDER_KEY
    config.install_dir = install_dir
    config.__file__ = __file__
    config.__proprep_runtime_key__ = True
    sys.modules["modeller.config"] = config
    return STATUS_CONFIGURED
