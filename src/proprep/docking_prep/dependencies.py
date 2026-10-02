"""
What the docking module needs, and whether this installation has it.

The conda ProPrep environment carries RDKit and gemmi through AmberTools
(``ambertools-dac``), but the AmberTools CMake build installs RDKit only in its
miniconda mode, and neither build installs Meeko or Vina. The module therefore
checks at start-up and names exactly what is missing and what it is needed for,
instead of failing on the first import.

Versions listed as ``tested`` are the ones the docking code was validated
against; a different version is reported, not refused.
"""

from __future__ import annotations

import importlib
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass
from typing import Callable, List, Optional


@dataclass
class Dependency:
    name: str
    needed_for: str
    tested: str
    install_hint: str = ""
    available: bool = False
    version: Optional[str] = None
    detail: str = ""                   # where it was found, or why it could not be used

    def describe(self) -> str:
        if not self.available:
            text = f"{self.name}: MISSING (needed for {self.needed_for}); {self.detail}."
            return f"{text} Install: {self.install_hint}." if self.install_hint else text
        if self.version is None:
            return f"{self.name} (version not reported; validated with {self.tested})"
        note = "" if self.version == self.tested else f" (validated with {self.tested})"
        return f"{self.name} {self.version}{note}"


def _python_module(module: str, version_attr: str = "__version__") -> Callable[[Dependency], None]:
    def check(dep: Dependency) -> None:
        try:
            imported = importlib.import_module(module)
        except Exception as error:                      # ImportError, or a broken binary extension
            dep.detail = f"import {module} failed: {error}"
            return
        dep.available = True
        dep.version = str(getattr(imported, version_attr, "") or "") or None
    return check


def find_executable(name: str) -> Optional[str]:
    """``name`` on PATH, else next to the running Python (a conda env's bin/ is not always on PATH)."""
    found = shutil.which(name)
    if found:
        return found
    candidate = os.path.join(os.path.dirname(os.path.abspath(sys.executable)), name)
    return candidate if os.access(candidate, os.X_OK) else None


def _executable(name: str, version_flag: str, version_prefix: str) -> Callable[[Dependency], None]:
    """Found on PATH or beside Python; version is the first output line after ``version_prefix``."""
    def check(dep: Dependency) -> None:
        path = find_executable(name)
        if path is None:
            dep.detail = f"no {name} on PATH or beside {sys.executable}"
            return
        dep.available = True
        dep.detail = path
        try:
            output = subprocess.run([path, version_flag], capture_output=True, text=True, timeout=30)
        except (OSError, subprocess.TimeoutExpired):
            return
        for line in (output.stdout + output.stderr).splitlines():
            if line.strip().startswith(version_prefix):
                dep.version = line.strip()[len(version_prefix):].strip() or None
                break
    return check


_CHECKS = [
    (Dependency("RDKit", "reading ligands and building 3D structures", "2025.03.6"), _python_module("rdkit")),
    (Dependency("gemmi", "reading Chemical Component Dictionary entries", "0.7.5"), _python_module("gemmi")),
    (Dependency("Meeko", "AutoDock atom typing, torsion trees and PDBQT files", "0.8.0",
                install_hint="conda install -c mjgplab meeko=0.8.0 (ProPrep's build; conda-forge's meeko "
                             "requires prody, which conda-forge cannot install on Apple Silicon with Python "
                             "3.12+, and Meeko does not need it), or "
                             "pip install --no-deps meeko==0.8.0"),
     _python_module("meeko")),
    (Dependency("AutoDock Vina", "docking with the Vina and AD4 scoring functions", "1.2.7",
                install_hint="conda install -c conda-forge vina=1.2.7"), _python_module("vina")),
    # conda-forge's autogrid package 4.2.9 reports itself as "AutoGrid 4.2.7.x"
    (Dependency("autogrid4", "AutoDock4 scoring maps (only for the AD4 scoring function)", "4.2.7.x",
                install_hint="conda install -c conda-forge autogrid"), _executable("autogrid4", "--version", "AutoGrid")),
]


def check_dependencies() -> List[Dependency]:
    """A fresh availability report for every docking dependency."""
    report = []
    for template, check in _CHECKS:
        dep = Dependency(template.name, template.needed_for, template.tested, install_hint=template.install_hint)
        check(dep)
        report.append(dep)
    return report


def missing(report: List[Dependency], *, include_ad4: bool) -> List[Dependency]:
    """Dependencies that block docking; autogrid4 counts only when AD4 scoring is wanted."""
    return [d for d in report if not d.available and (include_ad4 or d.name != "autogrid4")]
