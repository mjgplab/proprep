"""
Remove what AmberTools' bundled ProPrep leaves inside the installed package.

AmberTools 26 (``ambertools-dac``) ships its own, older ProPrep (1.0.0) in the
same ``site-packages/proprep`` folder as the standalone package. Installing
ProPrep over it (the install script's force-reinstall, the installers, the
AmberTools updater) replaces every file the two share, but leaves the files
only the old copy had: on 1.22.0 they were 35, among them a retired built-in
MD workflow the MD Manager still offered, old MD templates, a backup module
and old heme and Fe4S4 parameter files.

conda records in ``conda-meta/proprep-<version>-<build>.json`` every file the
installed ProPrep package owns. Anything else under the package's folder came
from elsewhere and is removed here. Nothing is removed when conda has no
record of ProPrep (a pip or editable install, the AmberTools build itself).

Run after installing: ``python -m proprep.utils.bundled_copy`` (``--dry-run``
lists without removing).
"""

from __future__ import annotations

import glob
import json
import os
import sys
from typing import List, Optional, Tuple


def _package_dir() -> str:
    import proprep
    return os.path.dirname(os.path.abspath(proprep.__file__))


def _conda_record(prefix: str) -> Optional[dict]:
    """conda's record of the installed proprep package in ``prefix``, or None."""
    records = [path for path in glob.glob(os.path.join(prefix, "conda-meta", "proprep-*.json"))
               if os.path.basename(path).split("-")[1][:1].isdigit()]       # not proprep-web etc.
    if len(records) != 1:
        return None
    with open(records[0]) as handle:
        return json.load(handle)


def files_not_in_package(prefix: Optional[str] = None,
                         package_dir: Optional[str] = None) -> Tuple[Optional[str], List[str]]:
    """(the package's conda record name, files under the package folder it does not own)."""
    prefix = os.path.abspath(prefix or sys.prefix)
    package_dir = os.path.abspath(package_dir or _package_dir())
    record = _conda_record(prefix)
    if record is None:
        return None, []
    owned = {os.path.normpath(os.path.join(prefix, f)) for f in record.get("files", [])}
    strangers = []
    for folder, dirs, files in os.walk(package_dir):
        dirs[:] = [d for d in dirs if d != "__pycache__"]
        for name in files:
            path = os.path.normpath(os.path.join(folder, name))
            if path not in owned:
                strangers.append(path)
    return f"{record.get('name')}-{record.get('version')}-{record.get('build')}", sorted(strangers)


def remove_files_not_in_package(prefix: Optional[str] = None, package_dir: Optional[str] = None,
                                dry_run: bool = False) -> Tuple[Optional[str], List[str]]:
    """Remove the files under the package folder that the installed package does not own.

    Folders left empty are removed too. Returns what ``files_not_in_package``
    found, which is what was (or, with ``dry_run``, would be) removed.
    """
    package_dir = os.path.abspath(package_dir or _package_dir())
    record, strangers = files_not_in_package(prefix, package_dir)
    if dry_run:
        return record, strangers
    for path in strangers:
        os.remove(path)
    for folder, dirs, files in os.walk(package_dir, topdown=False):
        if folder != package_dir and not os.listdir(folder):
            os.rmdir(folder)
    return record, strangers


def main(argv: Optional[List[str]] = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    dry_run = "--dry-run" in argv
    record, strangers = remove_files_not_in_package(dry_run=dry_run)
    if record is None:
        print("ProPrep is not installed as a conda package here; nothing to check.")
        return 0
    if not strangers:
        print(f"Every file under the ProPrep package belongs to {record}.")
        return 0
    verb = "Would remove" if dry_run else "Removed"
    print(f"{verb} {len(strangers)} file(s) left in the ProPrep package by another copy "
          f"(not part of {record}):")
    for path in strangers:
        print(f"  {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
