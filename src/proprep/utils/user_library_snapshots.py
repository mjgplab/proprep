"""Snapshot, reset and restore the user library, never destructively.

The user library is the part of ``~/.proprep`` that ProPrep fills as you
prepare systems: deposited force-field parameter sets (``forcefield_params/``)
and saved transformers (``transformers/``). Practice runs for a manuscript
leave it cluttered, and "start clean" used to mean moving directories by
hand. Here every operation that would remove anything first writes a
snapshot under ``~/.proprep/library_snapshots/<name>/``, so any state can be
brought back.

Settings, license keys, templates, cluster profiles, run plans and web
sessions are not part of the library and are never touched.

Command line: ``proprep-library {status,snapshot,list,reset,restore}``.
"""
from __future__ import annotations

import argparse
import datetime as _dt
import json
import os
import re
import shutil
import sys
from pathlib import Path
from typing import Dict, List, Optional

LIBRARY_PARTS = ("forcefield_params", "transformers")
SNAPSHOT_DIRNAME = "library_snapshots"
MANIFEST = "manifest.json"
_IGNORE = shutil.ignore_patterns("__pycache__", ".DS_Store", "*.pyc")


def default_root() -> Path:
    return Path(os.environ.get("PROPREP_USER_DIR") or (Path.home() / ".proprep"))


def _stamp() -> str:
    return _dt.datetime.now().strftime("%Y%m%d_%H%M%S")


def sanitize_name(raw: Optional[str], prefix: str = "snapshot") -> str:
    name = re.sub(r"[^A-Za-z0-9._-]+", "_", (raw or "").strip()).strip("._")
    return name or f"{prefix}_{_stamp()}"


def _count(path: Path) -> Dict[str, int]:
    files = size = 0
    if path.exists():
        for p in path.rglob("*"):
            if p.is_file() and "__pycache__" not in p.parts and p.suffix != ".pyc":
                files += 1
                size += p.stat().st_size
    return {"files": files, "bytes": size}


class UserLibrary:
    def __init__(self, root: Optional[Path] = None):
        self.root = Path(root) if root else default_root()
        self.snapshots_dir = self.root / SNAPSHOT_DIRNAME

    # ---- state ----------------------------------------------------------
    def status(self) -> Dict[str, Dict[str, int]]:
        return {part: _count(self.root / part) for part in LIBRARY_PARTS}

    def is_empty(self) -> bool:
        return all(v["files"] == 0 for v in self.status().values())

    # ---- snapshots ------------------------------------------------------
    def create_snapshot(self, name: Optional[str] = None, note: str = "") -> Path:
        name = sanitize_name(name)
        dest = self.snapshots_dir / name
        if dest.exists():
            raise FileExistsError(f"snapshot '{name}' already exists: {dest}")
        dest.mkdir(parents=True)
        parts: Dict[str, Dict[str, int]] = {}
        for part in LIBRARY_PARTS:
            src = self.root / part
            if src.exists():
                shutil.copytree(src, dest / part, ignore=_IGNORE, symlinks=True)
            parts[part] = _count(dest / part)
        manifest = {
            "name": name,
            "created": _dt.datetime.now().isoformat(timespec="seconds"),
            "note": note,
            "parts": parts,
        }
        try:
            from proprep import __version__ as _v
            manifest["proprep_version"] = _v
        except Exception:
            pass
        (dest / MANIFEST).write_text(json.dumps(manifest, indent=2))
        return dest

    def list_snapshots(self) -> List[Dict]:
        out = []
        if not self.snapshots_dir.exists():
            return out
        for d in sorted(self.snapshots_dir.iterdir()):
            m = d / MANIFEST
            if d.is_dir() and m.exists():
                try:
                    out.append(json.loads(m.read_text()))
                except ValueError:
                    out.append({"name": d.name, "created": "?", "note": "(unreadable manifest)", "parts": {}})
        return out

    def snapshot_path(self, name: str) -> Path:
        p = self.snapshots_dir / name
        if not (p / MANIFEST).exists():
            raise FileNotFoundError(f"no snapshot named '{name}' in {self.snapshots_dir}")
        return p

    # ---- destructive-looking operations, each guarded by a snapshot -------
    def _clear_parts(self) -> None:
        for part in LIBRARY_PARTS:
            target = self.root / part
            if target.exists():
                shutil.rmtree(target)
            target.mkdir(parents=True, exist_ok=True)   # code paths expect the directories

    def reset(self, note: str = "") -> Optional[Path]:
        """Empty the library. The current contents go to a snapshot first
        (unless there is nothing to keep). Returns the snapshot path or None."""
        saved = None
        if not self.is_empty():
            saved = self.create_snapshot(f"before_reset_{_stamp()}", note or "automatic, taken by reset")
        self._clear_parts()
        return saved

    def restore(self, name: str, note: str = "") -> Optional[Path]:
        """Replace the library with snapshot ``name``. The current contents go
        to a snapshot first (unless empty). Returns that safety snapshot."""
        src = self.snapshot_path(name)
        saved = None
        if not self.is_empty():
            saved = self.create_snapshot(f"before_restore_{_stamp()}", note or f"automatic, before restoring '{name}'")
        self._clear_parts()
        for part in LIBRARY_PARTS:
            if (src / part).exists():
                shutil.rmtree(self.root / part)
                shutil.copytree(src / part, self.root / part, ignore=_IGNORE, symlinks=True)
        return saved


# ---- command line -----------------------------------------------------------

def _fmt_parts(parts: Dict[str, Dict[str, int]]) -> str:
    return ", ".join(f"{k}: {v['files']} file{'s' if v['files'] != 1 else ''}" for k, v in parts.items())


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(
        prog="proprep-library",
        description="Snapshot, reset and restore the ProPrep user library "
                    "(~/.proprep/forcefield_params and ~/.proprep/transformers). "
                    "Nothing is ever deleted: reset and restore write a snapshot of the current library first.")
    ap.add_argument("--root", help="user directory (default: $PROPREP_USER_DIR or ~/.proprep)")
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("status", help="what the live library holds")
    p = sub.add_parser("snapshot", help="copy the current library to a named snapshot")
    p.add_argument("name", nargs="?", help="snapshot name (default: timestamp)")
    p.add_argument("-m", "--note", default="", help="a note stored in the manifest")
    sub.add_parser("list", help="list snapshots")
    p = sub.add_parser("reset", help="snapshot the current library, then empty it")
    p.add_argument("-m", "--note", default="")
    p.add_argument("-y", "--yes", action="store_true", help="do not ask for confirmation")
    p = sub.add_parser("restore", help="snapshot the current library, then replace it with a snapshot")
    p.add_argument("name")
    p.add_argument("-y", "--yes", action="store_true")
    args = ap.parse_args(argv)

    lib = UserLibrary(Path(args.root) if args.root else None)
    out = sys.stdout

    if args.cmd == "status":
        st = lib.status()
        print(f"Library root: {lib.root}", file=out)
        for part, c in st.items():
            print(f"  {part:18s} {c['files']:5d} files  {c['bytes'] / 1024:8.1f} KB", file=out)
        snaps = lib.list_snapshots()
        print(f"  snapshots: {len(snaps)} in {lib.snapshots_dir}", file=out)
        return 0

    if args.cmd == "snapshot":
        try:
            dest = lib.create_snapshot(args.name, args.note)
        except FileExistsError as e:
            print(f"error: {e}", file=sys.stderr); return 1
        print(f"Snapshot written: {dest}", file=out)
        return 0

    if args.cmd == "list":
        snaps = lib.list_snapshots()
        if not snaps:
            print(f"No snapshots in {lib.snapshots_dir}", file=out); return 0
        for s in snaps:
            note = f"  — {s['note']}" if s.get("note") else ""
            print(f"  {s['name']:40s} {s.get('created', '?'):20s} {_fmt_parts(s.get('parts', {}))}{note}", file=out)
        return 0

    if args.cmd in ("reset", "restore"):
        st = lib.status()
        what = ("empty the library" if args.cmd == "reset" else f"replace the library with snapshot '{args.name}'")
        if args.cmd == "restore":
            try:
                lib.snapshot_path(args.name)
            except FileNotFoundError as e:
                print(f"error: {e}", file=sys.stderr); return 1
        print(f"Live library: {_fmt_parts(st)}", file=out)
        if not args.yes:
            keep = "" if lib.is_empty() else " (its current contents are snapshotted first)"
            answer = input(f"About to {what}{keep}. Continue? [y/N] ").strip().lower()
            if answer not in ("y", "yes"):
                print("Nothing changed.", file=out); return 0
        saved = lib.reset(args.note) if args.cmd == "reset" else lib.restore(args.name)
        if saved:
            print(f"Previous library saved to: {saved}", file=out)
        print("Library is now empty." if args.cmd == "reset" else f"Library restored from '{args.name}'.", file=out)
        print("Restart ProPrep so open sessions see the change.", file=out)
        return 0
    return 2


if __name__ == "__main__":
    sys.exit(main())
