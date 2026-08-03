"""Cluster profile system for portable SLURM configuration.

Three-tier resolution: project (./.proprep/cluster_profiles/) → user
(~/.proprep/cluster_profiles/) → bundled (package data). User profiles layer
field-by-field over bundled profiles of the same name, so shipped templates
stay pristine while users fill in site-specific values (accounts, QoS, etc.).
"""

from __future__ import annotations

import json
import os
import shutil
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Any, Dict, List, Optional

from proprep.utils.paths import get_package_dir


SCHEMA_VERSION = "1"


def parse_slurm_time_to_seconds(value: str) -> Optional[int]:
    """Parse a SLURM wall-time string into seconds. Returns None if unparseable.

    Accepts `HH:MM:SS`, `D-HH:MM:SS`, `D-HH:MM`, and `D-HH`. Other SLURM formats
    (bare minutes, `MM:SS`) are ambiguous at the class-boundary check and are
    intentionally not supported — the prompt enforces `HH:MM:SS` already.
    """
    if not value:
        return None
    s = value.strip()
    days = 0
    if "-" in s:
        day_part, _, rest = s.partition("-")
        try:
            days = int(day_part)
        except ValueError:
            return None
        s = rest
    parts = s.split(":")
    try:
        nums = [int(p) for p in parts]
    except ValueError:
        return None
    if len(nums) == 3:
        h, m, sec = nums
    elif len(nums) == 2:
        h, m, sec = nums[0], nums[1], 0
    elif len(nums) == 1:
        h, m, sec = nums[0], 0, 0
    else:
        return None
    return days * 86400 + h * 3600 + m * 60 + sec


@dataclass
class ClusterProfile:
    """Hardware palette + cluster-wide SLURM configuration."""
    name: str
    display_name: str = ""
    description: str = ""
    schema_version: str = SCHEMA_VERSION
    required_user_fields: List[str] = field(default_factory=list)
    cluster: Dict[str, Any] = field(default_factory=dict)
    resource_classes: Dict[str, Dict[str, Any]] = field(default_factory=dict)
    defaults: Dict[str, str] = field(default_factory=dict)
    source: str = ""  # "bundled", "user", "project" — set at load time

    def to_json(self) -> Dict[str, Any]:
        data = asdict(self)
        data.pop("source", None)
        return data

    @classmethod
    def from_json(cls, data: Dict[str, Any], source: str = "") -> "ClusterProfile":
        known = {f for f in cls.__dataclass_fields__ if f != "source"}
        filtered = {k: v for k, v in data.items() if k in known}
        prof = cls(**filtered)
        prof.source = source
        return prof

    def resolve_class(self, class_name: str) -> Dict[str, Any]:
        """Return a copy of the resource class entry, with name embedded."""
        if class_name not in self.resource_classes:
            raise KeyError(
                f"Resource class '{class_name}' not in profile '{self.name}'. "
                f"Available: {list(self.resource_classes)}"
            )
        entry = dict(self.resource_classes[class_name])
        entry["class_name"] = class_name
        return entry

    def missing_required_fields(self) -> List[str]:
        missing = []
        for field_path in self.required_user_fields:
            if not _get_by_path(self.cluster, field_path):
                missing.append(field_path)
        return missing


def _get_by_path(data: Dict[str, Any], path: str) -> Any:
    """Dotted path lookup; returns falsy if absent or empty string."""
    cur: Any = data
    for part in path.split("."):
        if not isinstance(cur, dict) or part not in cur:
            return None
        cur = cur[part]
    return cur


def bundled_dir() -> Path:
    return get_package_dir() / "md_prep" / "cluster_profiles"


def user_dir() -> Path:
    return Path.home() / ".proprep" / "cluster_profiles"


def project_dir() -> Path:
    return Path.cwd() / ".proprep" / "cluster_profiles"


def _read_json(path: Path) -> Dict[str, Any]:
    with open(path, "r") as fh:
        return json.load(fh)


def _deep_merge(base: Dict[str, Any], override: Dict[str, Any]) -> Dict[str, Any]:
    """Field-level merge. Override wins for leaves; dicts recurse.

    Lists and non-dict values are replaced wholesale (user fully overrides).
    """
    result = dict(base)
    for key, value in override.items():
        if (
            key in result
            and isinstance(result[key], dict)
            and isinstance(value, dict)
        ):
            result[key] = _deep_merge(result[key], value)
        else:
            result[key] = value
    return result


def list_profiles() -> List[Dict[str, str]]:
    """Return metadata for every discoverable profile.

    Each entry: {"name", "display_name", "source", "path"}.
    If a profile exists in multiple scopes, all are returned so users can
    inspect the layering; load_profile merges them.
    """
    entries: List[Dict[str, str]] = []
    seen_names = set()
    for scope, directory in (
        ("project", project_dir()),
        ("user", user_dir()),
        ("bundled", bundled_dir()),
    ):
        if not directory.exists():
            continue
        for path in sorted(directory.glob("*.json")):
            if path.name.startswith("_") and scope != "bundled":
                continue
            try:
                data = _read_json(path)
            except Exception:
                continue
            name = data.get("name") or path.stem
            entries.append({
                "name": name,
                "display_name": data.get("display_name", name),
                "source": scope,
                "path": str(path),
            })
            seen_names.add(name)
    return entries


def load_profile(name: str) -> ClusterProfile:
    """Resolve a profile by name with project > user > bundled merge.

    Raises FileNotFoundError if the name is absent from all scopes.
    """
    layers: List[Dict[str, Any]] = []
    merged_source = ""
    for scope, directory in (
        ("bundled", bundled_dir()),
        ("user", user_dir()),
        ("project", project_dir()),
    ):
        path = directory / f"{name}.json"
        if path.exists():
            layers.append(_read_json(path))
            merged_source = scope  # last wins for display
    if not layers:
        raise FileNotFoundError(f"Cluster profile '{name}' not found")

    merged: Dict[str, Any] = {}
    for layer in layers:
        merged = _deep_merge(merged, layer)
    merged["name"] = name
    return ClusterProfile.from_json(merged, source=merged_source)


def save_profile(profile: ClusterProfile, scope: str = "user") -> Path:
    """Write a profile to user or project scope. Never writes to bundled."""
    if scope == "user":
        directory = user_dir()
    elif scope == "project":
        directory = project_dir()
    else:
        raise ValueError(f"Invalid scope '{scope}' (expected 'user' or 'project')")
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / f"{profile.name}.json"
    with open(path, "w") as fh:
        json.dump(profile.to_json(), fh, indent=2)
    return path


def delete_profile(name: str, scope: str = "user") -> bool:
    """Remove a profile from user or project scope. Returns True if removed."""
    if scope == "user":
        path = user_dir() / f"{name}.json"
    elif scope == "project":
        path = project_dir() / f"{name}.json"
    else:
        raise ValueError(f"Invalid scope '{scope}'")
    if path.exists():
        path.unlink()
        return True
    return False


def export_profile(name: str, dest: Path) -> Path:
    """Copy a resolved profile to an external path for sharing."""
    profile = load_profile(name)
    dest = Path(dest)
    if dest.is_dir():
        dest = dest / f"{profile.name}.json"
    with open(dest, "w") as fh:
        json.dump(profile.to_json(), fh, indent=2)
    return dest


def import_profile(path: Path, scope: str = "user", rename: Optional[str] = None) -> Path:
    """Copy an external profile JSON into user or project scope."""
    path = Path(path)
    data = _read_json(path)
    if rename:
        data["name"] = rename
    profile = ClusterProfile.from_json(data)
    return save_profile(profile, scope=scope)


def validate(profile: ClusterProfile) -> List[str]:
    """Return a list of human-readable errors; empty list means valid."""
    errors: List[str] = []
    if not profile.name or profile.name.startswith("_"):
        errors.append(f"Invalid profile name: '{profile.name}'")
    if profile.schema_version != SCHEMA_VERSION:
        errors.append(
            f"Unsupported schema_version '{profile.schema_version}' "
            f"(expected '{SCHEMA_VERSION}')"
        )
    if not profile.resource_classes:
        errors.append("Profile has no resource_classes defined")
    for cname, cdef in profile.resource_classes.items():
        mode = cdef.get("mode")
        if mode not in ("cpu", "gpu"):
            errors.append(f"Class '{cname}': mode must be 'cpu' or 'gpu'")
        if not cdef.get("partition"):
            errors.append(f"Class '{cname}': partition is required")
        if mode == "gpu" and not cdef.get("gres"):
            errors.append(f"Class '{cname}': gpu mode requires 'gres'")
        if not cdef.get("binary"):
            errors.append(f"Class '{cname}': binary is required")
    for role in ("cpu_class", "gpu_class"):
        cname = profile.defaults.get(role, "")
        if cname and cname not in profile.resource_classes:
            errors.append(f"defaults.{role} references unknown class '{cname}'")
    missing = profile.missing_required_fields()
    if missing:
        errors.append(f"Missing required user fields: {', '.join(missing)}")
    return errors
