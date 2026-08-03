"""Run plans: per-step resource assignments binding a protocol to a cluster profile.

A run plan is the third object in the SLURM configuration system. It stores,
for each step in a workflow preset, which resource class from a cluster
profile to use plus any overrides (wall time, memory, gpu count). Plans are
cluster-bound; loading a plan whose cluster no longer matches the active
profile triggers a remap wizard.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field, asdict
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

from proprep.md_prep.cluster_profile import ClusterProfile


SCHEMA_VERSION = "1"


@dataclass
class StepResource:
    """Per-step assignment: resource class + optional wall time + overrides.

    amber_flags is an explicit list of extra flags/values appended to the
    AMBER command line for this step (e.g. ``["-ref", "prior.rst7"]`` to
    pass reference coordinates only when there's a prior stage).
    """
    class_name: str
    time_limit: str = ""
    overrides: Dict[str, Any] = field(default_factory=dict)
    amber_flags: List[Any] = field(default_factory=list)

    def to_json(self) -> Dict[str, Any]:
        out: Dict[str, Any] = {"class": self.class_name}
        if self.time_limit:
            out["time_limit"] = self.time_limit
        if self.overrides:
            out["overrides"] = dict(self.overrides)
        if self.amber_flags:
            out["amber_flags"] = list(self.amber_flags)
        return out

    @classmethod
    def from_json(cls, data: Dict[str, Any]) -> "StepResource":
        return cls(
            class_name=data.get("class") or data.get("class_name", ""),
            time_limit=data.get("time_limit", ""),
            overrides=dict(data.get("overrides", {})),
            amber_flags=list(data.get("amber_flags", [])),
        )


@dataclass
class RunPlan:
    """Protocol + cluster + per-step assignments."""
    name: str
    protocol_name: str
    cluster_name: str
    protocol_version: str = ""
    description: str = ""
    schema_version: str = SCHEMA_VERSION
    created: str = ""
    step_resources: Dict[str, StepResource] = field(default_factory=dict)
    source: str = ""  # user / project, set at load

    def to_json(self) -> Dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "name": self.name,
            "protocol_name": self.protocol_name,
            "protocol_version": self.protocol_version,
            "cluster_name": self.cluster_name,
            "description": self.description,
            "created": self.created or datetime.now().isoformat(),
            "step_resources": {
                step_id: sr.to_json() for step_id, sr in self.step_resources.items()
            },
        }

    @classmethod
    def from_json(cls, data: Dict[str, Any], source: str = "") -> "RunPlan":
        step_resources = {
            step_id: StepResource.from_json(sr)
            for step_id, sr in data.get("step_resources", {}).items()
        }
        plan = cls(
            name=data.get("name", ""),
            protocol_name=data.get("protocol_name", ""),
            protocol_version=data.get("protocol_version", ""),
            cluster_name=data.get("cluster_name", ""),
            description=data.get("description", ""),
            schema_version=data.get("schema_version", SCHEMA_VERSION),
            created=data.get("created", ""),
            step_resources=step_resources,
        )
        plan.source = source
        return plan


def user_dir() -> Path:
    return Path.home() / ".proprep" / "run_plans"


def project_dir() -> Path:
    return Path.cwd() / ".proprep" / "run_plans"


def _read_json(path: Path) -> Dict[str, Any]:
    with open(path, "r") as fh:
        return json.load(fh)


def list_plans(protocol_name: Optional[str] = None) -> List[Dict[str, str]]:
    """Return metadata for all run plans, optionally filtered by protocol.

    Entry fields: name, protocol_name, cluster_name, source, path.
    """
    entries: List[Dict[str, str]] = []
    for scope, directory in (("project", project_dir()), ("user", user_dir())):
        if not directory.exists():
            continue
        for path in sorted(directory.glob("*.json")):
            try:
                data = _read_json(path)
            except Exception:
                continue
            if protocol_name and data.get("protocol_name") != protocol_name:
                continue
            entries.append({
                "name": data.get("name") or path.stem,
                "protocol_name": data.get("protocol_name", ""),
                "cluster_name": data.get("cluster_name", ""),
                "source": scope,
                "path": str(path),
            })
    return entries


def load_plan(name: str) -> RunPlan:
    """Load a run plan by name, project scope taking precedence over user."""
    for scope, directory in (("project", project_dir()), ("user", user_dir())):
        path = directory / f"{name}.json"
        if path.exists():
            data = _read_json(path)
            return RunPlan.from_json(data, source=scope)
    raise FileNotFoundError(f"Run plan '{name}' not found")


def save_plan(plan: RunPlan, scope: str = "user") -> Path:
    if scope == "user":
        directory = user_dir()
    elif scope == "project":
        directory = project_dir()
    else:
        raise ValueError(f"Invalid scope '{scope}'")
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / f"{plan.name}.json"
    with open(path, "w") as fh:
        json.dump(plan.to_json(), fh, indent=2)
    return path


def delete_plan(name: str, scope: str = "user") -> bool:
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


def export_plan(name: str, dest: Path) -> Path:
    plan = load_plan(name)
    dest = Path(dest)
    if dest.is_dir():
        dest = dest / f"{plan.name}.json"
    with open(dest, "w") as fh:
        json.dump(plan.to_json(), fh, indent=2)
    return dest


def import_plan(path: Path, scope: str = "user", rename: Optional[str] = None) -> Path:
    path = Path(path)
    data = _read_json(path)
    if rename:
        data["name"] = rename
    plan = RunPlan.from_json(data)
    return save_plan(plan, scope=scope)


def remap_for_cluster(
    plan: RunPlan,
    new_cluster: ClusterProfile,
    classifier_suggestions: Optional[Dict[str, str]] = None,
) -> List[Dict[str, Any]]:
    """Analyze a plan against a new cluster; return per-step remap decisions.

    Does not mutate; returns a list of decisions the wizard can present to
    the user. Each decision:
        {step_id, current_class, status, suggested_class, reason}
    status is one of: 'kept' (class name exists in new cluster),
    'auto-suggest' (defaults used), 'needs-pick' (no good default).
    """
    decisions: List[Dict[str, Any]] = []
    classifier_suggestions = classifier_suggestions or {}
    cpu_default = new_cluster.defaults.get("cpu_class", "")
    gpu_default = new_cluster.defaults.get("gpu_class", "")

    for step_id, sr in plan.step_resources.items():
        if sr.class_name in new_cluster.resource_classes:
            decisions.append({
                "step_id": step_id,
                "current_class": sr.class_name,
                "status": "kept",
                "suggested_class": sr.class_name,
                "reason": "class exists in new cluster",
            })
            continue
        suggested = classifier_suggestions.get(step_id, "")
        if not suggested:
            current_mode = _guess_mode_from_class_name(sr.class_name)
            suggested = gpu_default if current_mode == "gpu" else cpu_default
        status = "auto-suggest" if suggested else "needs-pick"
        decisions.append({
            "step_id": step_id,
            "current_class": sr.class_name,
            "status": status,
            "suggested_class": suggested,
            "reason": (
                "class missing; defaulting via mode"
                if status == "auto-suggest"
                else "class missing; no default available"
            ),
        })
    return decisions


def _guess_mode_from_class_name(class_name: str) -> str:
    lower = class_name.lower()
    if lower.startswith("gpu") or "gpu" in lower:
        return "gpu"
    if lower.startswith("cpu") or "cpu" in lower or "mpi" in lower:
        return "cpu"
    return "cpu"


def apply_remap(plan: RunPlan, new_cluster: ClusterProfile, decisions: List[Dict[str, Any]]) -> RunPlan:
    """Produce a new RunPlan bound to new_cluster with classes replaced per decisions."""
    new_steps: Dict[str, StepResource] = {}
    for d in decisions:
        step_id = d["step_id"]
        original = plan.step_resources[step_id]
        new_steps[step_id] = StepResource(
            class_name=d["suggested_class"],
            time_limit=original.time_limit,
            overrides=dict(original.overrides),
            amber_flags=list(original.amber_flags),
        )
    return RunPlan(
        name=f"{plan.name}__{new_cluster.name}",
        protocol_name=plan.protocol_name,
        protocol_version=plan.protocol_version,
        cluster_name=new_cluster.name,
        description=f"Remapped from '{plan.name}' onto cluster '{new_cluster.name}'",
        step_resources=new_steps,
    )
