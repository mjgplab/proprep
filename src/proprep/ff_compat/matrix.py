"""
Build a pairwise compatibility matrix across all discoverable FF sets.

Two FF sets are *incompatible* when they both declare an entry under the
same canonical key (same atom-type name, same canonical bond/angle/
dihedral/improper key) but with different numeric values. Shared keys
with identical values are harmless duplication and are not flagged.

Discovery walks bundled FF sets at
``forcefield_params/specialized_residues/`` and (optionally) user-derived FF
sets at ``~/.proprep/forcefield_params/specialized_residues/`` via the
existing loader API.

This module is intended to be invoked via the CLI in
``proprep.ff_compat.build_matrix``; the public function
``build_matrix`` is also importable for programmatic use (e.g., by the
later picker integration).
"""

from __future__ import annotations

import datetime as _dt
import json
import math
from dataclasses import asdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

from proprep.forcefield_params.loader import (
    discover_forcefield_files,
    get_available_cofactor_types,
    get_forcefield_base_path,
    get_user_forcefield_base_path,
    load_forcefield_metadata,
)
from proprep.ff_compat.parser import FFSignature, parse_set


SCHEMA_VERSION = 1


# ---------------------------------------------------------------------------
# Discovery
# ---------------------------------------------------------------------------

def discover_all_sets(include_user: bool = True) -> List[Dict]:
    """Walk metadata.json files; return one descriptor per forcefield_set.

    Each descriptor is a dict with:
        set_id          — "<cofactor_path>:<set_name>"
        cofactor_path   — e.g. "heme/bis_his_c_type"
        set_name        — e.g. "Henriques_Reduced"
        redox_state, spin_state, source
        frcmod          — absolute path
        libs            — list of absolute paths
    """
    descriptors: List[Dict] = []
    cofactors = get_available_cofactor_types()

    for cofactor_path, info in cofactors.items():
        if not info.get("valid"):
            continue
        # Filter user/bundled split if requested
        is_user_path = str(get_user_forcefield_base_path()) in info.get("path", "")
        if not include_user and is_user_path:
            continue

        metadata = load_forcefield_metadata(cofactor_path)
        redox_states = metadata.get("redox_states", {})
        for redox_state, redox_info in redox_states.items():
            for spin_state, _spin_info in redox_info.get("spin_states", {}).items():
                files_list = discover_forcefield_files(cofactor_path, redox_state, spin_state)
                for fs in files_list:
                    set_name = fs["name"]
                    libs = fs["lib"] if isinstance(fs["lib"], list) else [fs["lib"]]
                    # set_id must uniquely address the (cofactor, redox, spin, set)
                    # tuple. Same set_name (e.g., "AMBER_Library_Bonded") recurs
                    # across redox states within a cofactor; collapsing on
                    # cofactor_path:set_name alone loses real sets.
                    set_id = f"{cofactor_path}:{redox_state}:{spin_state}:{set_name}"
                    descriptors.append({
                        "set_id": set_id,
                        "cofactor_path": cofactor_path,
                        "set_name": set_name,
                        "redox_state": redox_state,
                        "spin_state": spin_state,
                        "source": "user" if is_user_path else "bundled",
                        "frcmod": fs["frcmod"],
                        "libs": libs,
                    })
    return descriptors


# ---------------------------------------------------------------------------
# Pairwise comparison
# ---------------------------------------------------------------------------

# Two emitters legitimately round the SAME parameter to different precision
# (e.g. antechamber writes a bond force constant as 354.25 while an
# MCPB-extracted frcmod writes 354.2). That is harmless duplication, not a
# conflict, so numeric values are compared with a tolerance. Real conflicts
# (different force constants, a 0-vs-180 phase, a 2-vs-3 periodicity) differ by
# far more than one rounded decimal and are still caught.
_VAL_REL_TOL = 1e-3
_VAL_ABS_TOL = 0.06   # covers rounding to one decimal (max error 0.05)


def _values_equivalent(av, bv) -> bool:
    """True if two parameter values are the same to within rounding precision."""
    # Non-numeric values (e.g. library unit-type names) must match exactly.
    if isinstance(av, str) or isinstance(bv, str):
        return av == bv
    # Scalars (MASS).
    if isinstance(av, (int, float)) and isinstance(bv, (int, float)):
        return math.isclose(float(av), float(bv),
                            rel_tol=_VAL_REL_TOL, abs_tol=_VAL_ABS_TOL)
    # Tuples/lists of scalars (NONB/BOND/ANGL/DIHE/IMPR).
    try:
        if len(av) != len(bv):
            return False
        return all(
            math.isclose(float(x), float(y), rel_tol=_VAL_REL_TOL, abs_tol=_VAL_ABS_TOL)
            for x, y in zip(av, bv)
        )
    except (TypeError, ValueError):
        return av == bv


def _dict_diffs(a: Dict, b: Dict) -> List[Tuple]:
    """Return (key, a_value, b_value) tuples for keys shared with different values."""
    diffs = []
    for k in set(a) & set(b):
        if not _values_equivalent(a[k], b[k]):
            diffs.append((k, a[k], b[k]))
    return diffs


def compare_signatures(a: FFSignature, b: FFSignature) -> Dict:
    """Compare two FFSignatures. Returns a conflict report dict (empty if compatible)."""
    entries: List[Dict] = []

    # MASS
    for k, av, bv in _dict_diffs(a.mass, b.mass):
        entries.append({"section": "MASS", "key": k, "a": av, "b": bv})

    # NONB
    for k, av, bv in _dict_diffs(a.nonb, b.nonb):
        entries.append({"section": "NONB", "key": list(k) if isinstance(k, tuple) else k,
                        "a": list(av), "b": list(bv)})

    # BOND
    for k, av, bv in _dict_diffs(a.bonds, b.bonds):
        entries.append({"section": "BOND", "key": list(k), "a": list(av), "b": list(bv)})

    # ANGL
    for k, av, bv in _dict_diffs(a.angles, b.angles):
        entries.append({"section": "ANGL", "key": list(k), "a": list(av), "b": list(bv)})

    # DIHE — keys are (canonical_quad, periodicity_int)
    for k, av, bv in _dict_diffs(a.dihedrals, b.dihedrals):
        quad, period = k
        entries.append({"section": "DIHE", "key": list(quad), "periodicity": period,
                        "a": list(av), "b": list(bv)})

    # IMPR
    for k, av, bv in _dict_diffs(a.impropers, b.impropers):
        quad, period = k
        entries.append({"section": "IMPR", "key": list(quad), "periodicity": period,
                        "a": list(av), "b": list(bv)})

    # Library unit-type assignments
    for k, av, bv in _dict_diffs(a.library_unit_types, b.library_unit_types):
        entries.append({"section": "LIB", "key": list(k), "a": av, "b": bv})

    if not entries:
        return {}

    shared_types = sorted(a.all_types() & b.all_types())
    return {
        "shared_types": shared_types,
        "entries": entries,
    }


# ---------------------------------------------------------------------------
# Matrix assembly
# ---------------------------------------------------------------------------

def build_matrix(include_user: bool = True,
                 verbose: bool = False,
                 detailed: bool = False) -> Dict:
    """Walk every discoverable FF set, parse each, do every pairwise compare.

    Args:
        include_user: also scan ~/.proprep/forcefield_params/.
        verbose: print per-set parse progress to stdout.
        detailed: include the full conflict entry list in each pair's
            report. Default False (ships only counts + shared types,
            keeping the matrix git-reviewable). Set True when called
            from Tool 2 (resolver) where the entry list is needed.

    Returns the matrix as a Python dict (call ``write_matrix`` to serialize).
    """
    descriptors = discover_all_sets(include_user=include_user)
    descriptors.sort(key=lambda d: d["set_id"])

    if verbose:
        print(f"Discovered {len(descriptors)} FF sets")

    signatures: Dict[str, FFSignature] = {}
    parse_errors: List[Tuple[str, str]] = []
    for d in descriptors:
        try:
            signatures[d["set_id"]] = parse_set(d["set_id"], d["frcmod"], d["libs"])
            if verbose:
                print(f"  parsed {d['set_id']}: {len(signatures[d['set_id']].all_types())} types")
        except Exception as e:
            parse_errors.append((d["set_id"], f"{type(e).__name__}: {e}"))
            if verbose:
                print(f"  FAILED to parse {d['set_id']}: {e}")

    set_ids = sorted(signatures)
    conflicts: Dict[str, Dict] = {}
    for i, a_id in enumerate(set_ids):
        for b_id in set_ids[i + 1:]:
            report = compare_signatures(signatures[a_id], signatures[b_id])
            if not report:
                continue
            if detailed:
                conflicts[f"{a_id}|{b_id}"] = report
            else:
                # Summary form: drop the full entry list, keep shape stats.
                # Tool 2 re-computes the entry list on demand.
                from collections import Counter
                by_section = Counter(e["section"] for e in report["entries"])
                conflicts[f"{a_id}|{b_id}"] = {
                    "shared_types": report["shared_types"],
                    "entry_count": len(report["entries"]),
                    "entries_by_section": dict(by_section),
                }

    return {
        "schema_version": SCHEMA_VERSION,
        "generated_at": _dt.datetime.now().isoformat(timespec="seconds"),
        "set_count": len(signatures),
        "detailed": detailed,
        "set_ids": set_ids,
        "parse_errors": [{"set_id": s, "error": e} for s, e in parse_errors],
        "conflicts": conflicts,
    }


def write_matrix(matrix: Dict, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(matrix, f, indent=2, sort_keys=False)


# ---------------------------------------------------------------------------
# Convenience formatters
# ---------------------------------------------------------------------------

def summarize_matrix(matrix: Dict) -> str:
    """Human-readable summary suitable for stdout."""
    lines = [
        f"FF compatibility matrix (schema v{matrix['schema_version']})",
        f"  Generated: {matrix['generated_at']}",
        f"  Sets:      {matrix['set_count']}",
        f"  Conflicts: {len(matrix['conflicts'])} pair(s)",
    ]
    if matrix.get("parse_errors"):
        lines.append(f"  Parse errors: {len(matrix['parse_errors'])}")
        for err in matrix["parse_errors"]:
            lines.append(f"    {err['set_id']}: {err['error']}")
    if matrix["conflicts"]:
        lines.append("")
        lines.append("Conflicting pairs:")
        for pair, report in matrix["conflicts"].items():
            shared = ", ".join(report["shared_types"]) or "(none in name space, conflict via shared keys)"
            # Support both detailed (entries list) and summary (entry_count) form.
            if "entries" in report:
                n_entries = len(report["entries"])
            else:
                n_entries = report.get("entry_count", 0)
            lines.append(f"  {pair}")
            lines.append(f"    shared types: {shared}")
            lines.append(f"    differing entries: {n_entries}")
    return "\n".join(lines)
