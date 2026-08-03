#!/usr/bin/env python3
"""Serialize a v3 transformer authoring session into a JSON transformer spec.

The v3 creator (:mod:`transformer_creator_v3`) captures a `TransformerSpecV3`
(roles, relationship phrases, editing passes, parameters) and today emits a
Python module via a code generator. This module instead translates that same
captured state into the executor-aligned JSON spec that :class:`SpecTransformer`
interprets, so authoring produces reusable *data* (no generated code) that the
Redox Site Preparer discovers from ``~/.proprep/transformers/*.json``.

The translation is pure (no I/O). It:
  * materializes each role's transferable identity — an element connectivity
    fingerprint (reused from :mod:`auto_rename`) and atom-connectivity
    discriminators read from the detected site's bonds;
  * maps each v3 ``EditCommand`` to an executor ``op``/``selector``/``action``;
  * re-keys the v3 ``ParameterSpec.name_mappings`` into per-operation
    ``name_by_state`` plus a ``parameters``/``parameter_order``/``reference_state``
    block.
"""

from __future__ import annotations

import itertools
import logging
from typing import Any, Dict, List, Optional, Tuple

from proprep.redoxsite_prep.transformation.auto_rename import connectivity_signature
from proprep.redoxsite_prep.transformation.spec_transformer import (
    residue_bond_discriminators,
)
from proprep.redoxsite_prep.transformation.transformer_creator_v3.data_models import (
    CommandType,
)

logger = logging.getLogger(__name__)

SCHEMA_VERSION = "0.3"


def _role_key(label: str, roles: Dict[str, Any],
              resolved: Dict[str, Tuple[str, int]]) -> Optional[Tuple[str, int]]:
    if label in resolved and resolved[label]:
        return resolved[label]
    role = roles.get(label)
    if role is not None and getattr(role, "chain", "") and getattr(role, "resid", 0):
        return (role.chain, role.resid)
    return None


def _created_labels(spec) -> set:
    created: set = set()
    for ep in spec.editing_passes:
        for cmd in ep.commands:
            if cmd.command_type == CommandType.MOVE_ATOMS_NEW and cmd.new_role_label:
                created.add(cmd.new_role_label)
    return created


def _state_combos(param_spec) -> Tuple[List[Tuple[str, str]], List[str],
                                       Dict[str, str], str]:
    """Return (combos, parameter_order, reference_state, reference_my_key).

    combos is a list of (my_key '|'-joined, v3_key '_'-joined) over the Cartesian
    product of the choice parameters' options, in declared order — regenerated
    identically to how the v3 creator built name_mappings' keys.
    """
    if param_spec is None:
        return [], [], {}, ""
    choice = [p for p in param_spec.parameters if p.param_type == "choice"]
    order = [p.name for p in choice]
    option_lists = [p.options for p in choice]
    combos: List[Tuple[str, str]] = []
    for combo in itertools.product(*option_lists):
        combos.append(("|".join(combo), "_".join(combo)))
    reference = {p.name: (p.default or (p.options[0] if p.options else ""))
                 for p in choice}
    ref_my_key = "|".join(reference[n] for n in order)
    return combos, order, reference, ref_my_key


def _name_by_state(role_label: str, param_spec,
                   combos: List[Tuple[str, str]]) -> Dict[str, str]:
    """Per-state target names for a state-dependent role, keyed '|'-joined."""
    if param_spec is None or role_label not in (param_spec.state_dependent_roles or []):
        return {}
    nbs: Dict[str, str] = {}
    for my_key, v3_key in combos:
        rn = (param_spec.name_mappings.get(v3_key, {}) or {}).get(role_label)
        if rn:
            nbs[my_key] = rn
    return nbs


def _translate_command(cmd, roles: Dict[str, Any],
                       resolved: Dict[str, Tuple[str, int]],
                       param_spec, combos, ref_my_key) -> Optional[Dict[str, Any]]:
    ct = cmd.command_type
    role = cmd.source_role

    if ct == CommandType.RENAME_RESIDUE:
        op: Dict[str, Any] = {"op": "rename_residue", "selector": {"role": role},
                              "action": {"change_residue_name": cmd.new_resname}}
        nbs = _name_by_state(role, param_spec, combos)
        if nbs:
            op["name_by_state"] = nbs
            op["action"]["change_residue_name"] = nbs.get(ref_my_key, cmd.new_resname)
        return op

    if ct == CommandType.RENAME_ATOMS:
        mapping = dict(zip(cmd.old_atom_names or [], cmd.new_atom_names or []))
        return {"op": "rename_atom", "selector": {"role": role},
                "action": {"rename_atoms": mapping}}

    if ct == CommandType.HETATM:
        return {"op": "set_record_type", "selector": {"role": role},
                "action": {"convert_to_hetatm": True}}

    if ct == CommandType.ATOM:
        return {"op": "set_record_type", "selector": {"role": role},
                "action": {"convert_to_atom": True}}

    if ct == CommandType.MOVE_ATOMS:
        target = cmd.target_role
        target_rn = getattr(roles.get(target), "resname", None)
        op = {"op": "move_to_existing_residue",
              "selector": {"role": role, "atom_names": list(cmd.atom_names or [])},
              "target_role": target, "action": {}}
        if target_rn:
            op["action"]["change_residue_name"] = target_rn
        return op

    if ct == CommandType.MOVE_ATOMS_NEW:
        new_label = cmd.new_role_label
        src_key = _role_key(role, roles, resolved)
        new_role = roles.get(new_label)
        offset = 1
        if new_role is not None and src_key and getattr(new_role, "resid", 0):
            offset = new_role.resid - src_key[1]
            if offset < 1:
                offset = 1
        return {"op": "move_to_new_residue",
                "selector": {"role": role, "atom_names": list(cmd.atom_names or [])},
                "creates_role": new_label,
                "relative": {"new_resid_offset": offset},
                "action": {"change_residue_name": cmd.new_residue_name}}

    logger.warning("Unknown v3 command type %r; skipped in JSON spec.", ct)
    return None


def build_json_spec(spec, roles: Dict[str, Any],
                    resolved: Dict[str, Tuple[str, int]], site) -> Dict[str, Any]:
    """Translate a captured `TransformerSpecV3` + roles + resolved + detected
    site into the JSON transformer spec consumed by :class:`SpecTransformer`."""
    created = _created_labels(spec)
    sig = connectivity_signature(site) if site is not None else {}

    role_entries: List[Dict[str, Any]] = []
    for label, role in roles.items():
        entry: Dict[str, Any] = {"label": label, "resname": role.resname}
        if label in created:
            entry["created"] = True
        else:
            key = _role_key(label, roles, resolved)
            entry["fingerprint"] = sig.get(key) if key else None
            entry["discriminators"] = [
                {"source_atom": s, "partner_resname": p, "target_atom": t}
                for (s, p, t) in sorted(residue_bond_discriminators(site, key))
            ] if (site is not None and key) else []
        role_entries.append(entry)

    param_spec = spec.parameters
    combos, order, reference, ref_my_key = _state_combos(param_spec)

    operations: List[Dict[str, Any]] = []
    for ep in spec.editing_passes:
        for cmd in ep.commands:
            op = _translate_command(cmd, roles, resolved, param_spec, combos, ref_my_key)
            if op is not None:
                operations.append(op)

    out: Dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "name": spec.name,
        "description": spec.description or "",
        "source": "transformer_creator_v3",
        "roles": role_entries,
        "operations": operations,
    }

    if order:
        out["parameters"] = [
            {"name": p.name, "options": list(p.options),
             "default": p.default or (p.options[0] if p.options else "")}
            for p in param_spec.parameters if p.param_type == "choice"
        ]
        out["parameter_order"] = order
        out["reference_state"] = reference

    # Forcefield block: path from the spec, plus any fixed redox/spin params.
    ff: Dict[str, Any] = {}
    if getattr(spec, "forcefield_path", None):
        ff["path"] = spec.forcefield_path
    if param_spec is not None:
        for p in param_spec.parameters:
            if p.param_type == "fixed" and p.name in ("redox_state", "spin_state"):
                ff[p.name] = p.fixed_value
    if ff:
        out["forcefield"] = ff

    return out
