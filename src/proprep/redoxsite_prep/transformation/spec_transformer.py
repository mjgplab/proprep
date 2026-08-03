#!/usr/bin/env python3
"""Data-driven, transferable RedoxSite transformer built from a JSON spec.

Where :class:`AutoRenameTransformerBase` interprets a rename-only ``RENAMES``
table, :class:`SpecTransformer` interprets a full recipe *spec* — a list of
role-tagged editing operations (residue/atom renames, chain/insertion-code
changes, HETATM conversion, moving atoms into an existing residue) plus,
optionally, state-dependent target names. The spec is pure data; all the
matching and replay logic lives here once, exactly as with the auto-rename base.

Role binding (which residue in a *new* protein fills each role) is a three-tier
ladder, deterministic and geometry-free:

  1. resname                 — a role whose residue name is unique in the site.
  2. element fingerprint      — the existing Weisfeiler-Lehman coordination label
                                (reused from :mod:`auto_rename`), for same-name
                                residues told apart by their metal environment.
  3. atom-connectivity        — the decisive tier for same-name residues that tie
     discriminator             under the (element-only) fingerprint: match by
                                *which atom bonds which partner atom*, read from
                                the site's user-specified bonds (e.g. a c-type
                                heme's two Cys, separated by SG->CAB vs SG->CAC).

If those cannot resolve a role, the transformer declares the same ``"_ambiguous"``
block the manager already knows how to resolve interactively — nothing new
downstream. Bonds are authoritative: an under-specified site (missing the bonds
that tell same-name residues apart) is reported as a missing/ambiguous role, not
guessed by distance.

Specs live under ``~/.proprep/transformers/*.json`` and are discovered by
:func:`load_user_spec_transformers`, called alongside the existing ``*.py``
loader. The collision policy never lets a user spec override an already
registered transformer (built-in or auto-emitted).
"""

from __future__ import annotations

import json
import logging
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from proprep.redoxsite_prep.transformation.auto_rename import (
    AutoRenameTransformerBase,
    connectivity_signature,
    DEFAULT_USER_TRANSFORMER_DIR,
    _class_name,
    _norm_resname,
    _sanitize_name,
)

logger = logging.getLogger(__name__)

ResidueKey = Tuple[str, int]


# ---------------------------------------------------------------------------
# Atom-connectivity discriminators (tier 3), read from the site's bonds.
# All bond kinds count here (a c-type heme thioether is covalent, not a metal
# coordinate bond, so it is invisible to connectivity_signature but essential
# to tell the two cysteines apart).
# ---------------------------------------------------------------------------
def residue_bond_discriminators(redox_site, key: ResidueKey) -> set:
    """Set of ``(source_atom, partner_resname_normalized, target_atom)`` for
    every site bond incident to residue ``key``."""
    chain, resid = key
    resname_by_res: Dict[ResidueKey, str] = {}
    for a in getattr(redox_site, "atoms", []) or []:
        resname_by_res[(a.chain, a.resid)] = a.resname

    out: set = set()
    for bond in getattr(redox_site, "bonds", []) or []:
        c1 = getattr(bond, "atom1_coords", None)
        c2 = getattr(bond, "atom2_coords", None)
        i1 = (redox_site.coord_to_pdb.get(c1, {}) if c1 is not None else {}) or {}
        i2 = (redox_site.coord_to_pdb.get(c2, {}) if c2 is not None else {}) or {}
        for near, far in ((i1, i2), (i2, i1)):
            if near.get("chain") == chain and near.get("resid") == resid:
                partner_key = (far.get("chain"), far.get("resid"))
                partner_rn = _norm_resname(
                    resname_by_res.get(partner_key, far.get("resname", "")))
                out.add((near.get("atom_name"), partner_rn, far.get("atom_name")))
    return out


def _role_discriminators(role: Dict[str, Any]) -> List[Tuple[str, str, str]]:
    return [(d["source_atom"], _norm_resname(d["partner_resname"]), d["target_atom"])
            for d in role.get("discriminators", []) or []]


# ---------------------------------------------------------------------------
class SpecTransformer(AutoRenameTransformerBase):
    """Transferable transformer whose recipe is a JSON ``SPEC``.

    Subclasses (built dynamically from a spec file) supply only class
    attributes: ``TRANSFORMER_NAME``, ``DESCRIPTION``, ``SUPPORTED_SITE_TYPES``,
    ``SPEC`` (the parsed recipe), and — synthesized from the spec's roles for the
    inherited evaluate/requirements logic — ``RENAMES``.
    """

    SPEC: Dict[str, Any] = {}

    # ---- roles -----------------------------------------------------------
    @classmethod
    def _roles(cls) -> List[Dict[str, Any]]:
        return [r for r in cls.SPEC.get("roles", []) if not r.get("created")]

    @classmethod
    def _param_order(cls) -> List[str]:
        order = cls.SPEC.get("parameter_order")
        if order:
            return list(order)
        return [p["name"] for p in cls.SPEC.get("parameters", []) or []]

    # ---- matching --------------------------------------------------------
    @classmethod
    def _satisfies(cls, redox_site, cand: ResidueKey, role: Dict[str, Any]) -> bool:
        required = _role_discriminators(role)
        if not required:
            return True  # vacuous: no discriminator constraint on this role
        have = residue_bond_discriminators(redox_site, cand)
        return all(req in have for req in required)

    @classmethod
    def _bijection_by_discriminator(cls, redox_site, roles: List[Dict[str, Any]],
                                    candidates: List[ResidueKey]
                                    ) -> Optional[Dict[str, ResidueKey]]:
        """Assign each role to a unique candidate using atom-connectivity
        discriminators. Constrained-first greedy: repeatedly assign any role that
        has exactly one still-available satisfying candidate. Returns
        {label -> (chain, resid)} on a clean full assignment, else None."""
        remaining = list(candidates)
        assignment: Dict[str, ResidueKey] = {}
        unassigned = list(roles)
        progressed = True
        while unassigned and progressed:
            progressed = False
            for role in list(unassigned):
                avail = [c for c in remaining if cls._satisfies(redox_site, c, role)]
                if len(avail) == 1:
                    assignment[role["label"]] = avail[0]
                    remaining.remove(avail[0])
                    unassigned.remove(role)
                    progressed = True
                elif not avail:
                    return None
        if unassigned:
            return None
        return assignment

    @classmethod
    def match_components(cls, redox_site) -> Tuple[Dict[str, Any], List[str]]:
        matched: Dict[str, Any] = {}
        missing: List[str] = []
        labels = connectivity_signature(redox_site)

        groups: Dict[str, List[Dict[str, Any]]] = defaultdict(list)
        for role in cls._roles():
            groups[_norm_resname(role["resname"])].append(role)

        for _resnorm, roles in groups.items():
            candidates = cls._site_residues_by_resname(redox_site, roles[0]["resname"])
            if len(candidates) < len(roles):
                missing.append(f"{roles[0]['resname']} "
                               f"(need {len(roles)}, found {len(candidates)})")
                continue

            # Tier 1: a single role with a single candidate -> bind directly.
            if len(roles) == 1 and len(candidates) == 1:
                cls._write_role(matched, roles[0]["label"], candidates)
                continue

            # Tier 3: atom-connectivity discriminators (decisive for the ties
            # the element fingerprint cannot break).
            assignment = cls._bijection_by_discriminator(redox_site, roles, candidates)

            # Tier 2 fallback: element-fingerprint bijection (reuses the base's
            # signature matcher) when discriminators don't resolve.
            if assignment is None:
                entries = [{"target": r["label"], "signature": r.get("fingerprint")}
                           for r in roles]
                assignment = cls._match_by_signature(labels, entries, candidates)

            if assignment is not None:
                for label, (chain, resid) in assignment.items():
                    cls._write_role(matched, label, [(chain, resid)])
            else:
                diagnostic = cls._signature_diagnostic(
                    labels,
                    [{"target": r["label"], "signature": r.get("fingerprint")} for r in roles],
                    candidates,
                )
                cls._declare_ambiguous(
                    matched, roles[0]["resname"],
                    [r["label"] for r in roles], candidates, diagnostic)

        # Splits: assign each created (new) residue its own reserved ID = the
        # source residue's ID + the authored offset, mirroring bis_his_c_type_heme.
        cls._assign_created_residues(matched)
        return matched, missing

    # ---- split / new-residue slots ---------------------------------------
    @classmethod
    def _split_ops(cls) -> List[Dict[str, Any]]:
        return [op for op in cls.SPEC.get("operations", [])
                if op.get("op") == "move_to_new_residue" and op.get("creates_role")]

    @classmethod
    def _assign_created_residues(cls, components: Dict[str, Any]) -> None:
        """Give each split-off residue an absolute ID = source ID + offset."""
        for op in cls._split_ops():
            created = op["creates_role"]
            src = op.get("selector", {}).get("role")
            src_key = cls._resolve_role_res(components, src) if src else None
            if src_key is None:
                continue
            offset = int(op.get("relative", {}).get("new_resid_offset", 1))
            components[f"{created}_chain"] = src_key[0]
            components[f"{created}_id"] = src_key[1] + offset

    @classmethod
    def get_required_residue_count(cls) -> int:
        # Parent residue(s) keep their IDs; each split reserves one new slot.
        created = {op["creates_role"] for op in cls._split_ops()}
        return 1 + len(created)

    @classmethod
    def get_residue_space_plan(cls, components: Dict[str, Any]) -> Dict[str, int]:
        plan: Dict[str, int] = {r["label"]: 0 for r in cls._roles()}
        for op in cls._split_ops():
            plan[op["creates_role"]] = int(op.get("relative", {}).get("new_resid_offset", 1))
        return plan if plan else {"center": 0}

    @classmethod
    def update_components_with_id_mapping(cls, components: Dict[str, Any],
                                          id_mapping: Dict[Tuple[str, int], int]
                                          ) -> Dict[str, Any]:
        # Base remaps every <role>_id present in id_mapping; then recompute the
        # split-off residue IDs from the (possibly relocated) source residue.
        updated = super().update_components_with_id_mapping(components, id_mapping)
        cls._assign_created_residues(updated)
        return updated

    # ---- replay ----------------------------------------------------------
    @classmethod
    def _resolved_name(cls, op: Dict[str, Any], parameters: Dict[str, Any],
                       default: str) -> str:
        nbs = op.get("name_by_state")
        if not nbs:
            return default
        order = cls._param_order()
        if not order:
            return default
        key = "|".join(str(parameters.get(p)) for p in order)
        return nbs.get(key, default)

    @classmethod
    def _resolve_role_res(cls, components: Dict[str, Any], label: str
                          ) -> Optional[ResidueKey]:
        for chain, resid in cls._iter_role_residues(components, label):
            if resid is not None:
                return (chain, resid)
        return None

    @classmethod
    def get_transformation_sequence(cls, components: Dict[str, Any],
                                    parameters: Dict[str, Any]) -> List[Dict[str, Any]]:
        sequence: List[Dict[str, Any]] = []
        for i, op in enumerate(cls.SPEC.get("operations", [])):
            optype = op["op"]
            sel = op.get("selector", {})
            role = sel.get("role")
            action_in = op.get("action", {})

            if role is None:
                logger.warning("SpecTransformer %s: op %d has no role selector; "
                               "skipping.", cls.TRANSFORMER_NAME, i)
                continue
            key = cls._resolve_role_res(components, role)
            if key is None:
                logger.warning("SpecTransformer %s: role %r unbound; skipping op %d.",
                               cls.TRANSFORMER_NAME, role, i)
                continue
            chain, resid = key

            selector: Dict[str, Any] = {"chain_id": chain, "residue_id": resid}
            if sel.get("atom_names"):
                selector["atom_names"] = list(sel["atom_names"])

            action: Dict[str, Any] = {}
            if optype == "move_to_new_residue":
                # Split: move the selected atom subset out of the source residue
                # into the reserved new residue ID (from components, absolute).
                created = op.get("creates_role")
                created_key = cls._resolve_role_res(components, created) if created else None
                if created_key is None:
                    logger.warning("SpecTransformer %s: created role %r unbound; "
                                   "skipping split op %d.",
                                   cls.TRANSFORMER_NAME, created, i)
                    continue
                action["change_residue_name"] = cls._resolved_name(
                    op, parameters, action_in.get("change_residue_name"))
                action["change_residue_id"] = created_key[1]
                action["change_chain_id"] = created_key[0]
                if "rename_atoms" in action_in:
                    action["rename_atoms"] = dict(action_in["rename_atoms"])
            elif optype == "move_to_existing_residue":
                target_key = None
                if op.get("target_role"):
                    target_key = cls._resolve_role_res(components, op["target_role"])
                if target_key is None:
                    logger.warning("SpecTransformer %s: move target role unbound; "
                                   "skipping op %d.", cls.TRANSFORMER_NAME, i)
                    continue
                action["change_chain_id"] = target_key[0]
                action["change_residue_id"] = target_key[1]
                if "change_residue_name" in action_in:
                    action["change_residue_name"] = action_in["change_residue_name"]
            else:
                if "change_residue_name" in action_in:
                    action["change_residue_name"] = cls._resolved_name(
                        op, parameters, action_in["change_residue_name"])
                if "rename_atoms" in action_in:
                    action["rename_atoms"] = dict(action_in["rename_atoms"])
                if "change_chain_id" in action_in:
                    action["change_chain_id"] = action_in["change_chain_id"]
                if "change_insertion_code" in action_in:
                    action["change_insertion_code"] = action_in["change_insertion_code"]
                if action_in.get("convert_to_hetatm"):
                    action["convert_to_hetatm"] = True

            if not action:
                continue
            sequence.append({
                "id": f"{cls.TRANSFORMER_NAME}_op{i}_{optype}",
                "description": f"{optype} on role {role} ({chain}:{resid})",
                "selector": selector,
                "action": action,
            })
        return sequence

    # ---- parameters ------------------------------------------------------
    @classmethod
    def get_parameter_definitions(cls) -> Dict[str, Any]:
        # Start from the inherited fixed redox_state / spin_state (baked from the
        # spec's forcefield block via REDOX_STATE / SPIN_STATE) so those are
        # ALWAYS present -- the topology generator's transformer_info reads
        # parameters['redox_state'/'spin_state'] to locate the deposited FF. A
        # spec parameter of the same name overrides the fixed one with a choice.
        defs: Dict[str, Any] = dict(super().get_parameter_definitions())
        for p in cls.SPEC.get("parameters", []) or []:
            defs[p["name"]] = {
                "description": f"{p['name']} (recipe parameter)",
                "type": "choice",
                "options": list(p["options"]),
                "default": p.get("default", p["options"][0]),
            }
        return defs


# ---------------------------------------------------------------------------
# Build a concrete transformer class from a spec dict.
# ---------------------------------------------------------------------------
def build_spec_transformer_class(spec: Dict[str, Any]) -> type:
    """Create a registerable :class:`SpecTransformer` subclass from a spec."""
    name = _sanitize_name(spec.get("name", "spec_transformer"))
    roles = [r for r in spec.get("roles", []) if not r.get("created")]
    # RENAMES lets the inherited evaluate_redox_site / get_site_requirements work:
    # both only read the resname column.
    renames = [{"resname": r["resname"], "target": r["label"],
                "signature": r.get("fingerprint")} for r in roles]
    ff = spec.get("forcefield") or {}
    attrs = {
        "TRANSFORMER_NAME": name,
        "DESCRIPTION": spec.get("description") or f"Spec transformer {name}",
        "SUPPORTED_SITE_TYPES": spec.get("site_types") or ["spec"],
        "SPEC": spec,
        "RENAMES": renames,
        "FORCEFIELD_PATH": ff.get("path"),
        "REDOX_STATE": ff.get("redox_state", "unknown"),
        "SPIN_STATE": ff.get("spin_state", "unknown"),
    }
    return type(_class_name(name), (SpecTransformer,), attrs)


# ---------------------------------------------------------------------------
# Loader — discover ~/.proprep/transformers/*.json with a never-override policy.
# ---------------------------------------------------------------------------
_LOADED_SPEC_PATHS: set = set()


def load_user_spec_transformers(user_dir: Optional[Path] = None) -> List[str]:
    """Register every ``*.json`` spec in the user transformer dir (idempotent).

    Collision policy: a spec never overrides an already-registered transformer
    (built-in or auto-emitted ``.py``). On a name clash the spec is registered
    under ``<name>__spec`` instead, and a warning is logged, so user content can
    only ADD transformers, never silently replace existing behavior. Returns the
    names newly registered this call.
    """
    user_dir = Path(user_dir) if user_dir else DEFAULT_USER_TRANSFORMER_DIR
    if not user_dir.is_dir():
        return []

    from proprep.redoxsite_prep.transformation.redox_transformer_framework import (
        redox_transformer_registry,
    )

    newly: List[str] = []
    for path in sorted(user_dir.glob("*.json")):
        if path.name.startswith("_"):
            continue
        resolved = str(path.resolve())
        if resolved in _LOADED_SPEC_PATHS:
            continue
        _LOADED_SPEC_PATHS.add(resolved)
        try:
            spec = json.loads(path.read_text())
            cls = build_spec_transformer_class(spec)
            name = cls.TRANSFORMER_NAME
            existing = set(redox_transformer_registry.list_transformer_names())
            if name in existing:
                qualified = f"{name}__spec"
                if qualified in existing:
                    logger.warning("Spec transformer %r collides and %r is also "
                                   "taken; skipping %s.", name, qualified, path.name)
                    continue
                logger.warning("Spec transformer %r collides with an existing "
                               "transformer; registering as %r instead "
                               "(existing transformers are never overridden).",
                               name, qualified)
                cls.TRANSFORMER_NAME = qualified
                name = qualified
            redox_transformer_registry.register(cls)
            newly.append(name)
        except Exception as exc:  # pragma: no cover - defensive
            logger.warning("Skipping spec transformer %s: %s", path, exc)

    return newly
