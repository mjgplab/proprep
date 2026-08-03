#!/usr/bin/env python3
"""Auto-emitted residue-rename transformers.

A parameterizer (MCPB metal sites today; small-molecule / modified-AA later)
already computes the exact residue-name changes it applies to produce a
force-field-ready "prepared" PDB. This module turns that rename map into a
reusable RedoxSite transformer so the same edits can be re-applied to another
instance of the site later — without re-running the parameterization — letting
the user reuse the parameters deposited in ``~/.proprep``.

Design (see the long design discussion): the emitted transformer is *pure
data* — a tiny subclass of :class:`AutoRenameTransformerBase` carrying only a
``RENAMES`` table. All matching intelligence is hand-written *once* here in the
base class, as a generic interpreter of that table:

  * Tier 1 — a residue name that occurs once maps to a single target: rename it.
  * Tier 2 — several residues of the same name map to the *same* target
    (e.g. 4 Cys -> CYZ): rename them all; order is irrelevant.
  * Tier 3 — several residues of the same name map to *different* targets
    (e.g. two Glu -> GL1/GL2): try to disambiguate automatically from each
    residue's coordination signature (which/how-many metals it bonds and via
    which atoms). If that uniquely resolves, assign silently; if the site is
    symmetric so it cannot, declare an ``"_ambiguous"`` block and let the
    transformation manager prompt the user (the existing resolver even
    highlights candidates in the 3D viewer). Either way the component-matching
    summary lets the user review before anything is applied.

Bespoke transformers (``bis_his_c``, ``zinc_cys4``, …) are untouched: they
subclass :class:`RedoxSiteTransformerBase` directly and keep their own
chemistry-specific matching. This is the automatic breadth for the long tail
of user-parameterized sites, alongside their curated depth.

Discovery: the built-in transformers register via static imports in
``transformers/__init__.py``. Auto-emitted transformers live under
``~/.proprep/transformers`` and are picked up by :func:`load_user_transformers`,
which the transformation manager calls at startup.
"""

from __future__ import annotations

import hashlib
import importlib.util
import logging
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from proprep.redoxsite_prep.transformation.redox_transformer_framework import (
    RedoxSiteTransformerBase,
    TransformerEvaluation,
    TransformerEvaluationDetail,
    register_redox_transformer,
)

logger = logging.getLogger(__name__)


# METALS is the element set the redox detector uses to classify coordinate
# bonds; we reuse it so "which partner is the metal" is decided identically
# at emit time and apply time.
try:
    from proprep.structure_prep.comprehensive_redox_detector import METALS
except Exception:  # pragma: no cover - defensive; detector should always import
    METALS = {
        "MN", "FE", "ZN", "CU", "NI", "CO", "MG", "CA", "MO", "W", "V", "CD",
        "NA", "K",
    }


DEFAULT_USER_TRANSFORMER_DIR = Path.home() / ".proprep" / "transformers"


# Protonation / bonding variants collapse to their standard base residue for
# matching. A parameterizer records names AS THEY WERE at MCPB time — i.e. after
# protonation assignment (HID/HIE/HIP, ASH, GLH, CYX…) — but a reuse run applies
# the transformer to a fresh structure that still has the standard names (the
# protonation analyzer runs AFTER the redox site preparer). So an "HID" entry
# must be satisfied by an "HIS" residue, etc.
_PROTONATION_FAMILY = {
    "HID": "HIS", "HIE": "HIS", "HIP": "HIS", "HSD": "HIS", "HSE": "HIS", "HSP": "HIS",
    "ASH": "ASP",
    "GLH": "GLU",
    "LYN": "LYS",
    "CYM": "CYS", "CYX": "CYS",
    "TYM": "TYR",
    # Water aliases collapse too: a metal-restrained water is WAT at MCPB time
    # but often HOH in a freshly-loaded structure, and it coordinates the metal
    # so it participates in the connectivity fingerprint.
    "HOH": "WAT", "H2O": "WAT", "SOL": "WAT", "TIP3": "WAT", "TIP4": "WAT", "TP3": "WAT",
}


def _norm_resname(name: str) -> str:
    """Collapse a protonation/bonding/water variant to a canonical base name,
    so a name recorded at MCPB time (post-protonation HID, water WAT) matches
    the same residue in a freshly-loaded structure (HIS, HOH)."""
    upper = (name or "").upper()
    return _PROTONATION_FAMILY.get(upper, upper)


# ---------------------------------------------------------------------------
# Coordination signature — shared by emit (parameterized site) and apply
# (brought site) so the two are directly comparable.
# ---------------------------------------------------------------------------

def coordination_signature(redox_site, chain: str, resid: int) -> Dict[str, Any]:
    """Summarize how residue ``(chain, resid)`` coordinates metals in a site.

    Returns ``{"n_metals": int, "coord_atoms": [str, ...]}`` where ``n_metals``
    is the number of *distinct* metal atoms the residue has coordinate bonds to
    (1 = terminal, >=2 = bridging) and ``coord_atoms`` are the residue's own
    atom names participating in those bonds. This is the discriminator used to
    tell same-name residues apart (e.g. a bridging Glu vs a terminal Glu).
    """
    coord_atoms: set = set()
    metal_partners: set = set()

    for bond in getattr(redox_site, "bonds", []) or []:
        if getattr(bond, "chemical_type", None) != "coordinate":
            continue
        c1 = getattr(bond, "atom1_coords", None)
        c2 = getattr(bond, "atom2_coords", None)
        info1 = (redox_site.coord_to_pdb.get(c1, {}) if c1 is not None else {}) or {}
        info2 = (redox_site.coord_to_pdb.get(c2, {}) if c2 is not None else {}) or {}
        e1 = (getattr(bond, "atom1_element", "") or info1.get("element") or "").upper()
        e2 = (getattr(bond, "atom2_element", "") or info2.get("element") or "").upper()

        # Consider each endpoint as possibly-this-residue, the other as metal.
        for info_self, info_metal, e_metal, coord_metal in (
            (info1, info2, e2, c2),
            (info2, info1, e1, c1),
        ):
            if (info_self.get("chain") == chain
                    and info_self.get("resid") == resid
                    and e_metal in METALS
                    and coord_metal is not None):
                name = info_self.get("atom_name")
                if name:
                    coord_atoms.add(name)
                metal_partners.add(coord_metal)

    return {"n_metals": len(metal_partners), "coord_atoms": sorted(coord_atoms)}


def _canon_signature(sig: Optional[Dict[str, Any]]) -> Optional[Tuple]:
    """Hashable canonical form for equality comparison, or None if absent."""
    if not sig:
        return None
    return (int(sig.get("n_metals", 0)), tuple(sorted(sig.get("coord_atoms", []))))


# ---------------------------------------------------------------------------
# Connectivity signature — a Weisfeiler-Lehman neighborhood fingerprint over
# the site's coordination graph. Because the connectivity of a given site type
# is invariant across proteins, so is this label: a residue's role is pinned by
# its coordination *environment* (which metals it bonds, and those metals' other
# ligands), not just its own atoms. This distinguishes same-name residues that
# a local signature cannot — e.g. two metals told apart by their differing
# ligand sets, or two Glu told apart by which metal each coordinates. It ties
# (leaving the site ambiguous) only for genuinely symmetric/automorphic sites,
# which is exactly when a human must choose.
# ---------------------------------------------------------------------------

_WL_ROUNDS = 4


def connectivity_signature(redox_site) -> Dict[Tuple[str, int], str]:
    """Canonical WL label per site residue, keyed by (chain, resid).

    Labels are protein-independent (built only from residue names, metal
    elements, coordinating atom names and graph structure — never resids or
    chains), so a label computed on the parameterized site matches the same
    role computed on another instance of the site.
    """
    # Nodes: every site residue EXCEPT waters. Waters are excluded on purpose:
    # they are never renamed, their naming varies (HOH/WAT/TIP3/...), and their
    # metal coordination is frequently replaced by an MD distance restraint that
    # the user reasonably omits from the reuse site's bonds. Including them would
    # make the fingerprint depend on whether water bonds happen to be defined,
    # which differs between MCPB time and a fresh reuse structure and shifts
    # EVERY label in the connected component (WL propagates). A residue is a
    # "metal" node if any of its atoms is a metal element.
    residues: Dict[Tuple[str, int], Dict[str, Any]] = {}
    for a in getattr(redox_site, "atoms", []) or []:
        if _norm_resname(a.resname) == "WAT":
            continue
        key = (a.chain, a.resid)
        elem = (getattr(a, "element", "") or "").upper()
        rec = residues.setdefault(
            key, {"resname": a.resname, "element": elem, "is_metal": False})
        if elem in METALS:
            rec["is_metal"] = True
            rec["element"] = elem

    # Edges: only metal<->ligand coordinate bonds. Metal-metal bonds are excluded
    # for the same reason as waters -- their definition is optional and often
    # skipped on the reuse side, and the two metals are already told apart by
    # their differing ligand sets. Each edge stores, from each endpoint's view,
    # the atom name that endpoint contributes. Bonds touching an excluded node
    # (e.g. a water) are dropped by the membership check below.
    edges: Dict[Tuple[str, int], List[Tuple[Tuple[str, int], str]]] = defaultdict(list)
    for bond in getattr(redox_site, "bonds", []) or []:
        if getattr(bond, "chemical_type", None) != "coordinate":
            continue
        c1 = getattr(bond, "atom1_coords", None)
        c2 = getattr(bond, "atom2_coords", None)
        i1 = (redox_site.coord_to_pdb.get(c1, {}) if c1 is not None else {}) or {}
        i2 = (redox_site.coord_to_pdb.get(c2, {}) if c2 is not None else {}) or {}
        k1 = (i1.get("chain"), i1.get("resid"))
        k2 = (i2.get("chain"), i2.get("resid"))
        if None in k1 or None in k2 or k1 not in residues or k2 not in residues:
            continue
        # Label the edge by the endpoint's bonding ELEMENT, never its atom NAME.
        # A residue's coordinating atom name is often symmetry-degenerate: a Glu
        # carboxylate's OE1/OE2 (or Asp OD1/OD2) are chemically equivalent and
        # which one is nearest the metal flips between structures and even
        # between snapshots. Keying on the name made the fingerprint depend on
        # that arbitrary choice, and because WL propagates, a single OE1<->OE2
        # flip cascaded to shift every label reachable from that metal. Element
        # ("O"/"N"/"S") keeps the real discriminator (what donor type bonds)
        # while being invariant to the degenerate naming.
        e_k1 = (i1.get("element") or "").upper()
        e_k2 = (i2.get("element") or "").upper()
        edges[k1].append((k2, e_k1))
        edges[k2].append((k1, e_k2))

    # Base labels: metals by element, ligands by NORMALIZED residue name (so the
    # fingerprint is invariant to protonation/water naming that differs between
    # MCPB time and a fresh reuse structure — HID vs HIS, WAT vs HOH).
    base: Dict[Tuple[str, int], str] = {}
    for key, rec in residues.items():
        base[key] = (f"M:{rec['element']}" if rec["is_metal"]
                     else f"L:{_norm_resname(rec['resname'])}")

    # WL refinement: relabel each node by a short stable hash of (its label,
    # sorted neighbour (label, my-bonding-atom) pairs). Hashing each round keeps
    # labels compact (naive repr-nesting explodes exponentially). Fixed rounds;
    # over-refinement is safe — WL never separates automorphic nodes.
    label: Dict[Tuple[str, int], str] = dict(base)
    for _ in range(_WL_ROUNDS):
        new_label: Dict[Tuple[str, int], str] = {}
        for key in residues:
            nbrs = sorted(
                (label[nk], atom) for (nk, atom) in edges.get(key, [])
            )
            new_label[key] = _wl_digest((label[key], nbrs))
        label = new_label

    # Final signature: readable base type + the refined fingerprint.
    return {key: f"{base[key]}#{label[key]}" for key in residues}


def _wl_digest(payload) -> str:
    """Short, deterministic hash of a canonical payload (not Python's salted
    hash, which varies per process)."""
    return hashlib.blake2b(repr(payload).encode(), digest_size=8).hexdigest()


# ---------------------------------------------------------------------------
# Generic base class — all the matching intelligence lives here.
# ---------------------------------------------------------------------------

class AutoRenameTransformerBase(RedoxSiteTransformerBase):
    """Data-driven residue-rename transformer.

    Subclasses (auto-emitted) supply only class attributes:

    * ``TRANSFORMER_NAME`` / ``DESCRIPTION`` / ``SUPPORTED_SITE_TYPES`` (str/list)
    * ``RENAMES`` — list of ``{"resname": str, "target": str,
      "signature": {...}?}``; one entry per residue that was renamed.
    * ``REDOX_STATE`` / ``SPIN_STATE`` — baked from the parameterization,
      threaded through as fixed parameters.
    * ``PROVENANCE`` — free-form dict (library path, source, …) for the record.
    """

    RENAMES: List[Dict[str, Any]] = []
    REDOX_STATE: str = "unknown"
    SPIN_STATE: str = "unknown"
    PROVENANCE: Dict[str, Any] = {}

    # ---- helpers ----------------------------------------------------------

    @classmethod
    def _rename_groups(cls) -> Dict[str, List[Dict[str, Any]]]:
        groups: Dict[str, List[Dict[str, Any]]] = defaultdict(list)
        for entry in cls.RENAMES:
            groups[entry["resname"]].append(entry)
        return groups

    @staticmethod
    def _site_residues_by_resname(redox_site, resname: str) -> List[Tuple[str, int]]:
        """(chain, resid) of every residue matching ``resname`` in the site.

        Matching is protonation-family aware: an ``HID`` request matches an
        ``HIS`` residue (and vice versa), etc.
        """
        target = _norm_resname(resname)
        residues = {
            (a.chain, a.resid)
            for a in getattr(redox_site, "atoms", [])
            if _norm_resname(a.resname) == target
        }
        return sorted(residues)

    @staticmethod
    def _write_role(matched: Dict[str, Any], target: str,
                    residues: List[Tuple[str, int]]) -> None:
        """Record a target's assigned residues using the resolver's convention.

        Single residue -> ``<target>_chain`` / ``<target>_id``; multiple ->
        ``<target>_chain_<n>`` / ``<target>_id_<n>`` (1-indexed). This mirrors
        exactly what :meth:`_resolve_ambiguities` writes back, so auto-resolved
        and user-resolved components look identical downstream.
        """
        if len(residues) == 1:
            chain, resid = residues[0]
            matched[f"{target}_chain"] = chain
            matched[f"{target}_id"] = resid
        else:
            for n, (chain, resid) in enumerate(residues, start=1):
                matched[f"{target}_chain_{n}"] = chain
                matched[f"{target}_id_{n}"] = resid

    @staticmethod
    def _iter_role_residues(components: Dict[str, Any], target: str):
        """Yield (chain, resid) for a target, single- or multi-slot form."""
        if f"{target}_id" in components:
            yield components[f"{target}_chain"], components[f"{target}_id"]
            return
        n = 1
        while f"{target}_id_{n}" in components:
            yield components.get(f"{target}_chain_{n}"), components[f"{target}_id_{n}"]
            n += 1

    @classmethod
    def _match_by_signature(cls, labels: Dict[Tuple[str, int], str],
                            entries: List[Dict[str, Any]],
                            candidates: List[Tuple[str, int]]
                            ) -> Optional[Dict[str, Tuple[str, int]]]:
        """Try to bijectively map targets -> candidate residues by WL label.

        ``labels`` is ``connectivity_signature(brought_site)``. Each entry's
        baked ``signature`` is the WL label it had in the parameterized site.
        Returns ``{target: (chain, resid)}`` on a clean unique match, else None
        (caller then defers to the user). Conservative: any signature missing,
        any tie among target labels, or any leftover candidate -> None.
        """
        target_labels: Dict[str, Optional[str]] = {
            e["target"]: e.get("signature") for e in entries
        }
        if any(s is None for s in target_labels.values()):
            return None
        # Target labels must be distinct, else the targets are indistinguishable
        # even in principle (a symmetric site) -> defer to the user.
        if len(set(target_labels.values())) != len(target_labels):
            return None

        cand_by_label: Dict[str, List[Tuple[str, int]]] = defaultdict(list)
        for chain, resid in candidates:
            cand_by_label[labels.get((chain, resid))].append((chain, resid))

        assignment: Dict[str, Tuple[str, int]] = {}
        used: set = set()
        for target, wl in target_labels.items():
            options = [cr for cr in cand_by_label.get(wl, []) if cr not in used]
            if len(options) != 1:
                return None
            assignment[target] = options[0]
            used.add(options[0])
        if len(used) != len(candidates):
            return None
        return assignment

    # ---- framework interface ---------------------------------------------

    @classmethod
    def evaluate_redox_site(cls, redox_site) -> TransformerEvaluation:
        need = Counter(e["resname"] for e in cls.RENAMES)
        details: List[TransformerEvaluationDetail] = []
        met = 0
        for resname, count in need.items():
            found = len(cls._site_residues_by_resname(redox_site, resname))
            ok = found >= count
            details.append(TransformerEvaluationDetail(
                description=f"{count}x {resname}",
                passed=ok, value_found=found, value_expected=count,
            ))
            if ok:
                met += 1
        total = len(need)
        return TransformerEvaluation(
            is_valid=(total > 0 and met == total),
            confidence=(met / total if total else 0.0),
            description=cls.DESCRIPTION,
            requirements_met=met,
            total_requirements=total,
            details=details,
        )

    @classmethod
    def get_site_requirements(cls) -> Dict[str, Any]:
        need = Counter(e["resname"] for e in cls.RENAMES)
        return {
            "centers": [],
            "atoms": {
                "required_residues": {
                    resname: {"min_count": count} for resname, count in need.items()
                }
            },
            "bonds": [],
        }

    @classmethod
    def match_components(cls, redox_site) -> Tuple[Dict[str, Any], List[str]]:
        matched: Dict[str, Any] = {}
        missing: List[str] = []

        # WL connectivity labels for the brought site, computed once.
        labels = connectivity_signature(redox_site)

        for resname, entries in cls._rename_groups().items():
            candidates = cls._site_residues_by_resname(redox_site, resname)
            distinct_targets = list(dict.fromkeys(e["target"] for e in entries))

            if not candidates:
                missing.append(resname)
                continue

            # Tier 1 / Tier 2: one target for this residue name -> rename all.
            if len(distinct_targets) == 1:
                cls._write_role(matched, distinct_targets[0], candidates)
                continue

            # Tier 3: same name, different targets -> must disambiguate.
            if len(candidates) != len(entries):
                missing.append(
                    f"{resname} (need {len(entries)}, found {len(candidates)})"
                )
                continue

            assignment = cls._match_by_signature(labels, entries, candidates)
            if assignment is not None:
                for target, (chain, resid) in assignment.items():
                    cls._write_role(matched, target, [(chain, resid)])
            else:
                diagnostic = cls._signature_diagnostic(labels, entries, candidates)
                cls._declare_ambiguous(matched, resname, distinct_targets,
                                       candidates, diagnostic)

        return matched, missing

    @classmethod
    def _signature_diagnostic(cls, labels: Dict[Tuple[str, int], str],
                              entries: List[Dict[str, Any]],
                              candidates: List[Tuple[str, int]]) -> str:
        """Explain, in one human-readable block, why auto-resolution failed.

        Shows each target's baked (parameterization-time) fingerprint and each
        candidate residue's freshly-computed fingerprint, marking where they do
        and don't line up, so a mismatch is self-diagnosing at the prompt rather
        than requiring an offline archaeology of the two coordination graphs.
        """
        baked = {e["target"]: e.get("signature") for e in entries}
        if any(s is None for s in baked.values()):
            return ("This transformer was emitted before connectivity "
                    "fingerprinting, so roles can't be auto-assigned. "
                    "Re-run the parameterization to re-emit it, or assign below.")

        # Reverse index: which target(s) does each computed label satisfy.
        target_by_sig: Dict[str, List[str]] = defaultdict(list)
        for tgt, sig in baked.items():
            target_by_sig[sig].append(tgt)

        if len(set(baked.values())) != len(baked):
            head = ("The parameterized fingerprints for these roles are "
                    "identical (a symmetric site), so they can't be told apart "
                    "automatically. Assign each below.")
        else:
            head = ("Auto-assignment failed: the coordination fingerprints in "
                    "this structure don't match the ones saved at "
                    "parameterization. Most often the bonds defined here differ "
                    "from those at parameterization (e.g. a symmetry-equivalent "
                    "donor atom, or a missing/extra ligand). Re-emit by re-running "
                    "the parameterization, or assign manually below.")

        lines = [head, "", "  expected roles (from parameterization):"]
        for tgt, sig in baked.items():
            lines.append(f"    {tgt:<6} {sig}")
        lines.append("  residues in this site:")
        for chain, resid in candidates:
            sig = labels.get((chain, resid), "<no signature>")
            hit = target_by_sig.get(sig)
            note = f"-> matches {hit[0]}" if hit and len(hit) == 1 else "-> no baked match"
            lines.append(f"    {chain}:{resid:<6} {sig}  {note}")
        return "\n".join(lines)

    @classmethod
    def _declare_ambiguous(cls, matched: Dict[str, Any], resname: str,
                           targets: List[str],
                           candidates: List[Tuple[str, int]],
                           diagnostic: Optional[str] = None) -> None:
        cand_dicts = [
            {"chain": chain, "resid": resid, "display": f"{resname} {chain}:{resid}"}
            for chain, resid in candidates
        ]
        block: Dict[str, Any] = {
            "label": f"{resname} role assignment",
            "description": (
                f"{len(candidates)} {resname} residues map to distinct names "
                f"({', '.join(targets)}); assign each."
            ),
            "candidates": cand_dicts,
            "roles": {target: 1 for target in targets},
        }
        if diagnostic:
            block["diagnostic"] = diagnostic
        matched.setdefault("_ambiguous", []).append(block)

    @classmethod
    def get_transformation_sequence(cls, components: Dict[str, Any],
                                    parameters: Dict[str, Any]) -> List[Dict[str, Any]]:
        sequence: List[Dict[str, Any]] = []
        seen_targets: set = set()
        for entry in cls.RENAMES:
            target = entry["target"]
            if target in seen_targets:
                continue
            seen_targets.add(target)
            resname = entry["resname"]
            for chain, resid in cls._iter_role_residues(components, target):
                if resid is None:
                    continue
                # Selector omits residue_name deliberately: at apply time the
                # residue may carry a protonation variant of the baked name
                # (HIS vs HID), and the apply engine filters atoms by exact
                # residue_name. chain+resid resolves to unique coordinates, so
                # it is the correct, protonation-agnostic key.
                action: Dict[str, Any] = {"change_residue_name": target}
                # A small-molecule ligand's library carries antechamber's atom
                # names (the mol2 the params were built from), but the reuse
                # structure has the original PDB/crystallographic names. Rename
                # the residue's atoms to the library naming in the SAME step so
                # tleap's loadpdb matches them; the apply engine applies both the
                # residue rename and the atom rename to each matched atom line.
                atom_renames = entry.get("atom_renames")
                if atom_renames:
                    action["rename_atoms"] = dict(atom_renames)
                sequence.append({
                    "id": f"rename_{resname}_{chain}_{resid}_to_{target}",
                    "description": f"Rename {resname} {chain}:{resid} -> {target}",
                    "selector": {
                        "chain_id": chain,
                        "residue_id": resid,
                    },
                    "action": action,
                })
        return sequence

    @classmethod
    def get_required_residue_count(cls) -> int:
        # Rename-only: every residue keeps its own ID, no new slots allocated.
        return 1

    @classmethod
    def get_residue_space_plan(cls, components: Dict[str, Any]) -> Dict[str, int]:
        return {"center": 0}

    @classmethod
    def update_components_with_id_mapping(cls, components: Dict[str, Any],
                                          id_mapping: Dict[Tuple[str, int], int]
                                          ) -> Dict[str, Any]:
        """Remap our ``<target>_id`` / ``<target>_id_<n>`` roles after the
        microstate ID-reallocation step (mirrors zinc_cys4's override for its
        custom keys). Standard keys are handled by super().
        """
        updated = super().update_components_with_id_mapping(components, id_mapping)
        for key in list(updated.keys()):
            if key.endswith("_id"):
                chain_key = key[:-3] + "_chain"
            elif "_id_" in key:
                base, suffix = key.split("_id_", 1)
                chain_key = f"{base}_chain_{suffix}"
            else:
                continue
            chain = updated.get(chain_key)
            resid = updated.get(key)
            if chain is not None and isinstance(resid, int) and (chain, resid) in id_mapping:
                updated[key] = id_mapping[(chain, resid)]
        return updated

    @classmethod
    def get_parameter_definitions(cls) -> Dict[str, Any]:
        # No user choices: targets are baked in. Declare redox/spin as fixed so
        # downstream FF-set lookups (which read parameters['redox_state'] etc.)
        # thread through, matching zinc_cys4's approach.
        return {
            "redox_state": {
                "description": "Redox state (baked from parameterization)",
                "type": "fixed", "value": cls.REDOX_STATE,
            },
            "spin_state": {
                "description": "Spin state (baked from parameterization)",
                "type": "fixed", "value": cls.SPIN_STATE,
            },
        }

    @classmethod
    def validate_parameters(cls, parameters: Dict[str, Any]) -> Tuple[bool, str]:
        return True, "No parameters required"

    @classmethod
    def get_parameter_mappings(cls, parameters: Dict[str, Any]) -> Dict[str, str]:
        return {}


# ---------------------------------------------------------------------------
# Emitter — serialize a rename table into a data-only transformer module.
# ---------------------------------------------------------------------------

_MODULE_TEMPLATE = '''#!/usr/bin/env python3
"""Auto-generated residue-rename transformer. DO NOT EDIT BY HAND.

{description}

Generated by ProPrep's parameterizer to make the residue renames it applied
reusable via the Redox Site Preparer. Regenerating overwrites this file.
"""

from proprep.redoxsite_prep.transformation.redox_transformer_framework import (
    register_redox_transformer,
)
from proprep.redoxsite_prep.transformation.auto_rename import (
    AutoRenameTransformerBase,
)


@register_redox_transformer
class {class_name}(AutoRenameTransformerBase):
    TRANSFORMER_NAME = {name!r}
    DESCRIPTION = {description!r}
    SUPPORTED_SITE_TYPES = {site_types!r}

    # Relative cofactor path (under specialized_residues) for the deposited FF.
    # The Topology Generator discovers frcmod/lib/atom-types from this + the
    # redox/spin state, the same way the built-in transformers do; leaving it
    # None makes topology generation skip FF loading for the renamed residues.
    FORCEFIELD_PATH = {forcefield_path!r}

    REDOX_STATE = {redox_state!r}
    SPIN_STATE = {spin_state!r}
    PROVENANCE = {provenance!r}

    RENAMES = {renames}
'''


def _sanitize_name(name: str) -> str:
    """Force a string to a safe transformer-name token (lowercase, id-safe)."""
    token = re.sub(r"[^0-9a-zA-Z]+", "_", name).strip("_").lower()
    if not token:
        token = "auto_rename"
    if token[0].isdigit():
        token = f"t_{token}"
    return token


def _class_name(name: str) -> str:
    parts = [p for p in re.split(r"[^0-9a-zA-Z]+", name) if p]
    camel = "".join(p[:1].upper() + p[1:] for p in parts) or "AutoRename"
    if camel[0].isdigit():
        camel = f"T{camel}"
    return f"{camel}Transformer"


def _format_renames(renames: List[Dict[str, Any]]) -> str:
    """Render RENAMES as readable, deterministic source (one entry per line)."""
    lines = ["["]
    for entry in renames:
        lines.append(f"        {entry!r},")
    lines.append("    ]")
    return "\n".join(lines)


def emit_rename_transformer(rename_table: List[Dict[str, Any]], *,
                            name: str,
                            description: str,
                            redox_state: str = "unknown",
                            spin_state: str = "unknown",
                            provenance: Optional[Dict[str, Any]] = None,
                            site_types: Optional[List[str]] = None,
                            forcefield_path: Optional[str] = None,
                            target_dir: Optional[Path] = None) -> Path:
    """Write a data-only rename transformer module and return its path.

    ``rename_table`` is a list of ``{"resname", "target", "signature"?}`` dicts
    (one per renamed residue) — exactly what the caller derives from its own
    rename map. The file is written to ``~/.proprep/transformers/<name>.py``;
    re-emitting the same ``name`` overwrites it, so re-parameterizing a site
    updates rather than duplicates.
    """
    if not rename_table:
        raise ValueError("rename_table is empty; nothing to emit")

    token = _sanitize_name(name)
    class_name = _class_name(token)
    target_dir = Path(target_dir) if target_dir else DEFAULT_USER_TRANSFORMER_DIR
    target_dir.mkdir(parents=True, exist_ok=True)

    source = _MODULE_TEMPLATE.format(
        class_name=class_name,
        name=token,
        description=description,
        site_types=(site_types or ["metal"]),
        forcefield_path=forcefield_path,
        redox_state=redox_state,
        spin_state=spin_state,
        provenance=(provenance or {}),
        renames=_format_renames(rename_table),
    )

    out_path = target_dir / f"{token}.py"
    out_path.write_text(source)
    # debug, not info: every caller already prints its own user-facing
    # "transformer emitted" line, so info-level here just adds a noisy
    # timestamped duplicate to the console (e.g. modAA step 10).
    logger.debug("Emitted rename transformer %r -> %s", token, out_path)
    return out_path


# ---------------------------------------------------------------------------
# Loader — import user transformer modules so their decorators register.
# ---------------------------------------------------------------------------

_LOADED_USER_PATHS: set = set()


def load_user_transformers(user_dir: Optional[Path] = None) -> List[str]:
    """Import every ``*.py`` in the user transformer dir (idempotent).

    The built-in transformers register via static imports in
    ``transformers/__init__.py``; nothing loads the user directory, so this
    fills that gap. Each file is imported at most once per process; import
    errors are logged and skipped so one bad file can't break the preparer.
    Returns the transformer names newly registered this call.
    """
    user_dir = Path(user_dir) if user_dir else DEFAULT_USER_TRANSFORMER_DIR
    if not user_dir.is_dir():
        return []

    from proprep.redoxsite_prep.transformation.redox_transformer_framework import (
        redox_transformer_registry,
    )

    before = set(redox_transformer_registry.list_transformer_names())
    for path in sorted(user_dir.glob("*.py")):
        if path.name.startswith("_"):
            continue
        resolved = str(path.resolve())
        if resolved in _LOADED_USER_PATHS:
            continue
        _LOADED_USER_PATHS.add(resolved)
        module_name = f"proprep_user_transformer_{path.stem}"
        try:
            spec = importlib.util.spec_from_file_location(module_name, resolved)
            if spec is None or spec.loader is None:
                logger.warning("Could not load user transformer %s (no spec)", path)
                continue
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
        except Exception as exc:  # pragma: no cover - defensive
            logger.warning("Skipping user transformer %s: %s", path, exc)

    # Also discover JSON spec transformers from the same directory. Lazy import
    # avoids a module cycle (spec_transformer imports this module); guarded so a
    # spec-loader problem can never break the established .py discovery above.
    try:
        from proprep.redoxsite_prep.transformation.spec_transformer import (
            load_user_spec_transformers,
        )
        load_user_spec_transformers(user_dir)
    except Exception as exc:  # pragma: no cover - defensive
        logger.warning("Spec transformer discovery failed: %s", exc)

    after = set(redox_transformer_registry.list_transformer_names())
    return sorted(after - before)
