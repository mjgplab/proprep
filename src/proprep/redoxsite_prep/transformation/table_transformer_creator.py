#!/usr/bin/env python3
"""Table-driven transformer creator (the user-facing authoring UX).

The user sees the atoms of a detected redox site grouped by residue, then applies
ONE PDB-editing operation at a time; the table redraws after each edit because
order matters. When a subset of atoms is split into a new residue, the tool
cascades the ID shift onto downstream residues. At the end the recipe is saved as
a reusable JSON transformer spec under ``~/.proprep/transformers`` — the same spec
:class:`SpecTransformer` interprets, so it transfers to any protein with that site
type.

This replaces the abstract phrase-based v3 creator as the recommended path; the
v3 creator is kept as a deprecated legacy option. All editing/role logic is a
faithful port of the standalone prototype, adapted to operate on the package
``RedoxSite`` and to compute role fingerprints/discriminators with the same
functions :class:`SpecTransformer` uses at apply time.
"""

from __future__ import annotations

import copy
import json
import logging
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Set, Tuple

from proprep.utils.prompts import prompt_with_context, confirm_with_context
from proprep.redoxsite_prep.transformation.auto_rename import (
    connectivity_signature, DEFAULT_USER_TRANSFORMER_DIR,
)
from proprep.redoxsite_prep.transformation.spec_transformer import (
    residue_bond_discriminators,
)

logger = logging.getLogger(__name__)

SCHEMA_VERSION = "0.3"
ResidueKey = Tuple[str, int, str]   # (chain, resid, icode)


class RecipeError(Exception):
    pass


# ---------------------------------------------------------------------------
# Internal editable structure, built from a RedoxSite. uid is the immutable key
# (coordinates never change in a RedoxSite; uid mirrors that stability).
# ---------------------------------------------------------------------------
@dataclass
class _Atom:
    uid: int
    record: str
    name: str
    resname: str
    chain: str
    resid: int
    icode: str
    element: str


class _Structure:
    def __init__(self, atoms: List[_Atom]):
        self.atoms = atoms

    @classmethod
    def from_redox_site(cls, redox_site) -> "_Structure":
        atoms: List[_Atom] = []
        for uid, a in enumerate(getattr(redox_site, "atoms", []) or []):
            rec = "ATOM"
            props = getattr(a, "properties", {}) or {}
            if props.get("record") in ("ATOM", "HETATM"):
                rec = props["record"]
            atoms.append(_Atom(
                uid=uid, record=rec, name=a.atom_name, resname=a.resname,
                chain=a.chain, resid=a.resid,
                icode=(getattr(a, "insertion_code", "") or ""),
                element=(getattr(a, "element", "") or ""),
            ))
        return cls(atoms)

    def residue_order(self) -> List[ResidueKey]:
        seen: Dict[ResidueKey, None] = {}
        for a in self.atoms:
            seen.setdefault((a.chain, a.resid, a.icode), None)
        return sorted(seen.keys(), key=lambda k: (k[0], k[1], k[2]))

    def residue_atoms(self, key: ResidueKey) -> List[_Atom]:
        return [a for a in self.atoms if (a.chain, a.resid, a.icode) == key]

    def resname_of(self, chain: str, resid: int, icode: str = "") -> Optional[str]:
        for a in self.atoms:
            if a.chain == chain and a.resid == resid and a.icode == icode:
                return a.resname
        return None

    def find(self, chain: str, resid: int, icode: str = "",
             atom_names: Optional[Set[str]] = None) -> List[_Atom]:
        out = []
        for a in self.atoms:
            if a.chain == chain and a.resid == resid and a.icode == icode:
                if atom_names is None or a.name in atom_names:
                    out.append(a)
        return out


# ---------------------------------------------------------------------------
# Recipe builder (ported from the prototype; role identity by immutable uid).
# ---------------------------------------------------------------------------
class RecipeBuilder:
    def __init__(self, structure: _Structure, redox_site):
        self.structure = structure
        self.redox_site = redox_site
        self.operations: List[dict] = []
        self._undo: List[tuple] = []
        self.parameters: List[dict] = []
        self.reference_state: Dict[str, str] = {}
        self.role_meta: Dict[str, dict] = {}
        self._role_uids: Dict[str, Set[int]] = {}
        self._addr_to_role: Dict[ResidueKey, str] = {}
        self._seed_roles_from_site()

    # ---- role tagging ----------------------------------------------------
    def _seed_roles_from_site(self):
        sig = connectivity_signature(self.redox_site)
        seen: Dict[str, int] = {}
        order = self.structure.residue_order()
        for key in order:
            chain, resid, icode = key
            rn = self.structure.resname_of(chain, resid, icode) or "UNK"
            n_same = sum(1 for k in order
                         if (self.structure.resname_of(*k) or "UNK") == rn)
            seen[rn] = seen.get(rn, 0) + 1
            label = rn.lower() if n_same == 1 else f"{rn.lower()}_{seen[rn]}"
            discs = residue_bond_discriminators(self.redox_site, (chain, resid))
            self.role_meta[label] = {
                "resname": rn,
                "fingerprint": sig.get((chain, resid)),
                "discriminators": [
                    {"source_atom": s, "partner_resname": p, "target_atom": t}
                    for (s, p, t) in sorted(discs)
                ],
            }
            self._role_uids[label] = {a.uid for a in self.structure.residue_atoms(key)}
        self._rebuild_addr_map()

    def _rebuild_addr_map(self):
        self._addr_to_role = {}
        by_uid = {a.uid: a for a in self.structure.atoms}
        for label, uids in self._role_uids.items():
            keys: Dict[ResidueKey, int] = {}
            for u in uids:
                a = by_uid.get(u)
                if a is None:
                    continue
                k = (a.chain, a.resid, a.icode)
                keys[k] = keys.get(k, 0) + 1
            if keys:
                self._addr_to_role[max(keys, key=keys.get)] = label

    def _role_label(self, key: ResidueKey) -> Optional[str]:
        if key in self._addr_to_role:
            return self._addr_to_role[key]
        rn = self.structure.resname_of(*key) or "UNK"
        label = rn.lower()
        i = 1
        while label in self.role_meta and self.role_meta[label]["resname"] != rn:
            i += 1
            label = f"{rn.lower()}_{i}"
        self.role_meta.setdefault(label, {"resname": rn, "fingerprint": None,
                                           "discriminators": []})
        self._addr_to_role[key] = label
        return label

    def _tag(self, selector: dict, key: ResidueKey) -> dict:
        role = self._role_label(key)
        return {"role": role, **selector} if role else selector

    # ---- undo ------------------------------------------------------------
    def _checkpoint(self):
        self._undo.append((copy.deepcopy(self.structure.atoms),
                           copy.deepcopy(self.operations),
                           copy.deepcopy(self.parameters),
                           copy.deepcopy(self.reference_state),
                           copy.deepcopy(self.role_meta),
                           copy.deepcopy(self._role_uids)))

    def undo(self) -> str:
        if not self._undo:
            raise RecipeError("Nothing to undo.")
        (self.structure.atoms, self.operations, self.parameters,
         self.reference_state, self.role_meta, self._role_uids) = self._undo.pop()
        self._rebuild_addr_map()
        return "Reverted last operation."

    def _require_residue(self, chain: str, resid: int, icode: str = "") -> str:
        rn = self.structure.resname_of(chain, resid, icode)
        if rn is None:
            raise RecipeError(f"No residue {chain}/{resid}{icode or ''} in the site.")
        return rn

    # ---- atomic operations ----------------------------------------------
    def rename_residue(self, chain, resid, new_name, icode="") -> str:
        old = self._require_residue(chain, resid, icode)
        self._checkpoint()
        selector = self._tag({"chain": chain, "resid": resid, "icode": icode,
                              "resname": old}, (chain, resid, icode))
        for a in self.structure.find(chain, resid, icode):
            a.resname = new_name
        self.operations.append({"op": "rename_residue", "selector": selector,
                                "action": {"change_residue_name": new_name}})
        self._rebuild_addr_map()
        return f"Renamed residue {chain}/{resid} {old} -> {new_name}"

    def rename_atom(self, chain, resid, old_atom, new_atom, icode="") -> str:
        rn = self._require_residue(chain, resid, icode)
        hits = self.structure.find(chain, resid, icode, {old_atom})
        if not hits:
            raise RecipeError(f"Atom {old_atom} not found in {chain}/{resid} ({rn}).")
        self._checkpoint()
        selector = self._tag({"chain": chain, "resid": resid, "icode": icode,
                              "resname": rn}, (chain, resid, icode))
        for a in hits:
            a.name = new_atom
        self.operations.append({"op": "rename_atom", "selector": selector,
                                "action": {"rename_atoms": {old_atom: new_atom}}})
        self._rebuild_addr_map()
        return f"Renamed atom {old_atom} -> {new_atom} in {chain}/{resid} ({rn})"

    def change_residue_id(self, chain, resid, new_id, icode="") -> str:
        rn = self._require_residue(chain, resid, icode)
        if self.structure.resname_of(chain, new_id, icode) is not None and new_id != resid:
            raise RecipeError(f"Residue {chain}/{new_id} already exists -- would collide. "
                              f"Use 'movenew' for insertion semantics.")
        self._checkpoint()
        selector = self._tag({"chain": chain, "resid": resid, "icode": icode,
                              "resname": rn}, (chain, resid, icode))
        for a in self.structure.find(chain, resid, icode):
            a.resid = new_id
        self.operations.append({"op": "change_residue_id", "selector": selector,
                                "action": {"change_residue_id": new_id}})
        self._rebuild_addr_map()
        return f"Changed residue ID {chain}/{resid} -> {chain}/{new_id} ({rn})"

    def change_chain(self, chain, resid, new_chain, icode="") -> str:
        rn = self._require_residue(chain, resid, icode)
        self._checkpoint()
        selector = self._tag({"chain": chain, "resid": resid, "icode": icode,
                              "resname": rn}, (chain, resid, icode))
        for a in self.structure.find(chain, resid, icode):
            a.chain = new_chain
        self.operations.append({"op": "change_chain", "selector": selector,
                                "action": {"change_chain_id": new_chain}})
        self._rebuild_addr_map()
        return f"Changed chain {chain}/{resid} -> {new_chain}/{resid} ({rn})"

    def change_icode(self, chain, resid, new_icode, icode="") -> str:
        rn = self._require_residue(chain, resid, icode)
        self._checkpoint()
        selector = self._tag({"chain": chain, "resid": resid, "icode": icode,
                              "resname": rn}, (chain, resid, icode))
        for a in self.structure.find(chain, resid, icode):
            a.icode = new_icode
        self.operations.append({"op": "change_insertion_code", "selector": selector,
                                "action": {"change_insertion_code": new_icode}})
        self._rebuild_addr_map()
        return f"Changed insertion code {chain}/{resid} -> '{new_icode}'"

    def set_record_type(self, chain, resid, record, icode="") -> str:
        rn = self._require_residue(chain, resid, icode)
        record = record.upper()
        if record not in ("ATOM", "HETATM"):
            raise RecipeError("Record type must be ATOM or HETATM.")
        self._checkpoint()
        selector = self._tag({"chain": chain, "resid": resid, "icode": icode,
                              "resname": rn}, (chain, resid, icode))
        for a in self.structure.find(chain, resid, icode):
            a.record = record
        action = {"convert_to_hetatm": True} if record == "HETATM" else {"convert_to_atom": True}
        self.operations.append({"op": "set_record_type", "selector": selector, "action": action})
        self._rebuild_addr_map()
        return f"Set {chain}/{resid} ({rn}) record type -> {record}"

    def move_to_new_residue(self, chain, resid, atom_names, new_name, id_offset, icode="") -> str:
        if id_offset < 1:
            raise RecipeError("ID offset for a NEW residue must be >= 1.")
        rn = self._require_residue(chain, resid, icode)
        want = set(atom_names)
        movers = self.structure.find(chain, resid, icode, want)
        found = {a.name for a in movers}
        missing = want - found
        if missing:
            raise RecipeError(f"Atoms not found in {chain}/{resid}: {', '.join(sorted(missing))}")
        if len(found) == len(self.structure.find(chain, resid, icode)):
            raise RecipeError("Cannot move ALL atoms to a new residue -- that is a rename.")
        target_id = resid + id_offset
        mover_uids = {a.uid for a in movers}
        self._checkpoint()
        source_role = self._role_label((chain, resid, icode))
        selector = self._tag({"chain": chain, "resid": resid, "icode": icode,
                              "resname": rn, "atom_names": sorted(found)}, (chain, resid, icode))
        for a in self.structure.atoms:
            if a.uid in mover_uids:
                continue
            if a.chain == chain and a.resid >= target_id:
                a.resid += 1
        for a in movers:
            a.resid = target_id
            a.resname = new_name
            a.icode = ""
        if source_role is not None:
            self._role_uids[source_role] -= mover_uids
        creates_role = new_name.lower()
        i = 1
        while creates_role in self.role_meta:
            i += 1
            creates_role = f"{new_name.lower()}_{i}"
        self.role_meta[creates_role] = {"resname": new_name, "fingerprint": None,
                                        "discriminators": [], "created": True}
        self._role_uids[creates_role] = set(mover_uids)
        self.operations.append({
            "op": "move_to_new_residue", "selector": selector,
            "action": {"change_residue_name": new_name, "change_residue_id": target_id},
            "relative": {"new_resid_offset": id_offset},
            "cascade": {"scope": "chain", "chain": chain, "shift": 1,
                        "applies_to": f"resid >= {target_id}"},
            "creates_role": creates_role,
        })
        self._rebuild_addr_map()
        return (f"Split {len(found)} atom(s) from {chain}/{resid} into new residue "
                f"{new_name} at {chain}/{target_id} (downstream residues cascaded +1).")

    def move_to_existing_residue(self, chain, resid, atom_names, target_chain,
                                 target_resid, icode="", target_icode="") -> str:
        rn = self._require_residue(chain, resid, icode)
        target_rn = self.structure.resname_of(target_chain, target_resid, target_icode)
        if target_rn is None:
            raise RecipeError(f"Target residue {target_chain}/{target_resid} does not exist.")
        want = set(atom_names)
        movers = self.structure.find(chain, resid, icode, want)
        found = {a.name for a in movers}
        missing = want - found
        if missing:
            raise RecipeError(f"Atoms not found: {', '.join(sorted(missing))}")
        target_names = {a.name for a in self.structure.find(target_chain, target_resid, target_icode)}
        clash = found & target_names
        if clash:
            raise RecipeError(f"Target already has atom name(s) {', '.join(sorted(clash))}.")
        self._checkpoint()
        source_role = self._role_label((chain, resid, icode))
        target_role = self._role_label((target_chain, target_resid, target_icode))
        selector = self._tag({"chain": chain, "resid": resid, "icode": icode,
                              "resname": rn, "atom_names": sorted(found)}, (chain, resid, icode))
        mover_uids = {a.uid for a in movers}
        for a in movers:
            a.chain, a.resid, a.icode, a.resname = target_chain, target_resid, target_icode, target_rn
        if source_role is not None:
            self._role_uids[source_role] -= mover_uids
        if target_role is not None:
            self._role_uids[target_role] |= mover_uids
        op = {"op": "move_to_existing_residue", "selector": selector,
              "action": {"change_chain_id": target_chain, "change_residue_id": target_resid,
                         "change_residue_name": target_rn}}
        if target_role is not None:
            op["target_role"] = target_role
        self.operations.append(op)
        self._rebuild_addr_map()
        return (f"Moved {len(found)} atom(s) into existing {target_chain}/{target_resid} ({target_rn}).")

    # ---- state parameters ------------------------------------------------
    def declare_parameter(self, name, options) -> str:
        if not options:
            raise RecipeError("A parameter needs at least one option.")
        if any(p["name"] == name for p in self.parameters):
            raise RecipeError(f"Parameter '{name}' already declared.")
        self._checkpoint()
        self.parameters.append({"name": name, "options": options, "default": options[0]})
        self.reference_state.setdefault(name, options[0])
        return f"Declared parameter {name} = {{{', '.join(options)}}} (reference: {options[0]})"

    def set_reference_state(self, mapping) -> str:
        for k, v in mapping.items():
            p = next((p for p in self.parameters if p["name"] == k), None)
            if p is None:
                raise RecipeError(f"Unknown parameter '{k}'.")
            if v not in p["options"]:
                raise RecipeError(f"'{v}' is not an option of {k}.")
        self._checkpoint()
        self.reference_state.update(mapping)
        return "Reference state: " + "|".join(self._ref_values())

    def _param_order(self):
        return [p["name"] for p in self.parameters]

    def _ref_values(self):
        return [self.reference_state[n] for n in self._param_order()]

    def vary_name(self, op_index, state_values, new_name) -> str:
        if not self.parameters:
            raise RecipeError("Declare parameters first with 'param'.")
        if op_index < 1 or op_index > len(self.operations):
            raise RecipeError(f"No operation #{op_index}.")
        op = self.operations[op_index - 1]
        if "change_residue_name" not in op.get("action", {}):
            raise RecipeError(f"Operation #{op_index} has no residue name to vary.")
        if len(state_values) != len(self.parameters):
            raise RecipeError(f"Expected {len(self.parameters)} state value(s) "
                              f"[{', '.join(self._param_order())}].")
        for val, p in zip(state_values, self.parameters):
            if val not in p["options"]:
                raise RecipeError(f"'{val}' is not an option of {p['name']}.")
        self._checkpoint()
        nbs = op.setdefault("name_by_state", {})
        nbs.setdefault("|".join(self._ref_values()), op["action"]["change_residue_name"])
        nbs["|".join(state_values)] = new_name
        return f"Op #{op_index}: state [{'|'.join(state_values)}] -> {new_name}"

    # ---- serialize -------------------------------------------------------
    def to_recipe(self, name, description="", forcefield=None) -> dict:
        recipe = {
            "schema_version": SCHEMA_VERSION, "name": name,
            "description": description, "source": "table_transformer_creator",
            "roles": [{"label": lbl, **meta} for lbl, meta in self.role_meta.items()],
        }
        if self.parameters:
            recipe["parameters"] = self.parameters
            recipe["parameter_order"] = self._param_order()
            recipe["reference_state"] = dict(self.reference_state)
        recipe["operations"] = self.operations
        if forcefield:
            recipe["forcefield"] = forcefield
        return recipe


# ---------------------------------------------------------------------------
# Rendering + command dispatch
# ---------------------------------------------------------------------------
def render_table(structure: _Structure, highlight: Optional[Set[ResidueKey]] = None) -> str:
    highlight = highlight or set()
    lines = ["  " + "-" * 62,
             f"    {'chain':<5} {'resid':>5} {'resname':<7} {'rec':<6} atoms",
             "  " + "-" * 62]
    for key in structure.residue_order():
        chain, resid, icode = key
        atoms = structure.residue_atoms(key)
        resname = atoms[0].resname
        rec = "HETATM" if any(a.record == "HETATM" for a in atoms) else "ATOM"
        mark = "*" if key in highlight else " "
        rid = f"{resid}{icode}" if icode else f"{resid}"
        names = " ".join(a.name for a in atoms)
        prefix = f"  {mark} {chain:<5} {rid:>5} {resname:<7} {rec:<6} "
        cont = "  " + " " * (len(prefix) - 2)
        first, cur = True, ""
        for w in names.split():
            if len(cur) + len(w) + 1 > 30:
                lines.append((prefix if first else cont) + cur.strip())
                first, cur = False, ""
            cur += w + " "
        lines.append((prefix if first else cont) + cur.strip())
    lines.append("  " + "-" * 62)
    return "\n".join(lines)


def render_summary(builder: RecipeBuilder) -> str:
    if not builder.operations:
        return "  (no operations yet)"
    out = []
    for i, op in enumerate(builder.operations, 1):
        sel = op["selector"]
        loc = f"{sel.get('chain')}/{sel.get('resid')} ({sel.get('resname', '?')})"
        atoms = sel.get("atom_names")
        ap = f" atoms[{','.join(atoms)}]" if atoms else ""
        act = ", ".join(f"{k}={v}" for k, v in op["action"].items())
        out.append(f"  {i:>2}. {op['op']:<22} {loc}{ap} -> {act}")
        for sk, nm in op.get("name_by_state", {}).items():
            out.append(f"        state [{sk}] -> {nm}")
    return "\n".join(out)


def apply_command(builder: RecipeBuilder, tokens: List[str]) -> str:
    verb = tokens[0].lower()

    def need(n, usage):
        if len(tokens) < n:
            raise RecipeError(f"Usage: {usage}")

    if verb in ("rename_res", "renameres"):
        need(4, "rename_res <chain> <resid> <NEWNAME>")
        return builder.rename_residue(tokens[1], int(tokens[2]), tokens[3].upper())
    if verb in ("rename_atom", "renameatom"):
        need(5, "rename_atom <chain> <resid> <OLD> <NEW>")
        return builder.rename_atom(tokens[1], int(tokens[2]), tokens[3].upper(), tokens[4].upper())
    if verb == "id":
        need(4, "id <chain> <resid> <NEWID>")
        return builder.change_residue_id(tokens[1], int(tokens[2]), int(tokens[3]))
    if verb == "chain":
        need(4, "chain <chain> <resid> <NEWCHAIN>")
        return builder.change_chain(tokens[1], int(tokens[2]), tokens[3])
    if verb == "icode":
        need(4, "icode <chain> <resid> <CODE>")
        return builder.change_icode(tokens[1], int(tokens[2]), tokens[3])
    if verb in ("hetatm", "atom"):
        need(3, f"{verb} <chain> <resid>")
        return builder.set_record_type(tokens[1], int(tokens[2]), verb)
    if verb == "movenew":
        need(6, "movenew <chain> <resid> <A,B,C> <NEWNAME> <ID_OFFSET>")
        atoms = [a.upper() for a in tokens[3].split(",") if a]
        return builder.move_to_new_residue(tokens[1], int(tokens[2]), atoms,
                                           tokens[4].upper(), int(tokens[5]))
    if verb == "move":
        need(6, "move <chain> <resid> <A,B,C> <TCHAIN> <TRESID>")
        atoms = [a.upper() for a in tokens[3].split(",") if a]
        return builder.move_to_existing_residue(tokens[1], int(tokens[2]), atoms,
                                                tokens[4], int(tokens[5]))
    if verb == "param":
        need(3, "param <NAME> <opt1,opt2,...>")
        return builder.declare_parameter(tokens[1], [o for o in tokens[2].split(",") if o])
    if verb == "refstate":
        need(2, "refstate <param=value> ...")
        mapping = {}
        for kv in tokens[1:]:
            if "=" not in kv:
                raise RecipeError("refstate args look like param=value")
            k, v = kv.split("=", 1)
            mapping[k] = v
        return builder.set_reference_state(mapping)
    if verb == "vary":
        need(4, "vary <op#> <val1,val2,...> <NEWNAME>")
        return builder.vary_name(int(tokens[1]), [v for v in tokens[2].split(",") if v],
                                 tokens[3].upper())
    if verb == "undo":
        return builder.undo()
    raise RecipeError(f"Unknown command '{verb}'. Type 'help'.")


_HELP = """
Commands (one operation at a time; the table redraws after each):
  rename_res  <chain> <resid> <NEWNAME>                 change residue name
  rename_atom <chain> <resid> <OLD> <NEW>               rename one atom
  id          <chain> <resid> <NEWID>                   change residue ID
  chain       <chain> <resid> <NEWCHAIN>                change chain ID
  icode       <chain> <resid> <CODE>                    set insertion code
  hetatm      <chain> <resid>                           mark residue as HETATM
  atom        <chain> <resid>                           mark residue as ATOM
  movenew     <chain> <resid> <A,B,C> <NEWNAME> <OFF>   split atoms into a NEW residue
  move        <chain> <resid> <A,B,C> <TCHAIN> <TRESID> move atoms into an EXISTING residue
  param       <NAME> <opt1,opt2,...>                    declare a state axis (redox/spin)
  refstate    <param=value> ...                         set the authoring state
  vary        <op#> <val1,...> <NEWNAME>                give a rename a per-state name
  undo | show | summary | save | quit
"""

_MODULE = "Table Transformer Creator"


class TableTransformerCreator:
    """Interactive, table-driven authoring of a JSON transformer spec."""

    def __init__(self, processor, redox_site):
        self.processor = processor
        self.redox_site = redox_site
        self.console = processor.console

    def create(self) -> Optional[str]:
        if self.redox_site is None or not getattr(self.redox_site, "atoms", None):
            self.console.print("[red]No redox site atoms to edit.[/red]")
            return None
        structure = _Structure.from_redox_site(self.redox_site)
        builder = RecipeBuilder(structure, self.redox_site)

        self.console.print("\n[bold cyan]Transformer Creator (interactive PDB editor)[/bold cyan]")
        self.console.print("[grey50]Apply one edit at a time; order matters. Type 'help' for "
                           "commands, 'save' to finish.[/grey50]")
        self.console.print(render_table(structure))

        while True:
            raw = prompt_with_context(
                self.processor, "edit", module=_MODULE,
                description="PDB editing command (help/show/summary/save/quit)",
            ).strip()
            if not raw:
                continue
            tokens = raw.split()
            verb = tokens[0].lower()
            if verb == "help":
                self.console.print(_HELP)
                continue
            if verb == "show":
                self.console.print(render_table(structure))
                continue
            if verb == "summary":
                self.console.print(render_summary(builder))
                continue
            if verb == "quit":
                if confirm_with_context(self.processor, "Discard this transformer?",
                                        module=_MODULE, default=False):
                    self.console.print("[yellow]Discarded.[/yellow]")
                    return None
                continue
            if verb == "save":
                return self._save(builder)
            try:
                msg = apply_command(builder, tokens)
                self.console.print(f"  [green]ok:[/green] {msg}")
                self.console.print(render_table(structure))
            except (RecipeError, ValueError) as e:
                self.console.print(f"  [red]! {e}[/red]")

    def _save(self, builder: RecipeBuilder) -> Optional[str]:
        if not builder.operations:
            self.console.print("[yellow]No operations recorded; nothing to save.[/yellow]")
            return None
        self.console.print("\n[bold]Operation sequence:[/bold]")
        self.console.print(render_summary(builder))
        name = prompt_with_context(self.processor, "Template name",
                                   module=_MODULE, description="transformer name").strip()
        if not name:
            self.console.print("[yellow]No name given; not saved.[/yellow]")
            return None
        token = _sanitize(name)
        ff = None
        if confirm_with_context(self.processor,
                                "Link a deposited force-field library (path under "
                                "specialized_residues)?", module=_MODULE, default=False):
            path = prompt_with_context(self.processor, "Forcefield path (e.g. heme/bis_his_c_type)",
                                       module=_MODULE).strip()
            if path:
                ff = {"path": path}
                redox = prompt_with_context(self.processor, "redox_state (blank to skip)",
                                            module=_MODULE, default="").strip()
                spin = prompt_with_context(self.processor, "spin_state (blank to skip)",
                                           module=_MODULE, default="").strip()
                if redox:
                    ff["redox_state"] = redox
                if spin:
                    ff["spin_state"] = spin

        recipe = builder.to_recipe(token, description=name, forcefield=ff)
        DEFAULT_USER_TRANSFORMER_DIR.mkdir(parents=True, exist_ok=True)
        out_path = DEFAULT_USER_TRANSFORMER_DIR / f"{token}.json"
        with open(out_path, "w") as f:
            json.dump(recipe, f, indent=2)

        # Register immediately so it is usable in this session.
        try:
            from proprep.redoxsite_prep.transformation.spec_transformer import (
                load_user_spec_transformers,
            )
            load_user_spec_transformers()
        except Exception as exc:  # pragma: no cover - defensive
            logger.warning("Could not register new spec immediately: %s", exc)

        self.console.print(f"\n  [green]✓ Transformer spec saved:[/green] {out_path}")
        self.console.print("  [grey50]Reusable on any protein with this site type; discovered "
                           "automatically by the Redox Site Preparer.[/grey50]")
        return str(out_path)


def _sanitize(name: str) -> str:
    import re
    token = re.sub(r"[^0-9a-zA-Z]+", "_", name).strip("_").lower()
    if not token:
        token = "custom_transformer"
    if token[0].isdigit():
        token = f"t_{token}"
    return token
