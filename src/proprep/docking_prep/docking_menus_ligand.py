"""
Molecular Docking menus 2 and 7-9: the ligand, its rotatable bonds, the search box, the settings.

6 ligand: a residue of the structure (CCD chemistry on the crystal pose), a
  SMILES string, or a file loaded with the Structure Loader; its chemistry is
  shown and can be edited with typed commands
7 rotatable bonds: every switchable bond listed with its kind; hold any rigid,
  make an amide rotatable, keep or open a macrocycle
8 search box: around the ligand, around residues, or typed
9 settings: Vina's documented defaults, shown and changeable; the seed must not
  be 0, which Vina reads as "random"
"""

from __future__ import annotations

import os
from typing import List

from rdkit import Chem

from . import docking_ui as ui
from .ccd_chemistry import is_metal
from .docking_run import DockingSettings, SCORING_FUNCTIONS, box_around
from .docking_state import LigandChoices
from .receptor_prep import read_atoms

_WATER = {"HOH", "WAT", "DOD"}
SOURCES = {"s": "A residue of the receptor structure", "m": "A SMILES string", "f": "A file from the Structure Loader"}
VINA_DEFAULTS = dict(scoring="vina", seed=42, exhaustiveness=8, n_poses=9, energy_range=3.0,
                     min_rmsd=1.0, max_evals=0, spacing=0.375, cpu=0)


# ---------------------------------------------------------------- 6 ligand ---

def choose_ligand(m) -> None:
    from . import docking_menus_receptor as receptor_menus
    from .docking_menus_receptor import review_ccd_chemistry
    state = m.state
    key = ui.ask(m.processor, "Ligand source: [s] residue of the structure, [m] SMILES, [f] file",
                 choices=list(SOURCES), default="s", description="Ligand source", options_map=SOURCES)
    choices = LigandChoices()
    if key == "s":
        pdb = state.receptor.pdb_path
        if not pdb:
            m.console.print("[dark_orange3]Choose the receptor structure first (1); the ligand is one of its "
                            "residues.[/dark_orange3]")
            return
        residues = {}
        for atom in read_atoms(pdb):
            if atom.record == "HETATM" and atom.resname not in _WATER:
                residues.setdefault(atom.res_id, {"resname": atom.resname, "atom": atom, "n": 0,
                                                  "altlocs": set()})
                residues[atom.res_id]["n"] += 1
                if atom.altloc:
                    residues[atom.res_id]["altlocs"].add(atom.altloc)
        candidates = [(rid, info) for rid, info in residues.items() if info["n"] > 1]
        if not candidates:
            m.console.print("[dark_orange3]The structure has no multi-atom HETATM residue.[/dark_orange3]")
            return
        m.console.print(ui.table("Residues that can be the ligand", ["#", "residue", "name", "atoms"],
                                 [(n, rid, info["resname"], info["n"]) for n, (rid, info) in enumerate(candidates, 1)]))
        number = ui.ask_int(m.processor, "Ligand number", default=1, min_value=1, max_value=len(candidates),
                            description="Ligand residue")
        rid, info = candidates[number - 1]
        atom = info["atom"]
        altloc, altloc_fills = None, {}
        if info["altlocs"]:
            from .docking_menus_receptor import alternates_by_residue, pick_alternates
            earlier = state.ligand.source if state.ligand.source.get("kind") == "structure_residue" else {}
            same = f"{earlier.get('chain')}:{earlier.get('resseq')}{earlier.get('icode', '')}" == rid
            altloc, altloc_fills = pick_alternates(m, alternates_by_residue(pdb, [rid]),
                                                   {rid: earlier.get("altloc")} if same else {})[rid]
        store = {}
        earlier = state.receptor.cofactor_edits.get(info["resname"])
        if earlier and not any(state.receptor.pdb_path and k != rid and v["resname"] == info["resname"]
                               for k, v in residues.items() if k in state.receptor.keep_hetero):
            store[info["resname"]] = state.receptor.cofactor_edits.pop(info["resname"])
            m.console.print(f"  The edits made to {info['resname']} as a cofactor (5) now apply to the ligand.",
                            highlight=False)
        # The CCD entry is reviewed here only when it cannot be used as it is
        # (an internally inconsistent entry cannot be built at all). A sound
        # entry goes straight to the built ligand, whose one edit prompt
        # covers protonation and bond orders.
        from .chemistry_edits import ChemistryEdit, apply_edits
        component = receptor_menus.fetch_component(info["resname"])
        _, problems = apply_edits(component.mol, [ChemistryEdit.from_dict(e) for e in store.get(info["resname"], [])])
        if problems:
            review_ccd_chemistry(m, info["resname"], store, rid, altloc, altloc_fills)
        else:
            m.console.print(f"\n  {info['resname']}: {ui.escape(component.name)}", highlight=False)
        choices.source = {"kind": "structure_residue", "path": os.path.abspath(pdb), "chain": atom.chain,
                          "resseq": atom.resseq, "icode": atom.icode, "resname": info["resname"],
                          "altloc": altloc, "altloc_fills": altloc_fills, "ccd_code": info["resname"]}
        choices.ccd_edits = store.get(info["resname"], [])
        if rid in state.receptor.keep_hetero:
            state.receptor.keep_hetero.remove(rid)
            m.console.print(f"  {rid} is no longer part of the receptor.", highlight=False)
    elif key == "m":
        smiles = ui.ask(m.processor, "SMILES", description="Ligand SMILES")
        seed = ui.ask_int(m.processor, "Seed for the 3D build (RDKit ETKDG)", default=42, min_value=1,
                          description="ETKDG seed")
        optimizer = ui.ask(m.processor, "Minimise the 3D build with MMFF94?",
                           choices=["mmff94", "none"], default="mmff94", description="3D build minimisation")
        choices.source = {"kind": "smiles", "smiles": smiles, "embed_seed": seed, "optimizer": optimizer,
                          "residue_name": "LIG"}
    else:
        path = m.processor._get_workspace().get("ligand_file")
        if not path:
            m.console.print("[dark_orange3]No ligand file is loaded. Load one with the Structure Loader "
                            "(option 4), then come back.[/dark_orange3]")
            return
        record = ui.ask_int(m.processor, "Molecule number in the file (0 for the first)", default=0,
                            min_value=0, description="Record in the ligand file")
        choices.source = {"kind": "file", "path": path, "record": record, "residue_name": "LIG"}
    state.ligand = choices
    state.box = {} if key != "s" else state.box
    try:
        build = m.build_ligand()
    except Exception as error:                          # the source is recorded; say what is wrong
        m.console.print(f"[red]The ligand could not be built: {error}[/red]")
        m.save_state()
        return
    edit_ligand(m, build)
    m.save_state()


def show_ligand(m, build) -> None:
    mol = build.mol
    heavy = sum(a.GetAtomicNum() > 1 for a in mol.GetAtoms())
    charged = [f"{a.GetProp('name')} {a.GetFormalCharge():+d}" for a in mol.GetAtoms() if a.GetFormalCharge()]
    m.console.print(f"\n  Ligand: {heavy} heavy atoms, {mol.GetNumAtoms() - heavy} hydrogens, net charge "
                    f"{Chem.GetFormalCharge(mol):+d}; charged atoms: {', '.join(charged) or 'none'}", highlight=False)
    m.console.print(f"  SMILES: {Chem.MolToSmiles(Chem.RemoveHs(mol))}", highlight=False)
    ui.print_notes(m.console, build.notes)
    m.show_ligand_in_viewer(build)


def edit_ligand(m, build) -> None:
    """Typed edits to the built ligand (protonation, bond orders), then metal charges."""
    from .chemistry_edits import apply_edits
    state = m.state
    while True:
        show_ligand(m, build)
        answer = ui.ask(m.processor, f"Edit for the ligand ({ui.EDIT_HELP}; 'p' to pick in the viewer; "
                        f"{ui.DROP_HELP}; Enter when done)", default="", description="Chemistry edit for the ligand")
        if not answer:
            break
        if answer.lower() == "p":
            answer = _picked_edit(m, build)
            if not answer:
                continue
        if answer.lower() == ui.DROP:
            if not state.ligand.edits:
                m.console.print("  No edit to drop.")
                continue
            from .chemistry_edits import ChemistryEdit
            dropped = state.ligand.edits.pop()
            m.console.print(f"  Dropped: {ui.escape(ChemistryEdit.from_dict(dropped).describe())}", highlight=False)
        else:
            try:
                edit = ui.parse_edit(answer)
                trial = state.ligand.edit_objects() + [edit]
                apply_edits(m.build_ligand(apply_edits_=False).mol, trial)
                state.ligand.edits.append(edit.to_dict())
            except (ValueError, KeyError) as error:
                m.console.print(f"  [red]{error}[/red]")
                continue
        build = m.build_ligand()
    metals = [a.GetProp("name") for a in build.mol.GetAtoms() if is_metal(a)]
    for name in metals:
        state.ligand.metal_charges[name] = ui.ask_required_int(
            m.processor, m.console, f"Formal charge of ligand metal {name}", description="Ligand metal formal charge")
    state.ligand.rigid_macrocycles = None             # torsions must be reviewed for this ligand (7)


PICK_KINDS = {"a": "An atom (change its charge, remove it, add a hydrogen)", "b": "A bond (change its order)"}
ATOM_ACTIONS = {"c": "Set its formal charge", "r": "Remove it", "h": "Add a hydrogen to it"}


def _picked_edit(m, build) -> str:
    """Build a typed edit command from a click on the ligand (structure 1 of the viewer)."""
    kind = ui.ask(m.processor, "Pick [a] an atom or [b] a bond", choices=list(PICK_KINDS), default="a",
                  description="What to pick in the viewer", options_map=PICK_KINDS)
    if kind == "b":
        result = ui.pick(m, "bond", "Click the ligand bond to change", 1)
        if not result:
            return ""
        names = [build.mol.GetAtomWithIdx(result[k]["index"]).GetProp("name") for k in ("atom1", "atom2")]
        order = ui.ask(m.processor, f"Bond order for {names[0]}-{names[1]}", choices=["1", "2", "3"],
                       description="Bond order for the picked bond")
        return f"bond {names[0]} {names[1]} {order}"
    result = ui.pick(m, "atom", "Click the ligand atom to change", 1)
    if not result:
        return ""
    name = build.mol.GetAtomWithIdx(result["atom"]["index"]).GetProp("name")
    action = ui.ask(m.processor, f"For {name}: [c] formal charge, [r] remove, [h] add a hydrogen",
                    choices=list(ATOM_ACTIONS), description="Change to the picked atom", options_map=ATOM_ACTIONS)
    if action == "c":
        value = ui.ask_required_int(m.processor, m.console, f"Formal charge of {name}",
                                    description="Formal charge for the picked atom")
        return f"charge {name} {value}"
    return f"remove {name}" if action == "r" else f"add-h {name}"


# --------------------------------------------------------------- 7 torsions ---

def choose_torsions(m) -> None:
    from .ligand_prep import torsion_bonds
    state = m.state
    if not state.ligand.source:
        m.console.print("[dark_orange3]Choose a ligand first (2).[/dark_orange3]")
        return
    from meeko.macrocycle import DEFAULT_MAX_RING_SIZE, DEFAULT_MIN_RING_SIZE
    build = m.build_ligand()
    ring_sizes = [len(ring) for ring in build.mol.GetRingInfo().AtomRings()]
    macro = [n for n in ring_sizes if DEFAULT_MIN_RING_SIZE <= n <= DEFAULT_MAX_RING_SIZE]
    if macro:
        m.console.print(f"  The ligand has a ring of {max(macro)} atoms; Meeko treats rings of "
                        f"{DEFAULT_MIN_RING_SIZE} to {DEFAULT_MAX_RING_SIZE} atoms as macrocycles it can open.",
                        highlight=False)
        opened = ui.confirm(m.processor, "Let Meeko open the macrocycle so the ring flexes during docking?",
                            default=True, description="Macrocycle flexibility")
        state.ligand.rigid_macrocycles = not opened
    else:
        state.ligand.rigid_macrocycles = False       # no macrocycle: the setting changes nothing
    metal_charges = dict(state.ligand.metal_charges) or None
    bonds = torsion_bonds(build.mol, rigid_macrocycles=state.ligand.rigid_macrocycles, metal_charges=metal_charges)
    m.show_ligand_in_viewer(build)
    rigid = {frozenset(b) for b in state.ligand.rigid_bonds}
    amides = {frozenset(b) for b in state.ligand.rotatable_amides}
    while True:
        rows = []
        for n, bond in enumerate(bonds, 1):
            key = frozenset(bond.atoms)
            rotates = (bond.rotatable_by_default and key not in rigid) or key in amides
            rows.append((n, bond.label(), bond.kind, "rotates" if rotates else "rigid"))
        m.console.print(ui.table("Bonds that can rotate", ["#", "bond", "kind", "now"], rows))
        count = sum(1 for r in rows if r[3] == "rotates")
        m.console.print(f"  {count} rotatable bonds.", highlight=False)
        answer = ui.ask(m.processor, "Numbers to switch between rotating and rigid (e.g. 2,5-7), "
                        "'c' to hold every conjugated single bond rigid, 'p' to pick a bond in the viewer, "
                        "Enter when done", default="", description="Rotatable bond choices")
        if not answer:
            break
        if answer.lower() == "p":
            result = ui.pick(m, "bond", "Click the ligand bond to switch between rotating and rigid", 1)
            if not result:
                continue
            picked = frozenset(build.mol.GetAtomWithIdx(result[k]["index"]).GetProp("name") for k in ("atom1", "atom2"))
            numbered = [n for n, bond in enumerate(bonds, 1) if frozenset(bond.atoms) == picked]
            if not numbered:
                m.console.print(f"  [dark_orange3]{'-'.join(sorted(picked))} cannot rotate (ring, double bond, "
                                "or it moves nothing).[/dark_orange3]", highlight=False)
                continue
            answer = str(numbered[0])
        if answer.lower() == "c":
            rigid |= {frozenset(b.atoms) for b in bonds if b.kind == "conjugated single"}
            continue
        try:
            numbers = ui.parse_numbers(answer, len(bonds))
        except ValueError as error:
            m.console.print(f"[red]{error}[/red]")
            continue
        for n in numbers:
            bond = bonds[n - 1]
            key = frozenset(bond.atoms)
            if bond.kind == "amide":
                amides ^= {key}
            else:
                rigid ^= {key}
    state.ligand.rigid_bonds = [sorted(b) for b in rigid]
    state.ligand.rotatable_amides = [sorted(b) for b in amides]
    m.save_state()


# -------------------------------------------------------------------- 8 box ---

BOX_MODES = {"l": "Around the ligand as placed", "r": "Around residues", "c": "Typed centre and size"}


def choose_box(m) -> None:
    import numpy as np
    state = m.state
    mode = ui.ask(m.processor, "Search box: [l] around the ligand, [r] around residues, [c] typed centre and size",
                  choices=list(BOX_MODES), default="l", description="Search box definition", options_map=BOX_MODES)
    if mode == "c":
        center = [ui.ask_float(m.processor, f"Box centre {axis} (A)", default=0.0, description=f"Box centre {axis}")
                  for axis in "xyz"]
        size = [ui.ask_float(m.processor, f"Box edge {axis} (A)", default=20.0, min_value=1.0,
                             description=f"Box edge {axis}") for axis in "xyz"]
        state.box = {"center": center, "size": size, "source": "typed"}
    else:
        if mode == "l":
            if not state.ligand.source:
                m.console.print("[dark_orange3]Choose a ligand first (2).[/dark_orange3]")
                return
            coords = Chem.RemoveHs(m.build_ligand().mol).GetConformer().GetPositions()
            where = "the ligand as placed"
            if state.ligand.source.get("kind") == "smiles":
                m.console.print("  [dark_orange3]A ligand built from SMILES sits wherever RDKit put it, not in "
                                "the receptor; a box around it is rarely what you want.[/dark_orange3]")
        else:
            if not state.receptor.pdb_path:
                m.console.print("[dark_orange3]Choose the receptor structure first (1).[/dark_orange3]")
                return
            answer = ui.ask(m.processor, "Residues to enclose (e.g. A:43,A:88,A:128)", description="Box residues")
            wanted = set(ui.parse_list(answer))
            atoms = [a for a in read_atoms(state.receptor.pdb_path) if a.res_id in wanted]
            missing = wanted - {a.res_id for a in atoms}
            if missing or not atoms:
                m.console.print(f"[red]Not in the structure: {', '.join(sorted(missing)) or answer}.[/red]")
                return
            coords = np.array([a.xyz for a in atoms])
            where = "residues " + ", ".join(sorted(wanted))
        padding = ui.ask_float(m.processor, "Padding around them on every side (A)", default=8.0, min_value=0.0,
                               description="Search box padding")
        state.box = box_around(coords, padding, where).to_dict()
    size = state.box["size"]
    m.console.print(f"  Box centre {tuple(round(x, 2) for x in state.box['center'])}, edges "
                    f"{' x '.join(f'{x:.1f}' for x in size)} A ({size[0] * size[1] * size[2]:.0f} A^3); "
                    f"{state.box['source']}.", highlight=False)
    m.show_box_in_viewer()
    m.save_state()


# --------------------------------------------------------------- 9 settings ---

SCORING = {"vina": "AutoDock Vina", "vinardo": "Vinardo", "ad4": "AutoDock4 (maps from autogrid4)"}


SCORING_GUIDE = (
    "  Scoring functions (the search is the same; only how a pose is scored differs):\n"
    "  vina     AutoDock Vina's empirical function (Trott & Olson, 2010): steric, hydrophobic and\n"
    "           hydrogen-bond terms. No partial charges and no electrostatics. Fast.\n"
    "  vinardo  A re-parameterisation of Vina's terms (Quiroga & Villarreal, 2016); same speed.\n"
    "  ad4      AutoDock 4's force-field function (Huey et al., 2007): van der Waals, hydrogen bonds,\n"
    "           desolvation and electrostatics from the Gasteiger charges, so charges and metal ions\n"
    "           count. Needs autogrid4 maps, built once per receptor and box, so it is slower. In\n"
    "           ProPrep's validation, acetazolamide on carbonic anhydrase's zinc redocked within\n"
    "           1.9-2.8 A with ad4 and the sulfonamide as its anion, against 5.2 A with vina; with the\n"
    "           neutral form ad4 did no better.\n"
    "  Scores are in kcal/mol for all three but are not comparable between functions."
)


def choose_settings(m) -> None:
    from .dependencies import find_executable
    current = dict(VINA_DEFAULTS)
    current.update(m.state.settings)
    s = {}
    m.console.print(ui.escape(SCORING_GUIDE), highlight=False)
    if find_executable("autogrid4") is None:
        m.console.print("  [dark_orange3]ad4 needs autogrid4, which is not installed.[/dark_orange3]")
    s["scoring"] = ui.ask(m.processor, "Scoring function: vina, vinardo or ad4", choices=list(SCORING),
                          default=VINA_DEFAULTS["scoring"], description="Scoring function", options_map=SCORING)
    s["seed"] = ui.ask_int(m.processor, "Random seed (not 0: Vina reads 0 as a random seed)", default=42, min_value=1,
                           description="Vina seed")
    s["exhaustiveness"] = ui.ask_int(m.processor, "Exhaustiveness", default=8, min_value=1,
                                     description="Vina exhaustiveness")
    s["n_poses"] = ui.ask_int(m.processor, "Number of poses", default=9, min_value=1, description="Poses returned")
    s["energy_range"] = ui.ask_float(m.processor, "Energy range above the best pose (kcal/mol)", default=3.0,
                                     min_value=0.1, description="Pose energy range")
    s["min_rmsd"] = ui.ask_float(m.processor, "Minimum RMSD between poses (A)", default=1.0, min_value=0.0,
                                 description="Minimum RMSD between poses")
    s["spacing"] = ui.ask_float(m.processor, "Grid spacing (A)", default=0.375, min_value=0.1, max_value=1.0,
                                description="Map grid spacing")
    s["cpu"] = ui.ask_int(m.processor, "CPU cores (0 for all)", default=0, min_value=0, description="CPU cores")
    s["max_evals"] = 0
    DockingSettings(**s).validate()
    m.state.settings = s
    m.save_state()
