"""
Molecular Docking menus 1 and 3-6: the receptor (2, the ligand, is in docking_menus_ligand).

1 structure and chains (and alternate locations)
3 which HETATM residues stay (cofactors, metals, waters)
4 protonation and chain ends (pdb2pqr/PROPKA proposals, metal overrides; REMARK 465, SEQRES, OXT)
5 metals and cofactor chemistry (charges, AutoDock types, CCD problems and edits)
6 flexible side chains

Each step changes the module's DockingState and saves it; a change that makes
later decisions stale clears them and says so.
"""

from __future__ import annotations

import os
from collections import OrderedDict
from typing import Dict, List, Optional

from rich.panel import Panel

from . import docking_ui as ui
from .ccd_chemistry import fetch_component, is_metal_symbol
from .docking_state import metal_key
from .receptor_prep import read_atoms, _incomplete_standard_residues, _meeko_default_keys

_WATER = {"HOH", "WAT", "DOD"}
TERMINUS_CHOICES = {"n": "N-terminus", "c": "C-terminus", "b": "chain break"}
# What each treatment builds (Meeko's N-/C-prefixed templates, or the plain one).
TERMINUS_MEANING = {"N-terminus": "charged NH3+ (the chain starts here)",
                    "C-terminus": "charged COO- (the chain ends here)",
                    "chain break": "neutral, as if the chain went on (residues missing beyond it)"}
END_CHOICES = {"N": {"n": "N-terminus: charged NH3+ (the chain starts here)",
                     "b": "chain break: neutral, as if the chain went on"},
               "C": {"c": "C-terminus: charged COO- (the chain ends here)",
                     "b": "chain break: neutral, as if the chain went on"}}


def _ask_end(m, end, *, required: bool) -> None:
    """One chain end's treatment; only the two that make sense for a first or a last residue."""
    r = m.state.receptor
    which = "first" if end.end == "N" else "last"
    chain = end.res_id.split(":", 1)[0]
    choices = END_CHOICES[end.end]
    letters = "/".join(f"[{k}] {v.split(':')[0]}" for k, v in choices.items())
    current = {v: k for k, v in TERMINUS_CHOICES.items()}.get(r.termini.get(end.res_id))
    key = ui.ask(m.processor, f"{end.res_id} {end.resname}, the {which} residue of chain {chain}: {letters}",
                 choices=list(choices), default=None if required else current,
                 description="Treatment of one chain end", options_map=choices)
    r.termini[end.res_id] = TERMINUS_CHOICES[key]
    r.termini_sources[end.res_id] = "chosen"


def _decide_open_ends(m) -> None:
    """Ask about every uncapped chain end the file gives no evidence for: Meeko must not guess."""
    from .receptor_decisions import chain_ends
    r = m.state.receptor
    open_ends = [e for e in chain_ends(r.pdb_path, r.chains) if not e.capped and e.res_id not in r.termini]
    if not open_ends:
        return
    m.console.print(f"\n  Nothing in the file (REMARK 465, SEQRES, or an OXT) says whether these chain ends are the "
                    f"protein's own ends or breaks: {', '.join(e.res_id for e in open_ends)}.", highlight=False)
    for end in open_ends:
        _ask_end(m, end, required=True)


def _residues(pdb_path: str) -> "OrderedDict[str, dict]":
    """res_id -> {resname, record, chain, atoms, altlocs, metals}."""
    out: "OrderedDict[str, dict]" = OrderedDict()
    for atom in read_atoms(pdb_path):
        entry = out.setdefault(atom.res_id, {"resname": atom.resname, "record": atom.record, "chain": atom.chain,
                                             "atoms": 0, "altlocs": set(), "metals": []})
        entry["atoms"] += 1
        if atom.altloc:
            entry["altlocs"].add(atom.altloc)
        if atom.element and is_metal_symbol(atom.element) and atom.name not in entry["metals"]:
            entry["metals"].append(atom.name)
    return out


# ------------------------------------------------------------ 1 structure ---

def choose_structure(m) -> None:
    from proprep.utils.structure_selector import StructureSelector
    console, state = m.console, m.state
    path = StructureSelector(m.processor._get_workspace(), console).get_structure(silent=True)
    if not path:
        console.print("[dark_orange3]No structure is loaded. Load one with the Structure Loader.[/dark_orange3]")
        return
    console.print(f"\n  Current structure: {path}")
    if state.receptor.pdb_path and os.path.abspath(state.receptor.pdb_path) != os.path.abspath(path):
        console.print(f"  [dark_orange3]The receptor was {state.receptor.pdb_path}; switching clears the "
                      "receptor, ligand and box choices.[/dark_orange3]")
    if not ui.confirm(m.processor, "Use this structure as the receptor?", default=True,
                      description="Receptor structure"):
        return
    if state.receptor.pdb_path != path:
        from .docking_state import LigandChoices, ReceptorChoices
        state.receptor = ReceptorChoices(pdb_path=path)
        state.ligand = LigandChoices()
        state.box = {}
    residues = _residues(path)
    chains = list(OrderedDict.fromkeys(r["chain"] for r in residues.values() if r["record"] == "ATOM"))
    counts = {c: sum(1 for r in residues.values() if r["record"] == "ATOM" and r["chain"] == c) for c in chains}
    console.print("  Chains: " + ", ".join(f"{c} ({counts[c]} residues)" for c in chains), highlight=False)
    while True:
        answer = ui.ask(m.processor, "Chains to include (e.g. A or A,B; Enter for all)", default="all",
                        description="Receptor chains")
        chosen = chains if answer.lower() in ("all", "") else [c.upper() for c in ui.parse_list(answer)]
        unknown = [c for c in chosen if c not in chains]
        if not unknown:
            break
        console.print(f"[red]No chain {', '.join(unknown)}.[/red]")
    state.receptor.chains = chosen
    _check_complete(m, residues)
    choose_altlocs(m, residues, revisit=True)
    m.save_state()


def _check_complete(m, residues) -> None:
    by_residue = OrderedDict()
    for atom in read_atoms(m.state.receptor.pdb_path):
        if atom.record == "ATOM" and atom.chain in m.state.receptor.chains:
            by_residue.setdefault(atom.res_id, []).append(atom)
    incomplete = _incomplete_standard_residues(by_residue)
    if incomplete:
        m.console.print(Panel("\n".join(incomplete), title="Residues missing heavy atoms", border_style="dark_orange3"))
        m.console.print("  [dark_orange3]Docking cannot type these. Rebuild them with the Structure Fixer "
                        "first, then come back.[/dark_orange3]")


def alternates_by_residue(pdb_path: str, res_ids) -> "OrderedDict[str, object]":
    """res_id -> altloc_picker.ResidueAlternates for the residues in ``res_ids`` that have alternates."""
    from proprep.structure_prep.altloc_picker import residues_from_records
    wanted = set(res_ids)
    found = residues_from_records((a.chain, a.resseq, a.icode, a.resname, a.name, a.altloc, a.occupancy)
                                  for a in read_atoms(pdb_path) if a.res_id in wanted)
    return OrderedDict((f"{c}:{n}{i}", alt) for (c, n, i), alt in found.items())


def pick_alternates(m, alternates, earlier: Dict[str, str]) -> Dict[str, tuple]:
    """Ask, residue by residue, which alternate to keep; res_id -> (letter, fill plan).

    The Structure Fixer's picker: each alternate's occupancy, how many atoms it
    models, a partial alternate flagged with the atoms it takes from another,
    and with the viewer open each residue shown with its alternates coloured.
    A residue with a single alternate letter has nothing to choose; it is
    kept and listed.
    """
    from proprep.structure_prep import altloc_picker
    chosen: Dict[str, tuple] = {}
    lone = [(rid, alt) for rid, alt in alternates.items() if len(alt.letters) == 1]
    for rid, alt in lone:
        chosen[rid] = (alt.letters[0], {})
    if lone:
        m.console.print("  Only one alternate location is present, so it is kept: "
                        + ", ".join(f"{rid} ({alt.letters[0]})" for rid, alt in lone) + ".", highlight=False)
    to_ask = [(rid, alt) for rid, alt in alternates.items() if len(alt.letters) > 1]
    if not to_ask:
        return chosen
    m.console.print(f"\n  {len(to_ask)} residue(s) have alternate locations; choose one for each.", highlight=False)
    viewer_state = None
    if altloc_picker.ask_to_open_viewer(m.processor, ui.MODULE_NAME):
        atoms = read_atoms(m.state.receptor.pdb_path)
        viewer_state = altloc_picker.open_viewer(m.console, m.state.receptor.pdb_path, [a.xyz for a in atoms],
                                                 [(a.chain, a.resseq, a.icode) for a in atoms])
    try:
        for rid, alt in to_ask:
            chosen[rid] = altloc_picker.pick(m.processor, m.console, alt, viewer_state, module=ui.MODULE_NAME,
                                             default_letter=earlier.get(rid))
    finally:
        altloc_picker.close_viewer(viewer_state)
    return chosen


def choose_altlocs(m, residues=None, *, revisit: bool = False) -> None:
    """Choose an alternate location for every receptor residue that has them.

    ``revisit`` (menu 1) asks again about every such residue, offering the
    earlier choice as the default; otherwise (after 2 changes which HETATM
    residues stay) only residues without a choice are asked about. A changed
    choice moves atoms, so protonation (4) is cleared and must be redone.
    """
    r = m.state.receptor
    residues = residues or _residues(r.pdb_path)
    relevant = [rid for rid, info in residues.items() if info["altlocs"]
                and ((info["record"] == "ATOM" and info["chain"] in r.chains) or rid in r.keep_hetero)]
    r.altlocs = {k: v for k, v in r.altlocs.items() if k in relevant}
    r.altloc_fills = {k: v for k, v in r.altloc_fills.items() if k in r.altlocs}
    to_choose = relevant if revisit else [rid for rid in relevant if rid not in r.altlocs]
    if not to_choose:
        return
    before = {rid: (r.altlocs.get(rid), r.altloc_fills.get(rid, {})) for rid in to_choose}
    for rid, (letter, plan) in pick_alternates(m, alternates_by_residue(r.pdb_path, to_choose), r.altlocs).items():
        r.altlocs[rid] = letter
        if plan:
            r.altloc_fills[rid] = plan
        else:
            r.altloc_fills.pop(rid, None)
    changed = [rid for rid, (letter, _) in before.items()
               if letter is not None and letter != r.altlocs.get(rid)]
    if changed and r.ph is not None:
        r.ph, r.protonation, r.protonation_sources, r.termini, r.termini_sources = None, {}, {}, {}, {}
        m.console.print(f"  [dark_orange3]Changed alternates ({', '.join(changed)}) move atoms, so the "
                        "protonation choices were cleared; redo 4.[/dark_orange3]", highlight=False)


# --------------------------------------------------------------- 2 hetero ---

# A hydrogen-bond donor-acceptor distance: waters this close to a ligand heavy
# atom are listed in menu 3 by default (any other water can be named by id).
WATER_LISTING_DISTANCE = 3.5


def _ligand_heavy_atoms(m) -> List:
    """The crystal ligand's heavy atoms (with its chosen alternate), when the ligand is a residue of the receptor."""
    from proprep.structure_prep.altloc_picker import keeps_atom
    s = m.state.ligand.source
    rid = m.state.ligand_residue()
    if not rid:
        return []
    return [a for a in read_atoms(m.state.receptor.pdb_path)
            if a.res_id == rid and a.resname == s.get("resname") and a.element not in ("H", "D")
            and keeps_atom(a.name, a.altloc, s.get("altloc"), s.get("altloc_fills"))]


def _hetero_rows(m, residues) -> List[dict]:
    """Menu 2's table: every non-water HETATM residue, then the waters near the ligand."""
    import numpy as np
    r = m.state.receptor
    ligand_res = m.state.ligand_residue()
    atoms = read_atoms(r.pdb_path)
    ligand = _ligand_heavy_atoms(m)
    ligand_xyz = np.array([a.xyz for a in ligand]) if ligand else None
    polar = [a for a in atoms if a.element in ("N", "O") and a.resname not in _WATER and a.res_id != ligand_res]
    polar_xyz = np.array([a.xyz for a in polar]) if polar else None
    by_residue: Dict[str, list] = {}
    for a in atoms:
        if a.record == "HETATM":
            by_residue.setdefault(a.res_id, []).append(a)

    def nearest(rid):
        if ligand_xyz is None or rid == ligand_res:
            return None
        xyz = np.array([a.xyz for a in by_residue[rid] if a.element not in ("H", "D")])
        return float(np.min(np.linalg.norm(xyz[:, None, :] - ligand_xyz[None, :, :], axis=2)))

    rows = []
    for rid, info in residues.items():
        if info["record"] != "HETATM" or info["resname"] in _WATER:
            continue
        kind = "metal" if info["metals"] and info["atoms"] == len(info["metals"]) else (
            "cofactor with metal" if info["metals"] else "other")
        rows.append({"rid": rid, "resname": info["resname"], "atoms": info["atoms"], "kind": kind,
                     "nearest": nearest(rid), "contacts": "", "water": False,
                     "first_atom": by_residue[rid][0].name})
    if ligand_xyz is not None:
        for rid, info in residues.items():
            if info["record"] != "HETATM" or info["resname"] not in _WATER:
                continue
            distance = nearest(rid)
            if distance is None or distance > WATER_LISTING_DISTANCE:
                continue
            oxygen = np.array(by_residue[rid][0].xyz)
            partners = [f"{info_l.resname} {info_l.name}" for info_l in ligand
                        if np.linalg.norm(np.array(info_l.xyz) - oxygen) <= WATER_LISTING_DISTANCE
                        and info_l.element in ("N", "O")]
            if polar_xyz is not None:
                close = np.linalg.norm(polar_xyz - oxygen, axis=1) <= WATER_LISTING_DISTANCE
                partners += [f"{a.res_id} {a.name}" for a, c in zip(polar, close) if c]
            rows.append({"rid": rid, "resname": info["resname"], "atoms": info["atoms"], "kind": "water",
                         "nearest": distance, "contacts": ", ".join(partners), "water": True,
                         "first_atom": by_residue[rid][0].name})
    return rows


def _status(m, rid) -> str:
    if rid == m.state.ligand_residue():
        return "ligand to dock"
    return "keep" if rid in m.state.receptor.keep_hetero else "leave out"


def _show_hetero_in_viewer(m, rows) -> None:
    """Each row of menu 3 in the viewer: its table number as a label, kept residues as balls and
    sticks, residues left out as thin lines, each with its own entry in the representation list."""
    try:
        from proprep.docking_prep.docking_module import _ngl_residue
        from proprep.structure_prep.viewer_coordinator import viewer
        r = m.state.receptor
        listed = " or ".join([f"({_ngl_residue(row['rid'])})" for row in rows] + ["water"])
        m._show_files([r.pdb_path], set_aside=None, aside_selection=listed or None,
                      aside_label="HETATM residues as found (menu 3)")
        entries = []
        for n, row in enumerate(rows, 1):
            status = _status(m, row["rid"])
            selection = _ngl_residue(row["rid"])
            entries.append({"label": f"dock_het_{n}", "selection": selection,
                            "style": "line" if status == "leave out" else "ball+stick", "color": "element",
                            "text": str(n),        # its table number, on the same panel row
                            "display_label": f"#{n} {row['resname']} {row['rid']}: {status}"})
        # Every other water, so any one can be seen and clicked ('p'); kept ones drawn larger.
        listed_waters = {row["rid"] for row in rows if row["water"]}
        waters = [rid for rid, info in _residues(r.pdb_path).items()
                  if info["record"] == "HETATM" and info["resname"] in _WATER and rid not in listed_waters]
        kept = [rid for rid in waters if rid in r.keep_hetero]
        others = [rid for rid in waters if rid not in r.keep_hetero]
        if others:
            entries.append({"label": "dock_het_waters", "selection": " or ".join(f"({_ngl_residue(w)})" for w in others),
                            "style": "ball+stick", "color": "element",
                            "display_label": f"Waters left out ({len(others)}): 'p' picks one"})
        if kept:
            entries.append({"label": "dock_het_waters_kept", "selection": " or ".join(f"({_ngl_residue(w)})" for w in kept),
                            "style": "spacefill", "color": "element", "scale": 0.6,
                            "display_label": f"Waters kept ({len(kept)}), drawn larger"})
        viewer.replace_annotations("dock_het_", entries)
        ligand = m.state.ligand_residue()
        if ligand:
            viewer.focus_on(_ngl_residue(ligand))
    except Exception as error:
        import logging
        logging.getLogger(__name__).debug("viewer: %s", error)


def _parse_switches(answer: str, rows, residues) -> List[str]:
    """Numbers ('1,3-5') and residue ids ('B:378') -> residue ids; ValueError for anything else."""
    ids, numbers = [], []
    for part in ui.parse_list(answer):
        if ":" in part:
            rid = part.upper()
            if rid not in residues or residues[rid]["record"] != "HETATM":
                raise ValueError(f"{part} is not a HETATM residue of this structure")
            ids.append(rid)
        else:
            numbers.append(part)
    if numbers:
        ids += [rows[n - 1]["rid"] for n in ui.parse_numbers(",".join(numbers), len(rows))]
    return ids


def choose_hetero(m) -> None:
    """Which HETATM residues stay in the receptor: the rigid body the ligand is docked against.

    The proposal is stated, not judged: every non-water HETATM residue except
    the ligand, kept as found. ProPrep cannot tell a crystallisation additive
    from a cofactor by itself; the distance to the ligand (when it comes from
    the structure) shows which ones sit in the binding site. Waters are left
    out, and those within a hydrogen-bond distance of the ligand are listed
    with their polar contacts, so a water that bridges ligand and protein can
    be kept; any other water can be named by its residue id.
    """
    r = m.state.receptor
    if not r.pdb_path:
        m.console.print("[dark_orange3]Choose the receptor structure first (1).[/dark_orange3]")
        return
    residues = _residues(r.pdb_path)
    ligand_res = m.state.ligand_residue()
    waters = [rid for rid, info in residues.items() if info["record"] == "HETATM" and info["resname"] in _WATER]
    rows = _hetero_rows(m, residues)
    if not r.hetero_reviewed:
        r.keep_hetero = [row["rid"] for row in rows if not row["water"] and row["rid"] != ligand_res]
    m.console.print("\n  Kept residues are part of the rigid receptor the ligand is docked against. Proposed: "
                    "every HETATM residue other than water and the ligand, as found; ProPrep does not judge "
                    "which are crystallisation additives. Waters are left out.", highlight=False)
    if ligand_res:
        m.console.print(f"  Waters within {WATER_LISTING_DISTANCE} A (a hydrogen-bond distance) of the ligand are "
                        "listed with their polar contacts; a water bridging ligand and protein may be worth "
                        "keeping.", highlight=False)
    else:
        m.console.print("  With the ligand chosen first (2), the waters near it are listed here with their "
                        "contacts.", highlight=False)
    while True:
        table_rows = []
        for n, row in enumerate(rows, 1):
            distance = f"{row['nearest']:.1f}" if row["nearest"] is not None else ""
            table_rows.append((n, row["rid"], row["resname"], row["atoms"], row["kind"], distance,
                               row["contacts"], _status(m, row["rid"])))
        columns = ["#", "residue", "name", "atoms", "kind", "to ligand (A)", "polar contacts", "receptor"]
        if not ligand_res:
            columns, table_rows = (columns[:5] + columns[7:], [t[:5] + t[7:] for t in table_rows])
        m.console.print(ui.table("HETATM residues", columns, table_rows))
        if ligand_res is None:
            m.console.print("  [dark_orange3]If the ligand to dock is one of these, leave it out here, or choose it "
                            "first (2), which takes it out of the receptor.[/dark_orange3]", highlight=False)
        kept_waters = [w for w in waters if w in r.keep_hetero]
        m.console.print(f"  Waters: {len(waters)} in the file, {len(kept_waters)} kept"
                        + (f" ({', '.join(kept_waters[:8])}{' ...' if len(kept_waters) > 8 else ''})"
                           if kept_waters else "")
                        + ". To keep one water, type its residue id (e.g. B:312) or 'p' and click it in the "
                        f"viewer; 'w' keeps or leaves out all {len(waters)} at once.", highlight=False)
        _show_hetero_in_viewer(m, rows)
        answer = ui.ask(m.processor, "Numbers or residue ids to switch between keep and leave out "
                        "(e.g. 1,3-5 or B:312), 'p' to click one in the viewer, 'w' to switch every water, "
                        "Enter when done", default="", description="HETATM residues kept in the receptor")
        if not answer:
            break
        if answer.lower() == "p":
            picked = ui.pick(m, "atom", "Click a water or other HETATM residue to keep or leave out", 0)
            if not picked:
                continue
            atom = picked["atom"]
            answer = f"{atom['chain']}:{atom['resno']}{atom.get('inscode', '')}"
            if answer not in residues or residues[answer]["record"] != "HETATM":
                m.console.print(f"  [dark_orange3]{ui.escape(answer)} ({ui.escape(atom['resname'])}) is part of "
                                "the protein, not a HETATM residue.[/dark_orange3]", highlight=False)
                continue
        if answer.lower() == "w":
            if kept_waters:
                r.keep_hetero = [k for k in r.keep_hetero if k not in waters]
            else:
                r.keep_hetero += [w for w in waters if w not in r.keep_hetero]
            continue
        try:
            switches = _parse_switches(answer, rows, residues)
        except ValueError as error:
            m.console.print(f"[red]{ui.escape(str(error))}[/red]")
            continue
        for rid in switches:
            if rid == ligand_res:
                m.console.print(f"[dark_orange3]{rid} is the ligand being docked.[/dark_orange3]")
            elif rid in r.keep_hetero:
                r.keep_hetero.remove(rid)
            else:
                r.keep_hetero.append(rid)
    r.hetero_reviewed = True
    try:
        from proprep.structure_prep.viewer_coordinator import viewer
        viewer.replace_annotations("dock_het_", [])
    except Exception:
        pass
    if any(w in r.keep_hetero for w in waters):
        m.console.print("  Kept waters are built with hydrogens in a fixed orientation (Meeko's water template). "
                        "Vina does not use hydrogen positions; AD4's hydrogen-bond term does.", highlight=False)
    r.metal_charges = {k: v for k, v in r.metal_charges.items() if k.split("|")[0] in r.keep_hetero}
    choose_altlocs(m, residues)
    m.save_state()


# ----------------------------------------------------------- 3 protonation ---

def _template_family(resname: str) -> List[str]:
    from meeko import ResidueChemTemplates
    chem = ResidueChemTemplates.create_from_defaults()
    return [k for k in chem.ambiguous.get(resname, [resname]) if not k.startswith(("N", "C")) or k in ("CYS", "CYX", "CYX-")]


def choose_protonation(m) -> None:
    from . import receptor_decisions as rd
    r = m.state.receptor
    if not r.pdb_path or not r.chains:
        m.console.print("[dark_orange3]Choose the receptor structure and chains first (1).[/dark_orange3]")
        return
    r.ph = ui.ask_float(m.processor, "pH for protonation states", default=7.0, min_value=0.0, max_value=14.0,
                        description="pH for pdb2pqr/PROPKA")
    m.console.print(f"  Running pdb2pqr with PROPKA at pH {r.ph} ...", highlight=False)
    result = rd.propose_protonation(r.pdb_path, r.chains, r.ph, r.altlocs, r.altloc_fills)
    m.console.print(f"  {result.message}", highlight=False)
    proposals = {k: (v.value, v.source) for k, v in result.proposals.items()}
    metal_residues = [rid for rid in r.keep_hetero if _residues(r.pdb_path).get(rid, {}).get("metals")]
    contacts = rd.metal_contacts(r.pdb_path, metal_residues, r.altlocs, r.altloc_fills) if metal_residues else []
    overrides, problems = rd.metal_overrides(contacts)
    for rid, proposal in overrides.items():
        before = proposals.get(rid, (None,))[0]
        source = proposal.source + (f"; pdb2pqr said {before}" if before and before != proposal.value else "")
        proposals[rid] = (proposal.value, source)
    for note in rd.symmetry_mate_contacts(r.pdb_path, metal_residues):
        m.console.print(f"  [dark_orange3]{note}[/dark_orange3]", highlight=False)
    for problem in problems:
        m.console.print(f"  [red]{problem}[/red]", highlight=False)
    residues = _residues(r.pdb_path)
    his = [rid for rid, info in residues.items() if info["resname"] == "HIS" and info["record"] == "ATOM"
           and info["chain"] in r.chains]
    missing_his = [rid for rid in his if rid not in proposals]
    if missing_his:
        m.console.print(f"  [dark_orange3]No proposal for {len(missing_his)} histidines "
                        f"({', '.join(missing_his[:10])}{' ...' if len(missing_his) > 10 else ''}).[/dark_orange3]",
                        highlight=False)
        choice = ui.ask(m.processor, "Tautomer for the histidines without a proposal", choices=["HID", "HIE", "HIP"],
                        description="Histidine tautomer where pdb2pqr gave none")
        for rid in missing_his:
            proposals[rid] = (choice, "chosen for all histidines without a proposal")
    r.protonation = {k: v[0] for k, v in proposals.items()}
    r.protonation_sources = {k: v[1] for k, v in proposals.items()}
    termini = rd.propose_termini(r.pdb_path, r.chains)
    r.termini = {k: v.value for k, v in termini.items()}
    r.termini_sources = {k: v.source for k, v in termini.items()}
    _decide_open_ends(m)
    _review_protonation(m, residues)
    m.save_state()


END_LABELS = {"N-terminus": "NH3+", "C-terminus": "COO-", "chain break": "break"}


def _show_protonation_in_viewer(m, listed, ends) -> None:
    """Menu 3 in the viewer: each listed residue labelled with its template (side chain), each
    chain end with what it is built as (its terminal atoms), one panel row each."""
    try:
        from proprep.docking_prep.docking_module import _ngl_residue
        from proprep.structure_prep.viewer_coordinator import viewer
        r = m.state.receptor
        m._show_files([r.pdb_path], set_aside=m.state.ligand_residue())
        entries = []
        for rid in listed:
            template = r.protonation[rid]
            entries.append({"label": f"dock_prot_{rid}", "selection": f"({_ngl_residue(rid)}) and sidechainAttached",
                            "style": "ball+stick", "color": "element", "text": template,
                            "display_label": f"{rid} {template} ({r.protonation_sources.get(rid, 'chosen')})"})
        for end in ends:
            treated = r.termini.get(end.res_id)
            atoms = ".N" if end.end == "N" else "(.C or .O or .OXT)"
            entries.append({"label": f"dock_prot_end_{end.res_id}",
                            "selection": f"({_ngl_residue(end.res_id)}) and {atoms}",
                            "style": "ball+stick", "color": "element",
                            "text": f"{end.resname}{end.res_id.split(':')[1]} {END_LABELS.get(treated, '?')}",
                            "display_label": f"{end.res_id} {'first' if end.end == 'N' else 'last'}: "
                                             f"{treated or 'not decided'}"})
        viewer.replace_annotations("dock_prot_", entries)
    except Exception as error:
        import logging
        logging.getLogger(__name__).debug("viewer: %s", error)


def _clear_protonation_in_viewer() -> None:
    try:
        from proprep.structure_prep.viewer_coordinator import viewer
        viewer.replace_annotations("dock_prot_", [])
    except Exception:
        pass


def _review_protonation(m, residues) -> None:
    """Every titratable residue's state and every chain end's treatment; any one can be changed."""
    from .receptor_decisions import chain_ends
    r = m.state.receptor
    standard = {"ASP", "GLU", "LYS", "CYS", "ARG", "TYR"}
    show_all = False
    while True:
        unusual = [rid for rid, key in r.protonation.items()
                   if not (key in standard and residues.get(rid, {}).get("resname") == key)]
        listed = list(r.protonation) if show_all else unusual
        rows = [(rid, residues[rid]["resname"] if rid in residues else "?", r.protonation[rid],
                 r.protonation_sources.get(rid, "chosen")) for rid in listed]
        usual = len(r.protonation) - len(unusual)
        m.console.print(f"\n  Of the {len(r.protonation)} titratable residues, {usual} are in their usual charged "
                        f"form at pH {r.ph} (Asp-, Glu-, Lys+, ...) and {len(unusual)} are not, or are His (which "
                        "has no single usual form). "
                        + ("All are listed ('a' again for the short list)." if show_all else
                           f"Only those {len(unusual)} are listed ('a' lists all).")
                        + " Type any residue to change it, listed or not.", highlight=False)
        title = "Protonation: every titratable residue" if show_all else \
            "Protonation: residues not in their usual charged form, and every His"
        m.console.print(ui.table(title, ["residue", "name", "template", "from"], rows))
        ends = [e for e in chain_ends(r.pdb_path, r.chains) if not e.capped]
        m.console.print(ui.table("Chain ends: a real terminus is charged; a chain break is left neutral",
                                 ["residue", "name", "end", "treated as", "from"],
                                 [(e.res_id, e.resname, "first" if e.end == "N" else "last",
                                   f"{r.termini.get(e.res_id, 'not decided')}: "
                                   f"{TERMINUS_MEANING.get(r.termini.get(e.res_id), 'asked before a run')}",
                                   r.termini_sources.get(e.res_id, "")) for e in ends]))
        _show_protonation_in_viewer(m, listed, ends)
        answer = ui.ask(m.processor, "Residue or chain end to change (e.g. A:94), 'p' to click one in the viewer, "
                        "'a' to list every titratable residue, or Enter when done", default="",
                        description="Change one residue's protonation or terminus")
        if not answer:
            _clear_protonation_in_viewer()
            return
        if answer.lower() == "p":
            picked = ui.pick(m, "atom", "Click a residue or chain end to change", 0)
            if not picked:
                continue
            atom = picked["atom"]
            answer = f"{atom['chain']}:{atom['resno']}{atom.get('inscode', '')}"
        if answer.lower() == "a":
            show_all = not show_all
            continue
        rid = answer.upper()
        if rid not in residues:
            m.console.print(f"[red]No residue {ui.escape(answer)}.[/red]")
            continue
        end = next((e for e in ends if e.res_id == rid), None)
        if end is not None:
            _ask_end(m, end, required=False)
            if residues[rid]["resname"] != "HIS" and len(_template_family(residues[rid]["resname"])) < 2:
                continue
        resname = residues[rid]["resname"]
        family = _template_family(resname)
        if len(family) < 2:
            m.console.print(f"[dark_orange3]{resname} has a single form.[/dark_orange3]")
            continue
        choice = ui.ask(m.processor, f"Template for {rid} {resname}", choices=family,
                        default=r.protonation.get(rid, family[0]), description="Protonation of one residue")
        r.protonation[rid] = choice
        r.protonation_sources[rid] = "chosen"


# ------------------------------------------------------- 4 metals/cofactors ---

def choose_metals_and_cofactors(m) -> None:
    from .receptor_prep import vina_accepts_type
    r = m.state.receptor
    if not r.pdb_path:
        m.console.print("[dark_orange3]Choose the receptor structure first (1).[/dark_orange3]")
        return
    residues = _residues(r.pdb_path)
    metals = [(rid, atom) for rid in r.keep_hetero for atom in residues.get(rid, {}).get("metals", [])]
    if not metals:
        m.console.print("\n  No metals are kept in the receptor; no charges to set.")
    if metals:
        m.console.print("\n  Every kept metal needs a formal charge (its oxidation state); none is assumed.")
    for rid, atom in metals:
        key = metal_key(rid, atom)
        label = f"{atom} of {residues[rid]['resname']} {rid}"
        current = r.metal_charges.get(key)
        if current is not None:
            m.console.print(f"  {label}: now {current:+d}", highlight=False)
        r.metal_charges[key] = ui.ask_required_int(m.processor, m.console, f"Formal charge of {label}",
                                                   description="Metal formal charge")
        element = next(a.element for a in read_atoms(r.pdb_path) if a.res_id == rid and a.name == atom)
        if element not in ("Mg", "Ca", "Mn", "Fe", "Zn") and not vina_accepts_type(element):
            m.console.print(f"  [dark_orange3]Vina has no {element} atom type.[/dark_orange3]")
            answer = ui.ask(m.processor, f"AutoDock type for {label}, or 'omit' to leave it out",
                            description="AutoDock type for a metal Vina cannot type")
            if answer.lower() == "omit" or vina_accepts_type(answer):
                r.metal_types[key] = "omit" if answer.lower() == "omit" else answer
            else:
                m.console.print(f"  [red]Vina does not accept {answer!r} either; nothing recorded.[/red]")
    defaults = _meeko_default_keys()
    cofactors = sorted({residues[rid]["resname"] for rid in r.keep_hetero
                        if rid in residues and residues[rid]["resname"] not in defaults
                        and residues[rid]["atoms"] > len(residues[rid]["metals"])})
    if not cofactors:
        m.console.print("  No cofactors are kept (only standard residues and ions); no chemistry to review.")
    if cofactors:
        m.console.print("  Kept HETATM residues other than metals are cofactors of the receptor: their bond orders, "
                        "formal charges and hydrogens come from their CCD entries, reviewed here once per residue "
                        "name.", highlight=False)
    for name in cofactors:
        kept = [rid for rid in r.keep_hetero if residues.get(rid, {}).get("resname") == name]
        shown_on = kept[0]
        m.console.print(f"\n  {name} is kept in the receptor ({', '.join(kept)}; option 3). If it is the ligand "
                        "you mean to dock, choose it in option 2 instead, which takes it out of the receptor.",
                        highlight=False)
        review_ccd_chemistry(m, name, r.cofactor_edits, shown_on, r.altlocs.get(shown_on),
                             r.altloc_fills.get(shown_on))
    m.save_state()


def review_ccd_chemistry(m, code: str, edit_store: Dict[str, List[dict]], residue_id: Optional[str] = None,
                         altloc: Optional[str] = None, altloc_fills: Optional[dict] = None) -> None:
    """Show a CCD entry's problems and charged atoms; collect edits until it is valid and accepted.

    With ``residue_id`` (the residue of the structure this entry describes),
    the entry is shown in the viewer on that residue after every edit, and
    'p' picks an atom or bond of it there instead of typing names.
    """
    from types import SimpleNamespace
    from .chemistry_edits import apply_edits, ChemistryEdit
    from .docking_menus_ligand import _picked_edit
    component = fetch_component(code)
    edits = [ChemistryEdit.from_dict(e) for e in edit_store.get(code, [])]
    while True:
        mol, problems = apply_edits(component.mol, edits)
        charged = [(a.GetProp("name"), f"{a.GetFormalCharge():+d}") for a in mol.GetAtoms() if a.GetFormalCharge()]
        from rdkit import Chem
        total = sum(a.GetFormalCharge() for a in mol.GetAtoms())
        m.console.print(f"\n  [bold]{code}[/bold] {component.name}: {mol.GetNumAtoms()} atoms, net formal charge "
                        f"{total:+d}; charged atoms: {', '.join(f'{n} {q}' for n, q in charged) or 'none'}",
                        highlight=False)
        if edits:
            m.console.print("  Edits so far: " + "; ".join(e.describe() for e in edits), highlight=False)
        if problems:
            m.console.print(f"  [red]The CCD entry for {code} is inconsistent:[/red]")
            ui.print_notes(m.console, [p.describe() for p in problems], style="red")
        else:
            m.console.print("  The entry is consistent. Press Enter to use it as written; edit it only for another "
                            "protonation or redox state (a reduced flavin, a deprotonated acid, ...).",
                            highlight=False)
        shown = m.show_ccd_in_viewer(mol, code, residue_id, altloc, altloc_fills)
        text = (f"Edit for {code} ({ui.EDIT_HELP}; 'p' to pick in the viewer; {ui.DROP_HELP}; Enter when done)")
        answer = ui.ask(m.processor, text, default="", description=f"Chemistry edit for {code}")
        if answer.lower() == "p":
            if not shown:
                m.console.print(f"  [dark_orange3]{code} is not in the viewer (no open viewer, or no CCD "
                                "coordinates to place on the structure). Type the edit instead.[/dark_orange3]",
                                highlight=False)
                continue
            answer = _picked_edit(m, SimpleNamespace(mol=mol))
            if not answer:
                continue
        if not answer:
            if problems:
                m.console.print(f"  [red]{code} still has problems; it cannot be used until they are resolved.[/red]")
                if not ui.confirm(m.processor, f"Leave {code} unresolved for now?", default=False,
                                  description="Leave a CCD entry unresolved"):
                    continue
            edit_store[code] = [e.to_dict() for e in edits]
            return
        if answer.lower() == ui.DROP:
            if edits:
                m.console.print(f"  Dropped: {ui.escape(edits.pop().describe())}", highlight=False)
            else:
                m.console.print("  No edit to drop.")
            continue
        try:
            edit = ui.parse_edit(answer)
            apply_edits(component.mol, edits + [edit])
            edits.append(edit)
        except (ValueError, KeyError) as error:
            m.console.print(f"  [red]{error}[/red]")


# --------------------------------------------------------------- 5 flexible ---

def choose_flexible(m) -> None:
    r = m.state.receptor
    if not r.pdb_path:
        m.console.print("[dark_orange3]Choose the receptor structure first (1).[/dark_orange3]")
        return
    residues = _residues(r.pdb_path)
    near = _residues_near_box(m, residues)
    if near:
        m.console.print("  Residues with an atom inside the search box: " + ", ".join(near[:40])
                        + (" ..." if len(near) > 40 else ""), highlight=False)
    m.console.print(f"  Flexible now: {', '.join(r.flexible) or 'none'}", highlight=False)
    answer = ui.ask(m.processor, "Residues whose side chains move during docking (e.g. A:45,A:88; "
                    "'none' for a rigid receptor; 'keep' to leave as is; 'p' to pick in the viewer)", default="keep",
                    description="Flexible receptor side chains")
    if answer.lower() == "keep":
        return
    if answer.lower() == "p":
        picked = []
        while True:
            result = ui.pick(m, "atom", "Click an atom of a residue whose side chain should move", 0)
            if not result:
                break
            a = result["atom"]
            rid = f"{a['chain']}:{a['resno']}{a.get('inscode', '')}"
            if rid not in picked:
                picked.append(rid)
            m.console.print(f"  Picked {a['resname']} {rid}.", highlight=False)
            if not ui.confirm(m.processor, "Pick another residue?", default=False, description="Pick another flexible residue"):
                break
        if not picked:
            return
        answer = ",".join(sorted(set(r.flexible) | set(picked)))
    chosen = [] if answer.lower() == "none" else ui.parse_list(answer)
    unknown = [rid for rid in chosen if rid not in residues or residues[rid]["record"] != "ATOM"]
    if unknown:
        m.console.print(f"[red]Not protein residues of this structure: {', '.join(unknown)}.[/red]")
        return
    r.flexible = chosen
    m.save_state()


def _residues_near_box(m, residues) -> List[str]:
    box = m.state.box
    if not box:
        return []
    import numpy as np
    low = np.subtract(box["center"], np.divide(box["size"], 2))
    high = np.add(box["center"], np.divide(box["size"], 2))
    inside = []
    for atom in read_atoms(m.state.receptor.pdb_path):
        if atom.record == "ATOM" and atom.chain in m.state.receptor.chains and np.all(atom.xyz >= low) \
                and np.all(atom.xyz <= high) and atom.res_id not in inside:
            inside.append(atom.res_id)
    return inside
