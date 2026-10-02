"""
Molecular Docking: dock a ligand into the loaded structure with AutoDock Vina
(Vina, Vinardo or AutoDock4 scoring), preparing both with Meeko.

The menu walks the decisions in order and shows each one on the dashboard:

    1 receptor structure and chains     6 flexible side chains
    2 ligand                            7 rotatable bonds
    3 cofactors, metals, waters kept    8 search box
    4 protonation and termini           9 docking settings
    5 metals and cofactor chemistry     r run, c campaign (every ligand of a library), v results

The ligand comes second: whether a HETATM residue is the ligand or part of
the receptor, which waters sit near it and which kept residues are cofactors
all depend on it.

No force-field parameters are needed. Everything decided is kept in the
workspace (key ``docking_state``); each run writes its own folder under
``<working directory>/docking`` (see ``docking_pipeline``).
"""

from __future__ import annotations

import logging
import os
from typing import Dict, List, Optional

from rich.console import Console
from rich.panel import Panel

from proprep.utils.module_registry import ProcessingModule, register_module

from . import docking_ui as ui
from .docking_state import DockingState

logger = logging.getLogger(__name__)
MODULE_NAME = ui.MODULE_NAME
STATE_KEY = "docking_state"

MAIN_OPTIONS = {
    "1": "Receptor structure and chains",
    "2": "Ligand",
    "3": "Cofactors, metals and waters kept",
    "4": "Protonation and termini",
    "5": "Metals and cofactor chemistry",
    "6": "Flexible side chains",
    "7": "Rotatable bonds",
    "8": "Search box",
    "9": "Docking settings",
    "r": "Run docking",
    "c": "Docking campaign (run docking for every ligand of a library)",
    "v": "View results (of the run or the campaign)",
    "b": "Back",
}


def _ngl_residue(residue_id: str) -> str:
    """NGL selection for a residue id such as A:40 or A:40B (insertion code after ^)."""
    chain, number = residue_id.split(":", 1)
    if number[-1].isalpha():
        number = f"{number[:-1]}^{number[-1]}"
    return f":{chain} and {number}"


def _fit(source, target, point):
    """Where ``point`` goes under the rotation and translation that best superpose ``source`` on ``target`` (Kabsch)."""
    import numpy as np
    sc, tc = source.mean(axis=0), target.mean(axis=0)
    u, _, vt = np.linalg.svd((source - sc).T @ (target - tc))
    d = np.sign(np.linalg.det(u @ vt))
    return (point - sc) @ (u @ np.diag([1.0, 1.0, d]) @ vt) + tc


@register_module
class MolecularDockingModule(ProcessingModule):
    NAME = MODULE_NAME
    DESCRIPTION = "Dock a ligand into the structure with AutoDock Vina (Vina, Vinardo or AD4 scoring)"
    VERSION = "1.0.0"
    CATEGORY = "analysis"
    PRIORITY = 60

    def __init__(self):
        super().__init__()
        self.state = DockingState()

    @property
    def console(self) -> Console:
        if self.processor and hasattr(self.processor, "console"):
            return self.processor.console
        return Console()

    # ── framework ────────────────────────────────────────────────────────
    def get_workspace_requirements(self) -> list:
        return []

    def get_workspace_outputs(self) -> list:
        return [STATE_KEY, "docking_output_dir", "docking_poses_sdf", "docking_results"]

    def can_process(self, workspace) -> bool:
        from proprep.utils.structure_selector import StructureSelector
        try:
            return bool(StructureSelector(workspace, self.console).get_structure(silent=True))
        except Exception:
            return False

    def availability_note(self, workspace):
        # the pattern module_registry recommends when can_process is custom: the note cannot
        # disagree with the menu indicator
        return None if self.can_process(workspace) else "Needs a loaded structure"

    def get_menu_options(self) -> Dict[str, str]:
        return {"dock": "Dock a ligand into the structure"}

    def get_enhanced_menu_options(self, workspace):
        from proprep.utils.enhanced_menu import MenuOption, OptionStatus
        done = workspace.get("docking_poses_sdf") is not None
        ready = self.can_process(workspace)
        return [MenuOption(key="1", description="Dock a ligand into the structure",
                           status=OptionStatus.COMPLETED if done else (OptionStatus.AVAILABLE if ready else OptionStatus.BLOCKED),
                           dependency_text=self.availability_note(workspace) or "")]

    def handle_menu_option(self, option: str) -> bool:
        if option == "dock":
            return self.process(self.processor._get_workspace())
        return False

    # ── state ────────────────────────────────────────────────────────────
    def save_state(self) -> None:
        self.update_workspace(self.processor._get_workspace(), STATE_KEY, self.state.to_dict())

    def base_folder(self) -> str:
        workspace = self.processor._get_workspace()
        return os.path.join(workspace.get("working_directory", os.getcwd()), "docking")

    def build_ligand(self, apply_edits_: bool = True):
        from .docking_pipeline import build_ligand
        if apply_edits_:
            return build_ligand(self.state)
        bare = DockingState.from_dict(self.state.to_dict())
        bare.ligand.edits = []
        return build_ligand(bare)

    def _hetero_summary(self, shown: int = 4) -> str:
        """Dashboard line 3: what menu 3 has kept in the receptor, and what it left out."""
        s, r = self.state, self.state.receptor
        if not r.pdb_path:
            return "not set"
        from .docking_menus_receptor import _WATER, _residues
        residues = _residues(r.pdb_path)
        ligand = s.ligand_residue()

        def listing(ids):
            text = ", ".join(f"{residues[i]['resname']} {i}" for i in ids[:shown])
            return text + (f" and {len(ids) - shown} more" if len(ids) > shown else "")

        others = [i for i, info in residues.items() if info["record"] == "HETATM"
                  and info["resname"] not in _WATER and i != ligand]
        waters = [i for i, info in residues.items() if info["record"] == "HETATM" and info["resname"] in _WATER]
        tail = f"; {residues[ligand]['resname']} {ligand} is the ligand (2)" if ligand in residues else ""
        if not others:
            text = "no other HETATM residues in the file"
        elif not r.hetero_reviewed:
            text = f"not reviewed; {len(others)} in the file: {listing(others)}"
        else:
            kept = [i for i in others if i in r.keep_hetero]
            left = [i for i in others if i not in r.keep_hetero]
            text = (f"{len(kept)} of {len(others)} kept" + (f": {listing(kept)}" if kept else "")
                    + (f"; left out: {listing(left)}" if left else ""))
        kept_waters = sum(1 for w in waters if w in r.keep_hetero)
        if kept_waters:
            text += f"; {kept_waters} of {len(waters)} waters kept"
        return text + tail

    # ── viewer (best effort, like every module: never interrupts the menus) ──
    def show_ligand_in_viewer(self, build) -> None:
        if not self.state.receptor.pdb_path:
            return
        try:
            from rdkit import Chem
            from proprep.structure_prep.viewer_coordinator import viewer
            folder = self.base_folder()
            os.makedirs(folder, exist_ok=True)
            path = os.path.join(folder, "ligand_preview.sdf")
            Chem.MolToMolFile(build.mol, path)
            self._show_files([self.state.receptor.pdb_path, path], set_aside=self.state.ligand_residue())
            self.show_box_in_viewer()
        except Exception as error:
            logger.debug("viewer: %s", error)

    def show_ccd_in_viewer(self, mol, code: str, residue_id: Optional[str], altloc: Optional[str] = None,
                           altloc_fills: Optional[dict] = None) -> bool:
        """Show a CCD entry being reviewed, superposed on its residue in the receptor structure.

        The entry (every hydrogen, its bond orders and charges, with the edits
        so far) is fitted onto the crystal residue's heavy atoms, matched by
        name, and shown as structure 1 so its atoms and bonds can be picked.
        Without coordinates to fit (no CCD coordinates, or fewer than three
        atoms matched), the crystal residue is highlighted instead. Returns
        whether the entry itself is shown.
        """
        pdb = self.state.receptor.pdb_path
        if not pdb or not residue_id:
            return False
        try:
            from proprep.structure_prep.viewer_coordinator import viewer
            selection = _ngl_residue(residue_id)
            path = self._ccd_overlay(mol, code, residue_id, altloc, altloc_fills)
            if path:
                self._show_files([pdb, path], set_aside=residue_id)
            else:
                self._show_files([pdb])
                viewer.highlight(selection, style="ball+stick", color="element", label="docking_ccd_residue",
                                 display_label=f"{code} {residue_id}")
            viewer.focus_on(selection)
            return path is not None
        except Exception as error:
            logger.debug("viewer: %s", error)
            return False

    def _ccd_overlay(self, mol, code, residue_id, altloc, altloc_fills) -> Optional[str]:
        import numpy as np
        from rdkit import Chem
        from rdkit.Geometry import Point3D
        from proprep.structure_prep.altloc_picker import keeps_atom
        from .receptor_prep import read_atoms
        if not mol.GetNumConformers():
            return None
        crystal = {a.name: np.array(a.xyz) for a in read_atoms(self.state.receptor.pdb_path)
                   if a.res_id == residue_id and a.resname == code and a.element not in ("H", "D")
                   and keeps_atom(a.name, a.altloc, altloc, altloc_fills)}
        positions = mol.GetConformer().GetPositions()
        matched = [(a.GetIdx(), crystal[a.GetProp("name")]) for a in mol.GetAtoms()
                   if a.GetAtomicNum() > 1 and a.HasProp("name") and a.GetProp("name") in crystal]
        if len(matched) < 3:
            return None
        # Matched atoms go to their crystal positions; every other atom (the
        # hydrogens, and any atom the crystal lacks) is placed by fitting the
        # matched atoms within two bonds of it (three, or all, if fewer than
        # three are that close), so it keeps the CCD's local
        # geometry around the crystal conformation. A rigid fit of the whole
        # ideal conformer would not do: on 1SDU, indinavir's ideal and crystal
        # conformations differ by 4 A RMSD.
        placed = dict(matched)
        distances = Chem.GetDistanceMatrix(mol)
        everything = list(placed)
        moved = positions.copy()
        for index, xyz in placed.items():
            moved[index] = xyz
        for atom in mol.GetAtoms():
            index = atom.GetIdx()
            if index in placed:
                continue
            for reach in (2, 3, None):
                near = [i for i in everything if reach is None or distances[index][i] <= reach]
                if len(near) >= 3:
                    break
            moved[index] = _fit(positions[near], np.array([placed[i] for i in near]), positions[index])
        shown = Chem.Mol(mol)
        conformer = shown.GetConformer()
        for index, xyz in enumerate(moved):
            conformer.SetAtomPosition(index, Point3D(*xyz))
        folder = self.base_folder()
        os.makedirs(folder, exist_ok=True)
        path = os.path.join(folder, f"ccd_preview_{code}.sdf")
        with open(path, "w") as handle:
            handle.write(Chem.MolToMolBlock(shown, kekulize=False))
            handle.write("$$$$\n")
        return path

    def _show_files(self, paths, set_aside: Optional[str] = None, aside_selection: Optional[str] = None,
                    aside_label: Optional[str] = None) -> None:
        """Show ``paths`` in the viewer, reloading it only when something shown has changed.

        A preview file rewritten in place (the ligand after an edit) has the
        same path, which the coordinator treats as already shown; it is
        re-read when its contents differ from the last time. An unchanged
        view is left alone rather than reloaded.

        ``set_aside`` is a residue id of the receptor (structure 0) whose
        crystal copy lies under a preview of the same molecule: it gets its
        own "Crystal ligand" representation, hidden at first, since without
        bond orders it covers the preview's. ``aside_selection`` (with
        ``aside_label``) sets aside any NGL selection of structure 0 instead,
        for a menu that draws those residues its own way.
        """
        import hashlib
        from proprep.structure_prep.viewer_coordinator import viewer
        paths = [str(p) for p in paths]
        digests = []
        for path in paths:
            with open(path, "rb") as handle:
                digests.append(hashlib.sha256(handle.read()).hexdigest())
        if aside_selection:
            aside = {0: {"selection": aside_selection, "label": aside_label or "Set aside"}}
        elif set_aside:
            aside = {0: {"selection": _ngl_residue(set_aside), "label": f"Crystal ligand ({set_aside})"}}
        else:
            aside = None
        last = getattr(self, "_shown_digests", None)
        if viewer.current_structures() == paths and last and last[2] == aside:
            if last[:2] != (paths, digests):
                viewer.refresh_structure()
        else:
            viewer.show_structures(paths, set_aside=aside)
        self._shown_digests = (paths, digests, aside)

    def show_box_in_viewer(self) -> None:
        try:
            from proprep.structure_prep.viewer_coordinator import viewer
            box = self.state.box
            viewer.show_box(box.get("center"), box.get("size"), label="docking_box")
        except Exception as error:
            logger.debug("viewer: %s", error)

    # ── main loop ────────────────────────────────────────────────────────
    def process(self, workspace) -> bool:
        from .dependencies import check_dependencies, missing
        report = check_dependencies()
        blocking = missing(report, include_ad4=False)
        if blocking:
            self.console.print(Panel("\n".join(d.describe() for d in report), title="Docking needs",
                                     border_style="dark_orange3"))
            return False
        self.state = DockingState.from_dict(workspace.get(STATE_KEY))
        from . import docking_menus_batch as batch, docking_menus_ligand as lig, docking_menus_receptor as rec
        actions = {"1": rec.choose_structure, "2": lig.choose_ligand, "3": rec.choose_hetero,
                   "4": rec.choose_protonation, "5": rec.choose_metals_and_cofactors, "6": rec.choose_flexible,
                   "7": lig.choose_torsions, "8": lig.choose_box, "9": lig.choose_settings}
        while True:
            self._dashboard()
            choice = ui.ask(self.processor, "\nEnter choice", choices=list(MAIN_OPTIONS), default=self._suggest(),
                            description="Main menu selection", options_map=MAIN_OPTIONS).lower()
            if choice == "b":
                return True
            try:
                if choice in actions:
                    actions[choice](self)
                elif choice == "r":
                    self.run()
                elif choice == "v":
                    self._view_results(batch)
                elif choice == "c":
                    batch.campaign(self)
            except Exception as error:             # a failed step is reported; the menu stays
                logger.exception("docking step %s failed", choice)
                self.console.print(f"[red]{MAIN_OPTIONS[choice]} failed: {error}[/red]")

    def _suggest(self) -> str:
        needs = self.state.missing()
        if not needs:
            return "r"
        return needs[0].split("(")[-1].rstrip(")")

    def _dashboard(self) -> None:
        s, r, lig = self.state, self.state.receptor, self.state.ligand
        unset = "[dark_orange3]not set[/dark_orange3]"
        self.console.print(Panel(
            "[bold]Dock a ligand into the loaded structure with AutoDock Vina.[/bold]\n"
            "Ligand and receptor are prepared with Meeko; no force-field parameters are needed.",
            title=MODULE_NAME, border_style="bright_blue", width=64, padding=(0, 1)))
        source = lig.source
        ligand_text = unset
        if source.get("kind") == "structure_residue":
            ligand_text = f"{source['resname']} {source['chain']}:{source['resseq']}{source.get('icode', '')} (crystal pose)"
        elif source.get("kind") == "smiles":
            ligand_text = f"SMILES {source['smiles']}"
        elif source.get("kind") == "file":
            ligand_text = f"{os.path.basename(str(source['path']))} record {source['record']}"
        rotatable = ("reviewed" if lig.rigid_macrocycles is not None else unset) if source else unset
        box = (f"{' x '.join(f'{x:.1f}' for x in s.box['size'])} A, {s.box['source']}" if s.box else unset)
        settings = (f"{s.settings['scoring']}, seed {s.settings['seed']}, exhaustiveness {s.settings['exhaustiveness']}, "
                    f"{s.settings['n_poses']} poses, max evaluations: Vina's heuristic"
                    if s.settings else unset)
        metals = sum(1 for _ in r.metal_charges)
        B = ui.BLUE
        lines = [
            "  [bold]Current setup[/bold]",
            f"  [{B}]1[/{B}] Receptor:     {os.path.basename(r.pdb_path) if r.pdb_path else unset}"
            + (f", chains {','.join(r.chains)}" if r.chains else ""),
            f"  [{B}]2[/{B}] Ligand:       {ligand_text}",
            f"  [{B}]3[/{B}] Kept HETATM:  {ui.escape(self._hetero_summary())}",
            f"  [{B}]4[/{B}] Protonation:  " + (f"pH {r.ph}, {len(r.protonation)} residue(s) set" if r.ph is not None else unset),
            f"  [{B}]5[/{B}] Metals:       {metals} charge(s) set; cofactor edits for "
            + (", ".join(r.cofactor_edits) or "none"),
            f"  [{B}]6[/{B}] Flexible:     {', '.join(r.flexible) or 'none (rigid receptor)'}",
            f"  [{B}]7[/{B}] Rotatable:    {rotatable}",
            f"  [{B}]8[/{B}] Search box:   {box}",
            f"  [{B}]9[/{B}] Settings:     {settings}",
        ]
        needs = s.missing()
        lines.append("")
        lines.append(f"  [bold #1a7f37]r[/bold #1a7f37] Run docking" + (f"  [dark_orange3](needs {', '.join(needs)})[/dark_orange3]" if needs else ""))
        library = self.processor._get_workspace().get("ligand_file")
        campaign_needs = [n for n in needs if "ligand" not in n and "rotatable" not in n]
        lines.append(f"  [{B}]c[/{B}] Docking campaign: run docking for every ligand of "
                     + (os.path.basename(library) if library else
                        "[dark_orange3]a library (load one with the Structure Loader, its option 4)[/dark_orange3]")
                     + (f"  [dark_orange3](needs {', '.join(campaign_needs)})[/dark_orange3]" if campaign_needs else ""))
        lines.append(f"  [{B}]v[/{B}] View results: " + self._results_summary())
        lines.append(f"  [{B}]b[/{B}] Back")
        self.console.print("\n".join(lines), highlight=False)

    # ── run and results ──────────────────────────────────────────────────
    def _results_folders(self):
        """(single-run folder, campaign folder), each None when there is none."""
        run = self.state.last_run.get("folder") if self.state.last_run else None
        campaign = self.state.last_campaign.get("folder") if self.state.last_campaign else None
        return run, campaign

    @staticmethod
    def _newer(run_folder: str, campaign_folder: str) -> str:
        """'s' or 'c': which folder was written last (each folder name ends in its time stamp)."""
        stamp = lambda folder: os.path.basename(folder.rstrip(os.sep))[-15:]
        return "c" if stamp(campaign_folder) > stamp(run_folder) else "s"

    def _results_summary(self) -> str:
        run, campaign = self._results_folders()
        if not run and not campaign:
            return "[grey50]no run or campaign yet[/grey50]"
        parts = []
        if run:
            parts.append("the run " + os.path.basename(run))
        if campaign:
            parts.append("the campaign " + os.path.basename(campaign))
        return ui.escape(" and ".join(parts))

    def _view_results(self, batch) -> None:
        """The single run's poses or the campaign's ranking; asks which when there are both."""
        run, campaign = self._results_folders()
        if not run and not campaign:
            self.console.print("[dark_orange3]No docking run or campaign yet (r or c).[/dark_orange3]")
            return
        if run and campaign:
            newer = self._newer(run, campaign)
            self.console.print(f"  The run: {os.path.basename(run)}" + ("  (newer)" if newer == "s" else "")
                               + f"\n  The campaign: {os.path.basename(campaign)}" + ("  (newer)" if newer == "c" else ""),
                               highlight=False)
            which = ui.ask(self.processor, "Results of: [s] the single run, [c] the campaign", choices=["s", "c"],
                           default=newer, description="Which results to view",
                           options_map={"s": "The single-ligand run", "c": "The campaign"})
        else:
            which = "c" if campaign else "s"
        if which == "c":
            batch.results(self)
        else:
            self.results()
    def run(self) -> None:
        from .dependencies import check_dependencies, missing
        from .docking_pipeline import run_state
        needs = self.state.missing()
        if needs:
            self.console.print(f"[dark_orange3]Still needed: {', '.join(needs)}.[/dark_orange3]")
            return
        if self.state.settings.get("scoring") == "ad4":
            blocking = missing(check_dependencies(), include_ad4=True)
            if blocking:
                self.console.print("[red]" + "; ".join(d.describe() for d in blocking) + "[/red]")
                return
        self.console.print("  Building the ligand and receptor and docking; this takes seconds to minutes ...")
        outcome = run_state(self.state, self.base_folder())
        self.console.print(f"  Done in {outcome.run.seconds:.1f} s. Output: {outcome.folder}", highlight=False)
        ui.print_notes(self.console, outcome.receptor.notes + outcome.run.notes)
        self._print_poses(outcome.poses, outcome.reference_rmsd)
        workspace = self.processor._get_workspace()
        self.state.last_run = {
            "folder": outcome.folder, "files": outcome.files,
            "energies": [p.energies for p in outcome.poses],
            "rmsd": [p.rmsd_to_reference for p in outcome.poses] if outcome.reference_rmsd else None,
        }
        self.save_state()
        self.update_workspace(workspace, "docking_output_dir", outcome.folder)
        self.update_workspace(workspace, "docking_poses_sdf", outcome.files["poses_sdf"])
        self.update_workspace(workspace, "docking_results", self.state.last_run)
        self._show_pose(1)

    def _print_poses(self, poses, with_rmsd: bool) -> None:
        columns = ["pose", "total (kcal/mol)", "inter", "intra", "torsions"] + (["RMSD to crystal (A)"] if with_rmsd else [])
        rows = []
        for pose in poses:
            e = pose.energies
            row = [pose.rank, f"{e['total']:.2f}", f"{e['inter']:.2f}", f"{e['intra']:.2f}", f"{e['torsions']:.2f}"]
            if with_rmsd:
                row.append(f"{pose.rmsd_to_reference:.2f}")
            rows.append(row)
        self.console.print(ui.table("Docked poses (best first)", columns, rows))

    def results(self) -> None:
        last = self.state.last_run
        if not last:
            self.console.print("[dark_orange3]No docking run yet (r).[/dark_orange3]")
            return
        rmsd = last.get("rmsd")
        rows = [[n, f"{e['total']:.2f}", f"{e['inter']:.2f}"] + ([f"{rmsd[n - 1]:.2f}"] if rmsd else [])
                for n, e in enumerate(last["energies"], 1)]
        self.console.print(f"  Run folder: {last['folder']}", highlight=False)
        self.console.print(ui.table(f"Poses of {os.path.basename(last['folder'])}", ["pose", "total (kcal/mol)", "inter"]
                                    + (["RMSD to crystal (A)"] if rmsd else []), rows))
        number = ui.ask_int(self.processor, "Pose to show in the viewer (0 for none)", default=0, min_value=0,
                            max_value=len(last["energies"]), description="Pose to view")
        if number:
            self._show_pose(number, asked=True)

    def _show_pose(self, number: int, *, asked: bool = False) -> None:
        """Show pose ``number`` with the receptor and the box, and say so.

        After a run (``asked`` False) the pose goes to a viewer that is already
        open; with none open the terminal says how to see it. Asked for from
        the results (``asked`` True), the viewer is opened if it is not.
        """
        files = self.state.last_run.get("files", {})
        try:
            from proprep.structure_prep.viewer_coordinator import viewer, _is_web_shell_mode
            open_already = viewer.is_running() or _is_web_shell_mode()
            if not open_already and not asked:
                self.console.print("  The viewer is not open; 'v' shows any pose in it.", highlight=False)
                return
            shown = [self.state.receptor.pdb_path, files[f"pose_{number}"]]
            if f"flex_{number}" in files:
                shown.append(files[f"flex_{number}"])
            viewer.show_structures(shown, force=not open_already)
            self.show_box_in_viewer()
            self.console.print(f"  Viewer: pose {number} with the receptor and the search box"
                               + (" (opened in your browser)." if not open_already else "."), highlight=False)
        except Exception as error:
            logger.debug("viewer: %s", error)
            self.console.print(f"  [dark_orange3]Could not show pose {number}: {ui.escape(str(error))}[/dark_orange3]",
                               highlight=False)
