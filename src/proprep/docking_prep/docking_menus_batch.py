"""
Molecular Docking menu c: a docking campaign (one receptor, one box, a ligand library).

The receptor (1, 3-6), the box (8) and the settings (9) are the single-ligand ones;
the library is the file loaded with the Structure Loader (option 4: SDF, mol2 or
SMILES). The questions the single-ligand menus ask about each ligand are asked
once, as the campaign's policies (see ``batch``), and recorded in its manifest.
"""

from __future__ import annotations

import os

from rich.panel import Panel

from . import docking_ui as ui

STEREO = {"s": "Skip them", "d": "Dock them as built (the configuration RDKit chose is reported)"}
FLAT = {"b": "Build 3D with RDKit", "s": "Skip them"}


def campaign(m) -> None:
    from .batch import BatchPolicies, find_resumable, ranked, read_library, run_campaign, _manifest
    from .docking_pipeline import build_receptor
    state = m.state
    needs = [n for n in state.missing() if "ligand" not in n and "rotatable" not in n]
    if needs:
        m.console.print(f"[dark_orange3]A campaign needs {', '.join(needs)} first.[/dark_orange3]")
        return
    library = m.processor._get_workspace().get("ligand_file")
    if not library or not os.path.exists(library):
        m.console.print("[dark_orange3]No ligand library is loaded. Load an SDF, mol2 or SMILES file with the "
                        "Structure Loader (option 4), then come back.[/dark_orange3]")
        return
    entries = read_library(library)
    unreadable = sum(1 for e in entries if e.error)
    m.console.print(f"\n  Library: {library}\n  {len(entries)} molecules"
                    + (f", {unreadable} unreadable (they will be listed as skipped)" if unreadable else "") + ".",
                    highlight=False)
    m.console.print("  Protonation and formal charges are taken as written in the library: ProPrep does not "
                    "change a ligand's protonation for pH. Each ligand's net charge is reported.", highlight=False)
    stereo = ui.ask(m.processor, "Ligands whose stereocentres are not specified: [s] skip, [d] dock as built",
                    choices=list(STEREO), default="s", description="Campaign: unspecified stereo",
                    options_map=STEREO)
    flat = ui.ask(m.processor, "Records with 2D coordinates: [b] build 3D with RDKit, [s] skip",
                  choices=list(FLAT), default="b", description="Campaign: 2D records", options_map=FLAT)
    seed = ui.ask_int(m.processor, "Seed for 3D builds (RDKit ETKDG)", default=42, min_value=1,
                      description="Campaign: ETKDG seed")
    optimizer = ui.ask(m.processor, "Minimise 3D builds with MMFF94?", choices=["mmff94", "none"],
                       default="mmff94", description="Campaign: 3D build minimisation")
    rigid_conjugated = ui.confirm(m.processor, "Hold every conjugated single bond rigid (ester C(=O)-O, polyene "
                                  "single bonds)?", default=False, description="Campaign: conjugated bonds rigid")
    open_macrocycles = ui.confirm(m.processor, "Let Meeko open macrocycles so the ring flexes?", default=True,
                                  description="Campaign: macrocycle flexibility")
    policies = BatchPolicies("dock" if stereo == "d" else "skip", "build" if flat == "b" else "skip", seed,
                             optimizer, rigid_conjugated, not open_macrocycles)

    m.console.print("  Building the receptor ...")
    receptor = build_receptor(state)
    resume = find_resumable(m.base_folder(), _manifest(state, library, policies, receptor))
    if resume:
        done = _count_rows(resume)
        m.console.print(f"  An unfinished campaign made with exactly these choices is in {resume} "
                        f"({done} of {len(entries)} done).", highlight=False)
        if not ui.confirm(m.processor, "Resume it?", default=True, description="Campaign: resume"):
            resume = None
    todo = len(entries) - (_count_rows(resume) if resume else 0)
    m.console.print(f"  {todo} ligand(s) to dock.", highlight=False)
    if not ui.confirm(m.processor, "Start docking? (Ctrl-C stops; finished ligands are kept)",
                      default=True, description="Campaign: start"):
        return
    result = run_campaign(state, library, policies, m.base_folder(), resume_folder=resume,
                          progress=lambda text: m.console.print(f"  {ui.escape(text)}", highlight=False))
    state.last_campaign = {"folder": result.folder, "library": library}
    m.save_state()
    m.update_workspace(m.processor._get_workspace(), "docking_campaign_dir", result.folder)
    summarise(m, result.folder, result.rows)


def _count_rows(folder) -> int:
    import csv
    path = os.path.join(folder, "results.csv")
    if not os.path.exists(path):
        return 0
    with open(path) as handle:
        return sum(1 for _ in csv.DictReader(handle))


def _rows(folder):
    import csv
    with open(os.path.join(folder, "results.csv")) as handle:
        return list(csv.DictReader(handle))


def summarise(m, folder, rows, top: int = 10) -> None:
    from .batch import ranked
    counts = {status: sum(1 for r in rows if r["status"] == status) for status in ("docked", "skipped", "failed")}
    m.console.print(Panel(f"{counts['docked']} docked, {counts['skipped']} skipped, {counts['failed']} failed.\n"
                          f"Folder: {folder}\nresults.csv: every ligand, with the reason for any skip or failure\n"
                          "top_poses.sdf: each docked ligand's best pose, best score first",
                          title="Campaign", border_style="bright_blue"), highlight=False)
    best = ranked(rows)[:top]
    if best:
        m.console.print(ui.table(f"Top {len(best)} by score", ["rank", "ligand", "score (kcal/mol)",
                                                              "per heavy atom", "heavy atoms", "net charge"],
                                 [(n, r["name"], r["best_score"], r["score_per_heavy_atom"], r["heavy_atoms"],
                                   r["net_charge"]) for n, r in enumerate(best, 1)]))


def results(m) -> None:
    """Rank a campaign's ligands and show any one's best pose in the viewer."""
    from .batch import ranked, _safe
    folder = m.state.last_campaign.get("folder")
    if not folder or not os.path.exists(os.path.join(folder, "results.csv")):
        m.console.print("[dark_orange3]No campaign yet (c).[/dark_orange3]")
        return
    rows = _rows(folder)
    summarise(m, folder, rows, top=20)
    order = ranked(rows)
    if not order:
        return
    rank = ui.ask_int(m.processor, "Rank of the ligand to show in the viewer (0 for none)", default=0, min_value=0,
                      max_value=len(order), description="Campaign ligand to view")
    if not rank:
        return
    row = order[rank - 1]
    poses = os.path.join(folder, "ligands", f"{int(row['index']):05d}_{_safe(row['name'])}", "poses.sdf")
    try:
        from rdkit import Chem
        from proprep.structure_prep.viewer_coordinator import viewer
        best = next(iter(Chem.SDMolSupplier(poses, removeHs=False)))
        shown = os.path.join(os.path.dirname(poses), "best_pose.sdf")
        Chem.MolToMolFile(best, shown)
        viewer.show_structures([m.state.receptor.pdb_path, shown])
        m.show_box_in_viewer()
        m.console.print(f"  Showing {row['name']} (rank {rank}, {row['best_score']} kcal/mol).", highlight=False)
    except Exception as error:
        m.console.print(f"  [dark_orange3]Could not show it: {error}[/dark_orange3]")
