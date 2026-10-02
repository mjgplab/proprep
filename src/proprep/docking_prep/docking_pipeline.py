"""
From a DockingState to docked poses on disk: the menus stay thin, and a saved
state is enough to rebuild the ligand, the receptor and the run.

Output folder (``state.output_dir``), one per run, named by time:
    receptor_rigid.pdbqt, receptor_flex.pdbqt   what Vina docked into
    ligand.pdbqt                                what Vina docked
    poses.sdf                                   every pose, best first, energies as properties
    pose_<n>.sdf                                one file per pose, for the viewer
    flexible_residues_pose_<n>.pdb              moved side chains of that pose (flexible runs)
    docking.json                                the state, settings, box, energies and notes
"""

from __future__ import annotations

import json
import os
import time
from dataclasses import dataclass, field
from typing import Dict, List, Optional

from rdkit import Chem

from . import ligand_sources as ls
from .ccd_chemistry import fetch_component
from .chemistry_edits import ChemistryEdit, apply_edits
from .docking_results import Pose, heavy_atom_rmsd, poses_from_run, write_poses_sdf
from .docking_run import DockingRun, DockingSettings, SearchBox, run_docking
from .docking_state import DockingState
from .ligand_prep import PreparedLigand, prepare_ligand
from .receptor_prep import PreparedReceptor, prepare_receptor


def build_ligand(state: DockingState, components: Optional[Dict] = None) -> ls.LigandBuild:
    """The ligand from its recorded source, with the recorded chemistry edits applied."""
    source = dict(state.ligand.source)
    kind = source.get("kind")
    if kind == "smiles":
        build = ls.from_smiles(source["smiles"], int(source["embed_seed"]), source["optimizer"],
                               residue_name=source.get("residue_name", "LIG"))
    elif kind == "file":
        build = ls.from_file(source["path"], record=source.get("record"),
                             residue_name=source.get("residue_name", "LIG"))
    elif kind == "structure_residue":
        code = source["ccd_code"]
        component = (components or {}).get(code) or fetch_component(code)
        component.mol, component.problems = apply_edits(
            component.mol, [ChemistryEdit.from_dict(e) for e in state.ligand.ccd_edits])
        build = ls.from_structure_residue(source["path"], source["chain"], int(source["resseq"]),
                                          source.get("icode", ""), source["resname"], component,
                                          altloc=source.get("altloc"),
                                          altloc_fills=source.get("altloc_fills") or None)
    else:
        raise ValueError("No ligand has been chosen.")
    if state.ligand.edits:
        mol, problems = apply_edits(build.mol, state.ligand.edit_objects())
        if problems:
            raise ValueError("The ligand edits leave invalid chemistry: "
                             + "; ".join(p.describe() for p in problems))
        build.mol = mol
        build.notes.append(f"{len(state.ligand.edits)} chemistry edit(s) applied: "
                           + "; ".join(e.describe() for e in state.ligand.edit_objects()))
    return build


def prepare_state_ligand(state: DockingState, build: ls.LigandBuild) -> PreparedLigand:
    return prepare_ligand(build.mol, rigid=[tuple(b) for b in state.ligand.rigid_bonds],
                          rotatable_amides=[tuple(b) for b in state.ligand.rotatable_amides],
                          rigid_macrocycles=bool(state.ligand.rigid_macrocycles),
                          metal_charges=dict(state.ligand.metal_charges) or None)


def build_receptor(state: DockingState, components: Optional[Dict] = None) -> PreparedReceptor:
    r = state.receptor
    if not r.pdb_path or not r.chains:
        raise ValueError("Choose the receptor structure and its chains first.")
    return prepare_receptor(r.pdb_path, r.chains, r.keep_hetero, components=components, **r.build_arguments())


@dataclass
class DockingOutcome:
    folder: str
    run: DockingRun
    poses: List[Pose]
    receptor: PreparedReceptor
    ligand: PreparedLigand
    build: ls.LigandBuild
    reference_rmsd: bool = False
    files: Dict[str, str] = field(default_factory=dict)


def run_state(state: DockingState, base_folder: str, components: Optional[Dict] = None) -> DockingOutcome:
    """Build everything from ``state``, dock, and write the outputs into a new folder under ``base_folder``."""
    missing = state.missing()
    if missing:
        raise ValueError("Still needed: " + ", ".join(missing) + ".")
    build = build_ligand(state, components)
    ligand = prepare_state_ligand(state, build)
    receptor = build_receptor(state, components)
    box = SearchBox(tuple(state.box["center"]), tuple(state.box["size"]), state.box["source"])
    settings = DockingSettings(**state.settings)
    folder = os.path.join(base_folder, time.strftime("docking_%Y%m%d_%H%M%S"))
    suffix = 1
    while os.path.exists(folder):
        suffix += 1
        folder = os.path.join(base_folder, time.strftime("docking_%Y%m%d_%H%M%S") + f"_{suffix}")
    os.makedirs(folder)
    run = run_docking(receptor.rigid_pdbqt, receptor.flex_pdbqt, ligand.pdbqt, box, settings, workdir=folder)
    poses = poses_from_run(run.poses_pdbqt, build.mol, run.energies)

    # a ligand taken from the structure is its own reference pose
    reference = build.source.get("kind") == "structure_residue"
    if reference:
        for pose in poses:
            pose.rmsd_to_reference = heavy_atom_rmsd(pose.ligand, build.mol)

    files = {"receptor_rigid": os.path.join(folder, "receptor_rigid.pdbqt"),
             "ligand_pdbqt": os.path.join(folder, "ligand.pdbqt")}
    if receptor.flex_pdbqt:
        files["receptor_flex"] = os.path.join(folder, "receptor_flex.pdbqt")
    files["poses_sdf"] = os.path.join(folder, "poses.sdf")
    write_poses_sdf(files["poses_sdf"], poses, {"scoring": settings.scoring, "seed": settings.seed,
                                                "receptor": os.path.basename(state.receptor.pdb_path)})
    for pose in poses:
        path = os.path.join(folder, f"pose_{pose.rank}.sdf")
        write_poses_sdf(path, [pose], {"scoring": settings.scoring, "seed": settings.seed})
        files[f"pose_{pose.rank}"] = path
        if pose.flexible_residues:
            flex_path = os.path.join(folder, f"flexible_residues_pose_{pose.rank}.pdb")
            with open(flex_path, "w") as handle:
                for side_chain in pose.flexible_residues:
                    handle.write(Chem.MolToPDBBlock(side_chain).replace("END\n", ""))
                handle.write("END\n")
            files[f"flex_{pose.rank}"] = flex_path
    record = {
        "state": state.to_dict(),
        "energies": run.energies,
        "rmsd_to_reference": [p.rmsd_to_reference for p in poses] if reference else None,
        "seconds": run.seconds,
        "notes": {"ligand": build.notes + ligand.notes, "receptor": receptor.notes, "run": run.notes},
        "receptor_residues": [vars(d) for d in receptor.residues],
        "metals": receptor.metals,
        "links": receptor.links,
        "removed": receptor.removed,
    }
    files["record"] = os.path.join(folder, "docking.json")
    with open(files["record"], "w") as handle:
        json.dump(record, handle, indent=1, default=str)
    return DockingOutcome(folder=folder, run=run, poses=poses, receptor=receptor, ligand=ligand,
                          build=build, reference_rmsd=reference, files=files)
