"""
Docking campaigns: one receptor, one box, one set of settings, a library of ligands.

A campaign cannot stop to ask about each ligand, so the questions the
single-ligand menus ask are answered once, as policies, and recorded with the
results:

* protonation and formal charges: as written in the library (RDKit has no
  pKa-based ligand protonation; prepare the library in the state intended);
  every ligand's net charge is reported;
* SMILES with unspecified stereocentres, and 2D records: dock as built or skip;
* 3D builds (SMILES, 2D records): RDKit ETKDG with one seed, MMFF94 or not;
* torsions: Meeko's defaults, optionally every conjugated single bond rigid,
  macrocycles opened or kept;
* a metal in a ligand: its formal charge as written.

A ligand that cannot be read or prepared is recorded with the reason and the
campaign goes on. Each ligand is docked by its own Vina instance with the
campaign's seed, so its result does not depend on where it sits in the
library, and results are written as each ligand finishes: an interrupted
campaign resumes where it stopped when rerun with the same receptor, box,
settings, policies and library (checked against the manifest's checksums).

Campaign folder:
    manifest.json        every decision, the library and receptor checksums
    receptor_rigid.pdbqt, receptor_flex.pdbqt
    results.csv          one row per ligand: status, scores, charge, notes or reason
    top_poses.sdf        each docked ligand's best pose, best score first
    ligands/<n>_<name>/  ligand.pdbqt, poses.sdf, Vina's output
    ad4_maps/            AutoGrid maps built once (AD4 scoring)
"""

from __future__ import annotations

import csv
import hashlib
import json
import os
import re
import time
from dataclasses import asdict, dataclass, field
from typing import Callable, Dict, List, Optional

from rdkit import Chem
from rdkit.Chem import AllChem

LIBRARY_SUFFIXES = (".sdf", ".sd", ".mol", ".mol2", ".smi", ".smiles")
RESULT_COLUMNS = ["index", "name", "status", "best_score", "best_inter", "score_per_heavy_atom", "poses",
                  "heavy_atoms", "rotatable_bonds", "net_charge", "seconds", "notes"]


@dataclass(frozen=True)
class BatchPolicies:
    unspecified_stereo: str      # "dock" (as built, reported) or "skip"
    flat_records: str            # 2D records: "build" (RDKit 3D) or "skip"
    embed_seed: int
    optimizer: str               # "mmff94" or "none", for 3D builds
    rigid_conjugated: bool
    rigid_macrocycles: bool

    def validate(self) -> None:
        if self.unspecified_stereo not in ("dock", "skip") or self.flat_records not in ("build", "skip"):
            raise ValueError("unspecified_stereo must be dock/skip and flat_records build/skip.")
        if self.optimizer not in ("mmff94", "none") or self.embed_seed < 1:
            raise ValueError("optimizer must be mmff94/none and the seed at least 1.")

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass
class LibraryEntry:
    index: int                      # 1-based position in the file
    name: str
    mol: Optional[Chem.Mol] = None  # SDF/mol2 records
    smiles: Optional[str] = None    # SMILES files
    error: Optional[str] = None


class LigandSkipped(Exception):
    """A ligand the campaign's policies or its chemistry keep from being docked; the message says why."""


# ------------------------------------------------------------------ library ---

def _unique(names: List[str]) -> List[str]:
    seen, out = {}, []
    for index, name in enumerate(names, 1):
        clean = re.sub(r"\s+", "_", name.strip()) or f"mol_{index}"
        if clean in seen:
            clean = f"{clean}_{index}"
        seen[clean] = True
        out.append(clean)
    return out


def read_library(path: str) -> List[LibraryEntry]:
    """Every molecule of an SDF, mol2 or SMILES file, in file order; unreadable ones carry an error
    (recorded in the results; RDKit's own log lines are silenced so they do not repeat it)."""
    from rdkit import RDLogger
    RDLogger.DisableLog("rdApp.error")
    try:
        return _read_library(path)
    finally:
        RDLogger.EnableLog("rdApp.error")


def _read_library(path: str) -> List[LibraryEntry]:
    suffix = os.path.splitext(path)[1].lower()
    entries: List[LibraryEntry] = []
    if suffix in (".sdf", ".sd", ".mol"):
        supplier = Chem.SDMolSupplier(path, removeHs=False, sanitize=True)
        mols = [supplier[i] for i in range(len(supplier))]
        names = _unique([m.GetProp("_Name") if m is not None and m.HasProp("_Name") else "" for m in mols])
        for index, (mol, name) in enumerate(zip(mols, names), 1):
            entries.append(LibraryEntry(index, name, mol=mol,
                                        error=None if mol is not None else "RDKit could not read the record"))
    elif suffix == ".mol2":
        with open(path) as handle:
            blocks = ["@<TRIPOS>MOLECULE" + b for b in handle.read().split("@<TRIPOS>MOLECULE")[1:]]
        mols = [Chem.MolFromMol2Block(b, removeHs=False, sanitize=True) for b in blocks]
        names = _unique([b.splitlines()[1].strip() if len(b.splitlines()) > 1 else "" for b in blocks])
        for index, (mol, name) in enumerate(zip(mols, names), 1):
            entries.append(LibraryEntry(index, name, mol=mol,
                                        error=None if mol is not None else "RDKit could not read the record"))
    elif suffix in (".smi", ".smiles"):
        rows = []
        with open(path) as handle:
            for line in handle:
                text = line.strip()
                if text and not text.startswith("#"):
                    parts = text.split(None, 1)
                    rows.append((parts[0], parts[1] if len(parts) > 1 else ""))
        for index, ((smiles, _), name) in enumerate(zip(rows, _unique([r[1] for r in rows])), 1):
            error = None if Chem.MolFromSmiles(smiles) is not None else f"RDKit could not read the SMILES {smiles!r}"
            entries.append(LibraryEntry(index, name, smiles=smiles, error=error))
    else:
        raise ValueError(f"Unsupported library file {suffix!r}; use {', '.join(LIBRARY_SUFFIXES)}.")
    return entries


# ---------------------------------------------------------------- preparing ---

def _is_3d(mol: Chem.Mol) -> bool:
    return mol.GetNumConformers() > 0 and mol.GetConformer().Is3D() and \
        any(abs(p[2]) > 1e-4 for p in mol.GetConformer().GetPositions())


def _unspecified_stereo(mol: Chem.Mol) -> int:
    return sum(1 for info in Chem.FindPotentialStereo(mol) if info.specified == Chem.StereoSpecified.Unspecified)


def build_entry(entry: LibraryEntry, policies: BatchPolicies):
    """The ligand of one library entry, as a LigandBuild; LigandSkipped when a policy or its chemistry says no."""
    from . import ligand_sources as ls
    if entry.error:
        raise LigandSkipped(entry.error)
    if entry.smiles is not None:
        mol = Chem.MolFromSmiles(entry.smiles)
        if _unspecified_stereo(mol) and policies.unspecified_stereo == "skip":
            raise LigandSkipped(f"{_unspecified_stereo(mol)} unspecified stereo element(s) (policy: skip)")
        try:
            return ls.from_smiles(entry.smiles, policies.embed_seed, policies.optimizer)
        except (ValueError, RuntimeError) as error:
            raise LigandSkipped(str(error)) from None
    mol = Chem.Mol(entry.mol)
    notes = []
    if not _is_3d(mol):
        if policies.flat_records == "skip":
            raise LigandSkipped("2D record (policy: skip)")
        Chem.AssignChiralTypesFromBondDirs(mol)
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
        if _unspecified_stereo(mol) and policies.unspecified_stereo == "skip":
            raise LigandSkipped(f"{_unspecified_stereo(mol)} unspecified stereo element(s) (policy: skip)")
        mol = Chem.AddHs(Chem.RemoveHs(mol))
        params = AllChem.ETKDGv3()
        params.randomSeed = policies.embed_seed
        if AllChem.EmbedMolecule(mol, params) != 0:
            raise LigandSkipped("RDKit could not build 3D coordinates from the 2D record")
        notes.append(f"2D record; 3D built with RDKit ETKDGv3, seed {policies.embed_seed}.")
        if policies.optimizer == "mmff94":
            if not AllChem.MMFFHasAllMoleculeParams(mol):
                raise LigandSkipped("MMFF94 has no parameters for this ligand (policy: minimise)")
            AllChem.MMFFOptimizeMolecule(mol, maxIters=2000)
            notes.append("Minimised with MMFF94.")
    elif not any(a.GetAtomicNum() == 1 for a in mol.GetAtoms()):
        mol = Chem.AddHs(mol, addCoords=True)
        notes.append("No hydrogens in the record; added from valence at the charges written.")
    generated = ls.assign_atom_names(mol)
    if generated:
        notes.append(f"Atom names generated for {generated} atoms.")
    ls._set_residue_info(mol, "LIG")
    return ls.LigandBuild(mol=mol, notes=notes, source={"kind": "library", "index": entry.index, "name": entry.name})


def prepare_entry(entry: LibraryEntry, policies: BatchPolicies):
    """(LigandBuild, PreparedLigand) for one entry; LigandSkipped with the reason otherwise."""
    from .ccd_chemistry import is_metal
    from .ligand_prep import prepare_ligand, torsion_bonds
    build = build_entry(entry, policies)
    metal_charges = {a.GetProp("name"): a.GetFormalCharge() for a in build.mol.GetAtoms() if is_metal(a)} or None
    try:
        rigid = []
        if policies.rigid_conjugated:
            rigid = [b.atoms for b in torsion_bonds(build.mol, rigid_macrocycles=policies.rigid_macrocycles,
                                                     metal_charges=metal_charges) if b.kind == "conjugated single"]
        prepared = prepare_ligand(build.mol, rigid=rigid, rigid_macrocycles=policies.rigid_macrocycles,
                                  metal_charges=metal_charges)
    except Exception as error:                      # Meeko refuses; the campaign goes on
        raise LigandSkipped(f"preparation failed: {error}") from None
    if metal_charges:
        build.notes.append("Metal charges as written: " + ", ".join(f"{k} {v:+d}" for k, v in metal_charges.items()))
    return build, prepared


# ----------------------------------------------------------------- campaign ---

def _sha256(path_or_text: str, is_text: bool = False) -> str:
    digest = hashlib.sha256()
    if is_text:
        digest.update(path_or_text.encode())
    else:
        with open(path_or_text, "rb") as handle:
            for chunk in iter(lambda: handle.read(1 << 20), b""):
                digest.update(chunk)
    return digest.hexdigest()


@dataclass
class CampaignResult:
    folder: str
    rows: List[dict]
    interrupted: bool = False
    resumed: int = 0
    notes: List[str] = field(default_factory=list)


def _manifest(state, library_path: str, policies: BatchPolicies, receptor) -> dict:
    return {"library": os.path.abspath(library_path), "library_sha256": _sha256(library_path),
            "receptor_sha256": _sha256(receptor.rigid_pdbqt + "\n#flex\n" + receptor.flex_pdbqt, is_text=True),
            "box": state.box, "settings": state.settings, "policies": policies.to_dict(),
            "receptor_choices": state.to_dict()["receptor"]}


def find_resumable(base_folder: str, manifest: dict) -> Optional[str]:
    """The newest unfinished campaign folder under ``base_folder`` made with exactly this manifest."""
    if not os.path.isdir(base_folder):
        return None
    candidates = []
    for name in os.listdir(base_folder):
        folder = os.path.join(base_folder, name)
        path = os.path.join(folder, "manifest.json")
        if name.startswith("campaign_") and os.path.exists(path):
            with open(path) as handle:
                saved = json.load(handle)
            if saved.get("manifest") == manifest and not saved.get("finished"):
                candidates.append(folder)
    return sorted(candidates)[-1] if candidates else None


def _write_results(folder: str, rows: List[dict]) -> None:
    path = os.path.join(folder, "results.csv")
    with open(path + ".tmp", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=RESULT_COLUMNS)
        writer.writeheader()
        for row in sorted(rows, key=lambda r: r["index"]):
            writer.writerow({k: row.get(k, "") for k in RESULT_COLUMNS})
    os.replace(path + ".tmp", path)


def _safe(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]", "_", name)[:60]


def run_campaign(state, library_path: str, policies: BatchPolicies, base_folder: str, *,
                 resume_folder: Optional[str] = None,
                 progress: Optional[Callable[[str], None]] = None) -> CampaignResult:
    """Dock every ligand of ``library_path`` into the receptor and box of ``state``."""
    from .docking_pipeline import build_receptor
    from .docking_results import poses_from_run, write_poses_sdf
    from .docking_run import DockingSettings, SearchBox, build_ad4_maps, run_docking
    from .ligand_prep import pdbqt_atom_type
    say = progress or (lambda text: None)
    policies.validate()
    needs = [n for n in state.missing() if "ligand" not in n and "rotatable" not in n]
    if needs:
        raise ValueError("Still needed: " + ", ".join(needs) + ".")
    settings = DockingSettings(**state.settings)
    settings.validate()
    box = SearchBox(tuple(state.box["center"]), tuple(state.box["size"]), state.box["source"])
    receptor = build_receptor(state)
    manifest = _manifest(state, library_path, policies, receptor)

    if resume_folder:
        folder = resume_folder
        with open(os.path.join(folder, "manifest.json")) as handle:
            saved = json.load(handle)
        if saved.get("manifest") != manifest:
            raise ValueError("That campaign was made with different choices or a different library; start a new one.")
        rows = []
        if os.path.exists(os.path.join(folder, "results.csv")):
            with open(os.path.join(folder, "results.csv")) as handle:
                rows = [dict(r, index=int(r["index"])) for r in csv.DictReader(handle)]
    else:
        stem = _safe(os.path.splitext(os.path.basename(library_path))[0])
        folder = os.path.join(base_folder, f"campaign_{stem}_{time.strftime('%Y%m%d_%H%M%S')}")
        os.makedirs(folder)
        rows = []
    with open(os.path.join(folder, "manifest.json"), "w") as handle:
        json.dump({"manifest": manifest, "finished": False}, handle, indent=1)
    with open(os.path.join(folder, "receptor_rigid.pdbqt"), "w") as handle:
        handle.write(receptor.rigid_pdbqt)
    if receptor.flex_pdbqt:
        with open(os.path.join(folder, "receptor_flex.pdbqt"), "w") as handle:
            handle.write(receptor.flex_pdbqt)

    done = {r["index"] for r in rows}
    entries = read_library(library_path)
    say(f"{len(entries)} molecules in {os.path.basename(library_path)}; {len(done)} already done.")

    # 1. prepare every remaining ligand (fast; failures are recorded, not raised)
    prepared: Dict[int, tuple] = {}
    for entry in entries:
        if entry.index in done:
            continue
        try:
            prepared[entry.index] = prepare_entry(entry, policies)
        except LigandSkipped as reason:
            rows.append({"index": entry.index, "name": entry.name, "status": "skipped", "notes": str(reason)})
    _write_results(folder, rows)
    say(f"Prepared {len(prepared)}; skipped {sum(1 for r in rows if r['status'] == 'skipped')} (reasons in results.csv).")

    # 2. AD4: one set of maps for every type the library needs
    ad4_maps = None
    if settings.scoring == "ad4" and prepared:
        types: List[str] = []
        for _, ligand in prepared.values():
            for line in ligand.pdbqt.splitlines():
                if line.startswith(("ATOM", "HETATM")) and pdbqt_atom_type(line) not in types:
                    types.append(pdbqt_atom_type(line))
        for line in receptor.flex_pdbqt.splitlines():
            if line.startswith(("ATOM", "HETATM")) and pdbqt_atom_type(line) not in types:
                types.append(pdbqt_atom_type(line))
        maps_dir = os.path.join(folder, "ad4_maps")
        ad4_maps, map_notes = build_ad4_maps(maps_dir, receptor.rigid_pdbqt, types, box, settings.spacing)
        say("; ".join(map_notes))

    # 3. dock, one ligand at a time, results written as each finishes
    interrupted = False
    started = time.time()
    total = len(prepared)
    by_index = {e.index: e for e in entries}
    try:
        for n, (index, (build, ligand)) in enumerate(sorted(prepared.items()), 1):
            entry = by_index[index]
            workdir = os.path.join(folder, "ligands", f"{index:05d}_{_safe(entry.name)}")
            t0 = time.time()
            try:
                run = run_docking(receptor.rigid_pdbqt, receptor.flex_pdbqt, ligand.pdbqt, box, settings,
                                  workdir=workdir, ad4_maps=ad4_maps)
                poses = poses_from_run(run.poses_pdbqt, build.mol, run.energies)
                write_poses_sdf(os.path.join(workdir, "poses.sdf"), poses,
                                {"ligand": entry.name, "library_index": index, "scoring": settings.scoring,
                                 "seed": settings.seed})
                heavy = sum(1 for a in build.mol.GetAtoms() if a.GetAtomicNum() > 1)
                best = poses[0].energies
                rows.append({"index": index, "name": entry.name, "status": "docked",
                             "best_score": f"{best['total']:.3f}", "best_inter": f"{best['inter']:.3f}",
                             "score_per_heavy_atom": f"{best['total'] / heavy:.3f}", "poses": len(poses),
                             "heavy_atoms": heavy, "rotatable_bonds": len(ligand.rotatable),
                             "net_charge": Chem.GetFormalCharge(build.mol), "seconds": f"{time.time() - t0:.1f}",
                             "notes": " ".join(build.notes)})
                _write_results(folder, rows)          # on disk before anything else can interrupt
                say(f"[{n}/{total}] {entry.name}: {best['total']:.2f} kcal/mol ({time.time() - t0:.0f} s)")
            except KeyboardInterrupt:
                raise
            except Exception as error:
                rows.append({"index": index, "name": entry.name, "status": "failed", "notes": f"docking failed: {error}"})
                _write_results(folder, rows)
                say(f"[{n}/{total}] {entry.name}: failed ({error})")
    except KeyboardInterrupt:
        interrupted = True
        _write_results(folder, rows)                  # anything finished but not yet written
        say("Stopped. Finished ligands are kept; run the campaign again with the same choices to resume.")

    _write_top_poses(folder, rows)
    if not interrupted:
        with open(os.path.join(folder, "manifest.json"), "w") as handle:
            json.dump({"manifest": manifest, "finished": True}, handle, indent=1)
    return CampaignResult(folder=folder, rows=sorted(rows, key=lambda r: r["index"]), interrupted=interrupted,
                          resumed=len(done), notes=[f"{time.time() - started:.0f} s docking"])


def _write_top_poses(folder: str, rows: List[dict]) -> None:
    """Each docked ligand's best pose in one SDF, best score first."""
    docked = sorted((r for r in rows if r.get("status") == "docked"), key=lambda r: float(r["best_score"]))
    writer = Chem.SDWriter(os.path.join(folder, "top_poses.sdf"))
    try:
        for rank, row in enumerate(docked, 1):
            path = os.path.join(folder, "ligands", f"{int(row['index']):05d}_{_safe(row['name'])}", "poses.sdf")
            if not os.path.exists(path):
                continue
            best = next(iter(Chem.SDMolSupplier(path, removeHs=False)), None)
            if best is None:
                continue
            best.SetProp("_Name", row["name"])
            best.SetProp("campaign rank", str(rank))
            writer.write(best)
    finally:
        writer.close()


def ranked(rows: List[dict]) -> List[dict]:
    return sorted((r for r in rows if r.get("status") == "docked"), key=lambda r: float(r["best_score"]))
