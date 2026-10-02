"""
Proposals for the receptor decisions ``prepare_receptor`` requires.

Nothing here builds anything; each function returns proposals with their
source, for the user to review and change before the receptor is built.

* ``propose_protonation``: pdb2pqr with PROPKA at a chosen pH (the same flags
  ProPrep's Protonation State Analyzer uses), on a copy holding only the
  selected chains' standard residues. Its residue names map onto Meeko's
  templates. pdb2pqr gets carbonic anhydrase's three zinc-bound histidines
  right (HID, HID, HIE) even though it drops the zinc, where Meeko's own
  default (HIE for every HIS) puts a hydrogen on the zinc-bound nitrogen of two.
* ``metal_contacts`` / ``metal_overrides``: a nitrogen or sulfur bound to a metal
  carries no hydrogen. Contacts come from the structure's LINK records and from
  the distance cutoff ProPrep's pb_titrate uses for inner-sphere coordination;
  each contact says which. A histidine bound through ND1 becomes HIE, through
  NE2 becomes HID; a bound cysteine becomes the thiolate CYX-.
* ``propose_termini``: pdb2pqr treats the first and last OBSERVED residues as
  charged termini, which is wrong wherever residues are missing (streptavidin
  1STP starts at residue 13). REMARK 465 decides: missing residues before the
  first observed one make it a chain break. A file without REMARK 465 gives no
  proposal; the user decides.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from rdkit import Chem

from .ccd_chemistry import is_metal, is_metal_symbol
from .dependencies import find_executable
from .receptor_prep import chosen_alternate, read_atoms

# pdb2pqr --ffout=AMBER residue name -> Meeko template key, for the residues whose state varies
_PDB2PQR_TO_MEEKO = {
    "HID": "HID", "HIE": "HIE", "HIP": "HIP",
    "ASP": "ASP", "ASH": "ASH", "GLU": "GLU", "GLH": "GLH",
    "LYS": "LYS", "LYN": "LYN",
    "CYS": "CYS", "CYM": "CYX-", "CYX": "CYX",
}
_TITRATABLE_INPUT = {"HIS", "HID", "HIE", "HIP", "ASP", "ASH", "GLU", "GLH", "LYS", "LYN", "CYS", "CYM", "CYX"}

# side-chain atoms that can bind a metal, and the template that atom's binding implies
_BINDING = {
    ("HIS", "ND1"): "HIE", ("HIS", "NE2"): "HID",
    ("CYS", "SG"): "CYX-",
    ("ASP", "OD1"): "ASP", ("ASP", "OD2"): "ASP",
    ("GLU", "OE1"): "GLU", ("GLU", "OE2"): "GLU",
    ("LYS", "NZ"): "LYN",
}
_UNSUPPORTED_BINDING = {("TYR", "OH"): "a metal-bound tyrosinate, for which Meeko has no template"}


@dataclass
class Proposal:
    res_id: str
    resname: str
    value: str
    source: str


@dataclass
class ProtonationResult:
    proposals: Dict[str, Proposal]
    ph: float
    ran: bool
    message: str = ""


def _protein_only_copy(pdb_path: str, chains: Sequence[str], altlocs: Dict[str, str], out_path: str,
                       altloc_fills: Optional[Dict[str, Dict[str, str]]] = None) -> int:
    """ATOM records of ``chains`` with the chosen alternate locations; returns the atom count."""
    written = 0
    with open(out_path, "w") as out:
        for atom in read_atoms(pdb_path):
            if atom.record != "ATOM" or atom.chain not in chains:
                continue
            if not chosen_alternate(atom, altlocs, altloc_fills):
                continue
            out.write(atom.line[:16] + " " + atom.line[17:] + "\n")
            written += 1
        out.write("END\n")
    return written


def propose_protonation(pdb_path: str, chains: Sequence[str], ph: float,
                        altlocs: Dict[str, str],
                        altloc_fills: Optional[Dict[str, Dict[str, str]]] = None) -> ProtonationResult:
    """Titration states from pdb2pqr + PROPKA at ``ph``; never raises on a pdb2pqr failure.

    A failure (pdb2pqr gives up on a structure it cannot assign charges to, for
    example an acetyl cap) returns ``ran=False`` with pdb2pqr's own reason, and
    no proposals: the user decides every titratable residue.
    """
    executable = find_executable("pdb2pqr") or find_executable("pdb2pqr30")
    if executable is None:
        return ProtonationResult({}, ph, False, "pdb2pqr is not installed.")
    workdir = tempfile.mkdtemp(prefix="proprep_pdb2pqr_")
    try:
        source = os.path.join(workdir, "protein.pdb")
        if _protein_only_copy(pdb_path, chains, altlocs, source, altloc_fills) == 0:
            return ProtonationResult({}, ph, False, "The selected chains have no ATOM records.")
        pqr = os.path.join(workdir, "protein.pqr")
        command = [executable, "--ff=AMBER", "--ffout=AMBER", "--titration-state-method=propka",
                   f"--with-ph={ph}", "--keep-chain", source, pqr]
        run = subprocess.run(command, capture_output=True, text=True, timeout=1800)
        if run.returncode != 0 or not os.path.exists(pqr):
            reason = [l for l in (run.stderr + run.stdout).splitlines() if "CRITICAL" in l or "Error" in l]
            return ProtonationResult({}, ph, False, "pdb2pqr failed: " + (" | ".join(reason[-3:]) or
                                                                          f"exit code {run.returncode}"))
        names: Dict[str, str] = {}
        original = {a.res_id: a.resname for a in read_atoms(pdb_path) if a.record == "ATOM"}
        with open(pqr) as handle:
            for line in handle:
                if line.startswith(("ATOM", "HETATM")):
                    res_id = f"{line[21]}:{int(line[22:26])}{line[26].strip()}"
                    names[res_id] = line[17:20].strip()
        proposals = {}
        for res_id, name in names.items():
            if original.get(res_id) in _TITRATABLE_INPUT and name in _PDB2PQR_TO_MEEKO:
                proposals[res_id] = Proposal(res_id, original[res_id], _PDB2PQR_TO_MEEKO[name],
                                             f"pdb2pqr/PROPKA at pH {ph}")
        return ProtonationResult(proposals, ph, True, f"pdb2pqr assigned {len(proposals)} titratable residues.")
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


# ---------------------------------------------------------------- metals ---

@dataclass
class LinkRecord:
    atom1: Tuple[str, str, str, str]      # (res_id, resname, atom name, altloc)
    atom2: Tuple[str, str, str, str]
    symmetry1: str
    symmetry2: str
    distance: Optional[float]

    @property
    def same_asymmetric_unit(self) -> bool:
        """Both symmetry operators are identity (1555), or blank, which the PDB format reads as identity."""
        return self.symmetry1 in ("", "1555") and self.symmetry2 in ("", "1555")


def read_link_records(pdb_path: str) -> List[LinkRecord]:
    """LINK records with their alternate locations and symmetry operators.

    ProPrep's PDBMetadataExtractor drops both; a LINK to a symmetry mate (e.g.
    4565) then reads as a bond inside the asymmetric unit, 12-27 A long in 2D5M.
    """
    records = []
    with open(pdb_path) as handle:
        for line in handle:
            if not line.startswith("LINK"):
                continue
            line = line.rstrip("\n").ljust(80)
            def atom(offset):
                name, alt, resname = line[offset:offset + 4].strip(), line[offset + 4].strip(), line[offset + 5:offset + 8].strip()
                chain, resseq, icode = line[offset + 9], line[offset + 10:offset + 14].strip(), line[offset + 14].strip()
                return (f"{chain}:{int(resseq)}{icode}", resname, name, alt)
            distance = line[73:78].strip()
            records.append(LinkRecord(atom(12), atom(42), line[59:65].strip(), line[66:72].strip(),
                                      float(distance) if distance else None))
    return records


@dataclass
class MetalContact:
    metal_res: str
    metal_resname: str
    metal_atom: str
    res_id: str
    resname: str
    atom: str
    distance: float
    source: str          # "LINK record", "distance", or both


def metal_contacts(pdb_path: str, metal_residues: Sequence[str], altlocs: Dict[str, str],
                   altloc_fills: Optional[Dict[str, Dict[str, str]]] = None) -> List[MetalContact]:
    """Side-chain atoms bound to the metals in ``metal_residues`` (by residue id).

    A contact counts when the structure's LINK records declare it, or when the
    atom lies within pb_titrate's inner-sphere cutoff of the metal.
    """
    from proprep.pb_titrate.metal_coordination import DEFAULT_COORDINATION_CUTOFF

    atoms = [a for a in read_atoms(pdb_path) if chosen_alternate(a, altlocs, altloc_fills)]
    metals = [a for a in atoms if a.res_id in metal_residues and a.element and is_metal_symbol(a.element)]
    candidates = [a for a in atoms if (a.resname, a.name) in _BINDING or (a.resname, a.name) in _UNSUPPORTED_BINDING]

    def chosen(res_id, name, alt):
        from proprep.structure_prep.altloc_picker import keeps_atom
        return keeps_atom(name, alt, altlocs.get(res_id), (altloc_fills or {}).get(res_id))

    declared = set()
    for link in read_link_records(pdb_path):
        if not link.same_asymmetric_unit or not chosen(link.atom1[0], link.atom1[2], link.atom1[3]) or not chosen(link.atom2[0], link.atom2[2], link.atom2[3]):
            continue
        a, b = (link.atom1[0], link.atom1[2]), (link.atom2[0], link.atom2[2])
        declared.add((a, b))
        declared.add((b, a))
    contacts = []
    for metal in metals:
        for atom in candidates:
            distance = float(np.linalg.norm(np.subtract(metal.xyz, atom.xyz)))
            linked = ((metal.res_id, metal.name), (atom.res_id, atom.name)) in declared
            near = distance <= DEFAULT_COORDINATION_CUTOFF
            if linked or near:
                source = " and ".join(s for s, ok in (("LINK record", linked),
                                                      (f"distance <= {DEFAULT_COORDINATION_CUTOFF} A", near)) if ok)
                contacts.append(MetalContact(metal.res_id, metal.resname, metal.name, atom.res_id, atom.resname,
                                             atom.name, round(distance, 2), source))
    return sorted(contacts, key=lambda c: (c.metal_res, c.distance))


def symmetry_mate_contacts(pdb_path: str, metal_residues: Sequence[str]) -> List[str]:
    """LINK records binding a kept metal to an atom of a symmetry mate: an open site in this copy."""
    notes = []
    for link in read_link_records(pdb_path):
        if link.same_asymmetric_unit:
            continue
        for metal_side, other_side, other_symmetry in ((link.atom1, link.atom2, link.symmetry2),
                                                       (link.atom2, link.atom1, link.symmetry1)):
            if metal_side[0] in metal_residues:
                symmetry = link.symmetry1 if other_side is link.atom1 else link.symmetry2
                notes.append(f"{metal_side[1]} {metal_side[0]} is bound to {other_side[1]} {other_side[0]} "
                             f"{other_side[2]} of a symmetry mate (operator {symmetry}); that site is open "
                             "in the single copy being docked into.")
    return notes


def metal_overrides(contacts: Sequence[MetalContact]) -> Tuple[Dict[str, Proposal], List[str]]:
    """Templates implied by metal binding, and problems no template can express."""
    bound: Dict[str, set] = {}
    info: Dict[str, MetalContact] = {}
    for contact in contacts:
        bound.setdefault(contact.res_id, set()).add(contact.atom)
        info.setdefault(contact.res_id, contact)
    proposals, problems = {}, []
    for res_id, atom_names in bound.items():
        contact = info[res_id]
        if contact.resname == "HIS" and atom_names >= {"ND1", "NE2"}:
            problems.append(f"HIS {res_id} binds metals through both ND1 and NE2 (an imidazolate); "
                            "Meeko has no template for it.")
            continue
        unsupported = [(contact.resname, a) for a in atom_names if (contact.resname, a) in _UNSUPPORTED_BINDING]
        if unsupported:
            problems.append(f"{contact.resname} {res_id} {unsupported[0][1]} binds a metal: "
                            f"{_UNSUPPORTED_BINDING[unsupported[0]]}.")
            continue
        templates = {_BINDING[(contact.resname, a)] for a in atom_names if (contact.resname, a) in _BINDING}
        if len(templates) == 1:
            atoms_text = "/".join(sorted(atom_names))
            proposals[res_id] = Proposal(res_id, contact.resname, templates.pop(),
                                         f"{atoms_text} bound to {contact.metal_resname} {contact.metal_res} "
                                         f"({contact.source})")
    return proposals, problems


# --------------------------------------------------------------- termini ---

_CAPS = ("ACE", "NME", "NH2")


@dataclass
class ChainEnd:
    res_id: str
    resname: str
    end: str          # "N" (first observed residue of the chain) or "C" (last)
    capped: bool      # next to an acetyl or amide cap: no terminus choice


def chain_ends(pdb_path: str, chains: Sequence[str]) -> List[ChainEnd]:
    """The first and last observed polymer residue of each chain, in chain order."""
    atoms = read_atoms(pdb_path)
    ends = []
    for chain in chains:
        polymer = [a for a in atoms if a.record == "ATOM" and a.chain == chain]
        if not polymer:
            continue
        caps = {a.resseq for a in atoms if a.chain == chain and a.resname in _CAPS}
        order = list(dict.fromkeys((a.resseq, a.icode, a.resname) for a in polymer))
        for end, (num, icode, name) in (("N", order[0]), ("C", order[-1])):
            capped = (end == "N" and num - 1 in caps) or (end == "C" and num + 1 in caps)
            ends.append(ChainEnd(f"{chain}:{num}{icode}", name, end, capped))
    return ends


def _seqres(pdb_path: str) -> Dict[str, List[str]]:
    """Chain -> residue names from the SEQRES records (empty when the file has none)."""
    out: Dict[str, List[str]] = {}
    with open(pdb_path) as handle:
        for line in handle:
            if line.startswith("SEQRES"):
                out.setdefault(line[11], []).extend(line[19:70].split())
    return out


def propose_termini(pdb_path: str, chains: Sequence[str]) -> Dict[str, Proposal]:
    """Each uncapped chain end: a true terminus or a chain break, with the evidence for it.

    * REMARK 465 (residues in the sequence but not modelled): residues
      missing beyond an end make it a chain break; none missing, a terminus.
    * Without REMARK 465, SEQRES: when every SEQRES residue of the chain is
      modelled, in order, nothing is missing and both ends are the chain's
      own termini (a file lists REMARK 465 only when something is missing).
    * Otherwise an OXT on the last residue: the model has a free carboxylate.
    An end with none of these is left without a proposal, for the user.
    A residue next to an acetyl or amide cap (ACE, NME, NH2) needs no choice.
    """
    from proprep.structure_prep.pdb_loader import PDBMetadataExtractor
    metadata = PDBMetadataExtractor()
    metadata.parse_pdb_file(pdb_path)
    with open(pdb_path) as handle:
        has_remark_465 = any(line.startswith("REMARK 465") for line in handle)
    seqres = _seqres(pdb_path)
    atoms = read_atoms(pdb_path)
    proposals: Dict[str, Proposal] = {}
    for end in chain_ends(pdb_path, chains):
        if end.capped:
            continue
        chain = end.res_id.split(":", 1)[0]
        num = int("".join(ch for ch in end.res_id.split(":", 1)[1] if ch.isdigit() or ch == "-"))
        terminus = "N-terminus" if end.end == "N" else "C-terminus"
        if has_remark_465:
            missing = [r["residue_number"] for r in metadata.missing_res_records.get(chain, [])]
            gap = any(m < num for m in missing) if end.end == "N" else any(m > num for m in missing)
            proposals[end.res_id] = (
                Proposal(end.res_id, end.resname, "chain break", "REMARK 465: residues missing beyond this end")
                if gap else
                Proposal(end.res_id, end.resname, terminus, "REMARK 465: no residues missing beyond this end"))
            continue
        names = seqres.get(chain)
        if names:
            wanted = set(names)
            observed = [n for (_, _, n) in dict.fromkeys(
                (a.resseq, a.icode, a.resname) for a in atoms if a.chain == chain and a.resname in wanted)]
            if observed == names:
                proposals[end.res_id] = Proposal(
                    end.res_id, end.resname, terminus,
                    f"SEQRES: all {len(names)} residues modelled, none missing")
                continue
        if end.end == "C" and any(a.res_id == end.res_id and a.name == "OXT" for a in atoms):
            proposals[end.res_id] = Proposal(end.res_id, end.resname, terminus, "OXT modelled on this residue")
    return proposals
