"""Detect titratable residues that directly coordinate a metal ion.

A titrating side chain in inner-sphere contact with a metal cation is not
pH-titratable in the holo state — its protonation is metal-locked — and a
fixed-background continuum PB treatment cannot handle the charged ligand atom
sitting ~2-2.5 Å from a +2 cation (it produces wildly shifted, meaningless
pKas). Such residues should be (a) dropped from the titratable set and (b) held
in the background at the charge state consistent with their coordination, which
also restores real charge to the cluster and improves the environment every
*other* site sees.

Detection here is purely geometric (closest side-chain atom to any metal ion),
independent of `detected_redox_sites`: structural cations like calmodulin's
EF-hand Ca²⁺ are usually not registered as redox sites, so a distance scan is
the reliable signal.
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple

# Coordinating metals of interest: structural / catalytic cations that a
# titratable side chain can ligate. Keyed by atomic number so detection is
# name-agnostic (handles CA / Ca2+ / CAL etc.). Bulk monovalent ions
# (Na/K/Cl) are intentionally excluded — they are stripped before PB and do
# not form the kind of inner-sphere bond that locks protonation.
METAL_ATOMIC_NUMBERS = {
    12: "Mg", 20: "Ca", 25: "Mn", 26: "Fe",
    27: "Co", 28: "Ni", 29: "Cu", 30: "Zn",
}

# Default inner-sphere cutoff (Å). In the calmodulin/NOS test, every direct
# carboxylate-O→Ca contact was ≤ 2.5 Å with a clean gap to 4.6 Å before the
# next residue, so 3.0 Å cleanly isolates true coordinators.
DEFAULT_COORDINATION_CUTOFF = 3.0

# Side-chain atoms that actually ligate a metal, per residue. Used to report
# the coordinating atom and to keep the distance test on chemically relevant
# atoms (not, e.g., a backbone atom that happens to be close).
_LIGATING_ATOMS = {
    "AS4": {"OD1", "OD2"}, "GL4": {"OE1", "OE2"},
    "ASP": {"OD1", "OD2"}, "GLU": {"OE1", "OE2"},
    "CYS": {"SG"}, "TYR": {"OH"},
    "HIP": {"ND1", "NE2"}, "HID": {"ND1", "NE2"}, "HIE": {"ND1", "NE2"},
    "LYS": {"NZ"},
}


def _is_metal_atom(atom) -> Optional[str]:
    """Return the metal symbol if ``atom`` is a lone metal cation, else None."""
    if atom.atomic_number in METAL_ATOMIC_NUMBERS and len(atom.residue.atoms) == 1:
        return METAL_ATOMIC_NUMBERS[atom.atomic_number]
    return None


def find_metal_coordinated_sites(
        structure,
        sites: List[Any],
        cutoff: float = DEFAULT_COORDINATION_CUTOFF,
) -> List[Dict[str, Any]]:
    """Find which ``sites`` directly coordinate a metal ion in ``structure``.

    Returns one record per coordinated site (sorted by residue number):
        {resname, resnum, n_metals, min_distance, coord_atom,
         metals: [{element, resname, resnum, distance, coord_atom}, ...]}
    A site qualifies if any of its ligating side-chain atoms is within
    ``cutoff`` of a lone metal cation.
    """
    import numpy as np

    xyz = structure.coordinates
    # Collect metal atoms once.
    metals = []
    for a in structure.atoms:
        sym = _is_metal_atom(a)
        if sym is not None:
            metals.append((a, sym))
    if not metals:
        return []

    records: List[Dict[str, Any]] = []
    for site in sites:
        res = site.titrating_residue
        ligators = _LIGATING_ATOMS.get(site.resname)
        # Atoms to test: the known ligating atoms if we recognize the residue,
        # else every side-chain heavy atom as a fallback.
        cand_atoms = [a for a in res.atoms
                      if (ligators is None and a.atomic_number > 1)
                      or (ligators is not None and a.name in ligators)]
        if not cand_atoms:
            continue

        per_metal: Dict[Tuple[str, int], Dict[str, Any]] = {}
        for a in cand_atoms:
            pa = np.asarray(xyz[a.idx], dtype=float)
            for m_atom, sym in metals:
                d = float(np.linalg.norm(pa - np.asarray(xyz[m_atom.idx],
                                                          dtype=float)))
                if d <= cutoff:
                    mkey = (m_atom.residue.name, m_atom.residue.number + 1)
                    prev = per_metal.get(mkey)
                    if prev is None or d < prev["distance"]:
                        per_metal[mkey] = {
                            "element": sym,
                            "resname": m_atom.residue.name,
                            "resnum": m_atom.residue.number + 1,
                            "distance": round(d, 2),
                            "coord_atom": a.name,
                        }
        if per_metal:
            metals_sorted = sorted(per_metal.values(),
                                   key=lambda m: m["distance"])
            best = metals_sorted[0]
            records.append({
                "resname": site.resname,
                "resnum": site.resnum,
                "n_metals": len(metals_sorted),
                "min_distance": best["distance"],
                "coord_atom": best["coord_atom"],
                "metals": metals_sorted,
            })

    records.sort(key=lambda r: r["resnum"])
    return records


# Fix-state classification --------------------------------------------------
#
# What charge state to hold a coordinator at depends on HOW it ligates, not on
# the metal count alone:
#   - carboxylate (AS4/GL4, and PRN — the heme propionate) / thiolate (CYS) /
#     phenolate (TYR): the donor atom is deprotonated → fix at the residue's
#     deprotonated state (net -1). This is unambiguous and assigned
#     automatically. (PB-Titrate only runs on conste+constph topologies, so the
#     propionate is always the switchable PRN, never a static PRD/PRP.)
#   - HIS via a ring N: a single-metal His ligates through a NEUTRAL imidazole
#     (HID/HIE); a bridging His (≥2 metals) is a histidinate (-1) — but the
#     constph HIP residue has no -1 state, so that can't be represented from the
#     standard state set and must fall back to the prmtop charges. So HIS is
#     never auto-assigned; the caller asks.
#   - LYS: rarely a true coordinator; ask.

_AUTO_DEPROTONATED = {"AS4", "GL4", "ASP", "GLU", "CYS", "TYR", "PRN"}


def classify_fix(resname: str, n_metals: int) -> Dict[str, Any]:
    """Suggest how to fix a coordinator's background charge.

    Returns {fix_mode, needs_prompt, note} where fix_mode is one of
    'deprotonated' | 'neutral' | None (None = no representable state; exclude
    only and leave prmtop charges).
    """
    if resname in _AUTO_DEPROTONATED:
        return {"fix_mode": "deprotonated", "needs_prompt": False,
                "note": "carboxylate/thiolate/phenolate → deprotonated (-1)"}
    if resname == "HIP":
        if n_metals >= 2:
            return {"fix_mode": None, "needs_prompt": True,
                    "note": ("bridging His = histidinate (-1); constph has no "
                             "-1 state — choose neutral or keep prmtop charges")}
        return {"fix_mode": "neutral", "needs_prompt": True,
                "note": "single-metal His → neutral imidazole (HID/HIE)"}
    return {"fix_mode": None, "needs_prompt": True,
            "note": "uncommon coordinator — choose how to fix"}
