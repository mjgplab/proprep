"""Write a standard-residue-named PDB for a modern-FF tleap rebuild.

PB-titrate computes states on a constant-pH topology (AS4/GL4/HIP names,
sourced from oldff/ff10 which lacks the ff14SB/ff19SB dihedral & CMAP
corrections). For *conventional* fixed-charge MD the better path is to bake the
recommended protonation states into a PDB as standard residue names and rebuild
with a modern protein FF.

This module renames each titratable residue in the (dry, minimized) structure
to the standard Amber residue that encodes its recommended protonation state,
and writes a PDB the Topology Generator can pick up (workspace key
``pb_rename_pdb_file``, highest discovery priority).

Source structure: the dry prmtop+rst7 (solvent-stripped, Ca²⁺ + cofactors
intact, coordinates = the pre-PB minimization). tleap re-solvates fresh with the
water model matching the chosen FF (OPC for ff19SB, TIP3P for ff14SB).

Coverage: EVERY titratable residue must lose its constph name or tleap will
reject the PDB. Three sources are unioned, in priority order:
  1. the solve/review state map (titrated sites' recommended states),
  2. metal-coordinating residues excluded in pbt-2 (held at their fixed
     metal-bound state — deprotonated carboxylate → ASP/GLU),
  3. a standard-pH fallback for any titratable residue covered by neither
     (emits a warning; should not occur in a complete run).

Tyrosinate gap: neither ff14SB nor ff19SB defines a deprotonated-Tyr residue.
A Tyr recommended deprotonated is written as neutral TYR with a loud warning.

Heme propionate (PRN): PB-Titrate titrates each heme propionate as the
constant-pH/redox residue PRN (defined in conste.lib). The modern-FF rebuild
does not load conste.lib, so PRN is renamed to its static fixed-pH counterpart
bundled in ProPrep's heme cofactor libraries — PRP (protonated, COOH, net 0)
or PRD (deprotonated, COO⁻, net −1). Each ring's propionate is a separate PRN
residue and is renamed independently. Because PRP/PRD live only in the heme
library, the renamed PDB MUST be rebuilt with the matching ``*_FixedpH`` heme
forcefield set; build_rename_map surfaces this requirement when any PRN is
renamed.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from .residues import get_residue

# Modern-FF standard residue names for each (constph residue, protonation
# class). "prot" = protonated form, "deprot" = deprotonated form; for HIS the
# class is the tautomer label.
#
# Protein residues — verified present in both amino12.lib and amino19.lib:
#   ASP ASH GLU GLH HID HIE HIP LYS LYN CYS CYM CYX
# NOT present (the tyrosinate gap): no deprotonated-Tyr residue.
#
# Heme propionate (PRN) — the constant-pH/redox titratable propionate lives in
# conste.lib, which the modern-FF rebuild deliberately does NOT load. Its static
# fixed-pH counterparts are bundled in ProPrep's heme cofactor libraries
# (the ``*_FixedpH`` forcefield sets): PRP = protonated propionate (COOH, net 0,
# syn) and PRD = deprotonated propionate (COO⁻, net −1). These use parm10 atom
# types and need only ff14SB/ff19SB, so the renamed PDB must be rebuilt with the
# matching ``*_FixedpH`` heme set (see the build note in build_rename_map).
_RENAME = {
    "AS4": {"prot": "ASH", "deprot": "ASP"},
    "GL4": {"prot": "GLH", "deprot": "GLU"},
    "GLU": {"prot": "GLH", "deprot": "GLU"},   # if already standard-named
    "ASP": {"prot": "ASH", "deprot": "ASP"},
    "LYS": {"prot": "LYS", "deprot": "LYN"},
    "CYS": {"prot": "CYS", "deprot": "CYM"},
    "TYR": {"prot": "TYR", "deprot": "TYR"},   # deprot has no modern residue
    "PRN": {"prot": "PRP", "deprot": "PRD"},   # heme propionate (FixedpH set)
    "HIP": {"HID": "HID", "HIE": "HIE", "HIP": "HIP"},
    "HID": {"HID": "HID", "HIE": "HIE", "HIP": "HIP"},
    "HIE": {"HID": "HID", "HIE": "HIE", "HIP": "HIP"},
}

# Residue names produced only by the heme cofactor (FixedpH) libraries, not by
# any modern protein FF. A PDB containing these must be rebuilt with the
# matching ``*_FixedpH`` heme set; surfaced as a build note, not a warning.
_FIXEDPH_HEME_RESIDUES = {"PRP", "PRD"}

# Standard-pH fallback class when a titratable residue has no recommendation
# (neither in the state map nor metal-fixed). Conventional protonation at ~pH7.
_FALLBACK_CLASS = {
    "AS4": "deprot", "ASP": "deprot",
    "GL4": "deprot", "GLU": "deprot",
    "LYS": "prot", "CYS": "prot", "TYR": "prot",
    "PRN": "deprot",   # heme propionate: carboxylate (COO⁻) at ~pH7
    "HIP": "HIE", "HID": "HIE", "HIE": "HIE",
}

# Residues with no modern-FF deprotonated form — flagged when the recommendation
# is the deprotonated class.
_NO_DEPROT_RESIDUE = {"TYR"}


def _classify(resname: str, state) -> str:
    """Return the protonation class of ``state`` for ``resname``.

    For HIS-family residues: 'HID' | 'HIE' | 'HIP' (via state_label).
    For two-state residues: 'prot' | 'deprot' by comparing prot_count to the
    residue's protonated/deprotonated extremes (robust to label conventions).
    """
    res = get_residue(resname)
    if resname in ("HIP", "HID", "HIE"):
        return res.state_label(state)
    # two-prot-count residue: protonated = the max-prot_count state
    return "prot" if state.prot_count >= res.prot_state.prot_count else "deprot"


def build_rename_map(
        structure,
        state_map_named: Dict[Tuple[str, int], str],
        metal_fixed: Optional[List[Dict[str, Any]]] = None,
) -> Tuple[Dict[Tuple[str, int], str], List[Dict[str, Any]]]:
    """Build {(resname,resnum): new_resname} covering every titratable residue.

    Returns (rename_map, report) where report is a list of per-residue dicts:
        {resname, resnum, new_resname, source, klass, warning, note}
    source ∈ {'recommended','metal-fixed','fallback'}; warning is a string or
    None (set for the tyrosinate gap); note is a non-alarming string or None
    (set when the new name needs a specific rebuild FF, e.g. PRP/PRD → the
    matching ``*_FixedpH`` heme set).
    """
    metal_fixed = metal_fixed or []
    # Resolve the recommended state objects.
    rec_state = {}
    for (rn, num), name in (state_map_named or {}).items():
        res = get_residue(rn)
        st = next((s for s in res.states if s.name == name), None)
        if st is not None:
            rec_state[(rn, num)] = st

    # Metal-fixed residues → their fixed state (deprotonated/neutral).
    metal_class = {}
    for entry in metal_fixed:
        key = (entry["resname"], entry["resnum"])
        metal_class[key] = ("deprot" if entry.get("fix_mode") == "deprotonated"
                            else "prot")

    rename: Dict[Tuple[str, int], str] = {}
    report: List[Dict[str, Any]] = []

    for res in structure.residues:
        rn = res.name
        if rn not in _RENAME:
            continue  # not a titratable constph residue — leave as-is
        key = (rn, res.number + 1)
        warning = None

        if key in rec_state:
            klass = _classify(rn, rec_state[key])
            source = "recommended"
        elif key in metal_class:
            klass = metal_class[key]
            source = "metal-fixed"
        else:
            klass = _FALLBACK_CLASS.get(rn, "deprot")
            source = "fallback"
            warning = (f"no PB recommendation for {rn}-{res.number+1}; "
                       f"using standard-pH default ({klass})")

        # Tyrosinate gap: deprotonated Tyr has no modern residue.
        if rn == "TYR" and klass == "deprot":
            warning = (f"TYR-{res.number+1} predicted DEPROTONATED but no "
                       f"ff14SB/ff19SB tyrosinate residue exists — writing "
                       f"neutral TYR. Supply custom TYM params for a true "
                       f"tyrosinate.")

        new = _RENAME[rn][klass]

        # Heme propionate (PRP/PRD) is defined only in the FixedpH heme library,
        # not in any modern protein FF — flag the rebuild requirement (a note,
        # not a warning: this is expected, not a problem).
        note = None
        if new in _FIXEDPH_HEME_RESIDUES:
            note = (f"{new} is a heme propionate — rebuild with the matching "
                    f"'*_FixedpH' heme forcefield set (PRP/PRD are not in "
                    f"ff14SB/ff19SB; they need no leaprc.constph/conste).")

        rename[key] = new
        report.append({
            "resname": rn, "resnum": res.number + 1,
            "new_resname": new, "source": source, "klass": klass,
            "warning": warning, "note": note,
        })

    return rename, report


def write_renamed_pdb(
        prmtop: Path, rst7: Path, out_pdb: Path,
        rename_map: Dict[Tuple[str, int], str],
        strip_hydrogens: bool = True,
) -> int:
    """Apply ``rename_map`` to (prmtop, rst7) and write ``out_pdb``.

    Returns the number of residues renamed. Heavy-atom coordinates are
    untouched — only residue names change — so the PDB is geometrically the
    minimized structure with modern-FF protonation-state names.

    ``strip_hydrogens`` (default True): remove ALL hydrogens before writing.
    The PB-titrate topology is built with one force field (constph/ff10) whose
    hydrogen set and H-atom names differ from the modern FF the structure is
    being rebuilt under (ff14SB/ff19SB). Carrying the old H into the rebuild
    causes tleap "atom not found in residue template" / extra-atom errors and,
    for the renamed protonation states, the wrong H would be present anyway
    (e.g. an AS4 with both syn/anti carboxyl H rebuilt as ASH). Stripping H and
    letting tleap add them back per the chosen FF is the correct cross-FF
    handoff — tleap reconstructs H deterministically from the heavy-atom
    skeleton + residue template. (Mirrors the H-strip in propka_compare, where
    ProPKA likewise wants a heavy-atom PDB.)
    """
    import parmed as pmd
    s = pmd.load_file(str(prmtop), str(rst7))
    n = 0
    for res in s.residues:
        key = (res.name, res.number + 1)
        new = rename_map.get(key)
        if new and new != res.name:
            res.name = new
            n += 1
    if strip_hydrogens:
        s.strip("@/H")  # remove all hydrogens by element; tleap re-adds per FF
    out_pdb.parent.mkdir(parents=True, exist_ok=True)
    s.save(str(out_pdb), overwrite=True)
    return n
