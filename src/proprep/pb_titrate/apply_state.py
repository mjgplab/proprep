"""Apply a converged state map to a prmtop and write production files.

After mean-field or MC has settled on a state per titratable site, this
module bakes those states into the topology by setting per-atom charges.
The output is a prmtop+rst7 ready for MD.

Coordinates are unchanged. For constph-prepared topologies (AS4, GL4,
HIP, LYS, TYR, CYS) all titrating H atoms are already present in the
PROT positions, so swapping protonation states is purely a charge
update — no atom add/remove, no parmed contortions.

Net charge of the output = sum(state.net_charge for all titrated sites)
+ original non-titrated net charge. The optional `rebalance` step (see
`titrate.rebalance.rebalance_ions`) restores neutrality by removing bulk
counter-ions, ranked by distance from the protein.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Set, Tuple

import parmed as pmd

from .residues import ProtonationState
from .intrinsic import SiteKey, apply_state_charges
from .rebalance import rebalance_ions


# A finalized state map is the same shape as a background state map —
# {(resname, resnum): ProtonationState} — with one entry per site that
# the solver chose a state for.
StateMap = Dict[SiteKey, ProtonationState]


def apply_state_map(structure: pmd.Structure, state_map: StateMap) -> None:
    """Mutate `structure` in place: each site in `state_map` gets the
    charges from its assigned ProtonationState applied to its residue's
    atoms.

    Sites in `state_map` that aren't found in `structure` are silently
    skipped (legitimate when the topology is a subset of the model).
    """
    for (resname, resnum), state in state_map.items():
        for r in structure.residues:
            if r.name == resname and (r.number + 1) == resnum:
                for a in r.atoms:
                    if a.name in state.charges:
                        a.charge = state.charges[a.name]
                break


def build_production_prmtop(input_prmtop: Path,
                              input_rst7:   Path,
                              state_map:    StateMap,
                              output_prefix: Path,
                              strip_solvent: bool = False,
                              rebalance: bool = False,
                              keep_ions: Optional[Set[Tuple[Any, Any]]] = None,
                              cation_preference: Optional[List[str]] = None,
                              anion_preference:  Optional[List[str]] = None,
                              target_net_q: float = 0.0,
                              ) -> Dict[str, Any]:
    """Apply state_map to input prmtop+rst7 and write {prefix}.{prmtop,rst7}.

    Optional steps (in order):
      1. strip_solvent: remove waters and common bulk ions by name mask.
      2. apply_state_map: bake the per-residue state charges.
      3. rebalance: remove bulk single-atom ions to restore `target_net_q`
         (default 0). Skipped by default; opt-in with `rebalance=True`.
         `keep_ions` specifies single-atom residues that must NEVER be
         removed (typically structural ions from CryoEM).

    Returns a dict with diagnostic info: net charges before/after, atom
    counts, per-site charge verification, and (if rebalance ran) the
    rebalance audit.
    """
    s = pmd.load_file(str(input_prmtop), str(input_rst7))
    n_atoms_in = len(s.atoms)
    q_in = sum(a.charge for a in s.atoms)

    if strip_solvent:
        s.strip(":WAT,Na+,Cl-,K+")

    apply_state_map(s, state_map)

    # Per-site verification: confirm each residue's net charge matches
    # the state's expected net charge (within rounding).
    site_check = {}
    for (resname, resnum), state in state_map.items():
        for r in s.residues:
            if r.name == resname and (r.number + 1) == resnum:
                actual = sum(a.charge for a in r.atoms)
                expected = state.net_charge
                site_check[(resname, resnum)] = {
                    "actual": actual, "expected": expected,
                    "delta": actual - expected,
                }
                break

    q_after_apply = sum(a.charge for a in s.atoms)

    rebalance_info = None
    if rebalance:
        rebalance_info = rebalance_ions(
            s, target_net_q=target_net_q,
            keep_ions=keep_ions,
            cation_preference=cation_preference,
            anion_preference=anion_preference,
            verbose=False)

    q_out = sum(a.charge for a in s.atoms)
    out_prmtop = output_prefix.with_suffix(".prmtop")
    out_rst7   = output_prefix.with_suffix(".rst7")
    s.save(str(out_prmtop), overwrite=True)
    s.save(str(out_rst7),   overwrite=True)

    return {
        "n_atoms_in":  n_atoms_in,
        "n_atoms_out": len(s.atoms),
        "net_q_in":             q_in,
        "net_q_after_apply":    q_after_apply,
        "net_q_out":            q_out,
        "delta_q":              q_out - q_in,
        "site_check":           site_check,
        "rebalance":            rebalance_info,
        "output_prmtop":        out_prmtop,
        "output_rst7":          out_rst7,
    }


def state_map_from_sites(sites: Iterable["TitratableSite"]) -> StateMap:
    """Convenience: extract a state_map from a list of TitratableSite objects.

    Avoids importing TitratableSite at module load (circular) — accepts
    any object with `.resname`, `.resnum`, `.state` attributes.
    """
    return {(s.resname, s.resnum): s.state for s in sites}
