"""Site abstraction for pb_titrate.

A *site* in pb_titrate is keyed by a single titrating residue but carries
an *envelope* — the set of structurally coupled residues whose partial
charges sum to integer only as a group (e.g., MCPB-parameterized hemes:
HCO + 2×PRN + DHO + PHO + CYO).

For standard titratable residues (AS4/GL4/HIP/LYS/CYS/TYR not embedded in
a redox site), the envelope is just the residue itself — a single-residue
site, identical in behavior to titrate's original code path.

For titratable residues that belong to a `detected_redox_sites` entry
(e.g., heme PRN), the envelope is *all* residues of that redox-site
entry, so PB calculations always see integer total charge over the group.

Discovery merges two sources:
  1. A topology scan for every titratable residue name.
  2. ProPrep's `detected_redox_sites` workspace key — a list of
     `RedoxSite` objects (or compatible mappings) from
     `comprehensive_redox_detector`.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, Iterable, List, Optional, Set, Tuple

import parmed as pmd


# Residue names that pb_titrate treats as titratable. Matches
# tleap_input_generator._scan_topology_for_titratable_residues so a
# constph prmtop produced by ProPrep's Topology Generator is the
# canonical input.
TITRATABLE_RESNAMES: Set[str] = {"AS4", "GL4", "HIP", "LYS", "CYS", "TYR", "PRN"}


# A residue identifier used to bridge ProPrep's PDB-centric metadata
# (chain + 1-based resid + insertion_code) and parmed's structure model.
ResKey = Tuple[str, int, str]   # (chain, resid_1based, insertion_code)


@dataclass
class Site:
    """A pb_titrate titration site.

    `titrating_residue` is the parmed residue whose partial charges flip
    across protonation states. `envelope_residues` is the list of all
    residues (including the titrating one) that must travel together
    through PB cluster cutouts to keep integer charge — typically just
    `[titrating_residue]` for standard residues, or the full MCPB group
    for hemes / FeS / similar.

    `other_redox_envelopes` is the list of *other* redox-site envelopes
    in the system (each as a list of residue indices). The cluster
    cutout uses this to enforce the all-or-nothing rule on overlap:
    if any residue from another redox site is inside the cluster sphere,
    its full envelope is included, and partial inclusion is forbidden.
    """
    titrating_residue: pmd.Residue
    envelope_residues: List[pmd.Residue]
    redox_site_id: Optional[str] = None
    other_redox_envelopes: List[List[int]] = field(default_factory=list)

    @property
    def resname(self) -> str:
        return self.titrating_residue.name

    @property
    def resnum(self) -> int:
        """1-based residue number (PDB convention)."""
        return self.titrating_residue.number + 1

    @property
    def envelope_idxs(self) -> Set[int]:
        """0-based parmed residue indices for the full envelope."""
        return {r.idx for r in self.envelope_residues}

    @property
    def envelope_atoms(self) -> List[pmd.Atom]:
        """All atoms across the envelope (heavy + H)."""
        return [a for r in self.envelope_residues for a in r.atoms]

    @property
    def is_multi_residue(self) -> bool:
        return len(self.envelope_residues) > 1

    def __repr__(self) -> str:
        n = len(self.envelope_residues)
        if n == 1:
            return f"Site({self.resname}-{self.resnum})"
        return (f"Site({self.resname}-{self.resnum}, envelope={n} residues, "
                f"redox_site={self.redox_site_id})")


# ---------------------------------------------------------------------------
# Bridge: detected_redox_sites entries → envelope of parmed residues
# ---------------------------------------------------------------------------

def _residue_keys_in_redox_site(redox_site: Any) -> Set[ResKey]:
    """Pull the set of (chain, resid_1based, insertion_code) keys identifying
    every residue that belongs to a `RedoxSite` (or a compatible mapping).

    Two paths:
      (a) `redox_site.residue_groups` — the canonical attribute on
          `RedoxSite` from comprehensive_redox_detector. Keys are
          `(chain, resid, insertion_code)` tuples already.
      (b) Fallback: walk `redox_site.atoms` and collect (chain, resid, ic)
          per atom.

    Either path is acceptable; (a) is preferred because it's already the
    de-duplicated residue-level view.
    """
    rg = getattr(redox_site, "residue_groups", None)
    if isinstance(rg, dict) and rg:
        return {(str(c), int(rid), str(ic) or "") for (c, rid, ic) in rg.keys()}

    out: Set[ResKey] = set()
    atoms = getattr(redox_site, "atoms", None) or []
    for a in atoms:
        chain = getattr(a, "chain", "")
        resid = getattr(a, "resid", None)
        ic    = getattr(a, "insertion_code", "") or ""
        if resid is None:
            continue
        out.add((str(chain), int(resid), str(ic)))
    return out


def _parmed_reskey(r: pmd.Residue) -> ResKey:
    """Build the (chain, resid_1based, insertion_code) key for a parmed residue.

    parmed exposes `r.chain` (often "" for prmtop without chain info) and
    `r.number` (0-based in some builds, but the value mirrors the PDB
    number that tleap wrote). For pb_titrate the "1-based" convention
    used everywhere is `r.number + 1`, matching cpinutil expectations.
    """
    return (
        str(getattr(r, "chain", "") or ""),
        int(r.number) + 1,
        str(getattr(r, "insertion_code", "") or ""),
    )


def resolve_envelope(structure: pmd.Structure,
                       redox_site: Any) -> List[pmd.Residue]:
    """Return the parmed residues that constitute the envelope of one
    `detected_redox_sites` entry.

    Matches by `(chain, resid_1based, insertion_code)`. Residues present
    in the redox site's residue_groups but not found in the parmed
    structure are silently skipped — typical for hydrogens-only or
    metadata-only entries that don't survive into the prmtop.
    """
    target_keys = _residue_keys_in_redox_site(redox_site)
    if not target_keys:
        return []
    out = []
    for r in structure.residues:
        if _parmed_reskey(r) in target_keys:
            out.append(r)
    return out


# ---------------------------------------------------------------------------
# Discovery
# ---------------------------------------------------------------------------

def _build_envelope_index(structure: pmd.Structure,
                            redox_sites: Iterable[Any],
                            ) -> Tuple[Dict[ResKey, Any],
                                       Dict[Any, List[pmd.Residue]]]:
    """For each redox site, resolve its envelope of parmed residues, then
    invert to a map from each residue's ResKey → owning redox site.

    Returns:
        (residue_to_owner, owner_to_envelope)
    """
    owner_to_envelope: Dict[Any, List[pmd.Residue]] = {}
    residue_to_owner:  Dict[ResKey, Any] = {}
    for rs in redox_sites or []:
        envelope = resolve_envelope(structure, rs)
        if not envelope:
            continue
        owner_to_envelope[id(rs)] = envelope  # use id() as a stable handle
        for r in envelope:
            residue_to_owner[_parmed_reskey(r)] = rs
    return residue_to_owner, owner_to_envelope


_NTERM_MARKERS = {"H1", "H2", "H3"}
_CTERM_MARKERS = {"OXT"}


def _is_terminal(residue: pmd.Residue) -> Optional[str]:
    """Return 'N-terminal', 'C-terminal', or None for an amino-acid residue.

    Detected via atom-name markers Amber writes for terminal residues:
      H1/H2/H3  → N-terminal (the three H atoms on the protonated N⁺)
      OXT       → C-terminal (the extra carboxylate oxygen)
    Returns the first match found; both is rare (single-residue peptide).
    """
    names = {a.name for a in residue.atoms}
    if names & _NTERM_MARKERS:
        return "N-terminal"
    if names & _CTERM_MARKERS:
        return "C-terminal"
    return None


def discover_sites(structure: pmd.Structure,
                    redox_sites: Optional[Iterable[Any]] = None,
                    *,
                    titratable_resnames: Set[str] = TITRATABLE_RESNAMES,
                    ) -> Tuple[List[Site], List[Tuple[Site, str]]]:
    """Build the unified list of titration sites for a constph prmtop.

    For each titratable residue (by name) in `structure`:
      - If it belongs to a `detected_redox_sites` entry, its envelope is
        the full residue group of that entry (multi-residue site).
      - Otherwise its envelope is just itself (single-residue site).

    Terminal residues (N- or C-terminus, detected via H1/H2/H3 or OXT
    atom markers) are EXCLUDED from the titration list because Amber's
    constph framework lacks multi-state libraries for terminal forms
    (constph.lib only has internal AS2/AS4/GL4; aminoct10.lib only has
    single-state CLYS/CASP/etc.). Excluded sites are returned in a
    second list for the caller to surface to the user.

    Returns:
        (titration_sites, excluded_terminal_sites)
          where excluded is a list of (Site, "N-terminal"|"C-terminal").

    Sites carry a list of *other* redox-site envelopes so cluster cutouts
    can enforce the all-or-nothing rule on overlap with non-target redox
    sites (see `cluster_for_site` in `cluster.py`).
    """
    redox_sites = list(redox_sites or [])
    residue_to_owner, owner_to_envelope = _build_envelope_index(
        structure, redox_sites)

    # Snapshot of every other redox-site's residue-index list, used for
    # the partial-overlap fix in cluster cutouts.
    all_envelopes_by_owner_id: Dict[int, List[int]] = {
        oid: [r.idx for r in env]
        for oid, env in owner_to_envelope.items()
    }

    sites: List[Site] = []
    excluded: List[Tuple[Site, str]] = []
    for r in structure.residues:
        if r.name not in titratable_resnames:
            continue
        owner = residue_to_owner.get(_parmed_reskey(r))
        if owner is not None:
            envelope = owner_to_envelope[id(owner)]
            redox_site_id = getattr(owner, "site_id", None)
            other_owner_id = id(owner)
        else:
            envelope = [r]
            redox_site_id = None
            other_owner_id = None
        other_envelopes = [
            idx_list
            for oid, idx_list in all_envelopes_by_owner_id.items()
            if oid != other_owner_id
        ]
        site = Site(
            titrating_residue=r,
            envelope_residues=envelope,
            redox_site_id=redox_site_id,
            other_redox_envelopes=other_envelopes,
        )
        # Skip residues at the chain termini — constph has no multi-state
        # terminal libraries, so titrating them would silently corrupt
        # backbone charges by overwriting C-terminal/N-terminal-specific
        # backbone params with cpinutil's internal-residue values, and
        # terminal-only atoms (OXT, H1/H2/H3) would be left orphaned.
        # PRN doesn't have N/C-termini so this never affects redox sites.
        terminal_kind = _is_terminal(r) if r.name != "PRN" else None
        if terminal_kind is not None:
            excluded.append((site, terminal_kind))
            continue
        sites.append(site)
    return sites, excluded
