"""Titratable-residue charge tables driven by `cpinutil.py --describe`.

cpinutil.py ships with AmberTools and is the canonical source for AMBER
constant-pH residue charges. This module wraps it, parses the output, and
exposes a small dataclass per residue.

Example
-------
    >>> from titrate.residues import get_residue
    >>> tyr = get_residue("TYR")
    >>> tyr.pka_model
    9.6
    >>> tyr.states[0].name        # the protonated reference
    'STATE 0'
    >>> tyr.states[0].charges["HH"]
    0.3992
    >>> tyr.deprot_state.charges["HH"]
    0.0

Residues with more than two states (HIP has 3, AS4 has 5) expose all
states; pick the right pair for a 2-state titration via `prot_state` /
`deprot_state` (which use `Prot Cnt` to identify them).
"""
from __future__ import annotations

import os
import re
import subprocess
import sys
from dataclasses import dataclass, field
from functools import lru_cache
from typing import Dict, List, Optional


@dataclass(frozen=True)
class ProtonationState:
    """One column from the cpinutil --describe table."""
    name: str                      # 'STATE 0', 'STATE 1', ...
    charges: Dict[str, float]      # atom_name -> partial charge
    prot_count: int                # number of titrating H atoms in this state
    pka_corr: float                # pH-dependent correction term

    @property
    def net_charge(self) -> float:
        return sum(self.charges.values())


@dataclass(frozen=True)
class TitratableResidue:
    """All AMBER constant-pH information for one residue type."""
    name: str                                  # e.g. 'TYR', 'AS4', 'HIP'
    pka_model: float                           # solution pKa from cpinutil header
    states: List[ProtonationState] = field(default_factory=list)

    @property
    def n_states(self) -> int:
        return len(self.states)

    @property
    def prot_state(self) -> ProtonationState:
        """The state with the most titrating protons (the protonated form)."""
        return max(self.states, key=lambda s: s.prot_count)

    @property
    def deprot_state(self) -> ProtonationState:
        """The state with the fewest titrating protons (the deprotonated form).

        For multi-state residues with tautomers (e.g. HIP has HID/HIE both
        sharing prot_count=1), this returns the lowest-prot_count state;
        if there are ties, it returns the first encountered. Use
        `states_with_prot_count(N)` to enumerate tautomers explicitly.
        """
        return min(self.states, key=lambda s: s.prot_count)

    def states_with_prot_count(self, n: int) -> List[ProtonationState]:
        return [s for s in self.states if s.prot_count == n]

    def state_by_name(self, name: str) -> ProtonationState:
        for s in self.states:
            if s.name == name:
                return s
        raise KeyError(f"No state named {name!r} in {self.name}")

    def state_label(self, state: ProtonationState) -> str:
        """Human-readable chemistry label for one of this residue's states.

            PROT  / DEPROT  for two-prot-count residues (AS4, GL4, LYS,
                            TYR, CYS, PRN — the dominant case)
            HIP / HID / HIE for HIP's three chemistries (HID is the
                            +0.6 kcal/mol pka_corr tautomer; HIE the 0.0
                            tautomer; HIP the doubly-protonated +1 form)
            N=k             fallback for any residue with intermediate
                            prot_count values

        Useful for human-readable state-map output alongside the literal
        cpinutil 'STATE k' name (which is what cpin actually consumes).
        """
        n = state.prot_count
        if self.name == "HIP":
            if n == 2:
                return "HIP"
            if state.pka_corr > 0.3:
                return "HID"
            return "HIE"
        n_max = max(s.prot_count for s in self.states)
        n_min = min(s.prot_count for s in self.states)
        if n == n_max:
            return "PROT"
        if n == n_min:
            return "DEPROT"
        return f"N={n}"

    @property
    def neutral_state(self) -> ProtonationState:
        """The state whose net charge is closest to 0.

        For Bashford-style intrinsic pKa calculations the canonical
        convention (MEAD, Karlsberg, H++) is to hold every non-target
        titratable in its neutral form so the protein's reference
        charge distribution is physically sensible. Concretely:

            AS4/GL4/PRN  → protonated form  (any STATE with prot_count > 0)
            LYS, HIP     → deprotonated form (neutral; HID for HIP)
            TYR, CYS     → protonated form  (neutral)

        Ties (e.g., HIP's two prot_count=1 tautomers) break by lowest
        index so the choice is reproducible run-to-run.
        """
        return min(self.states,
                    key=lambda s: (abs(s.net_charge), self.states.index(s)))


# ---------------------------------------------------------------------------
# parser
# ---------------------------------------------------------------------------

def _resolve_cpinutil() -> List[str]:
    """Return the argv prefix to run cpinutil.py with a parmed-having Python.

    Do NOT trust PATH: on a system with AmberTools installed separately (e.g.
    /usr/local/amberXX/bin ahead of the proprep conda env), a bare
    ``cpinutil.py`` resolves to that copy, whose shebang interpreter has no
    parmed, so every call dies — and because get_residue is lru_cached, the
    failure only surfaces in fork-spawned solver workers with a cold cache
    (which then masquerades as convergence). Prefer the cpinutil.py sitting
    next to THIS interpreter (the env that imported proprep), invoked with
    sys.executable so the shebang is irrelevant. Fall back to a bare name only
    if no sibling copy exists.
    """
    sibling = os.path.join(os.path.dirname(sys.executable), "cpinutil.py")
    if os.path.exists(sibling):
        return [sys.executable, sibling]
    return ["cpinutil.py"]


def _run_cpinutil(resname: str) -> str:
    """Returns stdout of `cpinutil.py --describe RESNAME`."""
    r = subprocess.run(
        _resolve_cpinutil() + ["--describe", resname],
        capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(
            f"cpinutil.py --describe {resname} failed (exit {r.returncode}):\n"
            f"{r.stderr}")
    return r.stdout


_HEADER_RE = re.compile(r"^(\S+)\s+pKa\s*=\s*(-?\d+\.?\d*)")
_ATOM_RE   = re.compile(r"^\s+([A-Z0-9]+)\s+((?:-?\d+\.\d+\s+)+-?\d+\.\d+)\s*$")
_PROTCNT_RE= re.compile(r"^Prot Cnt\s+((?:\d+\s+)+\d+)\s*$")
_PKACORR_RE= re.compile(r"^pKa Corr\s+((?:-?\d+\.\d+\s+)+-?\d+\.\d+)\s*$")


def parse_describe(text: str) -> TitratableResidue:
    """Parse cpinutil --describe output into a TitratableResidue."""
    name: Optional[str] = None
    pka_model: Optional[float] = None
    n_states: Optional[int] = None
    charges: Dict[str, List[float]] = {}
    prot_counts: List[int] = []
    pka_corrs: List[float] = []

    for line in text.splitlines():
        m = _HEADER_RE.match(line)
        if m:
            name, pka_str = m.group(1), m.group(2)
            pka_model = float(pka_str)
            continue
        m = _ATOM_RE.match(line)
        if m:
            atom = m.group(1)
            vals = [float(x) for x in m.group(2).split()]
            if n_states is None:
                n_states = len(vals)
            if atom not in charges:
                charges[atom] = vals
            continue
        m = _PROTCNT_RE.match(line)
        if m:
            prot_counts = [int(x) for x in m.group(1).split()]
            continue
        m = _PKACORR_RE.match(line)
        if m:
            pka_corrs = [float(x) for x in m.group(1).split()]
            continue
        if line.startswith("Reference Energies"):
            break  # rest is GB/PB self-energies; not used here

    if name is None or pka_model is None:
        raise ValueError("Could not find header line in cpinutil output")
    if not charges:
        raise ValueError(f"No atom charges parsed for {name}")
    if not prot_counts or not pka_corrs:
        raise ValueError(f"Missing Prot Cnt or pKa Corr for {name}")
    if not (n_states == len(prot_counts) == len(pka_corrs)):
        raise ValueError(
            f"State-count mismatch for {name}: charges={n_states}, "
            f"prot_counts={len(prot_counts)}, pka_corrs={len(pka_corrs)}")

    states: List[ProtonationState] = []
    for i in range(n_states):
        state_charges = {atom: vals[i] for atom, vals in charges.items()}
        states.append(ProtonationState(
            name=f"STATE {i}",
            charges=state_charges,
            prot_count=prot_counts[i],
            pka_corr=pka_corrs[i],
        ))
    return TitratableResidue(name=name, pka_model=pka_model, states=states)


@lru_cache(maxsize=None)
def get_residue(resname: str) -> TitratableResidue:
    """Return the TitratableResidue for a name, calling cpinutil and parsing.

    Cached: subsequent calls with the same name return the same object.
    """
    return parse_describe(_run_cpinutil(resname))


# Standard set of titratable residues in AMBER's constph framework.
STANDARD_TITRATABLE = ("AS4", "GL4", "HIP", "LYS", "TYR", "CYS")
# Bundle-specific addition:
HEME_TITRATABLE = ("PRN",)
