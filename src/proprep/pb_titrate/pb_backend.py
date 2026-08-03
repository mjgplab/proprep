"""AmberTools `pbsa` wrapper.

Single-call PBSA driver. Saves a (parmed) Structure as prmtop+inpcrd,
writes a `pbsa.in` namelist, invokes pbsa, and parses the energy block
out of the output file.

All PB parameters are explicit kwargs with sensible defaults — this is
the one place in the package where the underlying solver knobs live.

Returns a dict of energy components (kcal/mol):
    EPB, EELEC, EEL14, VDWAALS, ECAVITY, EDISPER, Etot
"""
from __future__ import annotations

import os
import re
import subprocess
import warnings
from pathlib import Path
from typing import Dict, Optional

import parmed as pmd


# Env var pointing at an append-only log for non-fatal pbsa exit-time crashes
# (the known rdparm2 teardown bug: pbsa writes a complete, correct result then
# crashes during cleanup). Parallel pb_titrate runs (mean-field, coupling) set
# this to a per-run file and tally it for a single summary line, instead of
# flooding the console with one RuntimeWarning per PB call across thousands of
# calls. When unset (standalone use), the warning is emitted as before.
PBSA_SALVAGE_LOG_ENV = "PROPREP_PBSA_SALVAGE_LOG"


def _record_pbsa_salvage(label: str, returncode: int, epb) -> None:
    """Record one salvaged pbsa exit-time crash (complete result, nonzero exit).

    Appends a line to the log named by ``$PROPREP_PBSA_SALVAGE_LOG`` if set
    (POSIX appends from concurrent workers are atomic for short lines); else
    falls back to a RuntimeWarning so the information is never silently lost.
    """
    msg = (f"pbsa exited {returncode} for {label} but wrote a complete result "
           f"(EPB={epb}); using it. Non-fatal exit-time crash in pbsa, not a "
           f"calculation failure.")
    log_path = os.environ.get(PBSA_SALVAGE_LOG_ENV)
    if log_path:
        try:
            with open(log_path, "a") as f:
                f.write(f"{label}\trc={returncode}\tEPB={epb}\n")
            return
        except OSError:
            pass  # fall through to the warning if the log can't be written
    warnings.warn(msg, RuntimeWarning)

# Physical constant: RT·ln10 at 298.15 K, units kcal/mol per pH unit.
RT_LN10 = 1.36358


def _recenter_for_format(structure: "pmd.Structure", margin: float = 10.0) -> None:
    """Translate `structure` in place so all coordinates fit Amber's F12.7.

    Amber inpcrd/rst7 files are fixed-column ``(6F12.7)`` with NO delimiter
    between fields. A coordinate whose magnitude reaches 1000 Å needs 13
    columns (``-1234.5678901``) but the field is only 12 wide, so ParmEd's
    ``%12.7f`` writer spills one column and shifts every following number on
    the line. pbsa then reads misaligned tokens under the F edit descriptor
    and gfortran's formatted-read machinery walks off into invalid memory
    (SIGSEGV in ``formatted_transfer_scalar_read``) — a *read*-time crash,
    before any solve, so it is fully deterministic and survives every retry.

    Large multi-chain assemblies (e.g. a multi-heme build whose late chains
    sit far from the deposited origin) can easily push a sliced cluster's
    atoms past ±1000 Å. PB energy is translation-invariant, so shifting the
    whole cluster so its minimum corner lands at ``+margin`` on each axis is
    a physical no-op that guarantees small, all-positive coordinates well
    inside the format width.

    Only safe on the energy-only path: callers that emit a phi map (coupling)
    interpolate that grid against coordinates from the *original* frame, so
    recentering there would desync the two frames. ``run_pbsa`` therefore
    skips this when ``phiout`` is set.
    """
    xyz = structure.coordinates
    if xyz is None or len(xyz) == 0:
        return
    # subtract per-axis minimum, then offset by a small positive margin
    structure.coordinates = xyz - xyz.min(axis=0) + margin


PBSA_TEMPLATE = """ titrate pb_backend
 &cntrl
   ipb={ipb}, inp={inp}
 /
 &pb
   npbverb=0, istrng={istrng}, ivalence=1, iprob={iprob}, dprob={dprob},
   epsout={epsout}, epsin={epsin},
   radiopt={radiopt}, smoothopt=1, sasopt=0,
   bcopt={bcopt}, eneopt={eneopt},
   nfocus={nfocus}, fscale=8,
   space={space}, fillratio={fillratio},
   accept={accept}, maxitn={maxitn},
   solvopt=1{phiout_block}
 /
"""


def run_pbsa(structure: pmd.Structure,
             work_dir: Path,
             label: str,
             *,
             # PB parameters with defaults that match common usage.
             # Override per-call when needed.
             epsin:    float = 4.0,
             epsout:   float = 80.0,
             istrng:   float = 150.0,    # mM
             space:    float = 0.5,      # Å
             nfocus:   int   = 2,
             fillratio:float = 2.0,
             bcopt:    int   = 5,        # 5 = Debye-Hückel boundary potentials from grid charges (pbsa default)
             eneopt:   int   = 2,        # 2 = report EELEC + EPB separately
             ipb:      int   = 2,        # 2 = level-set dielectric interface (pbsa default); npbopt=0 (linear PB) is left as pbsa default
             inp:      int   = 0,        # 0 = no nonpolar term
             iprob:    float = 2.0,      # ion exclusion radius (Å)
             dprob:    float = 1.4,      # solvent probe radius (Å)
             radiopt:  int   = 0,        # 0 = use radii from prmtop
             accept:   float = 1e-3,
             maxitn:   int   = 10000,
             phiout:   int   = 0,        # 1 = write phi map ('pbsa.phi')
             phiform:  int   = 1,        # 1 = ASCII (use with read_pbsa_phi)
             reuse:    bool  = True,     # reuse an existing parseable {label}.out
             ) -> Dict[str, float]:
    """Run pbsa on `structure` and return parsed energies (kcal/mol).

    Files written under `work_dir`:
        {label}.prmtop, {label}.inpcrd, {label}.in, {label}.out
    pbsa is invoked with cwd=work_dir so its side-effect files (mdinfo, etc.)
    land there, not in the caller's CWD.

    Resume: when ``reuse`` is set (default) and a parseable ``{label}.out``
    already exists, its energies are returned without re-running pbsa. The PB
    result is deterministic for fixed inputs, so this lets an interrupted pbt-3
    / pbt-4 run skip every leg it already completed (and salvages legs where
    pbsa wrote a full result but crashed on exit-time cleanup). Callers that
    change PB parameters must clear the work dir (or pass ``reuse=False``) so
    stale outputs aren't picked up.
    """
    out_path = work_dir / f"{label}.out"
    if reuse and not phiout and out_path.exists():
        try:
            return parse_pbsa_output(out_path)
        except (RuntimeError, OSError):
            pass  # absent/partial/corrupt result: fall through and (re)run

    work_dir.mkdir(parents=True, exist_ok=True)
    # Drop any periodic box: a PB pKa is dilute-solution continuum
    # electrostatics -- a single solute in an ionic continuum with the
    # potential decaying to zero far away -- which is NON-PERIODIC by design.
    # This is a deliberate methodological choice, not just a crash fix:
    #   * The model compound (a capped tripeptide, no box) is computed open/
    #     non-periodic, so the cluster leg of the ddG cycle must match it or
    #     the boundary condition fails to cancel in the subtraction.
    #   * Periodic BC would add spurious protein/self-image coupling absent in
    #     dilute solution, and is ill-defined for the net-charged per-state
    #     legs (needs a neutralizing-background finite-size correction).
    #   * The focused FD-PB solver (ipb>=1, bcopt=5 Debye-Huckel far field)
    #     never reads the prmtop box anyway -- pbsa builds its own grid -- so
    #     dropping it is result-neutral (verified: bit-identical EPB/EELEC).
    # Do NOT "restore" the box thinking it's a bug. If a crystal/membrane mode
    # is ever added it must be opt-in and carry its own finite-size handling.
    # Stripping it also dodges a real pbsa bug: a cluster sliced from a
    # solvated system keeps IFBOX=1 with many disconnected molecules (protein
    # fragment + ions + cofactors), but pbsa's locmem allocates only ONE slot
    # for the per-molecule ATOMS_PER_MOLECULE array (rdparm.F90 `lasti=i70`),
    # so rdparm2 overruns the integer array by (nspm-1) elements and crashes
    # in a formatted read (SIGSEGV in formatted_transfer_scalar_read) when the
    # overrun lands on live heap metadata. IFBOX=0 makes pbsa skip that read.
    if structure.box is not None:
        structure.box = None
    if not phiout:
        # Guarantee coordinates fit the inpcrd F12.7 fields (see helper).
        # Skipped on the phiout path so the grid frame stays aligned with the
        # original-frame coordinates that coupling.py interpolates against.
        _recenter_for_format(structure)
    structure.save(str(work_dir / f"{label}.prmtop"), overwrite=True)
    structure.save(str(work_dir / f"{label}.inpcrd"), overwrite=True)
    phiout_block = (f", phiout={phiout}, phiform={phiform}"
                     if phiout else "")
    (work_dir / f"{label}.in").write_text(PBSA_TEMPLATE.format(
        epsin=epsin, epsout=epsout, istrng=istrng, space=space,
        nfocus=nfocus, fillratio=fillratio, bcopt=bcopt, eneopt=eneopt,
        ipb=ipb, inp=inp, iprob=iprob, dprob=dprob, radiopt=radiopt,
        accept=accept, maxitn=maxitn, phiout_block=phiout_block,
    ))
    cmd = ["pbsa", "-O", "-i", f"{label}.in", "-o", f"{label}.out",
           "-p", f"{label}.prmtop", "-c", f"{label}.inpcrd"]
    # pbsa on some clusters (e.g. TYR-120 in the NOS calmodulin system) hits a
    # non-deterministic heap bug: for the SAME input it may exit 0, or compute a
    # complete FINAL RESULTS block and THEN abort during exit-time cleanup
    # (glibc "free(): invalid size" -> -6, or a SIGSEGV in teardown -> -11; the
    # signal varies run to run with memory layout). The computation itself is
    # deterministic (identical EPB whenever it completes). So:
    #   * if pbsa wrote a complete result, use it even on a nonzero exit (the
    #     crash was after the answer was written), and
    #   * otherwise re-run a few times — a fresh attempt usually completes.
    max_attempts = 4
    energies = None
    last_rc, last_err = 0, ""
    for attempt in range(1, max_attempts + 1):
        r = subprocess.run(cmd, capture_output=True, text=True,
                           cwd=str(work_dir))
        try:
            energies = parse_pbsa_output(out_path)
        except (RuntimeError, OSError):
            energies = None  # no/incomplete output on this attempt
        if energies is not None:
            if r.returncode != 0:
                _record_pbsa_salvage(label, r.returncode, energies.get('EPB'))
            break
        last_rc, last_err = r.returncode, (r.stderr or r.stdout or "")
    if energies is None:
        raise RuntimeError(
            f"pbsa produced no usable result for {label} after "
            f"{max_attempts} attempts (last exit {last_rc}):\n"
            f"{last_err[-2000:]}")
    if phiout:
        # pbsa writes phi to 'pbsa_phi.phi' (phiform=1) or 'pbsa_phi.dx'
        # (phiform=2). Rename to per-label so concurrent calls don't collide.
        ext = "phi" if phiform == 1 else "dx"
        default_phi = work_dir / f"pbsa_phi.{ext}"
        labeled_phi = work_dir / f"{label}.{ext}"
        if default_phi.exists():
            default_phi.rename(labeled_phi)
            energies["phi_path"] = str(labeled_phi)
    return energies


def parse_pbsa_output(path: Path) -> Dict[str, float]:
    """Extract energy fields from a pbsa output file.

    Uses the FINAL RESULTS block when present, else the last NSTEP block.
    Raises if EPB cannot be parsed (a strong signal pbsa silently failed).
    """
    text = path.read_text()
    finals = re.split(r"FINAL RESULTS\s*\n", text)
    block = finals[-1] if len(finals) > 1 else text
    energies: Dict[str, float] = {}
    patterns = {
        "Etot":    r"Etot\s*=\s*(-?\d+\.\d+)",
        "EPB":     r"EPB\s*=\s*(-?\d+\.\d+)",
        "EELEC":   r"EELEC\s*=\s*(-?\d+\.\d+)",
        "EEL14":   r"1-4 EEL\s*=\s*(-?\d+\.\d+)",
        "VDWAALS": r"VDWAALS\s*=\s*(-?\d+\.\d+)",
        "ECAVITY": r"ECAVITY\s*=\s*(-?\d+\.\d+)",
        "EDISPER": r"EDISPER\s*=\s*(-?\d+\.\d+)",
    }
    for key, pat in patterns.items():
        m = re.search(pat, block)
        if m:
            energies[key] = float(m.group(1))
    if "EPB" not in energies:
        raise RuntimeError(f"Could not parse EPB from {path}")
    return energies


def total_electrostatic(e: Dict[str, float]) -> float:
    """EPB + EELEC. The two electrostatic components from `eneopt=2`.

    EELEC includes long-range Coulomb (between point charges in vacuum).
    EPB is the reaction-field energy (the dielectric/ionic environment's
    response). Their sum is the total electrostatic energy under the
    PB picture. Other Etot components (vdW, internal) are independent
    of charge state and cancel in any state-difference calculation.
    """
    return e.get("EPB", 0.0) + e.get("EELEC", 0.0)
