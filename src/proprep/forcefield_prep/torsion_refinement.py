"""
Shared torsion (dihedral) refinement engine.

One engine serves every parameterizer that can refine a dihedral against a
relaxed QM scan: the small-molecule route (checklist step sm-7) and both
modified-amino-acid routes (step 9). Each route decides WHICH dihedrals to
scan, WHERE the scan files live and WHAT template Gaussian input to derive a
scan from; the engine does the rest:

  * derive a relaxed-scan input from an existing Gaussian input (same atoms,
    same order, same restraints, one added ``D i j k l S n step`` line);
  * collect finished scan logs into geometry/energy datasets;
  * pool every dataset into ONE mdcrd + ONE energy file and run ONE paramfit
    fit with every selected dihedral free at once (a joint fit), against a
    prmtop the caller built from the frcmod as it stands NOW;
  * splice every fitted DIHE term back into that frcmod, so the result is a
    complete parameter file, never paramfit's fitted-terms-only output.

Why joint: a relaxed scan of one torsion lets its neighbours drift, so the
energy along it depends on the other torsions' parameters too. Fitting the
selected torsions together against all scans at once accounts for that;
fitting them one at a time against stale parameters does not.

Why one K: paramfit fits absolute energies with a single QM-MM offset. Scans
of the same molecule at the same level of theory share one QM zero, so one K
aligns them all. Scans at different levels do NOT, which is why
:func:`check_compatible` refuses to pool them.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
from dataclasses import dataclass, field
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

from rich.console import Console
from rich.panel import Panel
from rich.table import Table

HARTREE_TO_KCAL = 627.5095

# Route-line tokens that belong to an optimization / ESP job, not to a scan.
_STRIP_ROUTE_TOKENS = (
    re.compile(r"^freq(\b|=|\()", re.IGNORECASE),
    re.compile(r"^iop\(", re.IGNORECASE),
    re.compile(r"^pop=", re.IGNORECASE),
)


# =============================================================================
# Data model
# =============================================================================

@dataclass
class TorsionScan:
    """One relaxed scan of one dihedral, in one Gaussian atom ordering.

    ``atom_indices`` are 1-based indices into the scan input's atom list.
    ``quad`` is the four AMBER/GAFF atom types in that same order (None until
    resolved). ``geometries``/``energies`` are filled by :func:`collect_scan`.
    """
    label: str
    gjf: str
    log: str
    atom_indices: Tuple[int, int, int, int]
    quad: Optional[Tuple[str, str, str, str]] = None
    source: str = ""
    n_steps: Optional[int] = None
    step_size: Optional[float] = None
    route: Optional[str] = None
    n_atoms: int = 0
    geometries: List[list] = field(default_factory=list)
    energies: List[float] = field(default_factory=list)

    @property
    def log_exists(self) -> bool:
        return bool(self.log) and os.path.exists(self.log)

    @property
    def ready(self) -> bool:
        return bool(self.geometries) and len(self.geometries) == len(self.energies)

    @property
    def quad_name(self) -> str:
        return quad_name(self.quad) if self.quad else "-".join(str(i) for i in self.atom_indices)


# =============================================================================
# Type-quad helpers
# =============================================================================

def quad_from_name(name: str) -> Optional[Tuple[str, ...]]:
    """``"c3-c3-os-c "`` → ``('c3','c3','os','c')``; None unless exactly four."""
    parts = [p.strip() for p in name.replace(" ", "").split("-")]
    parts = [p for p in parts if p]
    return tuple(parts) if len(parts) == 4 else None


def quad_name(quad: Sequence[str]) -> str:
    return "-".join(quad)


def quad_key(quad: Sequence[str]) -> Tuple[str, ...]:
    """Orientation-independent identity: a-b-c-d and d-c-b-a are one dihedral."""
    fwd, rev = tuple(quad), tuple(reversed(quad))
    return min(fwd, rev)


def quad_matches(quad: Sequence[str], other: Sequence[str]) -> bool:
    return quad_key(quad) == quad_key(other)


def frcmod_dihe_types(line: str) -> Optional[Tuple[str, ...]]:
    """The four types on a DIHE line, read from its fixed 11-character field.

    AMBER's format is A2,1X,A2,1X,A2,1X,A2, so ``line[0:11]`` splits on '-'
    into four 2-character fields. Comparing tuples of types (not the reversed
    character string) is what makes the reversed orientation match for
    two-character types like ``c3``.
    """
    if len(line) < 11:
        return None
    parts = [p.strip() for p in line[:11].split("-")]
    return tuple(parts) if len(parts) == 4 and all(parts) else None


# =============================================================================
# MOL2 ↔ scan mapping
# =============================================================================

def _mol2_types(mol2_file: str) -> List[str]:
    from proprep.forcefield_prep.seminario_refinement import parse_mol2_connectivity
    return list(parse_mol2_connectivity(mol2_file)["atom_types"])


def types_for_indices(mol2_file: str, idxs: Sequence[int]) -> Optional[Tuple[str, ...]]:
    """Atom types (mol2 order, 1-based indices) for a scanned quad, or None."""
    try:
        types = _mol2_types(mol2_file)
    except Exception:
        return None
    if not all(1 <= i <= len(types) for i in idxs):
        return None
    return tuple(types[i - 1] for i in idxs)


def atoms_for_quad(quad: Sequence[str], mol2_file: str) -> Optional[Tuple[int, int, int, int]]:
    """First atom quadruple (1-based) whose types match ``quad`` in either orientation.

    A type quad names a CLASS of dihedrals; the scan drives one instance of it.
    The fitted term applies to every instance, because AMBER dihedrals are
    typed, not atom-indexed.
    """
    from proprep.forcefield_prep.pes_scan_refinement import map_dihedral_to_atoms
    quiet = Console(quiet=True)
    hit = map_dihedral_to_atoms((quad_name(quad), 0.0, "", "DIHE"), mol2_file, quiet)
    if hit:
        return hit
    return map_dihedral_to_atoms((quad_name(tuple(reversed(quad))), 0.0, "", "DIHE"), mol2_file, quiet)


# =============================================================================
# Gaussian input handling
# =============================================================================

def read_scan_line(gjf: str) -> Optional[Tuple[Tuple[int, int, int, int], int, float]]:
    """``(indices, n_steps, step_size)`` from a ``D i j k l S n step`` line, or None."""
    try:
        with open(gjf) as f:
            for line in f:
                m = re.match(r"\s*D\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+S\s+(\d+)\s+([-\d.]+)", line)
                if m:
                    idxs = tuple(int(m.group(k)) for k in range(1, 5))
                    return idxs, int(m.group(5)), float(m.group(6))
    except OSError:
        return None
    return None


def read_route_line(path: str) -> Optional[str]:
    """The normalized route line of a Gaussian input or log.

    Inputs carry it as ``#p ...``; logs echo it after a dashed line, possibly
    wrapped. Normalized to lowercase with collapsed whitespace so two scans
    can be compared for the pooling check. None if not found.
    """
    try:
        with open(path, errors="replace") as f:
            lines = f.readlines()
    except OSError:
        return None
    route: List[str] = []
    for i, line in enumerate(lines[:400]):
        s = line.strip()
        if s.startswith("#"):
            route.append(s)
            # Logs wrap long routes onto following lines until a dashed line.
            for nxt in lines[i + 1:i + 6]:
                t = nxt.strip()
                if not t or t.startswith("-"):
                    break
                route.append(t)
            break
    if not route:
        return None
    text = " ".join(route)
    text = re.sub(r"^#\w*\s*", "", text)
    return " ".join(text.lower().split())


def _scan_route(route: str) -> str:
    """Turn an opt/freq/ESP route into a relaxed-scan route."""
    tokens = route.split()
    out: List[str] = []
    saw_opt = False
    for tok in tokens:
        if tok.startswith("#"):
            continue
        if any(p.match(tok) for p in _STRIP_ROUTE_TOKENS):
            continue
        low = tok.lower()
        if low.startswith("opt"):
            saw_opt = True
            if "modredundant" not in low:
                tok = "Opt=ModRedundant"
        out.append(tok)
    if not saw_opt:
        out.insert(0, "Opt=ModRedundant")
    if not any(t.lower().startswith("nosym") for t in out):
        out.append("NoSymm")
    return " ".join(out)


def _split_gjf(gjf: str) -> Dict[str, Any]:
    """Break a Gaussian input into link0 / route / title / charge / atoms / tail."""
    with open(gjf) as f:
        raw = f.read().splitlines()
    i = 0
    link0: List[str] = []
    while i < len(raw) and raw[i].startswith("%"):
        link0.append(raw[i])
        i += 1
    route: List[str] = []
    while i < len(raw) and raw[i].strip():
        route.append(raw[i].strip())
        i += 1
    while i < len(raw) and not raw[i].strip():
        i += 1
    title: List[str] = []
    while i < len(raw) and raw[i].strip():
        title.append(raw[i])
        i += 1
    while i < len(raw) and not raw[i].strip():
        i += 1
    charge_line = raw[i].strip() if i < len(raw) else "0 1"
    i += 1
    atoms: List[str] = []
    while i < len(raw) and raw[i].strip():
        atoms.append(raw[i])
        i += 1
    while i < len(raw) and not raw[i].strip():
        i += 1
    tail: List[str] = []
    while i < len(raw) and raw[i].strip():
        tail.append(raw[i].strip())
        i += 1
    return {"link0": link0, "route": " ".join(route), "title": title,
            "charge_line": charge_line, "atoms": atoms, "tail": tail}


def _atom_fields(line: str) -> Tuple[str, Optional[str], Tuple[float, float, float]]:
    """``El [flag] x y z`` → (element, flag-or-None, coords)."""
    parts = line.split()
    if len(parts) >= 5 and re.fullmatch(r"-?\d+", parts[1]):
        return parts[0], parts[1], (float(parts[2]), float(parts[3]), float(parts[4]))
    return parts[0], None, (float(parts[1]), float(parts[2]), float(parts[3]))


def gjf_elements(gjf: str) -> List[str]:
    return [_atom_fields(a)[0] for a in _split_gjf(gjf)["atoms"]]


def last_orientation(log: str, n_atoms: int) -> Optional[List[Tuple[float, float, float]]]:
    """Last printed geometry block (Input or Standard orientation) of a log."""
    try:
        with open(log, errors="replace") as f:
            lines = f.readlines()
    except OSError:
        return None
    last = None
    i = 0
    while i < len(lines):
        if "orientation:" in lines[i]:
            block: List[Tuple[float, float, float]] = []
            j = i + 5
            while j < len(lines) and not re.match(r"\s*-{40,}", lines[j]):
                p = lines[j].split()
                if len(p) >= 6:
                    try:
                        block.append((float(p[3]), float(p[4]), float(p[5])))
                    except ValueError:
                        pass
                j += 1
            if len(block) == n_atoms:
                last = block
            i = j
        i += 1
    return last


def derive_scan_input(template_gjf: str, output_gjf: str,
                      idxs: Sequence[int], n_steps: int, step_size: float, *,
                      coords: Optional[Sequence[Tuple[float, float, float]]] = None,
                      chk: Optional[str] = None,
                      title: Optional[str] = None) -> str:
    """Write a relaxed-scan input derived from an existing Gaussian input.

    Everything that defines the model is kept verbatim: atom symbols and
    order, freeze flags, charge and multiplicity, memory and processors, and
    any ``D ... F`` restraint line (except one on the scanned torsion itself,
    which cannot be both frozen and driven). The route loses ``Freq``,
    ``IOp(...)`` and ``Pop=...`` and gains ``Opt=ModRedundant`` and ``NoSymm``.
    Any previous ``S`` scan line is dropped and the new one appended.

    ``coords`` replaces the template coordinates (same atom count required),
    typically the optimized geometry read back from the template's own log.
    """
    parts = _split_gjf(template_gjf)
    atoms = parts["atoms"]
    if coords is not None and len(coords) != len(atoms):
        raise ValueError(
            f"{template_gjf}: {len(atoms)} atoms but {len(coords)} coordinates supplied")
    quad = tuple(int(i) for i in idxs)
    if not all(1 <= i <= len(atoms) for i in quad):
        raise ValueError(f"scan indices {quad} outside 1..{len(atoms)}")

    link0 = [ln for ln in parts["link0"] if not ln.lower().startswith("%chk")]
    if chk:
        link0.insert(0, f"%chk={chk}")
    route = _scan_route(parts["route"])

    atom_lines: List[str] = []
    for k, line in enumerate(atoms):
        el, flag, xyz = _atom_fields(line)
        if coords is not None:
            xyz = tuple(coords[k])
        if flag is not None:
            atom_lines.append(f"{el:<2s} {int(flag):2d}  {xyz[0]:12.6f} {xyz[1]:12.6f} {xyz[2]:12.6f}")
        else:
            atom_lines.append(f"{el:<2s}    {xyz[0]:12.6f} {xyz[1]:12.6f} {xyz[2]:12.6f}")

    kept_tail: List[str] = []
    for ln in parts["tail"]:
        m = re.match(r"\s*D\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+([FS])\b", ln)
        if m:
            old = tuple(int(m.group(k)) for k in range(1, 5))
            if m.group(5) == "S":
                continue
            if old == quad or old == tuple(reversed(quad)):
                continue
        kept_tail.append(ln)
    kept_tail.append("D " + " ".join(str(i) for i in quad)
                     + f" S {int(n_steps)} {float(step_size):.1f}")

    title_lines = [title] if title else (parts["title"] or ["relaxed dihedral scan"])
    with open(output_gjf, "w") as f:
        for ln in link0:
            f.write(ln + "\n")
        f.write(f"#p {route}\n\n")
        for ln in title_lines:
            f.write(ln + "\n")
        f.write("\n")
        f.write(parts["charge_line"] + "\n")
        for ln in atom_lines:
            f.write(ln + "\n")
        f.write("\n")
        for ln in kept_tail:
            f.write(ln + "\n")
        f.write("\n")
    return output_gjf


# =============================================================================
# Collecting scans
# =============================================================================

def scan_from_gjf(gjf: str, log: str, label: str, mol2_file: Optional[str] = None,
                  source: str = "") -> Optional[TorsionScan]:
    """A TorsionScan described by an existing scan input (indices, steps, route)."""
    sl = read_scan_line(gjf)
    if not sl:
        return None
    idxs, n_steps, step = sl
    quad = types_for_indices(mol2_file, idxs) if mol2_file else None
    return TorsionScan(label=label, gjf=gjf, log=log, atom_indices=idxs, quad=quad,
                       source=source, n_steps=n_steps, step_size=step,
                       route=read_route_line(gjf))


def collect_scan(scan: TorsionScan, console: Console, min_points: int = 5) -> bool:
    """Parse the scan's log into geometries/energies. False if unusable."""
    from proprep.forcefield_prep.pes_scan_refinement import parse_pes_scan_log
    if not scan.log_exists:
        return False
    data = parse_pes_scan_log(scan.log, console, step_size=scan.step_size or 15.0)
    if not data.get("success") or len(data.get("geometries") or []) < min_points:
        n = len(data.get("geometries") or [])
        console.print(f"[yellow]{scan.log}: {n} usable scan point(s); need at least {min_points}[/yellow]")
        return False
    scan.geometries = data["geometries"]
    scan.energies = data["energies"][:len(scan.geometries)]
    scan.n_atoms = data["n_atoms"]
    if scan.route is None:
        scan.route = read_route_line(scan.log)
    return True


def check_compatible(scans: Sequence[TorsionScan]) -> Tuple[bool, str]:
    """Every pooled scan must share one atom count and one level of theory."""
    ready = [s for s in scans if s.ready]
    if not ready:
        return False, "no scan has usable data"
    n_atoms = {s.n_atoms for s in ready}
    if len(n_atoms) > 1:
        return False, f"scans have different atom counts: {sorted(n_atoms)}"
    routes = {}
    for s in ready:
        r = s.route or read_route_line(s.log) or read_route_line(s.gjf)
        routes.setdefault(_method_signature(r), []).append(s.label)
    if len(routes) > 1:
        detail = "; ".join(f"{k or 'unknown'}: {', '.join(v)}" for k, v in routes.items())
        return False, ("scans were run at different levels of theory and cannot share "
                       f"one energy offset ({detail})")
    return True, ""


def _method_signature(route: Optional[str]) -> Optional[str]:
    """The route minus job-control tokens: what decides the QM energy zero."""
    if not route:
        return None
    keep = []
    for tok in route.split():
        low = tok.lower()
        if low.startswith(("opt", "nosym", "geom=", "integral=", "int=", "scf=", "guess=")):
            continue
        keep.append(low)
    return " ".join(sorted(keep))


# =============================================================================
# paramfit inputs, fit, merge
# =============================================================================

def write_joint_inputs(scans: Sequence[TorsionScan], mol_name: str, work_dir: str,
                       console: Console) -> Tuple[str, str, int, List[Tuple[str, int, int]]]:
    """One mdcrd and one energy file for every ready scan, in order.

    Returns (mdcrd, energies_file, n_structures, boundaries) where boundaries
    lists ``(label, first_frame, last_frame_exclusive)`` per scan.
    """
    os.makedirs(work_dir, exist_ok=True)
    mdcrd = os.path.join(work_dir, f"{mol_name}_torsions.mdcrd")
    energies = os.path.join(work_dir, f"{mol_name}_torsions_qm.dat")
    boundaries: List[Tuple[str, int, int]] = []
    n = 0
    with open(mdcrd, "w") as fm, open(energies, "w") as fe:
        fm.write(f"{mol_name} pooled relaxed scans for paramfit\n")
        for s in scans:
            if not s.ready:
                continue
            start = n
            for geom, e in zip(s.geometries, s.energies):
                coords: List[float] = []
                for atom in geom:
                    coords.extend([atom[1], atom[2], atom[3]])
                for j in range(0, len(coords), 10):
                    fm.write("".join(f"{c:8.3f}" for c in coords[j:j + 10]) + "\n")
                fe.write(f"{e}\n")
                n += 1
            boundaries.append((s.label, start, n))
    console.print(f"[green]✓ Pooled {n} scan point(s) from {len(boundaries)} scan(s) → {mdcrd}[/green]")
    return mdcrd, energies, n, boundaries


def merge_fitted_dihedrals(original_frcmod: str, fitted_frcmod: str,
                           quads: Iterable[Sequence[str]], output_file: str,
                           console: Console) -> str:
    """Splice paramfit's fitted DIHE terms into a complete frcmod.

    paramfit's ``WRITE_FRCMOD`` emits ONLY the terms it fitted (verified in
    AmberTools ``paramfit/write_input.c``), so its output can never be loaded
    on its own. Every DIHE line of the original whose four types match one of
    ``quads`` in either orientation is dropped, and every DIHE line paramfit
    wrote is inserted in its place. Multi-term dihedrals are replaced whole.
    """
    keys = {quad_key(q) for q in quads}
    fitted: List[str] = []
    if os.path.exists(fitted_frcmod):
        with open(fitted_frcmod) as f:
            in_dihe = False
            for line in f:
                s = line.strip()
                if s == "DIHE":
                    in_dihe = True
                    continue
                if s in ("BOND", "ANGLE", "ANGL", "IMPROPER", "NONBON", "END", ""):
                    in_dihe = False
                    continue
                if in_dihe:
                    t = frcmod_dihe_types(line)
                    if t and quad_key(t) in keys:
                        fitted.append(line if line.endswith("\n") else line + "\n")
    if not fitted:
        console.print(f"[yellow]No fitted dihedral lines found in {fitted_frcmod}; frcmod unchanged[/yellow]")
        shutil.copy(original_frcmod, output_file)
        return output_file

    with open(original_frcmod) as f:
        original = f.readlines()
    out: List[str] = []
    in_dihe = False
    replaced = 0
    has_dihe = any(ln.strip() == "DIHE" for ln in original)
    for line in original:
        s = line.strip()
        if s == "DIHE":
            in_dihe = True
            out.append(line)
            continue
        if in_dihe and s in ("IMPROPER", "NONBON", "END", "BOND", "ANGLE"):
            out.extend(fitted)
            in_dihe = False
            out.append(line)
            continue
        if in_dihe and s:
            t = frcmod_dihe_types(line)
            if t and quad_key(t) in keys:
                replaced += 1
                continue
        out.append(line)
    if in_dihe:  # DIHE was the last section
        out.extend(fitted)
    if not has_dihe:
        # Insert a DIHE section before IMPROPER/NONBON/END, or at the end.
        idx = next((i for i, ln in enumerate(out)
                    if ln.strip() in ("IMPROPER", "NONBON", "END")), len(out))
        out[idx:idx] = ["DIHE\n"] + fitted + ["\n"]
    with open(output_file, "w") as f:
        f.writelines(out)
    console.print(f"[green]✓ Merged {len(fitted)} fitted dihedral line(s) "
                  f"(replacing {replaced}) → {output_file}[/green]")
    return output_file


def per_scan_residuals(energy_output: str,
                       boundaries: Sequence[Tuple[str, int, int]]) -> List[Tuple[str, float, float]]:
    """RMS and max |MM+K − QM| per scan from paramfit's WRITE_ENERGY file."""
    rows: List[float] = []
    try:
        with open(energy_output) as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                p = line.split()
                if len(p) >= 3:
                    try:
                        rows.append(float(p[1]) - float(p[2]))
                    except ValueError:
                        continue
    except OSError:
        return []
    out: List[Tuple[str, float, float]] = []
    for label, a, b in boundaries:
        seg = rows[a:b]
        if not seg:
            continue
        rms = (sum(r * r for r in seg) / len(seg)) ** 0.5
        out.append((label, rms, max(abs(r) for r in seg)))
    return out


def run_joint_torsion_fit(scans: Sequence[TorsionScan], *, mol_name: str,
                          prmtop: str, frcmod: str, work_dir: str, console: Console,
                          fit_equilibrium: bool = True) -> Dict[str, Any]:
    """Fit every ready scan's dihedral together; return the merged frcmod.

    ``prmtop`` must have been built from ``frcmod`` (the CURRENT parameters,
    after any Seminario refinement): paramfit reads the MM side from it, so
    every term not being fitted is held at whatever that file says.
    """
    from proprep.forcefield_prep.paramfit_refinement import (
        run_paramfit_k_fitting, run_paramfit_parameter_fitting,
        run_paramfit_set_params_automated,
    )
    result: Dict[str, Any] = {"refinement_success": False, "fitted_frcmod": None,
                              "merged_frcmod": None, "message": "", "k_value": None,
                              "residuals": [], "quads": []}
    ready = [s for s in scans if s.ready and s.quad]
    if not ready:
        result["message"] = "no scan with both data and resolved atom types"
        return result
    ok, why = check_compatible(ready)
    if not ok:
        result["message"] = why
        console.print(f"[red]Cannot pool scans: {why}[/red]")
        return result

    # De-duplicate dihedrals: several scans may drive the same type quad
    # (e.g. one torsion in several backbone conformers) — one fitted term.
    quads: List[Tuple[str, ...]] = []
    for s in ready:
        if not any(quad_matches(s.quad, q) for q in quads):
            quads.append(tuple(s.quad))
    result["quads"] = [quad_name(q) for q in quads]

    work_dir = os.path.abspath(work_dir)
    prmtop = os.path.abspath(prmtop)
    frcmod = os.path.abspath(frcmod)
    mdcrd, energies, n_structures, boundaries = write_joint_inputs(ready, mol_name, work_dir, console)
    console.print(Panel(
        "[bold]Joint torsion fit[/bold]\n\n"
        f"Dihedral(s) fitted together: [cyan]{', '.join(result['quads'])}[/cyan]\n"
        f"Scan points pooled: {n_structures} from {len(boundaries)} scan(s)\n"
        f"MM reference topology: {os.path.basename(prmtop)} (built from the current frcmod)\n"
        f"Free per term: {'barrier and phase (periodicity kept)' if fit_equilibrium else 'barrier only'}",
        title="Torsion Refinement", border_style="cyan", expand=False))

    original_dir = os.getcwd()
    os.chdir(work_dir)
    try:
        k = run_paramfit_k_fitting(prmtop, mdcrd, energies, n_structures, console)
        if k is None:
            result["message"] = "K (energy offset) fit failed"
            return result
        result["k_value"] = k
        params_file = f"{mol_name}_torsions.params"
        # Both orientations, so paramfit's prompt matches whichever way it
        # prints the quad.
        selected = []
        for q in quads:
            selected.append((quad_name(q), float("inf"), "SCAN", "DIHE"))
            selected.append((quad_name(tuple(reversed(q))), float("inf"), "SCAN", "DIHE"))
        if not run_paramfit_set_params_automated(prmtop, params_file, selected, console,
                                                 force_constants_only=not fit_equilibrium):
            result["message"] = "paramfit parameter selection failed"
            return result
        fitted = f"{mol_name}_torsions_fitted.frcmod"
        if not run_paramfit_parameter_fitting(prmtop, mdcrd, energies, params_file,
                                              n_structures, k, fitted, console):
            result["message"] = "paramfit fit failed"
            return result
        result["fitted_frcmod"] = os.path.abspath(fitted)
        merged = os.path.join(work_dir, f"{mol_name}_torsions_merged.frcmod")
        merge_fitted_dihedrals(frcmod, fitted, quads, merged, console)
        result["merged_frcmod"] = merged
        result["residuals"] = per_scan_residuals(fitted.replace(".frcmod", "_energies.dat"), boundaries)
        if result["residuals"]:
            table = Table(title="Fit residuals per scan (MM+K − QM, kcal/mol)", expand=False)
            table.add_column("Scan", style="cyan")
            table.add_column("RMS", justify="right")
            table.add_column("Max", justify="right")
            for label, rms, mx in result["residuals"]:
                table.add_row(label, f"{rms:.2f}", f"{mx:.2f}")
            console.print(table)
        result["refinement_success"] = True
        result["message"] = f"Refined {len(quads)} dihedral(s) from {len(boundaries)} scan(s)"
        return result
    finally:
        os.chdir(original_dir)


# =============================================================================
# Presentation helpers shared by the routes
# =============================================================================

def dihedral_penalty_rows(penalties: Sequence[Tuple[str, float, str, str]]) -> List[Tuple[int, Tuple[str, ...], float, str]]:
    """``(table_number, quad, score, status)`` for every DIHE row with a real quad.

    ``table_number`` is the 1-based row shown by ``analyze_frcmod_penalties``,
    so the user can pick with the numbers already on screen.
    """
    rows = []
    for i, (name, score, status, section) in enumerate(penalties, 1):
        if section != "DIHE":
            continue
        q = quad_from_name(name)
        if q and "X" not in q:
            rows.append((i, q, score, status))
    return rows


def parse_selection(text: str, valid_numbers: Sequence[int]) -> Optional[List[int]]:
    """``'1,3-5'`` / ``'all'`` / ``'none'`` → table numbers; None on a parse error."""
    s = text.strip().lower()
    if s in ("", "none", "n"):
        return []
    if s in ("all", "a"):
        return list(valid_numbers)
    picked: List[int] = []
    try:
        for part in s.replace(" ", "").split(","):
            if not part:
                continue
            if "-" in part:
                a, b = part.split("-", 1)
                picked.extend(range(int(a), int(b) + 1))
            else:
                picked.append(int(part))
    except ValueError:
        return None
    return sorted({p for p in picked if p in valid_numbers})


def print_scan_action_panel(pending: Sequence[TorsionScan], console: Console,
                            resume_hint: str) -> None:
    lines = [f"  • {s.quad_name}: {s.gjf}" for s in pending]
    cmds = [f"  g16 {s.gjf} > {s.log}" for s in pending]
    console.print(Panel(
        f"[bold yellow]Gaussian relaxed scan(s) needed for {len(pending)} dihedral(s):[/bold yellow]\n\n"
        + "\n".join(lines)
        + "\n\n[bold]To run:[/bold]\n" + "\n".join(cmds)
        + f"\n\n[grey50]{resume_hint}[/grey50]",
        title="Action Required", border_style="yellow", expand=False))


def build_mol2_prmtop(mol2_file: str, frcmod_file: str, out_prefix: str,
                      console: Console,
                      sources: Sequence[str] = ("leaprc.gaff2",)) -> Optional[str]:
    """A throwaway prmtop of one mol2 (loaded as-is) with the given frcmod.

    ``loadmol2`` keeps the mol2's atom order, so when the mol2 was made from
    the same Gaussian atom list the scans use (antechamber on that output),
    the prmtop and the scan frames line up atom for atom with no renaming.
    The frcmod is the CURRENT one (after Seminario), so paramfit's MM side
    carries every refinement made so far.
    """
    if not shutil.which("tleap"):
        console.print("[yellow]tleap not found; cannot build the fitting topology[/yellow]")
        return None
    out_dir = os.path.dirname(os.path.abspath(out_prefix)) or "."
    os.makedirs(out_dir, exist_ok=True)
    prmtop = f"{out_prefix}.parm7"
    inpcrd = f"{out_prefix}.rst7"
    script = f"{out_prefix}_tleap.in"
    log = f"{out_prefix}_tleap.log"
    for stale in (prmtop, inpcrd):
        if os.path.exists(stale):
            os.remove(stale)
    with open(script, "w") as f:
        for src in sources:
            f.write(f"source {src}\n")
        f.write(f'loadamberparams "{os.path.abspath(frcmod_file)}"\n')
        f.write(f'm = loadmol2 "{os.path.abspath(mol2_file)}"\n')
        f.write(f'saveamberparm m "{prmtop}" "{inpcrd}"\nquit\n')
    try:
        run = subprocess.run(["tleap", "-f", script], check=False,
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout=300)
        with open(log, "w") as f:
            f.write(run.stdout or "")
            f.write(run.stderr or "")
    except (OSError, subprocess.TimeoutExpired) as e:
        console.print(f"[yellow]tleap failed while building the fitting topology: {e}[/yellow]")
        return None
    if not os.path.exists(prmtop):
        console.print(f"[yellow]tleap wrote no prmtop (see {log}); the frcmod may be missing terms[/yellow]")
        return None
    return prmtop


def build_gaff2_prmtop(mol2_file: str, frcmod_file: str, out_prefix: str,
                       console: Console) -> Optional[str]:
    """Small-molecule variant of :func:`build_mol2_prmtop` (GAFF2 only)."""
    return build_mol2_prmtop(mol2_file, frcmod_file, out_prefix, console, sources=("leaprc.gaff2",))


def check_atom_order(mol2_file: str, gjf: str) -> Tuple[bool, str]:
    """Cheap guard that a mol2 and a Gaussian input list the same atoms in order.

    Compares atom counts and the element letter of each position (the mol2's
    atom-type initial against the input's element symbol). It cannot prove
    the orders agree, but it catches a reordered or re-hydrogenated file
    before paramfit evaluates MM energies on scrambled coordinates.
    """
    try:
        types = _mol2_types(mol2_file)
        elements = gjf_elements(gjf)
    except Exception as e:  # noqa: BLE001
        return False, f"could not read atoms: {e}"
    if len(types) != len(elements):
        return False, f"{mol2_file} has {len(types)} atoms but {gjf} has {len(elements)}"
    for k, (t, el) in enumerate(zip(types, elements), 1):
        if t[:1].upper() != el[:1].upper():
            return False, f"atom {k}: mol2 type {t} vs Gaussian element {el}"
    return True, ""
