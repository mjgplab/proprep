"""Compare PB-titrate single-site pKas against ProPKA on the SAME structure.

The Protonation State Analyzer may have run ProPKA earlier, but typically on the
*unminimized* input — whereas PB-titrate runs on the minimized structure. A
direct comparison of those two would conflate the method difference (ProPKA vs
PB) with a coordinate difference (unminimized vs minimized). To make the
comparison legitimate, this module re-runs ProPKA on the exact minimized,
solvent-stripped structure PB titrated (the ``pb_titrate_prmtop`` /
``pb_titrate_rst7`` "dry" files), then lines the two predictions up per residue.

ProPKA reads a PDB with standard residue names (ASP/GLU/HIS/...), so the constph
topology (AS4/GL4/HIP/...) is renamed on the way out. Matching is by
(standard_resname, resnum); chains are ignored (PB sites carry no chain).
"""
from __future__ import annotations

import re
import shutil
import subprocess
import sys
import math
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

# constph / alternate-protonation names -> standard PDB names ProPKA understands.
_TO_STANDARD = {
    "AS4": "ASP", "ASH": "ASP",
    "GL4": "GLU", "GLH": "GLU",
    "HIP": "HIS", "HID": "HIS", "HIE": "HIS",
    "CYM": "CYS", "CYX": "CYS",
    "LYN": "LYS",
}

# PB site resname -> standard name, for matching PB sites to ProPKA output.
PB_TO_STANDARD = {
    "AS4": "ASP", "GL4": "GLU", "HIP": "HIS",
    "CYS": "CYS", "TYR": "TYR", "LYS": "LYS",
}


def to_standard_resname(name: str) -> str:
    return _TO_STANDARD.get(name, name)


def find_propka() -> Optional[str]:
    """Locate the propka3 executable, or None if unavailable."""
    exe = shutil.which("propka3")
    if exe:
        return exe
    # Fall back to the interpreter's own bin dir (conda env).
    cand = Path(sys.executable).parent / "propka3"
    return str(cand) if cand.exists() else None


def write_propka_pdb(prmtop: Path, rst7: Path, out_pdb: Path) -> int:
    """Write a standard-named, hydrogen-free PDB of (prmtop, rst7) for ProPKA.

    Renames constph/alt-protonation residues to standard names and assigns a
    chain id where missing (ProPKA's output keys on a chain column). All
    hydrogens are stripped: ProPKA reconstructs protons internally from
    heavy-atom geometry (a heavy-atom PDB is its native input), and the constph
    titratable residues carry extra carboxyl/imidazole H — including the
    charge-0 "ghost" protons on the deprotonated AS4/GL4 oxygens — that are not
    present in a standard ASP/GLU/HIS and could otherwise skew ProPKA's
    hydrogen-bond/desolvation terms. (Empirically the effect is small, ≤~0.2
    pKa units on a few residues, but stripping H is the correct, defensible
    input.) Returns the number of residues renamed.
    """
    import parmed as pmd
    s = pmd.load_file(str(prmtop), str(rst7))
    n = 0
    for r in s.residues:
        std = _TO_STANDARD.get(r.name)
        if std is not None:
            r.name = std
            n += 1
        if not (getattr(r, "chain", "") or "").strip():
            r.chain = "A"
    s.strip("@/H")  # remove all hydrogens by element (ghost + real)
    out_pdb.parent.mkdir(parents=True, exist_ok=True)
    s.save(str(out_pdb), overwrite=True)
    return n


def run_propka(pdb_path: Path, work_dir: Path) -> Path:
    """Run propka3 on ``pdb_path`` inside ``work_dir``; return the .pka path.

    Raises FileNotFoundError if propka3 isn't installed, or RuntimeError if it
    fails to produce a .pka file.
    """
    exe = find_propka()
    if exe is None:
        raise FileNotFoundError(
            "propka3 not found on PATH. Install ProPKA (it ships with the "
            "Protonation State Analyzer's dependencies) to enable this "
            "comparison.")
    work_dir.mkdir(parents=True, exist_ok=True)
    # propka3 writes <stem>.pka in its CWD; run with cwd=work_dir.
    pka_path = work_dir / (pdb_path.stem + ".pka")
    proc = subprocess.run(
        [exe, pdb_path.name, "--quiet"],
        cwd=str(work_dir), capture_output=True, text=True)
    if not pka_path.exists():
        raise RuntimeError(
            f"propka3 did not produce {pka_path.name} (exit {proc.returncode}).\n"
            f"stderr: {proc.stderr[-1500:]}\nstdout: {proc.stdout[-500:]}")
    return pka_path


# SUMMARY line, e.g.:  "   ASP  45 A     2.55      3.80"
_SUMMARY_RE = re.compile(
    r"^\s*([A-Z]{2,3})\s+(\d+)\s+(\S)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)")


def parse_propka_summary(pka_path: Path) -> Dict[Tuple[str, int], float]:
    """Parse the SUMMARY block of a .pka file -> {(std_resname, resnum): pKa}.

    Termini (N+/C-) and ligand groups are skipped (they don't match standard
    amino-acid PB sites). When a residue appears more than once (rare), the last
    value wins.
    """
    out: Dict[Tuple[str, int], float] = {}
    in_summary = False
    for ln in pka_path.read_text().splitlines():
        if "SUMMARY OF THIS PREDICTION" in ln:
            in_summary = True
            continue
        if in_summary:
            stripped = ln.strip()
            if stripped.startswith("---") or stripped.startswith("Free energy"):
                break
            m = _SUMMARY_RE.match(ln)
            if m:
                resname, resnum, _chain, pka = (
                    m.group(1), int(m.group(2)), m.group(3), float(m.group(4)))
                out[(resname, resnum)] = pka
    return out


def build_comparison(
        pb_pka: Dict[Tuple[str, int], Optional[float]],
        propka: Dict[Tuple[str, int], float],
        effective_pka: Optional[Dict[Tuple[str, int], Any]] = None,
) -> List[Dict[str, object]]:
    """Line up PB (single-site), effective (coupled) and ProPKA pKas per site.

    Returns one row dict per PB site:
      {resname, resnum, pb_pka, effective_pka, propka_pka, delta, delta_eff}
    - ``pb_pka`` is the single-site (intrinsic) PB pKa; ``delta`` = pb_pka − ProPKA.
    - ``effective_pka`` (when ``effective_pka`` map given) is the coupling-aware
      pKa from a coupled titration — a float, a sentinel string ('< x'/'> x') for
      locked sites, or None if not available for that site. ``delta_eff`` =
      effective − ProPKA, only when both are real floats.
    propka_pka / delta are None when ProPKA has no value (e.g. disulfide CYS).
    """
    effective_pka = effective_pka or {}
    rows: List[Dict[str, object]] = []
    for (rn, num) in sorted(pb_pka, key=lambda k: (k[1], k[0])):
        std = PB_TO_STANDARD.get(rn, rn)
        pp = propka.get((std, num))
        pbv = pb_pka[(rn, num)]
        delta = (pbv - pp) if (pp is not None and pbv is not None) else None
        eff = effective_pka.get((rn, num))
        delta_eff = ((eff - pp) if (pp is not None and isinstance(eff, (int, float)))
                     else None)
        rows.append({
            "resname": rn, "resnum": num,
            "pb_pka": pbv, "effective_pka": eff, "propka_pka": pp,
            "delta": delta, "delta_eff": delta_eff,
        })
    return rows


def summary_statistics(rows: List[Dict[str, Any]],
                        pH: float, use: str = "effective"
                        ) -> Optional[Dict[str, Any]]:
    """Aggregate PB-vs-ProPKA agreement statistics over the comparison rows.

    ``use``: which PB pKa to compare against ProPKA —
      'effective' (default): the coupling-aware pKa (delta_eff column). Falls
        back to single-site automatically if NO row has an effective value
        (e.g. coupling wasn't run), so the stats are never empty when single-
        site data exists.
      'single_site': the intrinsic single-site PB pKa (delta column).

    Considers only rows where the chosen PB pKa AND ProPKA both exist. Returns
    None if there are no such rows. Pure-Python (no numpy).

    Keys: basis ('effective'|'single_site'), n, mean_delta, median_delta,
    mae (mean|delta|), rmsd, pearson_r (None if degenerate), frac_within_2
    (|delta| <= 2 pH units), binary_agreement + binary_n_agree (PB and ProPKA
    agree on whether the site is protonated at `pH`; protonated <=> pKa > pH), pH.
    """
    def _paired(pb_key, delta_key):
        return [(r[pb_key], r["propka_pka"], r[delta_key]) for r in rows
                if r.get(delta_key) is not None
                and isinstance(r.get(pb_key), (int, float))
                and r.get("propka_pka") is not None]

    basis = use
    if use == "effective":
        paired = _paired("effective_pka", "delta_eff")
        if not paired:                      # coupling not run → use single-site
            basis = "single_site"
            paired = _paired("pb_pka", "delta")
    else:
        basis = "single_site"
        paired = _paired("pb_pka", "delta")

    n = len(paired)
    if n == 0:
        return None

    deltas = [d for _, _, d in paired]
    mean_delta = sum(deltas) / n
    sd = sorted(deltas)
    median_delta = (sd[n // 2] if n % 2 == 1
                    else 0.5 * (sd[n // 2 - 1] + sd[n // 2]))
    # Mean absolute error: average magnitude of PB−ProPKA disagreement,
    # independent of the systematic offset captured by mean_delta.
    mae = sum(abs(d) for d in deltas) / n
    rmsd = math.sqrt(sum(d * d for d in deltas) / n)
    frac_within_2 = sum(1 for d in deltas if abs(d) <= 2.0) / n

    xs = [pb for pb, _, _ in paired]
    ys = [pp for _, pp, _ in paired]
    mx = sum(xs) / n
    my = sum(ys) / n
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    pearson_r = (sxy / math.sqrt(sxx * syy)) if (sxx > 0 and syy > 0) else None

    # Binary protonation call at pH: protonated when pKa > pH (true for both
    # acids and bases — the titratable proton is on below the pKa).
    n_agree = sum(1 for pb, pp, _ in paired if (pb > pH) == (pp > pH))

    return {
        "basis": basis,   # 'effective' or 'single_site' (what PB pKa was used)
        "n": n,
        "mean_delta": mean_delta,
        "median_delta": median_delta,
        "mae": mae,
        "rmsd": rmsd,
        "pearson_r": pearson_r,
        "frac_within_2": frac_within_2,
        "binary_agreement": n_agree / n,
        "binary_n_agree": n_agree,
        "pH": pH,
    }

