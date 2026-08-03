"""Local, blocking energy minimization to relieve bad contacts before PB.

PB pKa shifts are exquisitely sensitive to steric clashes in the input
coordinates: an overlapping atom pair produces a spurious reaction-field
energy that contaminates every leg of the Bashford cycle (pbt-3/pbt-4). This
module runs a short, optionally restrained minimization with an Amber MD
engine so the PB calculations see a clash-free structure.

It is deliberately self-contained and *synchronous*, mirroring how
``pb_backend`` runs pbsa: the pbt-m checklist step blocks until the engine
exits and the minimized rst7 is on disk. The mdin defaults track
``md_templates/builtin/basic/00_minimization.mdin`` (the aqueous-metalloprotein
protocol); the pbt-m prompts let the user override the key variables.

Multi-state note: the public entry point operates on a single (prmtop, rst7)
pair and returns the minimized rst7. The pbt-m handler stores results in a
prmtop-keyed map (``pb_titrate_minimized``), so extending to per-microstate
minimization is just a loop over the microstate pairs — no change here.
"""
from __future__ import annotations

import os
import re
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional

# Defaults mirror md_templates/builtin/basic/00_minimization.mdin.
DEFAULT_MAXCYC = 5000
DEFAULT_NCYC = 1000
DEFAULT_CUT = 10.0
DEFAULT_NTR = True
DEFAULT_RESTRAINT_WT = 10.0
DEFAULT_RESTRAINTMASK = "@CA,C,O,N&!:WAT"
DEFAULT_ENGINE = "pmemd"
DEFAULT_MPI_TASKS = 4

# Engines exposed to the user. pmemd.cuda is intentionally omitted: a clash-
# relief minimization is short and CPU engines avoid GPU-availability surprises
# inside an interactive checklist step.
SUPPORTED_ENGINES = ("sander", "pmemd", "pmemd.MPI")


@dataclass
class MinimizationParams:
    """User-tunable knobs for the pre-PB minimization (see module docstring)."""
    maxcyc: int = DEFAULT_MAXCYC
    ncyc: int = DEFAULT_NCYC
    cut: float = DEFAULT_CUT
    ntr: bool = DEFAULT_NTR
    restraint_wt: float = DEFAULT_RESTRAINT_WT
    restraintmask: str = DEFAULT_RESTRAINTMASK
    engine: str = DEFAULT_ENGINE
    mpi_tasks: int = DEFAULT_MPI_TASKS


def detect_periodic(prmtop: Path, rst7: Path) -> bool:
    """True if the system carries a periodic box (-> ntb=1), else False (ntb=0).

    A solvated topology from the Topology Generator normally has a box and
    should be minimized under constant-volume PBC. A dry/gas-phase system has
    none, so we fall back to a non-periodic (vacuum) minimization, which is
    still fine for relieving steric clashes.
    """
    import parmed as pmd
    try:
        s = pmd.load_file(str(prmtop), str(rst7))
        return s.box is not None
    except Exception:
        # If parmed can't decide, assume periodic (the common post-tLEaP case).
        return True


def build_mdin(params: MinimizationParams, periodic: bool) -> str:
    """Render the minimization mdin from ``params`` and box presence."""
    ntb = 1 if periodic else 0
    if params.ntr:
        restraint_block = (
            "\n  ! Restraints\n"
            "  ntr=1,                            ! Positional restraints (see -ref)\n"
            f"  restraint_wt={params.restraint_wt},                ! kcal/mol-A2\n"
            f"  restraintmask='{params.restraintmask}',\n"
        )
    else:
        restraint_block = (
            "\n  ! Restraints\n"
            "  ntr=0,                            ! Unrestrained minimization\n"
        )
    return (
        "! pb_titrate pre-PB clash relief; defaults track\n"
        "! md_templates/builtin/basic/00_minimization.mdin\n"
        "Energy minimization before PB titration\n"
        "&cntrl\n"
        "  ! Minimization control\n"
        "  imin=1,                           ! Energy minimization\n"
        "  ntmin=1,                          ! Steepest descent for ncyc, then CG\n"
        f"  ncyc={params.ncyc},                        ! SD -> CG switch\n"
        f"  maxcyc={params.maxcyc},                      ! Maximum cycles\n\n"
        "  ! System control\n"
        f"  ntb={ntb},                            ! {'Periodic, constant volume' if ntb else 'Non-periodic (no box)'}\n"
        f"  cut={params.cut},                         ! Nonbonded cutoff (A)\n\n"
        "  ! Output control\n"
        "  ntwx=0,                           ! No trajectory for a pre-step\n"
        "  ntpr=100,                         ! Energy print frequency\n"
        "  ntwr=1000,                        ! Restart write frequency\n"
        f"{restraint_block}"
        "/\n"
    )


def resolve_mpirun() -> str:
    """Return an ``mpirun`` matched to the MPI that ``pmemd.MPI`` links against.

    A bare ``mpirun`` on PATH is unreliable: when a conda env is active its
    bundled mpirun (often a different MPI family/version than the AmberTools
    build) shadows the right one, and launching an MPICH-linked ``pmemd.MPI``
    with OpenMPI's mpirun (or vice versa) makes every rank abort
    (``write_line error; fd=-1`` / ``prterun`` / ``MPI_Abort``). Resolve the
    launcher from the actual library ``pmemd.MPI`` is linked against:

      1. find pmemd.MPI on PATH,
      2. ``ldd`` it, read the resolved path of its ``libmpi*`` entry,
      3. that library's install prefix is ``<libdir>/..``; use
         ``<prefix>/bin/mpirun`` if it exists.

    Falls back to the ``$PROPREP_MPIRUN`` override, then a bare ``mpirun``, if
    detection fails (non-Linux, static link, unusual layout).
    """
    override = os.environ.get("PROPREP_MPIRUN")
    if override:
        return override

    pmemd = shutil.which("pmemd.MPI")
    if pmemd:
        try:
            ldd = subprocess.run(["ldd", pmemd], capture_output=True,
                                 text=True, timeout=20)
            for line in ldd.stdout.splitlines():
                # e.g. "libmpi.so.12 => /usr/local/mpich-3.2.1/lib/libmpi.so.12 (0x..)"
                m = re.search(r"\blibmpi\S*\s*=>\s*(\S+)", line)
                if not m:
                    continue
                libpath = Path(m.group(1))
                if not libpath.is_absolute():
                    continue
                # <prefix>/lib/libmpi... -> <prefix>/bin/mpirun
                prefix = libpath.parent.parent
                cand = prefix / "bin" / "mpirun"
                if cand.exists():
                    return str(cand)
        except (OSError, subprocess.SubprocessError):
            pass

    return "mpirun"


def build_command(
    params: MinimizationParams,
    mdin: Path,
    prmtop: Path,
    rst7: Path,
    mdout: Path,
    restart: Path,
    ref: Optional[Path] = None,
    info: Optional[Path] = None,
) -> List[str]:
    """Map the engine choice to a full argv for the minimization run."""
    base = [
        "-O",
        "-i", str(mdin),
        "-p", str(prmtop),
        "-c", str(rst7),
        "-o", str(mdout),
        "-r", str(restart),
    ]
    if ref is not None:
        base += ["-ref", str(ref)]

    engine = params.engine
    if engine == "sander":
        return ["sander"] + base
    if engine == "pmemd":
        if info is not None:
            base += ["-inf", str(info)]
        return ["pmemd"] + base
    if engine == "pmemd.MPI":
        if info is not None:
            base += ["-inf", str(info)]
        return [resolve_mpirun(), "-np", str(params.mpi_tasks),
                "pmemd.MPI"] + base
    raise ValueError(
        f"Unsupported engine {engine!r}; choose one of {SUPPORTED_ENGINES}.")


def _required_binaries(engine: str) -> List[str]:
    return ["mpirun", "pmemd.MPI"] if engine == "pmemd.MPI" else [engine]


def check_engine_available(params: MinimizationParams) -> None:
    """Raise FileNotFoundError if the engine's binaries aren't on PATH."""
    missing = [b for b in _required_binaries(params.engine)
               if shutil.which(b) is None]
    if missing:
        raise FileNotFoundError(
            f"Required executable(s) {', '.join(missing)} not found on PATH. "
            f"Source your Amber environment (e.g. `source $AMBERHOME/amber.sh`) "
            f"before running the minimization step, or pick a different engine.")


def _tail(path: Path, n: int = 40) -> str:
    try:
        return "\n".join(Path(path).read_text().splitlines()[-n:])
    except OSError:
        return ""


def latest_step_block(mdout: Path, max_lines: int = 14,
                      require_complete: bool = False) -> Optional[str]:
    """Return the most recent ``NSTEP ENERGY ...`` report block from an mdout.

    Drives the live minimization progress display: parses the growing mdout and
    returns the latest energy report — the periodic ntpr block during the run,
    or the FINAL RESULTS block once minimization ends. Returns None before the
    engine has printed its first report.

    A report block spans several lines (the NSTEP header, then BOND/ANGLE/...,
    ending with the ``EAMBER`` line). When ``require_complete`` is set, return
    None for a block whose terminal ``EAMBER`` line hasn't been written yet, so
    the caller can keep showing the last complete block rather than a frame
    caught mid-write.
    """
    try:
        lines = Path(mdout).read_text().splitlines()
    except OSError:
        return None
    hdr = None
    for i in range(len(lines) - 1, -1, -1):
        if "NSTEP" in lines[i] and "ENERGY" in lines[i]:
            hdr = i
            break
    if hdr is None:
        return None
    block = [lines[hdr]]
    complete = False
    for k in range(hdr + 1, min(hdr + max_lines, len(lines))):
        block.append(lines[k])
        if "EAMBER" in lines[k]:
            complete = True
            break
    if require_complete and not complete:
        return None
    return "\n".join(block).rstrip("\n") or None


def _progress_renderable(block: Optional[str], params: "MinimizationParams"):
    """Build the rich Panel shown while the engine runs (latest NSTEP block)."""
    from rich.panel import Panel
    from rich.text import Text
    if block:
        body = Text(block, style="green")
    else:
        body = Text("waiting for the first energy report (ntpr=100)…",
                    style="grey50")
    return Panel(
        body,
        title=f"[cyan]Minimizing — {params.engine}[/cyan]",
        subtitle="[grey50]live min.out — latest step[/grey50]",
        border_style="cyan", expand=False)


def _run_engine_streaming(cmd, work_dir, mdout_path, log_path, params,
                          console) -> int:
    """Run the engine, streaming the latest energy block to ``console``.

    Engine stdout/stderr are redirected to ``log_path`` so the live display
    stays clean; the latest NSTEP report is polled from the growing mdout.
    Returns the process exit code. Falls back to a plain blocking wait when no
    console is given or rich's Live is unavailable.
    """
    import time
    with open(log_path, "w") as logf:
        proc = subprocess.Popen(
            cmd, cwd=str(work_dir), stdout=logf,
            stderr=subprocess.STDOUT, text=True)

        if console is None:
            proc.wait()
            return proc.returncode

        try:
            from rich.live import Live
        except Exception:
            Live = None

        if Live is None:
            # No rich Live: reprint the block whenever it advances. Only print
            # fully-written blocks so a mid-write frame isn't shown.
            last = None
            while proc.poll() is None:
                block = latest_step_block(mdout_path, require_complete=True)
                if block and block != last:
                    console.print(block)
                    last = block
                time.sleep(2.0)
            return proc.returncode

        with Live(console=console, refresh_per_second=4) as live:
            last_block = None
            while proc.poll() is None:
                # Steady state: advance only on a fully-written block; while the
                # next block is mid-write keep showing the last complete one (no
                # truncated frames). But until the FIRST complete block appears,
                # show whatever partial block exists — there's nothing to flicker
                # between yet, and pmemd.MPI buffers min.out, so requiring a
                # complete block here would leave a long, hung-looking "waiting"
                # banner through the first output flush.
                block = latest_step_block(
                    mdout_path, require_complete=(last_block is not None))
                if block is not None:
                    last_block = block
                live.update(_progress_renderable(last_block, params))
                time.sleep(0.5)
            # Final refresh so the FINAL RESULTS block is the last frame shown,
            # falling back to the last complete block if it isn't flushed yet.
            final_block = latest_step_block(mdout_path, require_complete=True)
            live.update(_progress_renderable(final_block or last_block, params))
        return proc.returncode


def final_energy(mdout: Path) -> Optional[float]:
    """Best-effort: pull the FINAL minimization energy from an Amber mdout."""
    try:
        lines = Path(mdout).read_text().splitlines()
    except OSError:
        return None
    for i, line in enumerate(lines):
        if "FINAL RESULTS" in line:
            # The block is: blank, 'NSTEP ENERGY ...' header, then the values.
            for j in range(i + 1, min(i + 8, len(lines))):
                if "NSTEP" in lines[j] and "ENERGY" in lines[j]:
                    parts = lines[j + 1].split() if j + 1 < len(lines) else []
                    if len(parts) >= 2:
                        try:
                            return float(parts[1])
                        except ValueError:
                            return None
    return None


def run_minimization(
    prmtop,
    rst7,
    work_dir,
    params: MinimizationParams,
    console=None,
) -> Path:
    """Minimize ``(prmtop, rst7)`` in ``work_dir``; return the minimized rst7.

    Synchronous: blocks until the engine exits. Raises FileNotFoundError if the
    engine isn't on PATH, or RuntimeError if it fails / produces no restart.
    """
    prmtop = Path(prmtop).resolve()
    rst7 = Path(rst7).resolve()
    work_dir = Path(work_dir).resolve()
    work_dir.mkdir(parents=True, exist_ok=True)

    check_engine_available(params)

    periodic = detect_periodic(prmtop, rst7)
    mdin_path = work_dir / "min.in"
    mdout_path = work_dir / "min.out"
    restart_path = work_dir / "min.rst7"
    info_path = work_dir / "min.mdinfo"

    mdin_path.write_text(build_mdin(params, periodic))

    ref = rst7 if params.ntr else None
    cmd = build_command(params, mdin_path, prmtop, rst7, mdout_path,
                        restart_path, ref=ref, info=info_path)

    log_path = work_dir / "engine.log"

    if console is not None:
        console.print(f"  engine:   {params.engine}")
        console.print(f"  ntb:      {1 if periodic else 0} "
                      f"({'periodic box' if periodic else 'no box / vacuum'})")
        console.print(f"  command:  {' '.join(cmd)}")
        console.print("  [grey50]Running minimization (live progress below)...[/grey50]")

    returncode = _run_engine_streaming(
        cmd, work_dir, mdout_path, log_path, params, console)
    if returncode != 0:
        detail = _tail(mdout_path) or _tail(log_path) or ""
        # pmemd writes its real startup error (e.g. an MPI-launch mismatch) to
        # fort.116, not to the captured stderr — surface it when present.
        fort116 = work_dir / "fort.116"
        if fort116.exists():
            ftxt = fort116.read_text().strip()
            if ftxt:
                detail = f"{ftxt}\n{detail}".strip()
        raise RuntimeError(
            f"{params.engine} minimization failed (exit {returncode}).\n"
            f"--- engine output ---\n{detail}")
    if not restart_path.exists() or restart_path.stat().st_size == 0:
        raise RuntimeError(
            f"{params.engine} exited 0 but produced no restart at "
            f"{restart_path}.\n--- {mdout_path.name} tail ---\n"
            f"{_tail(mdout_path)}")

    if console is not None:
        e = final_energy(mdout_path)
        if e is not None:
            console.print(f"  [green]Minimized.[/green] Final energy: "
                          f"{e:.2f} kcal/mol")
        else:
            console.print("  [green]Minimized.[/green]")
        console.print(f"  minimized rst7: {restart_path}")

    return restart_path
