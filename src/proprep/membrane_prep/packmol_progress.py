"""
Read packmol's own log as it is written and say how the packing is going.

packmol-memgen's log (packmol-memgen.log) goes quiet at "Running Packmol...".
Everything after that, which is nearly all of the run time, is written by
packmol itself to packmol.log: which component it is packing, and after every
GENCAN loop the objective function and the two violations it is driving to
zero. PackmolProgress turns those lines into one line per loop.

What the numbers are (packmol's own names are kept so they can be checked
against packmol.log):

* objective function ``f``: zero when nothing overlaps and every molecule is
  inside its region. It falls within a run of loops; when packmol gives up on
  the worst-placed molecules it moves them and ``f`` jumps back up. That
  sawtooth is the normal course of a packmol-memgen bilayer, not a symptom: a
  protein-free 40 A POPC patch does it for all 100 all-together loops, ends
  "without perfect packing", and packmol-memgen keeps packmol's best structure.
* "violation of target distance" and "violation of the constraints": packmol
  declares a component packed when both are below its precision
  (``precision``, 0.01 unless packmol.inp sets it).
"""

import re
import time
from typing import Callable, Dict, List, Optional

# At most one loop line per this many seconds. Display only: small systems run
# hundreds of loops a second and would bury the terminal. Phase changes and
# "solved" lines are always shown.
LOOP_LINE_INTERVAL_S = 5.0

ALL_TOGETHER = "all"

_STRUCTURE = re.compile(r"^\s*Structure\s+(\d+)\s*:\s*(\S+?)\(\s*(\d+)\s+atoms\)")
_LIMIT_ALL = re.compile(r"Maximum number of GENCAN loops for all molecule packing:\s*(\d+)")
_LIMIT_TYPE = re.compile(r"Maximum number of GENCAN loops for type:\s*(\d+)\s*:\s*(\d+)")
_PHASE_TYPE = re.compile(r"^\s*Packing molecules of type:\s*(\d+)\s*$")
_PHASE_ALL = re.compile(r"^\s*Packing all molecules together")
_LOOP = re.compile(r"Starting GENCAN loop:\s*(\d+)")
_FVALUE = re.compile(r"Function value from last GENCAN loop: f =\s*(\S+)")
_DIST = re.compile(r"Maximum violation of target distance:\s*(\S+)")
_CONSTR = re.compile(r"Maximum violation of the constraints:\s*(\S+)")
_SOLVED = re.compile(r"Packing solved for molecules of type\s+(\d+)")
_SUCCESS = re.compile(r"^\s*Success!")
_NOT_PERFECT = re.compile(r"ENDED WITHOUT PERFECT PACKING")


def _number(text: str) -> Optional[float]:
    """packmol writes Fortran reals such as ``.83818E+05``."""
    try:
        return float(text)
    except ValueError:
        return None


def _duration(seconds: float) -> str:
    if seconds < 10:
        return f"{seconds:.1f} s"
    seconds = int(round(seconds))
    if seconds < 60:
        return f"{seconds} s"
    minutes, seconds = divmod(seconds, 60)
    if minutes < 60:
        return f"{minutes} min {seconds:02d} s"
    hours, minutes = divmod(minutes, 60)
    return f"{hours} h {minutes:02d} min"


class PackmolProgress:
    """Feed it packmol.log line by line; it returns the lines worth showing."""

    def __init__(self, clock: Callable[[], float] = time.monotonic):
        self._clock = clock
        self.structures: Dict[int, str] = {}
        self.loop_limit: Dict[object, int] = {}
        self.phase: Optional[object] = None
        self.loop: Optional[int] = None
        self.f: Optional[float] = None
        self.best_f: Optional[float] = None
        self.distance_violation: Optional[float] = None
        self.constraint_violation: Optional[float] = None
        self.solved: List[int] = []
        self.finished: Optional[str] = None      # "success" | "not perfect"
        self._first_report: Optional[float] = None     # when the phase's first loop was reported
        self._last_report: Optional[float] = None
        self._loops_done = 0
        self._last_shown: Optional[float] = None

    # ── Names ────────────────────────────────────────────────────────────

    def _phase_name(self) -> str:
        if self.phase == ALL_TOGETHER:
            return "All components together"
        name = self.structures.get(self.phase, "")
        name = name[:-4] if name.lower().endswith(".pdb") else name
        return f"{name} (component {self.phase} of {len(self.structures)})" if name \
            else f"Component {self.phase}"

    # ── Reading ──────────────────────────────────────────────────────────

    def feed(self, line: str) -> Optional[str]:
        """One line at a time; at most one loop line per LOOP_LINE_INTERVAL_S."""
        kind, text = self._consume(line)
        if kind != "loop":
            return text
        return text if self._may_show_loop() else None

    def feed_many(self, lines) -> List[str]:
        """A batch of lines that arrived together, as the log is read in chunks.

        Of several loops in one batch the LATEST is shown, not the first: the
        first is already out of date. A loop line pending when the phase changes
        is shown before the change, so each phase ends on its last state.
        """
        shown: List[str] = []
        pending: Optional[str] = None
        for line in lines:
            kind, text = self._consume(line)
            if kind == "loop":
                pending = text
            elif text:
                if pending and kind in ("phase", "packed", "end"):
                    shown.append(pending)
                    self._last_shown = self._clock()
                    pending = None
                shown.extend(text.split("\n"))
        if pending and self._may_show_loop():
            shown.append(pending)
        return shown

    def _may_show_loop(self) -> bool:
        now = self._clock()
        limit = self.loop_limit.get(self.phase)
        last_loop = limit is not None and self.loop is not None and self.loop >= limit
        if self._last_shown is not None and now - self._last_shown < LOOP_LINE_INTERVAL_S and not last_loop:
            return False
        self._last_shown = now
        return True

    def _consume(self, line: str):
        """Parse one line. Returns (kind, text): kind is "loop", "phase", "packed", "end" or None."""
        m = _STRUCTURE.match(line)
        if m:
            self.structures[int(m.group(1))] = m.group(2)
            return None, None
        m = _LIMIT_ALL.search(line)
        if m:
            self.loop_limit[ALL_TOGETHER] = int(m.group(1))
            return None, None
        m = _LIMIT_TYPE.search(line)
        if m:
            self.loop_limit[int(m.group(1))] = int(m.group(2))
            return None, None

        m = _PHASE_TYPE.match(line)
        if m:
            return "phase", self._start_phase(int(m.group(1)))
        if _PHASE_ALL.match(line):
            return "phase", self._start_phase(ALL_TOGETHER)

        m = _LOOP.search(line)
        if m:
            self.loop = int(m.group(1))
            return None, None
        m = _FVALUE.search(line)
        if m:
            self.f = _number(m.group(1))
            if self.f is not None and (self.best_f is None or self.f < self.best_f):
                self.best_f = self.f
            return None, None
        m = _DIST.search(line)
        if m and self.loop is not None:
            self.distance_violation = _number(m.group(1))
            return None, None
        m = _CONSTR.search(line)
        if m and self.loop is not None:
            # Last line of a loop's report.
            self.constraint_violation = _number(m.group(1))
            self._loops_done += 1
            now = self._clock()
            if self._first_report is None:
                self._first_report = now
            self._last_report = now
            return "loop", "    " + self._state()

        m = _SOLVED.search(line)
        if m:
            self.solved.append(int(m.group(1)))
            return "packed", f"    packed: {self._phase_name()}"
        if _SUCCESS.match(line):
            self.finished = "success"
            return "end", "    packmol: Success, every molecule placed without overlap"
        if _NOT_PERFECT.search(line):
            self.finished = "not perfect"
            # Not a failure: packmol-memgen accepts this and goes on. A protein-free
            # 40 A POPC patch ends the same way, after all 100 loops.
            return "end", ("    packmol used all its loops and stopped short of a perfect packing. "
                    "packmol-memgen accepts that and keeps packmol's best structure; the close "
                    "contacts left in it are for the minimization to relax")
        return None, None

    def _start_phase(self, phase) -> str:
        self.phase = phase
        self.loop = self.f = self.best_f = None
        self.distance_violation = self.constraint_violation = None
        self._first_report = self._last_report = None
        self._loops_done = 0
        self._last_shown = None
        limit = self.loop_limit.get(phase)
        upto = f", up to {limit} loops" if limit else ""
        text = f"  Packing: {self._phase_name()}{upto}"
        if phase == ALL_TOGETHER:
            # Said before an hour of numbers that never reach zero: a user watching
            # them took the run for one that was failing to converge.
            text += ("\n    packmol-memgen does not wait for a perfect packing: it lets packmol use these "
                     "loops and keeps the best structure it reached. The objective falls, jumps when "
                     "packmol moves its worst-placed molecules, and falls again; \"lowest so far\" is "
                     "the number to watch. What overlap is left, the minimization relaxes.")
        return text

    def _seconds_per_loop(self) -> Optional[float]:
        """Mean time from one loop report to the next. None until two reports have been seen:
        the first report of a phase says nothing about how long a loop takes."""
        if self._loops_done < 2 or self._first_report is None or self._last_report is None:
            return None
        if self._last_report <= self._first_report:
            return None             # every report so far arrived in one batch: no interval to go by
        return (self._last_report - self._first_report) / (self._loops_done - 1)

    def _state(self) -> str:
        # packmol's own loop number (it counts from 0), so the line can be found in
        # packmol.log.
        limit = self.loop_limit.get(self.phase)
        parts = [f"loop {self.loop}" + (f" of {limit}" if limit else "")]
        if self.f is not None:
            best = f" (lowest so far {self.best_f:.3g})" if self.best_f is not None \
                and self.best_f < self.f else ""
            parts.append(f"objective {self.f:.3g}{best}")
        if self.distance_violation is not None:
            parts.append(f"distance violation {self.distance_violation:.3g}")
        if self.constraint_violation is not None:
            parts.append(f"constraint violation {self.constraint_violation:.3g}")
        per_loop = self._seconds_per_loop()
        if per_loop is not None:
            parts.append(f"{_duration(per_loop)} per loop")
            if limit and self.loop is not None and self.loop < limit:
                parts.append(f"at most {_duration(per_loop * (limit - self.loop))} "
                             f"more if every loop is needed")
        return " | ".join(parts)

    # ── Where it stood ───────────────────────────────────────────────────

    def summary(self) -> str:
        """For a run that was stopped or failed: how far packmol had got."""
        if self.phase is None:
            return "packmol had not started packing."
        if self.finished == "success":
            return "packmol had finished packing successfully."
        text = f"packmol was packing: {self._phase_name()}, {self._state()}."
        unsolved = [i for i in sorted(self.structures) if i not in self.solved
                    and i in self.loop_limit and (self.phase == ALL_TOGETHER or i < self.phase)]
        if unsolved:
            names = ", ".join(f"{self.structures[i]} (component {i})" for i in unsolved)
            text += f" Not packed on their own within their loops: {names}."
        return text
