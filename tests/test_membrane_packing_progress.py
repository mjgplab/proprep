"""Membrane Builder reports packmol's progress, and sets no hidden time limit.

Seen on 6R2Q (about 550,000 atoms). After "Initiating all-together packing.
This might take a while!" ProPrep printed nothing for 51 minutes and then
"packmol-memgen timed out after 1 hour". It followed packmol-memgen.log, which
goes quiet once packmol starts; the loop-by-loop progress is in packmol.log.
The run was at all-together loop 51 of 100, 53 s a loop: it needed about
1 h 40 min, and the limit was fixed, unannounced and unchangeable.

``tests/data/packmol_6r2q_excerpt.txt`` is cut from that run's packmol.log
(lines untouched). Checked live against real packmol-memgen when written:
lines stream during a build, a time limit stops packmol too (no process left),
and a protein-free POPC patch ends "without perfect packing" after all 100
loops while packmol-memgen reports success.
"""

from pathlib import Path

from proprep.membrane_prep.membrane_config import MembraneConfig
from proprep.membrane_prep.packmol_progress import LOOP_LINE_INTERVAL_S, PackmolProgress
from proprep.membrane_prep.packmol_runner import _LeafletVolumes

# .txt, not .log: the repository ignores *.log, and a fixture that is not
# committed passes only on the machine that made it.
LOG = Path(__file__).parent / "data" / "packmol_6r2q_excerpt.txt"


def _replay(seconds_per_loop):
    now = [0.0]
    progress = PackmolProgress(clock=lambda: now[0])
    shown = []
    for line in LOG.read_text().splitlines():
        if "Starting GENCAN loop" in line:
            now[0] += seconds_per_loop
        out = progress.feed(line)
        if out:
            shown.append(out)
    return progress, shown


def test_each_phase_is_named_with_its_loop_budget():
    _, shown = _replay(53.0)
    phases = [l.splitlines()[0].strip() for l in shown if l.startswith("  Packing:")]
    assert phases == [
        "Packing: POPC (component 2 of 9), up to 20 loops",
        "Packing: POPC (component 3 of 9), up to 20 loops",
        "Packing: WAT (component 4 of 9), up to 20 loops",
        "Packing: All components together, up to 100 loops",
    ]


def test_a_loop_line_carries_packmols_numbers_and_the_pace():
    _, shown = _replay(53.0)
    first = next(l for l in shown if "loop 0 of 100" in l)
    # f = .83818E+05, violations 3.994228 and .30974E+03 in packmol.log
    assert "objective 8.38e+04" in first
    assert "distance violation 3.99" in first
    assert "constraint violation 310" in first
    # One report says nothing about how long a loop takes: it used to read
    # "0.0 s per loop | at most 0.0 s more".
    assert "per loop" not in first and "at most" not in first
    second = next(l for l in shown if "loop 1 of 100" in l)
    assert "53 s per loop" in second
    assert "at most 1 h 27 min more if every loop is needed" in second


def test_a_component_packed_on_its_own_is_said_so():
    _, shown = _replay(53.0)
    assert "    packed: WAT (component 4 of 9)" in shown


def test_summary_says_where_a_stopped_run_stood():
    progress, _ = _replay(53.0)
    summary = progress.summary()
    # packmol counts loops from 0: loop 50 is the 51st.
    assert summary.startswith("packmol was packing: All components together, loop 50 of 100")
    assert "lowest so far" in summary
    # The leaflets ran out of loops; the water packed. (Components 5-9 are cut
    # from the excerpt.)
    assert "POPC.pdb (component 2), POPC.pdb (component 3)" in summary
    assert "WAT.pdb (component 4)" not in summary


def test_fast_loops_do_not_flood_the_terminal():
    _, slow = _replay(LOOP_LINE_INTERVAL_S)
    _, fast = _replay(LOOP_LINE_INTERVAL_S / 50)
    count = lambda shown: sum("loop " in l for l in shown)
    assert count(fast) < count(slow)
    # Phase changes and "packed" are never dropped.
    assert [l for l in fast if "loop " not in l] == [l for l in slow if "loop " not in l]


def test_nothing_started_has_a_plain_summary():
    assert PackmolProgress().summary() == "packmol had not started packing."


# ── A protein that does not cross the membrane ──────────────────────────

def _leaflet_warning(upper, lower):
    check = _LeafletVolumes()
    lines = [f"         in upper leaflet     = {upper}          A^3",
             f"         in lower leaflet     = {lower}      A^3",
             "         in upper water box   = 0            A^3"]
    return [w for w in (check.feed(l) for l in lines) if w]


def test_protein_with_nothing_in_one_leaflet_is_flagged_once():
    """The 6R2Q values: MEMEMBED left the barrel complex under the membrane."""
    warnings = _leaflet_warning("0.0", "7694.67")
    assert len(warnings) == 1
    assert "no volume in the upper leaflet" in warnings[0] and "7695 A^3" in warnings[0]
    assert "PPM3" in warnings[0] and "peripheral" in warnings[0]


def test_protein_in_both_leaflets_or_in_neither_is_not_flagged():
    assert _leaflet_warning("41250.5", "39877.1") == []
    assert _leaflet_warning("0.0", "0.0") == []        # protein-free bilayer


# ── The time limit ──────────────────────────────────────────────────────

def test_there_is_no_time_limit_unless_one_is_set():
    assert MembraneConfig().time_limit_hours is None


# ── Arriving in batches, and what is said before the long phase ─────────

def _all_together_lines():
    text = LOG.read_text().splitlines()
    start = next(i for i, l in enumerate(text) if "Packing all molecules together" in l)
    return text[:60], text[start:]          # the header (names, loop budgets), then the phase


def test_of_several_loops_that_arrive_together_the_latest_is_shown():
    """packmol-memgen flushes packmol.log 8 KB at a time: five or six loops at once."""
    header, phase = _all_together_lines()
    now = [0.0]
    progress = PackmolProgress(clock=lambda: now[0])
    progress.feed_many(header)
    now[0] = 400.0
    shown = progress.feed_many(phase)
    loops = [l for l in shown if l.strip().startswith("loop ")]
    assert len(loops) == 1 and loops[0].strip().startswith("loop 50 of 100")      # not loop 0


def test_loops_that_all_arrived_at_once_give_no_pace_rather_than_zero():
    header, phase = _all_together_lines()
    progress = PackmolProgress(clock=lambda: 0.0)
    progress.feed_many(header)
    shown = progress.feed_many(phase)
    assert not any("per loop" in l for l in shown)


def test_the_long_phase_says_first_that_convergence_is_not_expected():
    _, shown = _replay(53.0)
    text = "\n".join(shown)
    assert text.count("does not wait for a perfect packing") == 1
    assert text.index("does not wait for a perfect packing") < text.index("loop 0 of 100")
    assert "lowest so far" in text.split("does not wait")[1].splitlines()[0]


# ── Getting the log written as packmol runs ─────────────────────────────

def test_a_python_packmol_memgen_is_started_through_the_live_log_launcher(tmp_path):
    import sys
    from proprep.membrane_prep import packmol_runner as pr
    script = tmp_path / "packmol-memgen"
    script.write_text(f"#!{sys.executable}\nprint('x')\n")
    command = pr._live_log_command(str(script))
    assert command[0] == sys.executable and command[-1] == str(script)
    assert command[1].endswith("packmol_memgen_live_log.py")


def test_anything_else_is_started_directly_as_before(tmp_path):
    from proprep.membrane_prep import packmol_runner as pr
    binary = tmp_path / "packmol-memgen"
    binary.write_bytes(b"\x7fELF\x02\x01")
    shell = tmp_path / "wrapper.sh"
    shell.write_text("#!/bin/sh\nexec packmol-memgen \"$@\"\n")
    gone = tmp_path / "packmol-memgen-gone"
    gone.write_text("#!/nonexistent/python\n")
    for path in (binary, shell, gone, tmp_path / "missing"):
        assert pr._live_log_command(str(path)) == [str(path)]


def test_the_launcher_line_buffers_log_files_and_nothing_else(tmp_path):
    from proprep.membrane_prep import packmol_memgen_live_log as launcher
    log = launcher._line_buffered_logs(str(tmp_path / "packmol.log"), "w")
    other = launcher._line_buffered_logs(str(tmp_path / "bilayer.pdb"), "w")
    binary = launcher._line_buffered_logs(str(tmp_path / "x.log"), "wb")
    try:
        assert log.line_buffering is True and other.line_buffering is False
        log.write("Starting GENCAN loop: 0\n")
        assert (tmp_path / "packmol.log").read_text() == "Starting GENCAN loop: 0\n"      # on disk, unflushed by us
    finally:
        for handle in (log, other, binary):
            handle.close()
