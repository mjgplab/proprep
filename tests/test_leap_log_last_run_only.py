"""tLEaP APPENDS to leap.log, so a directory that has seen more than one run
holds every run's messages. _parse_leap_log used to read the whole file and
pair the LAST run's summary line with EVERY run's error blocks, so a clean
build printed "Errors = 0" above a list of errors left by an earlier attempt.
Regression for that: the parser now reports one run, either from the byte
offset recorded before tLEaP was launched or, with no offset, the last
complete run in the file."""

from proprep.tleap_prep.tleap_input_generator import TLeapInputGenerator


_FAILED = (
    "> loadpdb old.pdb\n"
    "FATAL: Atom .R<MN 1>.A<MN 1> does not have a type.\n"
    "Bad: Error! bond not found\n"
    "\n"
    "Exiting LEaP: Errors = 1; Warnings = 0; Notes = 0.\n"
)
_CLEAN = (
    "> loadpdb new.pdb\n"
    "Ok: Warning! renaming residue\n"
    "\n"
    "Noted: Note. ignoring unknown\n"
    "\n"
    "Exiting LEaP: Errors = 0; Warnings = 1; Notes = 1.\n"
)


def _parse(path, offset=None):
    gen = object.__new__(TLeapInputGenerator)
    return gen._parse_leap_log(str(path), start_offset=offset)


def test_earlier_runs_do_not_leak_into_a_clean_build(tmp_path):
    log = tmp_path / "leap.log"
    log.write_text(_FAILED + _CLEAN)

    warnings, errors, notes, summary = _parse(log)

    assert errors == [], "an earlier run's errors were reported against this run"
    assert len(warnings) == 1
    assert len(notes) == 1
    assert summary == "Exiting LEaP: Errors = 0; Warnings = 1; Notes = 1."


def test_offset_recorded_before_the_run_wins(tmp_path):
    log = tmp_path / "leap.log"
    log.write_text(_FAILED)
    offset = TLeapInputGenerator._leap_log_size(str(log))
    with open(log, "a") as handle:
        handle.write(_CLEAN)

    warnings, errors, notes, summary = _parse(log, offset=offset)

    assert errors == []
    assert len(warnings) == 1
    assert summary == "Exiting LEaP: Errors = 0; Warnings = 1; Notes = 1."


def test_a_single_run_is_reported_whole(tmp_path):
    log = tmp_path / "leap.log"
    log.write_text(_FAILED)

    warnings, errors, notes, summary = _parse(log)

    assert len(errors) == 1
    assert summary == "Exiting LEaP: Errors = 1; Warnings = 0; Notes = 0."


def test_a_run_that_dies_before_its_summary_still_reports_its_errors(tmp_path):
    log = tmp_path / "leap.log"
    log.write_text("> loadpdb x.pdb\nBad: Error! died here\n\n")

    warnings, errors, notes, summary = _parse(log)

    assert len(errors) == 1
    assert summary is None


def test_size_of_a_log_that_does_not_exist_yet_is_zero(tmp_path):
    assert TLeapInputGenerator._leap_log_size(str(tmp_path / "absent.log")) == 0
