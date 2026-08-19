"""
Triage categories must survive a resume, or four steps silently do nothing.

triage_results is populated by step 3 and lives on the preprocessor instance.
A resumed run builds a NEW preprocessor with step 3 already marked complete, so
it never runs and the attribute stays empty. Every step that reads it then
concludes the structure has no such residues:

    step 6  organic          v == 'B'
    step 7  organometallic   v == 'C'
    step 8  metal clusters   v == 'F'
    step 9  isolated metals  v == 'D'

Each prints "No ... in structure - skipping" and completes SUCCESSFULLY, so the
checklist shows a tick for work that never happened. Reported on a 4UHX run
where step 3 had recorded "1309 protein, 2 organic, 2 metal cluster" and step 6
then found nothing.

The data was already being persisted -- _run_triage_only writes
``preprocessing_triage`` -- it was simply never read back.
"""

from collections import Counter

import pytest
from rich.console import Console

from proprep.forcefield_prep.structure_preprocessor import StructurePreprocessor


TRIAGE = {
    "A:1:MET": "A", "A:2:ASP": "A",
    "A:1310:FES": "F", "A:1312:MOS": "F",
    "A:1311:MTE": "B", "A:1313:FAD": "B",
    "A:900:ZN": "D", "A:2000:HOH": "E",
}


class _Workspace:
    def __init__(self, data=None):
        self.data = dict(data or {})

    def get(self, key, default=None):
        return self.data.get(key, default)

    def set(self, key, value):
        self.data[key] = value


def _resumed(workspace_data=None, pdb_file=None):
    """A preprocessor as a resume creates it: no triage on the instance."""
    preprocessor = StructurePreprocessor.__new__(StructurePreprocessor)
    preprocessor.triage_results = {}
    preprocessor.workspace = _Workspace(workspace_data)
    preprocessor.console = Console(quiet=True)
    preprocessor._pdb_file = pdb_file
    return preprocessor


# --------------------------------------------------------------------------- #
# the reported failure
# --------------------------------------------------------------------------- #

def test_triage_is_restored_from_the_workspace():
    preprocessor = _resumed({"preprocessing_triage": TRIAGE})

    assert preprocessor._ensure_triage_results() == TRIAGE


def test_the_organic_residues_are_found_again():
    """Step 6 reported "No organic small molecules" for a structure with two."""
    triage = _resumed({"preprocessing_triage": TRIAGE})._ensure_triage_results()

    assert sorted(k for k, v in triage.items() if v == "B") == [
        "A:1311:MTE", "A:1313:FAD"]


@pytest.mark.parametrize("category,expected", [
    ("B", 2),   # organic          -> step 6
    ("C", 0),   # organometallic   -> step 7
    ("F", 2),   # metal clusters   -> step 8
    ("D", 1),   # isolated metals  -> step 9
])
def test_every_dependent_step_sees_its_residues(category, expected):
    """All four read the same attribute, so all four were affected."""
    triage = _resumed({"preprocessing_triage": TRIAGE})._ensure_triage_results()

    assert sum(1 for v in triage.values() if v == category) == expected


def test_it_reads_the_key_that_is_already_written():
    """
    _run_triage_only persists `preprocessing_triage`. The data was there all
    along; nothing read it. A second key would have been redundant.
    """
    import inspect

    source = inspect.getsource(StructurePreprocessor._ensure_triage_results)

    assert "preprocessing_triage" in source


# --------------------------------------------------------------------------- #
# not disturbing the normal path
# --------------------------------------------------------------------------- #

def test_an_instance_that_already_has_triage_is_untouched():
    """No workspace round trip mid-run, and no chance of a stale overwrite."""
    preprocessor = _resumed({"preprocessing_triage": {"A:1:MET": "A"}})
    preprocessor.triage_results = dict(TRIAGE)

    assert preprocessor._ensure_triage_results() == TRIAGE


def test_no_workspace_and_no_structure_yields_empty():
    """Callers already handle an empty result; it must not raise."""
    preprocessor = _resumed()

    assert preprocessor._ensure_triage_results() == {}


def test_an_empty_saved_value_is_not_used():
    preprocessor = _resumed({"preprocessing_triage": {}})

    assert preprocessor._ensure_triage_results() == {}


def test_a_non_dict_saved_value_is_ignored():
    preprocessor = _resumed({"preprocessing_triage": "not a dict"})

    assert preprocessor._ensure_triage_results() == {}


def test_re_running_triage_is_the_fallback(monkeypatch):
    """Triage is deterministic and needs only the structure."""
    preprocessor = _resumed(pdb_file="/tmp/x.pdb")
    monkeypatch.setattr(StructurePreprocessor, "_run_triage_only",
                        lambda self, path: dict(TRIAGE))

    assert preprocessor._ensure_triage_results() == TRIAGE


def test_a_failing_re_run_does_not_raise(monkeypatch):
    def boom(self, path):
        raise RuntimeError("no structure")

    preprocessor = _resumed(pdb_file="/tmp/x.pdb")
    monkeypatch.setattr(StructurePreprocessor, "_run_triage_only", boom)

    assert preprocessor._ensure_triage_results() == {}
