"""Step 8 says which artifact is missing, and the RESP selection survives a resume.

"Missing ESP or AC file from previous steps" named both and identified neither,
directly after step 7 had announced the AC file -- so the natural reading was
the wrong one. And the step-5 selection lived only in memory, so a resumed
session re-selected and could demand QM the user had not run.
"""

from proprep.forcefield_prep.modified_amino_acid_parameterizer import (
    ModifiedAAWorkflowManager,
)


class _Workspace:
    def __init__(self):
        self._d = {}

    def get(self, key, default=None):
        return self._d.get(key, default)

    def set(self, key, value):
        self._d[key] = value


class _Processor:
    def __init__(self, workspace):
        self._workspace = workspace

    def _get_workspace(self):
        return self._workspace


def _manager(**attrs):
    m = ModifiedAAWorkflowManager.__new__(ModifiedAAWorkflowManager)
    m.step_results = {}
    m.amino_acid = "CS1"
    m.conformer_mode = "de_novo"      # skips the from-structure disk recovery
    m.processor = None
    for key, value in attrs.items():
        setattr(m, key, value)
    return m


# --------------------------------------------------------------------------
# what the failure message says
# --------------------------------------------------------------------------

def test_the_message_names_both_artifacts_when_both_are_missing():
    result = _manager()._run_step_8()
    assert result["success"] is False
    assert "combined ESP from step 6" in result["message"]
    assert "AC file from step 7" in result["message"]


def test_it_names_the_expected_esp_filename_and_where_it_looked():
    message = _manager()._run_step_8()["message"]
    assert "cs1_combined.esp" in message
    assert "/" in message                      # the directory it searched


def test_a_present_ac_file_is_not_blamed():
    m = _manager(step_results={"step_7": {"ac_file": "CS1.ac", "charge": 0}})
    message = m._run_step_8()["message"]
    assert "combined ESP" in message
    assert "AC file" not in message


def test_it_explains_why_step_6_may_have_produced_nothing():
    assert "pauses" in _manager()._run_step_8()["message"]


# --------------------------------------------------------------------------
# the selection survives the session
# --------------------------------------------------------------------------

def test_selection_round_trips_through_the_workspace():
    m = _manager(processor=_Processor(_Workspace()))
    m._persist_selection([("xtal", 8), ("xtal", 9)])
    assert m._persisted_selection() == [("xtal", 8), ("xtal", 9)]


def test_an_unscanned_conformer_key_round_trips():
    """Its scan-point slot is None, which has to survive the JSON round trip."""
    m = _manager(processor=_Processor(_Workspace()))
    m._persist_selection([("xtal", None)])
    assert m._persisted_selection() == [("xtal", None)]


def test_another_residues_selection_is_not_inherited():
    workspace = _Workspace()
    _manager(processor=_Processor(workspace))._persist_selection([("xtal", 3)])
    other = _manager(amino_acid="SEP", processor=_Processor(workspace))
    assert other._persisted_selection() is None


def test_a_persisted_selection_drives_the_candidate_list():
    m = _manager(processor=_Processor(_Workspace()))
    m._persist_selection([("xtal", 9)])
    m._structure_candidates = lambda: [{"key": ("xtal", 8)}, {"key": ("xtal", 9)}]
    assert [c["key"] for c in m._selected_candidates()] == [("xtal", 9)]


def test_persisting_without_a_workspace_is_harmless():
    m = _manager(processor=None)
    m._persist_selection([("xtal", 1)])        # must not raise
    assert m._persisted_selection() is None
