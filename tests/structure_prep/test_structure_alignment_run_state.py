"""The Structure Aligner forgets one run's state before the next.

The module instance lives for the whole session. ``hetatm_to_add`` and
``final_transformation_matrices`` are keyed by a structure's position in
``self.structures``, so entries carried over from an earlier run name
whatever structure sits at that position in the later run. A second
redox-site alignment in one session used to re-add the first run's ions and
waters, and it put them into the second run's reference structure
(calmodulin 1CLL: four spurious Ca2+/HOH at A:293-296).
"""

from unittest.mock import Mock, patch

from proprep.structure_prep.structure_alignment import StructureAlignmentModule


def _module_with_stale_state():
    module = StructureAlignmentModule()
    module.initialize()
    module.set_processor(Mock(_get_workspace=Mock(return_value={})))
    module.structures = [("a.pdb", object()), ("b.pdb", object())]
    module.reference_idx = 0
    module.hetatm_to_add = [(object(), object(), [1])]
    module.final_transformation_matrices = {1: object()}
    module.aligned_structures = {0: object(), 1: object()}
    module.alignment_results = [{"step": 1}]
    module.residue_mappings = [object()]
    return module


def _assert_clean(module):
    assert module.hetatm_to_add == []
    assert module.final_transformation_matrices == {}
    assert module.structures == []
    assert module.aligned_structures == {}
    assert module.alignment_results == []
    assert module.residue_mappings == []
    assert module.reference_idx is None


def test_plain_alignment_starts_from_a_clean_slate():
    module = _module_with_stale_state()
    # Fewer than two structures makes the run return right after loading,
    # which is after the reset under test.
    with patch.object(module, "_load_structures_from_input", return_value=[]):
        module._align_structures_interactive({})
    _assert_clean(module)


def test_redox_site_alignment_starts_from_a_clean_slate():
    module = _module_with_stale_state()
    # No detected sites in the workspace: the run returns before any prompt.
    module._align_on_redox_sites_interactive({"detected_redox_sites": []})
    _assert_clean(module)


def test_cleanup_forgets_pending_hetatm_additions():
    module = _module_with_stale_state()
    module.cleanup()
    _assert_clean(module)
