"""Redox Site Detector: the viewer shows the structure detection runs on.

Detection receives the structure the user picked in the structure selector
(e.g. "H-Stripped"), but the viewer helpers chose their own file from a
fixed key list that cannot return ``hstripped_pdb_file``. On 6R2Q, loaded
from the PDB with its 11567 hydrogens stripped, the sites were detected on
the stripped atoms and displayed on the deposited file, hydrogens and all.

Run with: pytest tests/test_redox_detector_viewer_structure.py
"""

import io
import sys
from pathlib import Path

from rich.console import Console

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from proprep.structure_prep.comprehensive_redox_detector import (  # noqa: E402
    ComprehensiveRedoxDetector,
    _find_workspace_structure,
)


class _Workspace(dict):
    def set(self, key, value):
        self[key] = value


class _Processor:
    def __init__(self, workspace):
        self._workspace = workspace

    def _get_workspace(self):
        return self._workspace


def _loaded_and_stripped(tmp_path):
    deposited = tmp_path / "6R2Q.pdb"
    stripped = tmp_path / "6R2Q_noH.pdb"
    deposited.write_text("END\n")
    stripped.write_text("END\n")
    workspace = _Workspace(rcsb_pdb_file=str(deposited), hstripped_pdb_file=str(stripped))
    return _Processor(workspace), deposited, stripped


def _detector(processor):
    return ComprehensiveRedoxDetector(console=Console(file=io.StringIO()), processor=processor)


def test_viewer_shows_the_structure_detection_runs_on(tmp_path):
    processor, _, stripped = _loaded_and_stripped(tmp_path)
    detector = _detector(processor)
    detector.source_pdb_file = str(stripped)
    detector._record_viewer_structure()
    assert _find_workspace_structure(processor) == str(stripped)


def test_a_new_detection_replaces_the_previous_choice(tmp_path):
    processor, deposited, stripped = _loaded_and_stripped(tmp_path)
    first = _detector(processor)
    first.source_pdb_file = str(stripped)
    first._record_viewer_structure()
    second = _detector(processor)
    second.source_pdb_file = str(deposited)
    second._record_viewer_structure()
    assert _find_workspace_structure(processor) == str(deposited)


def test_detection_without_a_file_falls_back_to_the_loaded_structure(tmp_path):
    """``source_pdb_file`` can be a placeholder ("structure.pdb") or unset."""
    processor, deposited, stripped = _loaded_and_stripped(tmp_path)
    first = _detector(processor)
    first.source_pdb_file = str(stripped)
    first._record_viewer_structure()
    for source in (None, "structure.pdb"):
        detector = _detector(processor)
        detector.source_pdb_file = source
        detector._record_viewer_structure()
        assert _find_workspace_structure(processor) == str(deposited), source
