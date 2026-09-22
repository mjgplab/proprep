"""The session log says which of its numbers is ProPrep's version.

A log began with "version": "1.1" and ended with "proprep_version": "1.21.0".
The first is the version of the file format; under a bare "version" it read as
ProPrep's, and was wrong. Nothing reads the key, so it is renamed; a log from
before the rename still replays.
"""

import json

from proprep.utils.session_recorder import SessionRecorder


def test_a_new_log_names_the_format_version_and_has_no_bare_version(tmp_path):
    recorder = SessionRecorder(str(tmp_path / "s.json"))
    recorder.start_recording({"proprep_version": "1.22.0"})
    recorder.stop_recording()
    data = json.loads((tmp_path / "s.json").read_text())
    assert data["session_format_version"] == "1.1"
    assert "version" not in data
    assert data["metadata"]["proprep_version"] == "1.22.0"
    assert list(data)[0] == "session_format_version"          # still the first block


def test_a_log_from_before_the_rename_still_replays(tmp_path):
    from proprep.utils.session_recorder import SessionReplayer
    old = {"version": "1.1", "start_time": "t", "end_time": "t", "metadata": {},
           "interactions": [{"index": 0, "timestamp": "t", "type": "prompt", "prompt": "Enter your choice",
                             "response": "1", "choices": ["1", "2"], "context": {}}]}
    path = tmp_path / "old.json"
    path.write_text(json.dumps(old))
    replayer = SessionReplayer(str(path))          # loads in the constructor
    replayer.start_replay()
    assert replayer.session_data["interactions"]
    assert replayer.get_next_response("prompt", "Enter your choice") == "1"
