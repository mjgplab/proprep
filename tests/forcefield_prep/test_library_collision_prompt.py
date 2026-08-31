"""The library-collision prompt has to show what it is asking about.

``options_map`` is session-recorder metadata: it never renders, and passing it
suppresses Rich's inline choice list. A prompt that relies on it alone reaches
the user as a bare "How would you like to proceed? (3):".
"""

from unittest.mock import patch

import pytest
from rich.console import Console

from proprep.forcefield_prep import library_promotion as lp
from proprep.forcefield_params import LibraryCollisionError, PromotionRequest
from proprep.forcefield_params import user_library

COLLISION_MSG = (
    "Set 'meagher_gdp_ff94' already exists for small_molecules/gdp "
    "[single_state/default]. Choose overwrite or version_bump."
)


def _request():
    return PromotionRequest(
        family="small_molecules", type_name="gdp",
        set_name="meagher_gdp_ff94", frcmod_src="x.frcmod", lib_srcs=["x.lib"],
    )


def _run(answer, *, raise_times=1):
    """Drive the retry loop with a canned answer; return (result, output, request)."""
    console = Console(record=True, width=100, force_terminal=False)
    request = _request()
    calls = {"n": 0}

    def promote(req):
        calls["n"] += 1
        if calls["n"] <= raise_times:
            raise LibraryCollisionError(COLLISION_MSG)
        return {"library_path": "/tmp/lib", "metadata_path": "/tmp/lib/metadata.json",
                "state_dir": "/tmp/lib/single_state/default",
                "set_name": req.set_name, "copied_files": []}

    with patch.object(lp, "promote_state", side_effect=promote), \
         patch.object(lp, "prompt_with_context", return_value=answer):
        result = lp._promote_with_collision_retry(console, None, request)
    return result, console.export_text(), request


def test_every_option_is_printed_with_its_key():
    _, out, _ = _run("3")
    assert "1. Save as a new version" in out
    assert "2. Overwrite the existing set" in out
    assert "3. Cancel" in out


def test_the_version_option_names_the_set_it_would_write():
    _, out, _ = _run("3")
    assert "meagher_gdp_ff94_v2" in out


def test_state_coordinates_survive_rich_markup():
    """'[single_state/default]' must not be parsed as a style tag and eaten."""
    _, out, _ = _run("3")
    assert "[single_state/default]" in out
    assert "gdp ." not in out          # the stray space left when it vanished


@pytest.mark.parametrize("answer,expected", [
    ("1", user_library.ON_COLLISION_VERSION),
    ("2", user_library.ON_COLLISION_OVERWRITE),
])
def test_choice_sets_the_collision_policy_and_retries(answer, expected):
    result, _, request = _run(answer)
    assert request.on_collision == expected
    assert result is not None


def test_cancel_writes_nothing():
    result, out, request = _run("3")
    assert result is None
    assert "cancelled" in out.lower()
    assert request.on_collision == user_library.ON_COLLISION_ERROR
