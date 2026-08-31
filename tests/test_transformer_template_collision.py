"""Saving a transformer must not silently replace another one.

The save is a plain open(path, "w"), and _sanitize collapses case and
punctuation -- "GDP ff94" and "gdp-ff94" become the same file -- so reusing a
name without noticing destroyed the earlier recipe.
"""

import json
from unittest.mock import patch

from rich.console import Console

import proprep.redoxsite_prep.transformation.table_transformer_creator as ttc
from proprep.redoxsite_prep.transformation.table_transformer_creator import (
    TableTransformerCreator,
)


def _creator(tmp_path, monkeypatch, console=None):
    monkeypatch.setattr(ttc, "DEFAULT_USER_TRANSFORMER_DIR", tmp_path)
    creator = TableTransformerCreator.__new__(TableTransformerCreator)
    creator.processor = None
    creator.console = console or Console(record=True, width=100, force_terminal=False)
    return creator


def _existing(tmp_path, token, description="an earlier recipe"):
    path = tmp_path / f"{token}.json"
    path.write_text(json.dumps({"name": token, "description": description}))
    return path


def test_a_free_name_is_taken_as_is(tmp_path, monkeypatch):
    creator = _creator(tmp_path, monkeypatch)
    assert creator._resolve_template_collision("GDP ff94") == ("GDP ff94", "gdp_ff94")


def test_a_taken_name_is_reported_before_anything_is_written(tmp_path, monkeypatch):
    console = Console(record=True, width=100, force_terminal=False)
    creator = _creator(tmp_path, monkeypatch, console)
    _existing(tmp_path, "gdp_ff94")
    with patch.object(ttc, "confirm_with_context", return_value=True):
        creator._resolve_template_collision("GDP ff94")
    out = console.export_text()
    assert "already" in out
    assert "an earlier recipe" in out


def test_punctuation_variants_are_recognised_as_the_same_file(tmp_path, monkeypatch):
    """'GDP ff94' and 'gdp-ff94' both sanitize to gdp_ff94."""
    creator = _creator(tmp_path, monkeypatch)
    _existing(tmp_path, "gdp_ff94")
    with patch.object(ttc, "confirm_with_context", return_value=True) as confirm:
        creator._resolve_template_collision("gdp-ff94")
    assert confirm.called


def test_declining_the_replacement_lets_a_new_name_be_given(tmp_path, monkeypatch):
    creator = _creator(tmp_path, monkeypatch)
    _existing(tmp_path, "gdp_ff94")
    with patch.object(ttc, "confirm_with_context", return_value=False), \
         patch.object(ttc, "prompt_with_context", return_value="gdp ff94 v2"):
        name, token = creator._resolve_template_collision("GDP ff94")
    assert (name, token) == ("gdp ff94 v2", "gdp_ff94_v2")


def test_an_empty_new_name_aborts_the_save(tmp_path, monkeypatch):
    creator = _creator(tmp_path, monkeypatch)
    _existing(tmp_path, "gdp_ff94")
    with patch.object(ttc, "confirm_with_context", return_value=False), \
         patch.object(ttc, "prompt_with_context", return_value="  "):
        _, token = creator._resolve_template_collision("GDP ff94")
    assert token is None


def test_a_malformed_existing_file_still_blocks_the_write(tmp_path, monkeypatch):
    creator = _creator(tmp_path, monkeypatch)
    (tmp_path / "gdp_ff94.json").write_text("{ not json")
    with patch.object(ttc, "confirm_with_context", return_value=True) as confirm:
        creator._resolve_template_collision("GDP ff94")
    assert confirm.called
