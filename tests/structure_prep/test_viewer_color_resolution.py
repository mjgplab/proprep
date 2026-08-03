"""Regression tests for ViewerCoordinator colour resolution.

NGL reads a representation ``color`` string first as a colour-*scheme* id
(``element``, ``resname``, ...). A bare CSS colour word matches no scheme and
is silently dropped, so the representation never renders. This bit the PROPKA
determinant viewer: the target residue, highlighted with the named "magenta",
was invisible while the hex-coloured driver residues showed fine.

``_resolve_color`` must therefore translate named colours to hex while leaving
scheme ids and hex values untouched.
"""

from proprep.structure_prep.viewer_coordinator import _resolve_color


def test_named_colors_become_hex():
    assert _resolve_color("magenta") == "#ff00ff"
    assert _resolve_color("orange") == "#ff8c00"
    assert _resolve_color("red") == "#ff0000"


def test_named_color_is_case_insensitive():
    assert _resolve_color("MAGENTA") == "#ff00ff"
    assert _resolve_color("Orange") == "#ff8c00"


def test_ngl_scheme_ids_pass_through():
    # These are NGL colour schemes, NOT colours — must not be rewritten.
    for scheme in ("element", "resname", "chainid", "bfactor", "sstruc"):
        assert _resolve_color(scheme) == scheme


def test_hex_passes_through():
    assert _resolve_color("#15803d") == "#15803d"
    assert _resolve_color("#ff00ff") == "#ff00ff"


def test_palette_syntax_resolves_to_hex():
    out = _resolve_color("palette:3")
    assert out.startswith("#") and len(out) == 7


def test_malformed_palette_falls_through():
    assert _resolve_color("palette:notanint") == "palette:notanint"
