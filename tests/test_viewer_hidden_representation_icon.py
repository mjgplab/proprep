"""A hidden representation keeps an icon on its visibility button.

Turning a representation off replaced the eye with an em dash in #444 on the
row's #353535, about 1.3:1: the button was still there, and there was nothing
to see of it. Both states are now drawn, an open eye and a crossed-out eye, so
they differ in shape as well as in colour and neither depends on how a system
renders an emoji. Checked on the real page in headless Chrome when written:
white when shown, amber with a slash when hidden, the button 26 px wide in
both, and a representation that starts hidden (ions) drawn the same way.

As for the other viewer features, the shipped functions are extracted from the
template and run under Node.
"""

import json
import re
import shutil
import subprocess
from pathlib import Path

import pytest

TEMPLATE = (Path(__file__).parent.parent / "src" / "proprep" / "structure_prep" / "templates" / "ngl_viewer.html")
needs_node = pytest.mark.skipif(shutil.which("node") is None, reason="Node not installed")


@pytest.fixture(scope="module")
def template_text():
    return TEMPLATE.read_text()


def _luminance(hex_colour):
    channels = [int(hex_colour[i:i + 2], 16) / 255 for i in (1, 3, 5)]
    linear = [c / 12.92 if c <= 0.03928 else ((c + 0.055) / 1.055) ** 2.4 for c in channels]
    return 0.2126 * linear[0] + 0.7152 * linear[1] + 0.0722 * linear[2]


def _contrast(a, b):
    hi, lo = sorted((_luminance(a), _luminance(b)), reverse=True)
    return (hi + 0.05) / (lo + 0.05)


def test_no_state_of_the_button_is_an_em_dash_or_an_emoji(template_text):
    row = template_text[template_text.index("function buildRepRow("):template_text.index("function buildStyleOptions(")]
    toggle = template_text[template_text.index("function toggleRepVisibility("):template_text.index("function eyeIcon(")]
    for source in (row, toggle):
        assert "eyeIcon(" in source
        assert "\\u2014" not in source and "\\u{1F441}" not in source


def test_both_states_stand_out_from_the_row(template_text):
    row_background = re.search(r"\.rep-row \{\s*background: (#[0-9a-fA-F]{6})", template_text).group(1)
    shown = re.search(r"\.icon-btn\.eye-btn \{ color: (#[0-9a-fA-F]{6})", template_text).group(1)
    hidden = re.search(r"\.icon-btn\.eye-btn\.eye-hidden \{ color: (#[0-9a-fA-F]{6})", template_text).group(1)
    # WCAG asks 3:1 of a graphical control; the old hidden state had 1.3:1.
    assert _contrast(shown, row_background) > 7 and _contrast(hidden, row_background) > 7
    assert shown.lower() != hidden.lower()


@needs_node
def test_the_hidden_icon_is_the_shown_one_with_a_slash_through_it(template_text):
    functions = re.search(r"( {8}function eyeIcon\(.*?\n {8}\}\n).*?( {8}function eyeTitle\(.*?\n {8}\}\n)", template_text, re.S)
    script = functions.group(1) + functions.group(2) + \
        "console.log(JSON.stringify({shown: eyeIcon(true), hidden: eyeIcon(false), t1: eyeTitle(true), t0: eyeTitle(false)}));"
    out = json.loads(subprocess.run(["node", "-e", script], capture_output=True, text=True, check=True).stdout)
    assert out["shown"].startswith("<svg") and out["hidden"].startswith("<svg")
    assert "<line" in out["hidden"] and "<line" not in out["shown"]
    assert out["hidden"].replace(re.search(r"<line[^>]*/>", out["hidden"]).group(0), "") == out["shown"]
    assert "currentColor" in out["shown"]                       # takes the button's colour
    assert out["t1"].startswith("Shown") and out["t0"].startswith("Hidden")
