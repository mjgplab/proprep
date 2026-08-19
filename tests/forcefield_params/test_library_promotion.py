"""Tests for the interactive promotion wizard (forcefield_prep/library_promotion.py).

The wizard's two prompt helpers are imported into its namespace, so we replace
them with scripted stubs and drive the flow without a TTY. The library base is
redirected to a tmp dir (as in test_user_library) so nothing touches ~/.proprep.
"""

from pathlib import Path

import pytest
from rich.console import Console

from proprep.forcefield_params import loader
from proprep.forcefield_prep import library_promotion


@pytest.fixture
def user_lib(tmp_path, monkeypatch):
    base = tmp_path / "userlib"
    monkeypatch.setattr(loader, "get_user_forcefield_base_path", lambda: base)
    return base


@pytest.fixture
def params(tmp_path):
    ws = tmp_path / "ws"; ws.mkdir()
    frcmod = ws / "lig.frcmod"; frcmod.write_text("MASS\n")
    lib = ws / "lig.lib"; lib.write_text('"LIG"\n')
    return {"frcmod_file": str(frcmod), "lib_file": str(lib), "prep_file": None}


def _script(monkeypatch, *, confirms, prompts):
    """Install scripted stubs for the two prompt helpers.

    confirms / prompts are lists consumed in call order. A prompt value of None
    means 'return the supplied default' (simulating the user pressing Enter).
    """
    cq = list(confirms)
    pq = list(prompts)

    def fake_confirm(processor, prompt, default=None, **kw):
        return cq.pop(0)

    def fake_prompt(processor, prompt, choices=None, default=None, **kw):
        val = pq.pop(0)
        return default if val is None else val

    monkeypatch.setattr(library_promotion, "confirm_with_context", fake_confirm)
    monkeypatch.setattr(library_promotion, "prompt_with_context", fake_prompt)


def _browse(monkeypatch, selections):
    """Script the shared file browser.

    _prompt_existing_file imports file_browser inside the function, so the
    patch has to land on the source module. Values are returned in call order;
    None stands for the user skipping an optional file.
    """
    from proprep.utils import file_browser as fb

    queue = list(selections)

    def fake_browser(*args, **kwargs):
        value = queue.pop(0)
        return fb.SKIP if value is None and kwargs.get("optional") else value

    monkeypatch.setattr(fb, "file_browser", fake_browser)


def test_offer_promotion_small_molecule_accept(user_lib, params, monkeypatch):
    # confirm: yes to "add to library"
    # prompts: type_name(default), description(default), set_name(default),
    #          references(default empty)
    _script(monkeypatch, confirms=[True],
            prompts=[None, None, None, None])

    result = library_promotion.offer_promotion(
        Console(quiet=True), processor=None,
        category="small_molecule", residue_name="LIG",
        parameter_files=params,
    )
    assert result is not None
    assert result["set_name"] == "user_LIG"

    sets = loader.discover_forcefield_files("small_molecules/LIG",
                                            "single_state", "default")
    assert [s["name"] for s in sets] == ["user_LIG"]
    # GAFF → no new atom types declared.
    meta = loader.load_forcefield_metadata("small_molecules/LIG")
    spin = meta["redox_states"]["single_state"]["spin_states"]["default"]
    assert spin["atom_types"] == []
    assert spin["residue_name"] == "LIG"


def test_offer_promotion_decline_writes_nothing(user_lib, params, monkeypatch):
    _script(monkeypatch, confirms=[False], prompts=[])
    result = library_promotion.offer_promotion(
        Console(quiet=True), processor=None,
        category="small_molecule", residue_name="LIG",
        parameter_files=params,
    )
    assert result is None
    assert not (user_lib / "small_molecules").exists()


def test_offer_promotion_no_frcmod_returns_none(user_lib, monkeypatch):
    _script(monkeypatch, confirms=[], prompts=[])  # never prompted
    result = library_promotion.offer_promotion(
        Console(quiet=True), processor=None,
        category="small_molecule", residue_name="LIG",
        parameter_files={"frcmod_file": None, "lib_file": None},
    )
    assert result is None


def test_offer_promotion_collision_version_bump(user_lib, params, monkeypatch):
    # First promotion accepted with defaults.
    _script(monkeypatch, confirms=[True], prompts=[None, None, None, None])
    library_promotion.offer_promotion(
        Console(quiet=True), processor=None,
        category="small_molecule", residue_name="LIG", parameter_files=params)

    # Second promotion: same defaults → collision → choose "1" (new version).
    _script(monkeypatch, confirms=[True],
            prompts=[None, None, None, None, "1"])
    result = library_promotion.offer_promotion(
        Console(quiet=True), processor=None,
        category="small_molecule", residue_name="LIG", parameter_files=params)
    assert result["set_name"] == "user_LIG_v2"

    sets = loader.discover_forcefield_files("small_molecules/LIG",
                                            "single_state", "default")
    assert {s["name"] for s in sets} == {"user_LIG", "user_LIG_v2"}


def test_metal_site_records_atom_types_and_state(user_lib, tmp_path, monkeypatch):
    ws = tmp_path / "ws2"; ws.mkdir()
    frcmod = ws / "site.frcmod"; frcmod.write_text("MASS\n")
    lib = ws / "site.lib"; lib.write_text('"FE1"\n')

    # type_name, description, set_name, references, redox_state, spin_state,
    # formal_charge → all defaults except we provide an explicit redox label.
    _script(monkeypatch, confirms=[True],
            prompts=[None, None, None, None, "oxidized", None, "-2"])

    result = library_promotion.offer_promotion(
        Console(quiet=True), processor=None,
        category="metal_site", residue_name="FE1",
        parameter_files={"frcmod_file": str(frcmod), "lib_file": str(lib)},
        atom_types=['{ "FE" "Fe" "sp3" }'],
    )
    assert result is not None
    meta = loader.load_forcefield_metadata("metal_sites/FE1")
    spin = meta["redox_states"]["oxidized"]["spin_states"]["default"]
    assert spin["atom_types"] == ['{ "FE" "Fe" "sp3" }']
    assert meta["redox_states"]["oxidized"]["formal_charge"] == -2


# --------------------------------------------------------------------------- #
# lib locator
# --------------------------------------------------------------------------- #

def test_find_lib_artifact_prefers_explicit(tmp_path):
    assert library_promotion.find_lib_artifact(
        tmp_path, parameter_files={"lib_file": "/x/y.lib"}) == "/x/y.lib"


def test_find_lib_artifact_stem_match_then_first(tmp_path):
    (tmp_path / "AAA.lib").write_text("a")
    (tmp_path / "LIG.lib").write_text("b")
    # stem match wins
    assert library_promotion.find_lib_artifact(
        tmp_path, residue_name="LIG").endswith("LIG.lib")
    # no match → first alphabetically
    assert library_promotion.find_lib_artifact(
        tmp_path, residue_name="ZZZ").endswith("AAA.lib")


def test_find_lib_artifact_none_when_empty(tmp_path):
    assert library_promotion.find_lib_artifact(tmp_path) is None


# --------------------------------------------------------------------------- #
# standalone import wizard
# --------------------------------------------------------------------------- #

def test_import_wizard_registers_loose_files(user_lib, tmp_path, monkeypatch):
    src = tmp_path / "loose"; src.mkdir()
    frcmod = src / "myligand.frcmod"; frcmod.write_text("MASS\n")
    lib = src / "myligand.lib"; lib.write_text('"MYL"\n')

    # File selection goes through the shared browser now, so the paths are
    # scripted there rather than as typed prompt answers.
    _browse(monkeypatch, [str(frcmod), str(lib), None])

    # prompts: category "1", type_name(default=lib stem), description,
    #          set_name, references
    _script(monkeypatch, confirms=[], prompts=["1", None, None, None, None])

    result = library_promotion.run_import_wizard(Console(quiet=True), processor=None)
    assert result is not None
    # default type_name is the lib stem "myligand"
    sets = loader.discover_forcefield_files("small_molecules/myligand",
                                            "single_state", "default")
    assert len(sets) == 1
    assert sets[0]["frcmod"].endswith("myligand.frcmod")


# --------------------------------------------------------------------------- #
# file selection goes through the shared browser
# --------------------------------------------------------------------------- #

def _capture_browser(monkeypatch, selections):
    """Script the browser and record the kwargs each call received."""
    from proprep.utils import file_browser as fb

    calls = []
    queue = list(selections)

    def fake_browser(*args, **kwargs):
        calls.append(kwargs)
        value = queue.pop(0)
        return fb.SKIP if value is None and kwargs.get("optional") else value

    monkeypatch.setattr(fb, "file_browser", fake_browser)
    return calls


def test_the_import_wizard_browses_rather_than_asking_for_a_path(
        user_lib, tmp_path, monkeypatch):
    """
    Importing means pointing at files produced elsewhere -- a prior project, a
    collaborator, a published set -- which are the paths a user is least likely
    to have memorised. It used to ask for each as typed text.
    """
    src = tmp_path / "loose"; src.mkdir()
    (src / "fad.frcmod").write_text("MASS\n")
    (src / "fad.lib").write_text('"FAD"\n')

    calls = _capture_browser(monkeypatch,
                             [str(src / "fad.frcmod"), str(src / "fad.lib"), None])
    _script(monkeypatch, confirms=[], prompts=["1", None, None, None, None])

    assert library_promotion.run_import_wizard(
        Console(quiet=True), processor=None) is not None
    assert len(calls) == 3


def test_each_file_type_is_filtered_to_its_extensions(user_lib, tmp_path, monkeypatch):
    src = tmp_path / "loose"; src.mkdir()
    (src / "fad.frcmod").write_text("MASS\n")
    (src / "fad.lib").write_text('"FAD"\n')

    calls = _capture_browser(monkeypatch,
                             [str(src / "fad.frcmod"), str(src / "fad.lib"), None])
    _script(monkeypatch, confirms=[], prompts=["1", None, None, None, None])
    library_promotion.run_import_wizard(Console(quiet=True), processor=None)

    assert calls[0]["extensions"] == [".frcmod"]
    assert calls[1]["extensions"] == [".lib", ".off"]
    assert calls[2]["extensions"] == [".prep", ".prepi", ".prepin"]


def test_the_companion_files_start_beside_the_frcmod(user_lib, tmp_path, monkeypatch):
    """A .lib almost always sits next to its .frcmod; don't send the user back."""
    src = tmp_path / "elsewhere"; src.mkdir()
    (src / "fad.frcmod").write_text("MASS\n")
    (src / "fad.lib").write_text('"FAD"\n')

    calls = _capture_browser(monkeypatch,
                             [str(src / "fad.frcmod"), str(src / "fad.lib"), None])
    _script(monkeypatch, confirms=[], prompts=["1", None, None, None, None])
    library_promotion.run_import_wizard(Console(quiet=True), processor=None)

    assert calls[1]["directory"] == str(src)
    assert calls[2]["directory"] == str(src)


def test_only_the_prep_file_is_optional(user_lib, tmp_path, monkeypatch):
    src = tmp_path / "loose"; src.mkdir()
    (src / "fad.frcmod").write_text("MASS\n")
    (src / "fad.lib").write_text('"FAD"\n')

    calls = _capture_browser(monkeypatch,
                             [str(src / "fad.frcmod"), str(src / "fad.lib"), None])
    _script(monkeypatch, confirms=[], prompts=["1", None, None, None, None])
    library_promotion.run_import_wizard(Console(quiet=True), processor=None)

    assert calls[0]["optional"] is False
    assert calls[1]["optional"] is False
    assert calls[2]["optional"] is True


def test_a_typed_path_is_still_accepted(user_lib, tmp_path, monkeypatch):
    """allow_path_jump keeps the old behaviour available for scripted runs."""
    src = tmp_path / "loose"; src.mkdir()
    (src / "fad.frcmod").write_text("MASS\n")
    (src / "fad.lib").write_text('"FAD"\n')

    calls = _capture_browser(monkeypatch,
                             [str(src / "fad.frcmod"), str(src / "fad.lib"), None])
    _script(monkeypatch, confirms=[], prompts=["1", None, None, None, None])
    library_promotion.run_import_wizard(Console(quiet=True), processor=None)

    assert all(c["allow_path_jump"] for c in calls)


def test_cancelling_the_frcmod_abandons_the_import(user_lib, monkeypatch):
    """Returning None from a required selection must not proceed."""
    _capture_browser(monkeypatch, [None])
    _script(monkeypatch, confirms=[False], prompts=[])

    assert library_promotion.run_import_wizard(
        Console(quiet=True), processor=None) is None


def test_skipping_the_prep_file_is_not_an_error(user_lib, tmp_path, monkeypatch):
    src = tmp_path / "loose"; src.mkdir()
    (src / "fad.frcmod").write_text("MASS\n")
    (src / "fad.lib").write_text('"FAD"\n')

    _capture_browser(monkeypatch,
                     [str(src / "fad.frcmod"), str(src / "fad.lib"), None])
    _script(monkeypatch, confirms=[], prompts=["1", None, None, None, None])

    result = library_promotion.run_import_wizard(Console(quiet=True), processor=None)

    assert result is not None


# --------------------------------------------------------------------------- #
# the category prompt has to show its options
# --------------------------------------------------------------------------- #

def test_the_category_options_are_displayed(monkeypatch):
    """
    It read as a bare "What kind of parameters are these? (1):" with nothing to
    go on. options_map only feeds the session recorder, and passing it
    SUPPRESSES the inline choice list -- show_choices defaults to False on the
    assumption the options are printed separately, which they were not.
    """
    console = Console(record=True, width=100)
    monkeypatch.setattr(library_promotion, "prompt_with_context",
                        lambda *a, **k: "1")

    library_promotion._prompt_category(console, None)
    text = console.export_text()

    for label in ("Small molecule", "Modified amino acid", "Metal site"):
        assert label in text


def test_each_option_says_where_it_lands():
    """The choice decides the library family, so make that visible."""
    families = {v[2] for v in library_promotion._CATEGORY_CHOICES.values()}

    assert families == {"small_molecules", "modified_aa", "metal_sites"}


def test_the_displayed_families_match_the_real_mapping(monkeypatch):
    """A drifting legend would be worse than none."""
    for key, (_label, _meaning, family) in library_promotion._CATEGORY_CHOICES.items():
        monkeypatch.setattr(library_promotion, "prompt_with_context",
                            lambda *a, _k=key, **kw: _k)
        category = library_promotion._prompt_category(Console(quiet=True), None)
        assert library_promotion._CATEGORY_FAMILY[category] == family


@pytest.mark.parametrize("choice,expected", [
    ("1", "small_molecule"), ("2", "modified_amino_acid"), ("3", "metal_site"),
])
def test_each_choice_maps_to_its_category(monkeypatch, choice, expected):
    monkeypatch.setattr(library_promotion, "prompt_with_context",
                        lambda *a, **k: choice)

    assert library_promotion._prompt_category(Console(quiet=True), None) == expected
