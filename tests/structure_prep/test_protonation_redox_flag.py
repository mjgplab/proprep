"""
Protonation tables mark the redox-site residues that were kept in the analysis.

Declining to exclude the redox-site residues is how you get to choose their
protonation state — a cysteine ligating an Fe-S cluster has to be set to CYM.
But the CYS table then lists it among every other cysteine with nothing to say
which ones are the ligands, so you have to carry the residue numbers in your
head from the exclusion prompt.

The flag is deliberately limited to residues still under analysis. An excluded
residue is not analyzed and never reaches a table, so there is nothing to mark;
that is what makes the column mean "you are being asked about this ligand".
"""

import pytest

from proprep.structure_prep.protonation_worker import ProtonationStateAnalyzer


class _WS:
    def __init__(self, d=None):
        self.d = dict(d or {})

    def get(self, k, default=None):
        return self.d.get(k, default)

    def set(self, k, v):
        self.d[k] = v


class _Processor:
    def __init__(self, workspace):
        self.workspace = workspace


def _worker(entries):
    w = ProtonationStateAnalyzer.__new__(ProtonationStateAnalyzer)
    w.processor = _Processor(_WS({"protonation_redox_site_residues": entries}))
    return w


def test_flags_are_keyed_by_chain_and_resid():
    w = _worker([
        {"chain": "A", "resid": 114, "site_id": "site_1"},
        {"chain": "A", "resid": 117, "site_id": "site_1"},
        {"chain": "A", "resid": 44, "site_id": "site_3"},
    ])

    flags = w._redox_site_flags()

    assert flags[("A", 114)] == "site_1"
    assert flags[("A", 44)] == "site_3"
    assert ("A", 999) not in flags


def test_string_resids_are_normalised():
    """The workspace round-trips through JSON for checklist resume."""
    w = _worker([{"chain": "A", "resid": "114", "site_id": "site_1"}])

    assert w._redox_site_flags() == {("A", 114): "site_1"}


@pytest.mark.parametrize("entries", [
    None,
    [],
    [{"chain": "A"}],                       # missing resid
    [{"chain": "A", "resid": "x"}],         # unparseable resid
    "not-a-list-of-dicts",
])
def test_malformed_state_yields_no_flags(entries):
    """A label must never break the analysis it annotates."""
    assert _worker(entries)._redox_site_flags() == {}


def test_absent_workspace_entry_is_empty():
    w = ProtonationStateAnalyzer.__new__(ProtonationStateAnalyzer)
    w.processor = _Processor(_WS({}))

    assert w._redox_site_flags() == {}


def test_no_processor_is_empty():
    w = ProtonationStateAnalyzer.__new__(ProtonationStateAnalyzer)
    w.processor = None

    assert w._redox_site_flags() == {}


# --------------------------------------------------------------------------- #
# what the module records
# --------------------------------------------------------------------------- #

from proprep.structure_prep.protonation_state_analyzer import ProtonationStateModule  # noqa: E402


class _Site:
    def __init__(self, site_id, residues):
        self.site_id = site_id
        self.site_type = "metal"
        # (chain, resid, icode) -> coords list
        self.residue_groups = {(c, r, ""): [(0.0, 0.0, 0.0)] for c, r in residues}
        self.coord_to_pdb = {(0.0, 0.0, 0.0): {"resname": "CYS"}}


def _module(sites):
    m = ProtonationStateModule.__new__(ProtonationStateModule)
    ws = _WS({"detected_redox_sites": sites})
    m.processor = _Processor(ws)
    m.get_from_workspace = lambda k, d=None: ws.get(k, d)
    m.update_workspace = ws.set
    return m, ws


def test_kept_residues_are_recorded_and_excluded_ones_are_not():
    """The 4UHX case: one cluster's cysteines kept, another's excluded."""
    sites = [
        _Site("site_1", [("A", 114), ("A", 117), ("A", 149), ("A", 151)]),
        _Site("site_3", [("A", 44), ("A", 49)]),
    ]
    m, ws = _module(sites)

    # The user excluded site_3's residues only.
    excluded = {("A", 44, ""), ("A", 49, "")}
    m._store_redox_site_flags(excluded)

    recorded = {(e["chain"], e["resid"]): e["site_id"]
                for e in ws.get("protonation_redox_site_residues")}

    assert recorded == {
        ("A", 114): "site_1", ("A", 117): "site_1",
        ("A", 149): "site_1", ("A", 151): "site_1",
    }
    assert ("A", 44) not in recorded, "an excluded residue is not analyzed, so not flagged"


def test_excluding_everything_records_nothing():
    sites = [_Site("site_1", [("A", 114), ("A", 117)])]
    m, ws = _module(sites)

    m._store_redox_site_flags({("A", 114, ""), ("A", 117, "")})

    assert ws.get("protonation_redox_site_residues") == []


def test_excluding_nothing_records_every_site_residue():
    """Pressing Enter at the exclusion prompt — the reported workflow."""
    sites = [_Site("site_1", [("A", 114), ("A", 117)])]
    m, ws = _module(sites)

    m._store_redox_site_flags(set())

    assert len(ws.get("protonation_redox_site_residues")) == 2


def test_no_redox_sites_records_nothing():
    m, ws = _module([])

    m._store_redox_site_flags(set())

    assert ws.get("protonation_redox_site_residues") == []


def test_round_trip_module_to_worker():
    """What the module writes is what the worker can read."""
    sites = [_Site("site_1", [("A", 114)])]
    m, ws = _module(sites)
    m._store_redox_site_flags(set())

    w = _worker(ws.get("protonation_redox_site_residues"))

    assert w._redox_site_flags() == {("A", 114): "site_1"}
