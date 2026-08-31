"""
Titration inputs: constant pH (cpin), constant Redox potential (cein), and the
joint proton/electron format (cpein).

AMBER ships three generators for the "which residues titrate, in what states,
with what reference energies" input file::

    cpinutil.py   -> cpin    residues with typ == "ph"
    ceinutil.py   -> cein    residues with typ == "redox"      (HEH only)
    cpeinutil.py  -> cpein   ALL titratable residues, unfiltered

They share an argparse skeleton, the same terminus filter, and the same
``TitratableResidueList.write_cpin()`` writer, differing in a mode string --
``"ph"``, ``"redox"``, ``"phredox"`` -- that doubles as the residue filter
above.

**cpin and cein are complementary, not alternatives.** A protein whose heme
titrates in redox while its propionates and carboxylates titrate in pH gets one
of each, partitioned by residue, and runs with both flag triples at once. That
is the documented workflow -- Amber tutorial 33 (C(pH,E)MD on microperoxidase
8) generates::

    ceinutil.py -resnames HEH     -p mp8_is.prmtop -igb 2 -o mp8_is.cein
    cpinutil.py -resnames PRN GL4 -p mp8_is.prmtop -igb 2 -o mp8_is.cpin

and runs production with::

    -cpin mp8_is.cpin -cpout ... -cprestrt ...
    -cein mp8_is.cein -ceout ... -cerestrt ...

setting ``icnstph``/``ntcnstph``/``solvph`` and ``icnste``/``ntcnste``/``solve``
together in one &cntrl. sander compares ``cpin_specified`` and
``cein_specified`` nowhere; the only exclusion in ``mdfil.F90`` (line 586) is
one-directional and guards ``cpein`` alone:

    cpein + cpin  -> error
    cpein + cein  -> error
    cpin  + cein  -> fine, and is the intended combination

**What cpein is for** follows from the filters rather than from any prose:
``cpinutil`` drops everything that is not ``typ == "ph"`` and ``ceinutil``
drops everything that is not ``typ == "redox"``, so a residue typed
``"phredox"`` fits in neither file. Such a residue titrates in proton *and*
electron on the same site, with a joint state manifold (``TYX``: four states,
protonated/deprotonated x oxidized/reduced). Proton-coupled electron transfer
on one site cannot be split across two files, so cpein carries everything in a
single list. ``TYX`` is currently the only ``phredox`` residue, and ProPrep
does not emit it.

This module is the one place these differences are written down, so the
topology generator (which builds the files) and the MD manager (which hands
them to the engine) cannot disagree about a flag name.
"""

from dataclasses import dataclass
from typing import Dict, FrozenSet, Iterable, List, Optional, Sequence, Tuple

__all__ = [
    "TitrationMode",
    "PH", "REDOX", "PHREDOX",
    "TITRATION_MODES",
    "PH_RESIDUE_PKA",
    "REDOX_RESIDUE_EO",
    "PHREDOX_RESIDUE_NAMES",
    "get_mode",
    "modes_for_residue_names",
    "partition_residues_by_mode",
    "engine_flags_for_modes",
    "mdin_keyword_sets",
    "parmed_titratable_atom_count",
    "parmed_residue_type",
]


# ---------------------------------------------------------------------------
# Mode keys
# ---------------------------------------------------------------------------

PH = "ph"
REDOX = "redox"
PHREDOX = "phredox"


# ---------------------------------------------------------------------------
# Residue reference data
# ---------------------------------------------------------------------------

# pKa values as printed by `cpinutil.py --describe`. These are the residue-TYPE
# defaults; per-site values computed by PB Titrate override them for display.
PH_RESIDUE_PKA: Dict[str, float] = {
    'AS4': 4.0,   # Titratable aspartate (NOT ASP!)
    'GL4': 4.4,   # Titratable glutamate (NOT GLU!)
    'HIP': 6.6,   # Titratable histidine (NOT HIE/HID!)
    'LYS': 10.4,  # Lysine
    'CYS': 8.5,   # Cysteine (free thiol)
    'TYR': 9.6,   # Tyrosine
    'PRN': 4.8,   # Propionate (for heme groups)
}

# Standard reduction potentials in Volts, as printed by `ceinutil.py --describe`.
# HEH is currently the ONLY residue in ParmEd with typ == "redox": the bis-His
# c-type heme of conste.lib, an 87-atom construct spanning the porphyrin plus
# both thioether Cys and both axial His side chains.
REDOX_RESIDUE_EO: Dict[str, float] = {
    'HEH': -0.203,
}

# Residues that titrate in proton AND electron on the same site (ParmEd
# typ == "phredox"). Neither cpinutil nor ceinutil will accept them, so their
# presence forces the joint cpein format.
PHREDOX_RESIDUE_NAMES: Tuple[str, ...] = ('TYX',)


@dataclass(frozen=True)
class TitrationMode:
    """Everything that differs between cpin / cein / cpein generation."""

    key: str
    label: str                  # human-readable, e.g. "Constant pH"
    short_label: str            # e.g. "CpHMD"
    utility: str                # e.g. "cpinutil.py"
    file_ext: str               # e.g. "cpin"

    # sander/pmemd command-line flag triple (see mdfil.F90:290-341).
    flag_input: str             # e.g. "-cpin"
    flag_output: str            # e.g. "-cpout"
    flag_restart: str           # e.g. "-cprestrt"
    ext_output: str             # e.g. "cpout"
    ext_restart: str            # e.g. "cprestrt"

    # `-op`: write a radii-modified prmtop for explicit solvent. ceinutil.py
    # imports changeRadii/change but never calls them and exposes no -op flag.
    # In a combined pH+redox run this does not matter -- cpinutil writes the
    # modified prmtop and both files share it, exactly as tutorial 33 does
    # (mp8_es.new.prmtop). It matters only for a redox-only explicit run.
    supports_output_prmtop: bool

    # Range-filter flag pair, or None if the utility has no such filter.
    range_flag_min: Optional[str]
    range_flag_max: Optional[str]
    range_label: str            # "pKa" / "Eo (V)"

    # mdin keywords this file drives. In a combined run BOTH sets are set in
    # the same &cntrl (tutorial 33 section 3).
    mdin_flag_keyword: str      # "icnstph" / "icnste"
    mdin_setpoint_keyword: str  # "solvph" / "solve"
    mdin_freq_keyword: str      # "ntcnstph" / "ntcnste"
    mdin_relax_keyword: str     # "ntrelax" / "ntrelaxe"

    # Residue names this utility will accept.
    residue_names: FrozenSet[str]

    # igb values for which the shipped reference energies are defined, split by
    # solvent treatment. Selecting an igb outside these gives `None` reference
    # energies and a meaningless output file.
    allowed_igb_implicit: FrozenSet[int]
    allowed_igb_explicit: FrozenSet[int]

    # Internal dielectric values with defined reference energies. All three
    # utilities hard-reject anything other than 1 or 2; on top of that, HEH's
    # dielc2 energies are registered with no arguments (all None), so intdiel=2
    # is unusable for any file that can contain HEH.
    allowed_intdiel: FrozenSet[float]

    def output_names(self, prefix: str) -> Tuple[str, str]:
        """(output_file, restart_file) basenames for a step named `prefix`."""
        return f"{prefix}.{self.ext_output}", f"{prefix}.{self.ext_restart}"

    def engine_flags(self, input_name: str, prefix: str) -> List[str]:
        """This file's three sander/pmemd flags."""
        out_name, restrt_name = self.output_names(prefix)
        return [
            self.flag_input, input_name,
            self.flag_output, out_name,
            self.flag_restart, restrt_name,
        ]

    def allowed_igb(self, sim_type: str) -> FrozenSet[int]:
        return (self.allowed_igb_explicit if sim_type == 'explicit'
                else self.allowed_igb_implicit)

    @property
    def titrates_ph(self) -> bool:
        return self.key in (PH, PHREDOX)

    @property
    def titrates_redox(self) -> bool:
        return self.key in (REDOX, PHREDOX)


# All five igb models cpinutil accepts. The shipped pH reference energies cover
# every one of them for the standard residues.
_ALL_IGB = frozenset({1, 2, 5, 7, 8})

# HEH's reference energies (parmed.amber.titratable_residues) are registered
# for igb 2/5/7/8 implicit but only 2/5/7 explicit, and NOTHING for igb=1 --
# the reduced state's _ReferenceEnergy is built without an igb1 argument.
_HEH_IGB_IMPLICIT = frozenset({2, 5, 7, 8})
_HEH_IGB_EXPLICIT = frozenset({2, 5, 7})


_PH_MODE = TitrationMode(
    key=PH,
    label="Constant pH",
    short_label="CpHMD",
    utility="cpinutil.py",
    file_ext="cpin",
    flag_input="-cpin",
    flag_output="-cpout",
    flag_restart="-cprestrt",
    ext_output="cpout",
    ext_restart="cprestrt",
    supports_output_prmtop=True,
    range_flag_min="-minpKa",
    range_flag_max="-maxpKa",
    range_label="pKa",
    mdin_flag_keyword="icnstph",
    mdin_setpoint_keyword="solvph",
    mdin_freq_keyword="ntcnstph",
    mdin_relax_keyword="ntrelax",
    residue_names=frozenset(PH_RESIDUE_PKA),
    allowed_igb_implicit=_ALL_IGB,
    allowed_igb_explicit=_ALL_IGB,
    allowed_intdiel=frozenset({1.0, 2.0}),
)

_REDOX_MODE = TitrationMode(
    key=REDOX,
    label="Constant Redox potential",
    short_label="CEMD",
    utility="ceinutil.py",
    file_ext="cein",
    flag_input="-cein",
    flag_output="-ceout",
    flag_restart="-cerestrt",
    ext_output="ceout",
    ext_restart="cerestrt",
    supports_output_prmtop=False,
    range_flag_min="-mineo",
    range_flag_max="-maxeo",
    range_label="Eo (V)",
    mdin_flag_keyword="icnste",
    mdin_setpoint_keyword="solve",
    mdin_freq_keyword="ntcnste",
    mdin_relax_keyword="ntrelaxe",
    residue_names=frozenset(REDOX_RESIDUE_EO),
    allowed_igb_implicit=_HEH_IGB_IMPLICIT,
    allowed_igb_explicit=_HEH_IGB_EXPLICIT,
    allowed_intdiel=frozenset({1.0}),
)

_PHREDOX_MODE = TitrationMode(
    key=PHREDOX,
    label="Joint pH/Redox (proton-coupled electron transfer)",
    short_label="C(pH,E)MD",
    utility="cpeinutil.py",
    file_ext="cpein",
    flag_input="-cpein",
    flag_output="-cpeout",
    flag_restart="-cperestrt",
    ext_output="cpeout",
    ext_restart="cperestrt",
    supports_output_prmtop=True,
    # cpeinutil exposes no range filter at all.
    range_flag_min=None,
    range_flag_max=None,
    range_label="pKa / Eo",
    # A cpein run still sets both keyword families; these two fields name the
    # primary pair only, so consumers should use `mdin_keyword_sets` instead.
    mdin_flag_keyword="icnstph",
    mdin_setpoint_keyword="solvph",
    mdin_freq_keyword="ntcnstph",
    mdin_relax_keyword="ntrelax",
    # cpeinutil applies no typ filter, so a cpein may hold anything.
    residue_names=frozenset(PH_RESIDUE_PKA)
                  | frozenset(REDOX_RESIDUE_EO)
                  | frozenset(PHREDOX_RESIDUE_NAMES),
    # A cpein containing HEH inherits HEH's narrower coverage. (TYX is
    # narrower still -- only igb=2 has all four of its states set.)
    allowed_igb_implicit=_HEH_IGB_IMPLICIT,
    allowed_igb_explicit=_HEH_IGB_EXPLICIT,
    allowed_intdiel=frozenset({1.0}),
)


TITRATION_MODES: Dict[str, TitrationMode] = {
    PH: _PH_MODE,
    REDOX: _REDOX_MODE,
    PHREDOX: _PHREDOX_MODE,
}


def get_mode(key: Optional[str]) -> TitrationMode:
    """Resolve a mode key, defaulting to constant pH.

    ``None`` maps to :data:`PH` so that workspace records and WorkflowConfigs
    written before the redox modes existed (which carry only ``cpin_file``)
    keep behaving exactly as they did.
    """
    if not key:
        return _PH_MODE
    try:
        return TITRATION_MODES[key]
    except KeyError:
        raise ValueError(
            f"Unknown titration mode '{key}'. "
            f"Expected one of: {', '.join(sorted(TITRATION_MODES))}"
        ) from None


def modes_for_residue_names(resnames: Iterable[str]) -> List[str]:
    """Which titration files a topology needs, in generation order.

    Returns a LIST because a structure that titrates in both pH and redox needs
    a cpin *and* a cein -- they are complementary files run together, not
    competing options (tutorial 33).

        ['ph']            only pH-titratable residues
        ['redox']         only redox-titratable residues (HEH)
        ['ph', 'redox']   both families present -- generate both files
        ['phredox']       a proton-coupled site (TYX) is present, which fits
                          in neither cpin nor cein; the joint file carries
                          everything, and sander refuses it alongside either
                          of the others

    Empty list when nothing titratable was found.
    """
    present = {str(r).strip().upper() for r in resnames}

    # A phredox residue cannot go in a cpin or a cein, and a cpein cannot be
    # combined with either -- so its presence collapses everything into one
    # joint file.
    if present & frozenset(PHREDOX_RESIDUE_NAMES):
        return [PHREDOX]

    modes: List[str] = []
    if present & frozenset(PH_RESIDUE_PKA):
        modes.append(PH)
    if present & frozenset(REDOX_RESIDUE_EO):
        modes.append(REDOX)
    return modes


def partition_residues_by_mode(residues: Sequence[dict],
                               modes: Sequence[str]) -> Dict[str, List[dict]]:
    """Split a scanned residue list into the file each residue belongs in.

    Partitioning matters beyond tidiness: handing ceinutil a `-resnums` list
    that includes an AS4 is a hard error ("Residue number N [AS4] is not
    titratable"), and handing cpinutil a HEH is the same. Tutorial 33
    partitions with `-resnames`; ProPrep partitions by residue so individual
    sites can still be deselected.
    """
    out: Dict[str, List[dict]] = {m: [] for m in modes}
    for mode_key in modes:
        allowed = get_mode(mode_key).residue_names
        out[mode_key] = [r for r in residues if r.get('resname') in allowed]
    return out


def engine_flags_for_modes(files_by_mode: Dict[str, str],
                           prefix: str) -> List[str]:
    """Flags for every generated titration file, concatenated.

    A combined run gets all six: ``-cpin ... -cpout ... -cprestrt ... -cein ...
    -ceout ... -cerestrt ...``, matching tutorial 33 section 3.
    """
    flags: List[str] = []
    for mode_key in sorted(files_by_mode, key=_mode_order):
        mode = get_mode(mode_key)
        flags.extend(mode.engine_flags(files_by_mode[mode_key], prefix))
    return flags


def mdin_keyword_sets(modes: Sequence[str]) -> List[TitrationMode]:
    """The keyword families a set of modes drives, deduplicated.

    A cpein run sets both families from one file, so it maps to the pH and
    redox descriptors rather than to itself.
    """
    wanted: List[str] = []
    for mode_key in modes:
        mode = get_mode(mode_key)
        if mode.titrates_ph and PH not in wanted:
            wanted.append(PH)
        if mode.titrates_redox and REDOX not in wanted:
            wanted.append(REDOX)
    return [get_mode(k) for k in wanted]


def _mode_order(key: str) -> int:
    """Stable ordering: pH first, then redox, then joint (tutorial order)."""
    return {PH: 0, REDOX: 1, PHREDOX: 2}.get(key, 99)


# ---------------------------------------------------------------------------
# ParmEd cross-checks
# ---------------------------------------------------------------------------

def _parmed_residue(resname: str):
    """Fetch ParmEd's TitratableResidue definition, or None."""
    try:
        from parmed.amber import titratable_residues as _res
    except ImportError:
        return None
    if resname not in getattr(_res, 'titratable_residues', ()):
        return None
    return getattr(_res, resname, None)


def parmed_titratable_atom_count(resname: str) -> Optional[int]:
    """Number of atoms ParmEd expects in a titratable residue, or None.

    This matters because cpinutil/ceinutil/cpeinutil all filter termini by
    comparing the topology residue's atom count against ``len(res.atom_list)``
    and **silently skip** any mismatch (ceinutil.py:236) -- no message, and the
    run still exits 0 as long as something survived.
    """
    res = _parmed_residue(resname)
    if res is None:
        return None
    try:
        return len(res.atom_list)
    except (AttributeError, TypeError):
        return None


def parmed_residue_type(resname: str) -> Optional[str]:
    """ParmEd's ``typ`` for a residue ('ph' / 'redox' / 'phredox'), or None."""
    res = _parmed_residue(resname)
    if res is None:
        return None
    return getattr(res, 'typ', None)
