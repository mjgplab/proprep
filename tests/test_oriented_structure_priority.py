"""The Structure Orientator wrote `<stem>_oriented.pdb` and registered it, but
nothing downstream ranked it: the Topology Generator's priority list had no
entry for it and the shared registry put it below "Filtered", so tLEaP always
loaded the unoriented coordinates and the orientation step had no effect on the
built system. Orientation is not cosmetic for a solvated build, because the
bounding box sets the box size and the number of waters.

Order required: the membrane build outranks the oriented structure, and the
oriented structure outranks every protein-stage structure below it."""

import re
from pathlib import Path

import proprep.tleap_prep.tleap_input_generator as tig
from proprep.utils.structure_selector import StructureRegistry


def _topology_priority_keys():
    """The ordered workspace keys from the Topology Generator's own list."""
    source = Path(tig.__file__).read_text()
    block = source[source.index("priority_keys = ["):]
    block = block[: block.index("]")]
    return [m.group(1) for m in re.finditer(r'\("([a-z_]+)"', block)]


def test_topology_generator_prefers_the_oriented_structure():
    keys = _topology_priority_keys()
    assert "oriented_pdb_file" in keys, (
        "the Topology Generator would build from unoriented coordinates"
    )
    for lower in (
        "protonation_pdb_file",
        "transformed_pdb_file",
        "repaired_pdb_file",
        "filtered_pdb_file",
    ):
        assert keys.index("oriented_pdb_file") < keys.index(lower), (
            f"oriented must outrank {lower}"
        )


def test_membrane_build_outranks_the_oriented_structure():
    keys = _topology_priority_keys()
    assert keys.index("membrane_packed_pdb") < keys.index("oriented_pdb_file"), (
        "a packed bilayer already contains an oriented protein and must win"
    )


def test_registry_ranks_oriented_above_protonation_and_filtered():
    by_key = {s.workspace_key: s.priority for s in StructureRegistry.get_all()}
    assert by_key["oriented_pdb_file"] < by_key["protonation_pdb_file"]
    assert by_key["oriented_pdb_file"] < by_key["filtered_pdb_file"]
