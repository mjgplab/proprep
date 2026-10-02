"""Which channel built this copy of ProPrep.

None, the value in the source tree, means no channel of ours built it: the
AmberTools build (CMake runs setup.py on this tree), a source checkout, or an
install from the public source. The mjgplab conda recipe (recipe/meta.yaml)
rewrites the line below to CHANNEL = "mjgplab" when it builds the package,
and the standalone installers are made from that package.

Anything that points users at ProPrep's own GitHub releases (the startup
notice of a newer release, proprep.utils.update_check) runs only when
CHANNEL is "mjgplab". AmberTools users get ProPrep updates with AmberTools,
so their copy must never send them to ours, and a copy that was not built
by the recipe says nothing rather than guessing.
"""

CHANNEL = None
