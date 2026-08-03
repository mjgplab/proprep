"""
FF compatibility tooling.

Two-phase design (see docs/ff_collision_plan.md):

- ``proprep.ff_compat.parser``: parse a bundled or user FF set
  (frcmod + lib(s)) into a normalized ``FFSignature`` capturing every
  parameter that could collide with another set's parameters.

- ``proprep.ff_compat.matrix``: pairwise compare every bundled +
  user FF set, producing a JSON compatibility matrix. Standalone CLI:
  ``python -m proprep.ff_compat.build_matrix``.

- ``proprep.ff_compat.resolver`` (later): given a set of selected FFs
  with known conflicts, rewrite one or more sets into a workspace-local
  directory with renamed atom types so the load is conflict-free.

This package is intentionally not user-facing in its current form; it's
a maintainer / CI tool. Tool 2's resolver will be hooked into the FF
picker once it's stable.
"""

from proprep.ff_compat.parser import FFSignature, parse_set

__all__ = ["FFSignature", "parse_set"]
