"""ProPrep test suite.

This file exists to make ``tests`` a package. Without it, pytest takes
``tests/structure_prep`` (which has its own ``__init__.py``) as the package
root and imports its modules as ``structure_prep.test_*``, which collides with
the real ``proprep.structure_prep`` package and fails collection for the whole
directory. With it, the modules import as ``tests.structure_prep.test_*``.
"""
