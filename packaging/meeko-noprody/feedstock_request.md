# Request to conda-forge/meeko-feedstock

Opened 2026-10-02 as https://github.com/conda-forge/meeko-feedstock/issues/13
(by mgp's account). If conda-forge's meeko stops requiring prody, ProPrep's recipe can take it unchanged and the mjgplab build (this folder) can be retired.

---

**Title:** Make prody optional: meeko 0.8.0 cannot be installed on osx-arm64 with Python >= 3.12

`meeko` lists `prody` as a run requirement (added in #7), and conda-forge has no `prody` that installs on osx-arm64 with Python 3.12 or newer: prody 2.6.x is built for linux-64 and osx-64 only, the noarch 2.4.1 and 2.5.0 builds require Python < 3.12, and one noarch 2.5.0 build also requires `gcc_linux-64` (the attempts at osx-arm64 builds in conda-forge/prody-feedstock#54, #56, #57 and #58 were closed unmerged). So on Apple Silicon with Python >= 3.12, meeko 0.8.0 cannot be installed:

```
conda create -n t -c conda-forge python=3.12 meeko=0.8.0
```

fails to solve on osx-arm64, and without the version pin the solver silently picks meeko 0.5.0, from before prody was required.

Meeko itself declares no dependencies in its PyPI metadata, and imports prody unconditionally only in `meeko/covalentbuilder.py` and `meeko/utils/prodyutils.py`; in `meeko/__init__.py` and `meeko/polymer.py` the import is inside `try` blocks. Ligand and receptor preparation (`MoleculePreparation`, `Polymer.from_pdb_string`, PDBQT writing) work without it: we use meeko 0.8.0 without prody for docking on osx-arm64 and linux-64.

Would you consider moving `prody` from `run` to `run_constrained` (or otherwise making it optional), so that meeko installs on every platform and users who need the ProDy-based features install prody themselves? We'd be glad to open a PR if that helps.

Context: we package docking in ProPrep (https://github.com/mjgplab/proprep), which depends on meeko; for now we ship meeko 0.8.0 without the prody requirement on our own channel, and would rather depend on conda-forge's.
