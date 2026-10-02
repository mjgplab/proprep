# Meeko without prody, for the mjgplab channel

ProPrep's Molecular Docking prepares ligands and receptors with
[Meeko](https://github.com/forlilab/Meeko). ProPrep's conda package declares it
as a dependency, built from this recipe.

## Why not conda-forge's meeko

conda-forge's `meeko` declares `prody`, and conda-forge has no `prody` that
installs on Apple Silicon with Python 3.12 or newer: 2.6.x is built for linux-64
and osx-64 only, and the noarch 2.4.1 and 2.5.0 builds need Python < 3.12 (one
2.5.0 build also requires `gcc_linux-64`). So on Apple Silicon conda-forge's
meeko 0.8.0 cannot be installed next to ProPrep (Python >= 3.12); left unpinned,
the solver silently picks meeko 0.5.0, from before prody was required. Meeko itself declares
no dependencies (its PyPI metadata has none). It imports prody only in
`covalentbuilder.py` and `utils/prodyutils.py` (covalent docking, ProDy object
input), and inside guarded `try` blocks elsewhere; ProPrep's docking uses
neither. This recipe builds the unmodified PyPI source with prody left out of
the run requirements (build string `noprody_py_0`).

Meeko is LGPL-2.1: the package ships Meeko's own LICENSE and unmodified source.

## Build and upload (mgp)

```bash
conda build packaging/meeko-noprody -c conda-forge
anaconda upload -u mjgplab "$(conda build packaging/meeko-noprody -c conda-forge --output)"
```

The build is `noarch: python`, so one upload serves osx-arm64, osx-64 and
linux-64. Its test prepares a ligand end to end (atom types, torsion tree,
PDBQT) and checks that the prody-only module reports what it needs.

ProPrep's install commands list mjgplab first, so under strict channel priority
this build is the one installed.

## Retiring it

If conda-forge's meeko-feedstock makes prody optional (requested in conda-forge/meeko-feedstock#13; see
`feedstock_request.md`), drop this package from the channel, and ProPrep's
recipe takes conda-forge's meeko unchanged. When ProPrep moves to a newer Meeko,
update `version` and `sha256` here (both are on PyPI) and rebuild, or retire
this build if conda-forge's no longer needs prody.
