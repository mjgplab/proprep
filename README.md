# ProPrep Installer

ProPrep is a protein preparation assistant that guides users through the
complete workflow of preparing PDB structures for molecular dynamics
simulations with the Amber suite — from PDB retrieval through force field
parameterization, tLEaP input generation, and MD input preparation.

This repository hosts the installer only. The source code is distributed via
the [`mjgplab` conda channel](https://anaconda.org/mjgplab/proprep) and will
also ship with a future release of AmberTools.

## Install

One-liner (pipes to bash):

```
curl -fsSL https://raw.githubusercontent.com/mjgplab/proprep/main/install_proprep.sh | bash
```

Or, to inspect before running:

```
curl -fsSLO https://raw.githubusercontent.com/mjgplab/proprep/main/install_proprep.sh
less install_proprep.sh
bash install_proprep.sh
```

The script creates a conda environment named `ProPrep`, installs the package
from the `mjgplab` channel, and pip-installs `tmtools` (structure alignment).
It pins and verifies the exact ProPrep version, so a partial conda solve can no
longer silently leave you on an older release.

## Run

Activate the environment and launch ProPrep:

```
conda activate ProPrep
proprep              # interactive command-line interface
proprep-web          # browser-based UI (web shell)
```

To use the bundled Amber tools in the same shell:

```
source $CONDA_PREFIX/amber.sh
```

## Update

Rerun the installer — it now installs the pinned version and aborts if conda
resolves anything different:

```
curl -fsSL https://raw.githubusercontent.com/mjgplab/proprep/main/install_proprep.sh | bash
```

It detects the existing `ProPrep` environment and offers to update or rebuild
it. To update manually instead:

```
conda install -n ProPrep -c mjgplab -c dacase -c salilab -c bioconda -c conda-forge proprep -y
conda run -n ProPrep pip install --upgrade tmtools
```

## License

MIT — see [LICENSE](LICENSE).
