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

## Update

Rerun the script, or inside the environment:

```
conda update -n ProPrep -c mjgplab -c dacase -c salilab -c conda-forge proprep -y
conda run -n ProPrep pip install --upgrade tmtools
```

## License

MIT — see [LICENSE](LICENSE).
