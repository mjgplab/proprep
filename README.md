# ProPrep

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21764540.svg)](https://doi.org/10.5281/zenodo.21764540)

ProPrep is a prompt-driven workflow manager for preparing protein and
metalloprotein structures for molecular dynamics (MD) and QM/MM simulations
with the Amber suite — from RCSB PDB retrieval through structure curation,
force field parameterization, tLEaP topology generation, simulation setup,
and analysis.

ProPrep automates the labor of preparation, not the decisions. Each step is a
small, mostly single-keystroke choice presented with the what, why, and how of
that step, so a user stays informed and in control rather than accepting
hidden defaults. The interface is textual throughout, which keeps it usable
over a remote connection and by visually impaired researchers; a graphical
structure viewer is available but optional.

## What it does

- **Structure loading and curation** — RCSB search and retrieval, metadata
  browsing, chain and component filtering, textual descriptions of composition
  and quaternary structure.
- **Repair and mutagenesis** — missing residues (from `REMARK 465`, SEQRES vs
  ATOM records, a deposited FASTA, or numbering discontinuities), missing
  atoms, alternate locations, loop modeling via MODELLER, and substitutions to
  standard or modified amino acids.
- **Protonation states** — PROPKA-based assignment with per-residue rationale,
  including constant-pH and fixed-pH treatments.
- **Force fields and parameterization** — browsing existing AMBER parameters
  and generating new ones for small organic molecules, modified amino acids,
  and metal coordination sites (a reimplementation of the MCPB.py methodology).
  Finished parameters are deposited into a local library so they can be
  discovered and reused in later structures.
- **Redox sites** — detection, grouping into sites, and template-driven
  transformation to force-field-compatible residue names and states.
- **Topology generation** — tLEaP input assembled with custom atom types,
  parameter load commands, and bond directives, viewable and editable before
  it runs; multi-microstate batches supported.
- **Simulation setup and analysis** — an MD manager for protocols, restraints,
  hardware configuration, and queueing, with an interface to pytraj for
  structural, energetic, and geometric analyses, plus a plugin interface for
  third-party analyses.

## Session recording and replay

Every user input is recorded, with timestamps, to a session log that doubles
as replayable history. A log can be used to resume an interrupted session,
document or demonstrate a protocol, undo a mistaken choice through the
interactive session editor, or become a template for high-throughput batch
processing.

A template promotes one or more recorded answers to variables; a CSV then
supplies a set of values per run. ProPrep validates that every run supplies
every variable before any of them start.

Session logs are tied to the ProPrep version that recorded them: a later
release that adds or reworks a prompt can leave an older log without an answer
for it. Record and replay with the same version.

## Requirements

- Python 3.12 or later
- macOS or Linux natively; Windows through the Windows Subsystem for Linux,
  which supplies the POSIX environment AmberTools requires
- `proprep-web` (the browser-based shell) depends on POSIX pseudoterminals and
  is therefore macOS and Linux only
- AmberTools, for tLEaP and the Amber programs ProPrep drives
- A MODELLER license key for structure repair and loop modeling
  (free for academics — https://salilab.org/modeller/registration.html)
- Gaussian, only if you generate new QM-derived parameters

ProPrep also ships with AmberTools (from AmberTools26), which is released
annually. The conda channel below carries current releases between those.

## Install

The installer creates a self-contained conda environment named `ProPrep`,
installs the package from the [`mjgplab` conda
channel](https://anaconda.org/mjgplab/proprep), and pip-installs `tmtools` for
structure alignment. It pins and verifies the exact version, so a partial
conda solve cannot silently leave you on an older release.

Requires Miniforge or Anaconda/Miniconda
([download](https://github.com/conda-forge/miniforge#download)).

One-liner:

```
curl -fsSL https://raw.githubusercontent.com/mjgplab/proprep/main/install_proprep.sh | bash
```

Or inspect it before running:

```
curl -fsSLO https://raw.githubusercontent.com/mjgplab/proprep/main/install_proprep.sh
less install_proprep.sh
bash install_proprep.sh
```

## Run

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

Rerun the installer. It installs the pinned version and aborts if conda
resolves anything different:

```
curl -fsSL https://raw.githubusercontent.com/mjgplab/proprep/main/install_proprep.sh | bash
```

It detects an existing `ProPrep` environment and offers to update or rebuild
it. To update manually instead:

```
conda install -n ProPrep -c mjgplab -c dacase -c salilab -c bioconda -c conda-forge proprep -y
conda run -n ProPrep pip install --upgrade tmtools
```

### Updating inside an existing AmberTools environment

If you installed AmberTools through conda (`dacase::ambertools-dac`) and want
the current ProPrep in that same environment rather than a separate one, use
the in-place updater. AmberTools bundles an older ProPrep at the same paths,
so the helper installs the current package, forces it to take precedence,
clears the stale metadata that would make `proprep --version` misreport, and
verifies the result:

```
curl -fsSL https://raw.githubusercontent.com/mjgplab/proprep/main/update_proprep_in_ambertools.sh | bash -s -- <env-name>
```

Replace `<env-name>` with the conda environment holding your AmberTools; omit
it to use the active environment or be prompted. The helper only runs against
a conda AmberTools — on a source build it stops rather than pulling in a
second copy. Rerun it if you later reinstall or update `ambertools-dac`.

## Repository layout

```
src/proprep/     source
recipe/          conda build recipe
tests/           test suite
docs/            documentation
examples/        example inputs
scripts/         maintenance and audit scripts
```

## Reporting issues

ProPrep has integrated issue reporting: bugs, unexpected behavior, and feature
requests can be submitted from within the program. Submissions include the
Python version and operating system, and optionally an encrypted copy of the
session history. Issues can also be opened directly on this repository.

## Citing ProPrep

See [CITATION.cff](CITATION.cff). Please cite both the software (archived
release) and the accompanying article.

## Changelog

Release-by-release notes are in [CHANGELOG.md](CHANGELOG.md). Versions
correspond to the `proprep` package on the `mjgplab` conda channel and to
tagged releases here.

## License

MIT — see [LICENSE](LICENSE).
