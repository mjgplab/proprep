# Glossary terms introduced by w1 (Quick start, Ch. 1-3, App. C)

- **topology** — The file (`.prmtop`) that lists every atom of the simulated system with its bonds and force-field parameters. Amber's simulation programs read it together with a coordinate file.
- **coordinates** / **coordinate file** — The file (`.rst7`) that gives the position of every atom. It is paired with the topology; after a simulation step it also holds the atom velocities.
- **Protein Data Bank (PDB)** — The public archive of experimentally solved structures of proteins and nucleic acids. Each entry has a four-character code such as 1UBQ.
- **pH** — A measure of how acidic the surroundings are. It decides whether groups such as histidine, aspartate or glutamate carry a proton (and a charge).
- **MD (molecular dynamics)** — A calculation that moves every atom a tiny step at a time using the forces between atoms, producing a short movie of how the molecule moves.
- **trajectory** — The file (usually `.nc`) in which an MD simulation saves snapshots of the atom positions over time.
- **force field** — The set of numbers (bond lengths, angles, charges, and so on) that describes the forces between atoms in a simulation, for example ff19SB for proteins.
- **session file** — The file (`proprep_session_<date>_<time>.json`) in which ProPrep records every answer you give, so the work can be resumed, corrected, replayed or used as a template.
- **DOI** — Digital Object Identifier: a permanent web address for a publication, software release or data set.
- **AmberTools** — The free part of the Amber simulation suite (tLEaP, sander, cpptraj, antechamber, ParmEd, MCPB.py, ...). ProPrep's installer includes it.
- **Windows Subsystem for Linux (WSL)** — A feature of Windows 10/11 that runs a small Linux system alongside Windows; ProPrep runs inside it on Windows.
- **Ubuntu** — The Linux system that WSL installs by default; its window is the Terminal you use for ProPrep on Windows.
- **conda** — A program that installs scientific software and its dependencies into separate, self-contained folders.
- **conda environment** — One such self-contained folder of software. Option B installs ProPrep into an environment named `ProPrep`, switched on with `conda activate ProPrep`.
- **MODELLER** — A program from the Sali lab that builds missing parts of protein structures (loops, residues, atoms) and homology models. Installed with ProPrep; needs a free academic license key.
- **license key** — A personal code you receive after registering for a program (here MODELLER) that lets it run. ProPrep reads the MODELLER key from `~/.proprep/modeller_key` or the `KEY_MODELLER` setting.
- **AMBERHOME** — A setting in the Terminal that names the folder of an Amber installation. If it points to another Amber, ProPrep uses that one instead of its own.
