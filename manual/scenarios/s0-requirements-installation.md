# s0 — Requirements and Installation (recorded on this Mac)

Recorded 2026-11-21 on macOS 27 (Darwin 27.0.0), Apple silicon (arm64), with
Miniforge conda 26.7.2 already installed. ProPrep version installed: **1.23.0**.

## Feature: Requirements
**What it is for (plain English):** what a computer needs before ProPrep will run.

| Need | Required? | Notes from this machine |
|---|---|---|
| macOS or Linux (Windows only via WSL2) | Yes | macOS arm64 worked |
| Python 3.12+ | Yes (installer brings its own) | System Python 3.15 was NOT used; the env has 3.12 |
| AmberTools (tleap, sander, antechamber, parmchk2, MCPB.py, cpptraj, pdb4amber) | Yes | Installer brought AmberTools 26 |
| MODELLER + free academic licence key | Only for structure repair / loop modelling | MODELLER installed by the installer, but without a key it refuses to start ("Invalid license key: XXXX"); ProPrep still runs and says repair is unavailable |
| Gaussian (g16) | Only for new QM-derived parameters / ONIOM runs | Not installed here |
| ORCA | Only to run ORCA QM/MM inputs | Not installed here |
| BLAST+ | Possibly for homology search (see s3) | Not installed here |
| pmemd / pmemd.cuda (Amber, licensed) | Optional, faster MD on CPU/GPU | Not installed; `sander` (free) is |
| AutoDock Vina, packmol, packmol-memgen, propka, pdb2pqr, reduce | Used by docking, membrane, protonation | All brought by the installer |
| Disk space | ~5 GB for the conda env (measured 4.8 GB); INSTALL.md says ~8 GB free for the self-contained installer | |
| Internet | Needed to fetch structures from RCSB/AlphaFold, and for the conda install route | |
| Web browser | Only for `proprep-web` and the structure viewer | |

## Feature: Installation

### Scenario 1: Run from the source folder in an existing AmberTools env — RUN
- Starting point: an older conda env (`etanalyze`) with AmberTools 25 and Python 3.12, no ProPrep.
- Command: `PYTHONPATH=<repo>/src <env>/bin/python -m proprep.main`
- Result: ProPrep 1.23.0 started and showed its menu-mode prompt. It warned
  "MODELLER not found. Structure repair (missing loops/residues) will be unavailable."
  pydantic/psutil/web packages were missing in that env. **Not a recommended route for users.**

### Scenario 2: Option B — conda installer script `install_proprep.sh` — RUN
- Starting point: Miniforge installed, no `ProPrep` env.
- Command: `bash install_proprep.sh` (run from the repository; the README one-liner downloads the same script with curl).
- What the user sees:
  1. Asks: "Enter your MODELLER license key (or press Enter to skip):" (typing shows `*`). We skipped → "NOTE: Skipping MODELLER key. Homology modeling will not be available until you configure the key manually."
  2. If a `ProPrep` env already exists it offers: 1 Update, 2 Full reinstall, 3 Cancel.
  3. Creates env `ProPrep` with Python 3.12, installs `proprep=1.23.0` from channels mjgplab, dacase, conda-forge, salilab, bioconda, then `pip install tmtools`, then force-reinstalls ProPrep so it wins over the older copy bundled in AmberTools.
  4. Verifies: "ProPrep installed successfully: ProPrep v1.23.0 - Proper Protein Preparation workflow manager", then prints "Installation complete!" with how to start.
- Time: several minutes (download-dominated). Env size 4.8 GB.
- Result: `/opt/homebrew/Caskroom/miniforge/base/envs/ProPrep` with `proprep`, `proprep-web`, `proprep-library`, `mpsa`.
- To use afterwards: `conda activate ProPrep`, `source $CONDA_PREFIX/amber.sh`, `proprep` or `proprep-web`.
- Problems: none. Note the installer's skip message says "Homology modeling" — it means MODELLER-based repair.

### Scenario 3: Option A — self-contained installer `ProPrep-<version>-<platform>.sh` — RUN (macOS arm64, 2026-11-27)
- Downloaded `ProPrep-1.23.0-MacOSX-arm64.sh` (915 MB) from the v1.23.0 GitHub release
  (other assets: MacOSX-x86_64 978 MB, Linux-x86_64 1410 MB). Header: "Created by constructor 3.16.1".
- Ran in a clean environment: `bash ~/Downloads/ProPrep-1.23.0-MacOSX-arm64.sh -b -p ~/ProPrep` → exit 0 in 126 s.
  Output: `PREFIX=...`, "Unpacking bootstrapper...", "Unpacking payload...", "Installing base environment...",
  a MODELLER post-link note ("Edit .../modeller/config.py and replace XXXX with your Modeller license key"),
  a list of removed leftover files from the AmberTools-bundled ProPrep, then
  "ProPrep v1.23.0 - Proper Protein Preparation workflow manager" / "ProPrep is installed." / "To run:" with
  proprep, proprep-web, proprep-library, and the `~/.proprep/modeller_key` instructions.
- Size: 4.8 GB. Contains tleap, sander, MCPB.py, vina, packmol-memgen, propka3, MODELLER 10.8 (no key);
  no pmemd, xtb, crest.
- `~/ProPrep/bin/proprep --version` → "ProPrep v1.23.0 - Proper Protein Preparation workflow manager".
- `~/ProPrep/bin/proprep-web --no-browser --port 8477` started (uvicorn lines) and `/healthz` answered `ok`.
- With nothing activated, ProPrep sets AMBERHOME to `~/ProPrep` itself (amber_env bootstrap).
- FAILED part: `~/ProPrep/bin/proprep-library` → "No such file or directory" (shebang `#!.../bin/python3.13`, env has 3.12;
  cause: missing from recipe entry_points, fixed in mjgplab/proprep#2). Workaround works:
  `~/ProPrep/bin/python -m proprep.utils.user_library_snapshots status`.
- Update: same command with `-u`. Remove: delete `~/ProPrep`. Windows (WSL) and Linux not tried.

### Scenario 4: Updating inside an existing AmberTools env — SOURCE-ONLY
- `update_proprep_in_ambertools.sh <env-name>`; only for conda AmberTools (`dacase::ambertools-dac`).

### Scenario 5: MODELLER licence key — SOURCE-ONLY (key not available)
- `mkdir -p ~/.proprep && echo 'YOUR_KEY' > ~/.proprep/modeller_key`, or `export KEY_MODELLER=...` before installing.
- Without it: start-up message says structure repair is unavailable; everything else works.

## Notes for the manual writer
- For non-technical users recommend Option A (one file, no conda). Option B suits people who already use conda.
- Explain "terminal", "conda environment", "licence key", "AmberTools" in one sentence each.
- The `etanalyze` plugin was auto-loaded when ProPrep ran in the user's older env — plugins are discovered from the environment.

## Feature: Using the right environment (found during recording)

### Scenario 6: Another AmberTools env is active in the shell — RUN (observed by s5, s6, s7, s8)
- On this Mac the shell had a different conda env (`etanalyze_env`, AmberTools 25) activated, so `AMBERHOME`
  and `PATH` pointed there. ProPrep respects an existing `AMBERHOME`, so tLEaP, sander, ParmEd, pytraj,
  packmol-memgen and the force-field files all came from that other env, not ProPrep's own AmberTools 26.
  ProPrep gives no warning about this.
- Advice for the manual: open a fresh terminal and run `conda activate ProPrep` (Option B), or use
  `~/ProPrep/bin/proprep` (Option A); if `echo $AMBERHOME` shows a different folder, run `unset AMBERHOME`
  (or `source $CONDA_PREFIX/amber.sh` after activating ProPrep).

### Scenario 7: Clean check with only the ProPrep env — RUN
- Re-ran s7's ProPrep-generated ubiquitin input (`1UBQ_tleap.in`: ff19SB, TIP3P truncated octahedron 10 Å,
  13 Na+ / 13 Cl-) with `PATH` and `AMBERHOME` set to the ProPrep env only.
- Result: tLEaP exit 0, "Errors = 0; Warnings = 13", 5065 residues, 3.1 MB `clean.prmtop` + `clean.rst7`,
  force-field files read from `envs/ProPrep/dat/leap/...`. Files: `runs/s0-cleanenv/`.
- Conclusion: the topology results recorded by the agents are not an artefact of the wrong env; the
  standard path works with ProPrep's own AmberTools.
