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

### Scenario 3: Option A — self-contained installer `ProPrep-<version>-<platform>.sh` — SOURCE-ONLY (documented in INSTALL.md, not run here)
- `bash ~/Downloads/ProPrep-<version>-MacOSX-arm64.sh -b -p ~/ProPrep` (Intel Mac: `MacOSX-x86_64`; Linux and Windows/WSL: `Linux-x86_64`).
- 1–1.5 GB download, ~8 GB free disk, 2–5 min install, no internet needed during install. Prints "Unpacking payload..." then "ProPrep is installed."
- Test: `~/ProPrep/bin/proprep --version`. Start: `~/ProPrep/bin/proprep-web` (browser) or `~/ProPrep/bin/proprep` (terminal).
- Update: same command with `-u`. Remove: delete `~/ProPrep`.
- Windows: `wsl --install` in an admin PowerShell, restart, create Ubuntu user, run the Linux installer from `/mnt/c/Users/<name>/Downloads/`, open http://127.0.0.1:8000/ in the Windows browser.
- INSTALL.md already has a good plain-language troubleshooting list (Permission denied → add `bash`; already exists → `-u`; no space; port 8000 busy → 8001).

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
