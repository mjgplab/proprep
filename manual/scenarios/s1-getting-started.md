# Getting started & the application shell (scope s1)

All scenarios were driven against ProPrep 1.23.0 (`/opt/homebrew/Caskroom/miniforge/base/envs/ProPrep/bin/proprep`)
on macOS, terminal 140 columns wide, on 2026-11-21. Preferences and `--set-default` runs used an isolated
home folder (`HOME=manual/runs/s1/home`), so the real `~/.proprep` was not changed (checked afterwards: it still
has no `settings.json`).

Transcripts are in `manual/runs/_sessions/<session>/transcript.txt`. Run folders are in `manual/runs/s1/`.

| Session name | What it covered |
|---|---|
| `s1-first` | First launch with no flags, Full Menu, unavailable tool, invalid input, `undo`, Help, Workspace Options (all 12), Feedback, Exit |
| `s1-home1` | Isolated HOME: first-run update notice, Workflow Mode chosen and saved as default, Project = current directory, all workflow stages, Utilities submenu, all Preferences options, load 1UBQ, switch `a`/`g` |
| `s1-home2` | Isolated HOME: "ask" mode, "Previous session files found!" prompt, `--menu-grid` |
| `s1-home3` | Isolated HOME: `--full-menu --menu-list --no-session` (list overrides saved grid), `undo` without recording, Ctrl-C |
| `s1-home4` | Isolated HOME: saved default `workflow` skips the menu-mode question |
| `s1-debug` | `--debug` |
| `s1-verbose` | `--verbose` |
| `s1-resume`, `s1-resume2` | Launching inside an existing project folder, session browser, project name with `/` (crash) and with a space, Exit from a submenu |
| `s1-outdir` | `--output-dir` with a relative nested path |
| `s1-wf2` | `--workflow` flag, `b` on the first stage, Reset workspace inside Workflow Mode |

---

## Feature: Starting ProPrep for the first time (no flags)

**What it is for (plain English):** Typing `proprep` opens the program. On the first run it asks which kind
of menu you want and where to keep your files. After that you are at the main menu.
**Needs before you start:** The ProPrep environment is installed and you have a terminal open. An internet
connection is needed only for downloads and the once-a-day update check.
**Where to find it:** Type `proprep` in a terminal. The folder you are in when you type it matters (see Scenario 6).

### Scenario 1: Very first launch, choose Full Menu, do not save, named project — RUN
- Starting point: real HOME with no `~/.proprep/settings.json`. Folder `runs/s1/first`. Command `proprep`.
- Steps:
  1. Start-up messages appear before anything else (see "What the user sees").
  2. Prompt `Select menu mode  [1/2]:` -> pressed Enter with no answer. There is **no default** here:
     `Received '', which is not one of the available options`, and the question is asked again.
  3. Same prompt -> typed `3` -> `Received '3', which is not one of the available options`.
  4. Same prompt -> typed `2` (Full Menu Mode).
  5. Prompt `Save this as your default preference?  [y/n] (n):` -> typed `maybe` -> rejected the same way.
  6. Same prompt -> typed `n`. Message: `Tip: Use --workflow or --full-menu flags to skip this prompt`.
  7. "Project Setup" screen -> prompt `Project name (.):` -> typed `myproject`.
- What the user sees:
  - Before the banner, every launch prints a MODELLER warning:
    `MODELLER not found. Structure repair (missing loops/residues) will be unavailable.` followed by install and
    licence hints (`conda install -c salilab modeller`, `export KEY_MODELLER='your_license_key_here'`), then
    `Successfully imported core modules`. (On this machine MODELLER is installed but has no valid licence key.)
  - The ProPrep banner box: ASCII logo, "Perform Protein Preparation like a Pro  v1.23.0", "Learn the process,
    skip the pain", "Guberman-Pfeffer Lab, New Mexico State University".
  - The "Choose Your Experience" box describes `1. Workflow Mode (Recommended for new users)` (guided,
    step-by-step, shows only the tools for the current stage) and `2. Full Menu Mode (For experienced users)`
    (all tools visible, no guided progression).
  - Project Setup text: "Enter a name for your project. All structures and outputs will be saved in a directory
    with this name. Press Enter to use the current directory."
  - After the name: `✓ Project directory: .../runs/s1/first/myproject`, then
    `Session recording to: ./proprep_session_20261121_123119.json`,
    `Tip: Use --no-session to disable automatic session recording`,
    `Tip: Type 'undo' at any prompt to rewind to an earlier answer`, then the Full Menu.
- Result / files produced: the folder `runs/s1/first/myproject/` was created and ProPrep moved into it.
  `proprep_session_20261121_123119.json` holds a record of every answer you give after this point. It can be
  replayed.
- Problems / confusing bits:
  - The MODELLER warning appears on every start, even when you will never repair loops. It can be ignored
    unless you need Structure Fixer loop repair.
  - The menu-mode question has no default, so Enter alone does nothing. You must type 1 or 2.
  - The menu-mode, save-default, previous-session and project-name answers are given **before** recording
    starts, so they are not in the session file and `undo` cannot go back to them.
- Session log: `runs/_sessions/s1-first/transcript.txt`

### Scenario 2: First launch (isolated HOME), choose Workflow, save as default, press Enter for project — RUN
- Starting point: empty isolated HOME, folder `runs/s1/wf`, command `HOME=.../runs/s1/home proprep`.
- Steps:
  1. Prompt `Select menu mode  [1/2]:` -> `1`.
  2. Prompt `Save this as your default preference?  [y/n] (n):` -> `y`.
     Messages: `✓ Saved 'workflow' as default` and `Tip: Use --set-default to change this later`.
  3. Prompt `Project name (.):` -> pressed Enter (the default `.` means "the current folder").
- What the user sees:
  - Under the banner, on the first run only:
    `ProPrep checks GitHub once a day for a newer release. To turn this off: Preferences (p) > Check for new ProPrep releases at startup.`
  - `✓ Project directory: .../runs/s1/wf`. No new folder is made; files go straight into the current folder.
  - Session-recording lines as in Scenario 1, then the Workflow Mode screen (see the Workflow Mode feature).
- Result / files produced: `runs/s1/home/.proprep/settings.json` was created with `"default_menu_mode": "workflow"`.
  `runs/s1/home/.proprep/update_check.json` stores when GitHub was last asked.
  `runs/s1/wf/proprep_session_20261121_123926.json`.
- Session log: `runs/_sessions/s1-home1/transcript.txt`

### Scenario 3: Saved default skips the menu-mode question — RUN
- Starting point: isolated HOME after `proprep --set-default workflow`. Command `proprep --no-session`.
- What the user sees: banner, then straight to "Project Setup". No "Choose Your Experience" box. After Enter,
  `Navigation: Main Menu (Workflow Mode)`.
- With the saved default set to `ask` (session `s1-home2`), the question comes back on every launch.
- Session log: `runs/_sessions/s1-home4/transcript.txt`

### Scenario 4: Update-check message — RUN (partly SOURCE-ONLY)
- Seen: the one-time explanation line in Scenario 2. Later launches with the same HOME printed nothing,
  because GitHub had been asked less than 24 hours before and 1.23.0 is the newest release.
- From source (`src/proprep/utils/update_check.py`, not seen because no newer release exists):
  - When a newer release exists, ProPrep prints a bold line `ProPrep X.Y.Z is available (you have 1.23.0).`,
    then the update command for your type of install and a link to the release page. Nothing is downloaded.
  - GitHub is asked at most once every 24 hours, in the background. Start-up waits for it for at most
    1.5 seconds. With no internet connection nothing is printed.
  - Only copies from the mjgplab conda channel or the ProPrep installers check. This copy is one
    (`CHANNEL = "mjgplab"`). Copies that came with AmberTools never check. Their Preferences option 4 says
    "not available in this copy".
  - The check can be turned off in Preferences, option 4 (Scenario P4).

### Scenario 5: Project name variants — RUN
- `myproject` -> a new folder is created inside the current folder (Scenario 1).
- Enter -> current folder is used (Scenario 2).
- `my project` (with a space) -> accepted: `✓ Project directory: .../myproject/my project`. Spaces in folder
  names may cause trouble for the AMBER command-line tools later. Advise users to avoid them.
- `sub/deeper` (a folder inside a folder that does not exist yet) -> **ProPrep crashes and quits**:
  `Error: [Errno 2] No such file or directory: 'sub/deeper'` / `Run with --debug flag for more details`, then
  the goodbye box. Only one folder level can be created. (Session `s1-resume`.)
- An existing folder name is fine (`mkdir(exist_ok=True)`), and existing files in it are kept.

### Scenario 6: Launching inside a folder that already has session files — RUN
- Starting point: `cd runs/s1/first/myproject` (contains a `proprep_session_*.json`), `proprep --full-menu`.
- What the user sees, after the banner (and after the menu-mode question when one is asked), **before** Project Setup:
  ```
  Previous session files found!
  What would you like to do?
    1. Start fresh session (record new interactions)
    2. Resume a previous session from any point
    3. Continue without session recording
  Select option  [1/2/3] (1):
  ```
  - `2` -> a table "Available Session Files" (#, Filename, Date/Time, Interactions, Notes), then
    `Select session to resume (0 to cancel)  [1/0] (0):`. `0` -> goes on as a fresh session.
  - Picking `1` -> "Session Info" box (number of interactions) and the options `1. Continue from end (append new
    interactions)`, `2. Open session editor (edit/delete multiple interactions)`, `3. Cancel`. `3` -> fresh session.
    (Replaying or editing a session is covered in another scope.)
  - `3` (Continue without session recording) -> Project Setup, then the menu with **no** "Session recording to"
    lines. No session file is written.
- ProPrep looks for session files only in the folder you launch from (`--session-dir`, default `.`). Session
  files are saved inside the project folder, so to be offered "Resume" you must `cd` into the project folder
  first. Launching from the parent folder and typing the project name again does not offer to resume.
- Session logs: `runs/_sessions/s1-resume/`, `s1-resume2/`, `s1-home2/`

### Options reference (start-up prompts)
| Prompt | Choices | Default | Effect |
|---|---|---|---|
| `Select menu mode [1/2]` | 1 = Workflow Mode, 2 = Full Menu Mode | none (must type) | Mode for this run. Skipped if `--workflow`/`--full-menu` is given or a default is saved |
| `Save this as your default preference? [y/n]` | y / n | n | y writes `default_menu_mode` to `~/.proprep/settings.json`, so the question is not asked again |
| `Previous session files found! Select option [1/2/3]` | 1 fresh, 2 resume, 3 no recording | 1 | Shown only when `proprep_session_*.json` exists in the launch folder |
| `Project name (.)` | any folder name, or Enter | `.` (current folder) | Creates one folder level and moves into it. All outputs go there |

---

## Feature: Full Menu Mode (main menu)

**What it is for (plain English):** One screen listing every ProPrep tool, grouped by stage, with a mark
showing whether each tool can be used yet.
**Needs before you start:** Nothing.
**Where to find it:** Choose 2 at start-up, run `proprep --full-menu`, or press `a` in Workflow Mode.

### Scenario 1: Reading the main menu — RUN
- Layout (list style, the default): sections `1. STRUCTURE INPUT` (1a Structure Loader, 1b Biological Assembly
  Generator), `2. REDOX SITE DETECTION` (2a Redox Site Detector, 2b Force Field Explorer, 2c Force Field
  Parameterizer), `3. STRUCTURE COMPARISON` (3a Homology Searcher, 3b Structure Aligner), `4. STRUCTURE PREPARATION`
  (4a PDB Filter, 4b Amino Acid Mutator, 4c Structure Fixer, 4d Redox Site Preparer, 4e Protonation State Analyzer,
  4f Structure Orientator, 4g Membrane Builder), `5. SIMULATION SETUP & EXECUTION` (5a Topology Generator, 5b
  Molecular Dynamics Manager, 5c QM/MM Preparator, 5d Molecular Docking), `UTILITIES` (ua Structure Viewer), and
  `WORKSPACE & SETTINGS` (w, p, h, f, g, x). Each tool has a one-line description.
- Symbols:
  - `✓` (green) = the tool can be used now.
  - `○` (yellow) = the tool is not ready yet because something it needs is missing.
  - `⚠` line under a `○` tool = why it is not ready, e.g. `⚠ Needs a loaded structure`,
    `⚠ Needs detected redox sites — run Redox Site Detector first`,
    `⚠ Needs a topology — run Topology Generator first`,
    `⚠ Needs an experimental structure (RCSB or local PDB; not AlphaFold)`,
    `⚠ Needs a structure, a packed membrane, or a parm7/rst7 pair`.
  - `[exp]` (orange, after a name) = experimental tool (from source; none showed in this build).
  - Plugin tools would show `[plugin vX]` (from source; none were installed).
- With nothing loaded, the ✓ tools were 1a, 2a, 2b, 3a, 4g. After loading 1UBQ, everything became ✓ except
  4d Redox Site Preparer, 5b MD Manager and 5c QM/MM Preparator.
- At the top of every screen there is a breadcrumb such as `Navigation: Main Menu > Structure Loader`.
- Prompt: `Enter your choice (1a):`. Enter alone opens 1a Structure Loader.
- Session log: `runs/_sessions/s1-first/transcript.txt`

### Scenario 2: Choosing a tool that is not yet available — RUN
- Typed `1b` (○ Biological Assembly Generator) with no structure loaded.
- What the user sees: the tool's own menu still opens:
  ```
  ══ Biological Assembly Generator Menu ══
  ⚠ Needs a loaded structure
    1. ○ Generate biological assembly from asymmetric unit
         ⚠ Load a structure first
  [b]ack · [m]main · [x]exit
  ```
  Typing `1` only repeats `⚠ Load a structure first` and shows the same menu again. Nothing breaks.
  You leave with `b` (back), `m` (main menu) or `x` (exit).
- Invalid input inside a tool menu (`zz`): `Received 'zz', which is not one of the available options`.

### Scenario 3: Invalid input at the main menu — RUN
- `9z` -> `Invalid choice: 9z. Please try again.` / `Hint: Use format like '1a', '3c', 'w', 'h', or 'x'`.
- `W` (capital) -> also rejected. **The Full Menu is case-sensitive** (Workflow Mode is not).

### Scenario 4: `undo` — RUN
- At `Enter your choice (1a):` typed `undo`.
- What the user sees: a table "Recent answers in this session" (columns #, module, question, answer), here
  0 = `1b Biological Assembly Generator`, 1 = `1 Generate biological assembly...`, 2 = `b Back`, 3 = `9z`, 4 = `W`.
  Prompt: `Rewind to which answer? number to rewind to, [m]ore, or [x] cancel (4):`.
  - `m` showed the list again (there were no older answers).
  - Typed `1` -> `[r]e-answer it live, or [c]hange it and replay what followed? ([x] cancel) (r):` -> Enter (r).
  - Messages: `[Rewinding: question 1 will be asked again; 4 answers from there on are set aside in a backup of
    the session log]`, `Rewinding session: proprep_session_..json`,
    `[Session backup: proprep_session_20261121_123119.json.bak-20261121_123224]`, "Hybrid mode started ...".
    ProPrep replays answer 0 (`Enter your choice [REPLAY: 1b]`) and stops at answer 1 in the Biological Assembly
    menu: `[Replay complete - now recording new interactions]`.
- Files: a backup `proprep_session_..json.bak-<time>` is written next to the session file.
- Notes: invalid main-menu answers (`9z`, `W`) are recorded and appear in the undo list. Invalid answers to a
  fixed-choice question (`zz`) are not. The rewind screen also prints a technical log line
  (`... proprep.main - INFO - main.py:768 - Rewind to interaction 1 ...`) and a stray `[0m` code in the table.
  With `--no-session` (or "Continue without session recording") `undo` is not available:
  `Invalid choice: undo`.

### Scenario 5: Returning to the main menu — RUN
- In a tool menu: `m` = straight to main menu, `b` = back one level (to the main menu from a tool's first screen).
  Help pages end with `Press Enter to return to main menu...`. Workspace Options/Preferences have a numbered
  "Back to main menu" item (12 and 6).

### Scenario 6: `x` Exit and the exit screen — RUN
- `x` -> `Are you sure you want to exit? [y/n] (n):`. `n` (or Enter) stays in the menu. `y` ->
  `Goodbye!`, `[Session recording stopped]`, then the box "Thank you for using ProPrep": please cite Walker, A. &
  Guberman-Pfeffer, M.J. (2026) ... bioRxiv https://doi.org/10.64898/2026.02.26.708365, and "Feedback (positive
  or negative) is welcome at: mjgp@nmsu.edu".
- `x` from inside a tool menu asks the same question. Answering `n` there takes you to the **main menu**, not
  back to the tool you were in.
- Ctrl-C at any prompt: `Application interrupted by user`, then the same thank-you box. The session file is still
  saved up to that point (session `s1-home3`).
- The code is meant to print `Session saved to: <file>` at exit, but this never appeared in any run (bug, see notes).

### Scenario 7: `--menu-grid` vs `--menu-list` — RUN
- `proprep --menu-grid`: the same tools drawn as bordered panels in two columns (1|2, 3|4, 5 on the left), with
  UTILITIES and WORKSPACE & SETTINGS side by side at the bottom. Descriptions wrap inside the panels. Keys,
  ✓/○/⚠ are the same. Described as "good for screenshots". Needs a wide terminal (140 columns was fine).
- `proprep --menu-list`: single column (the default).
- The flags apply to this run only and override the saved Preferences layout. With grid saved in Preferences,
  `--menu-list` still showed the list (session `s1-home3`).
- Using both together: `proprep: error: argument --menu-list: not allowed with argument --menu-grid`.
- Layout affects Full Menu only. Workflow Mode always looks the same.

### Options reference (Full Menu)
| Key | Action |
|---|---|
| `1a`–`5d`, `ua` | Open that tool's menu (opens even when ○, but its actions are blocked) |
| `w` | Workspace Options |
| `p` | Preferences |
| `h` | Help & Workflows |
| `f` | Feedback & Support |
| `g` | Switch to Workflow Mode (`Switching to Workflow Mode...`) |
| `x` | Exit (asks to confirm) |
| `undo` | Rewind to an earlier answer (needs session recording) |
| Enter | Same as `1a` |
| Inside tool menus | `b` back · `m` main menu · `x` exit |

---

## Feature: Workflow Mode (guided)

**What it is for (plain English):** A step-by-step version of the menu. It shows only the tools for the current
stage (Load → Detect → Compare → Fix → Simulate) and suggests what to do next.
**Needs before you start:** Nothing.
**Where to find it:** Choose 1 at start-up, run `proprep --workflow`, or press `g` in the Full Menu.

### Scenario 1: The Workflow screen — RUN
- What the user sees:
  ```
  ╭──── WORKFLOW GUIDE ────╮
  │ Load:▶ → Detect:○ → Compare:○ → Fix:○ → Simulate:○
  │ Start by loading a structure from databases or local files and then optionally applying symmetry
  │ transformations to generate a biological assembly.
  ══ 1. STRUCTURE INPUT ══
    1. ✓ Structure Loader
    2. ○ Biological Assembly Generator   ⚠ Needs a loaded structure
    v. Structure Viewer  |  u. Utilities...
  [n]ext | [l]oad [d]etect [c]ompare [f]ix [s]imulate | [a]ll tools | e[x]it
  Enter your choice (1):
  ```
- Stage marks in the guide: `▶` current, `✓` done, `⊘` skipped, `○` not reached yet.
- The suggestion text changes per stage:
  - Detect: "Detect redox-active sites, explore the parameter space of AMBER FFs, and if needed, parameterize
    small organic molecules, modified amino acids, or metal sites."
  - Compare: "Compare sequences or align structures".
  - Fix: "Filter, mutate, repair, prepare redox sites, assign protonation states, orient your structure, or
    embed it in a membrane".
  - Simulate: "Build topologies, configure, execute, or analyze MD simulations, or setup QM/MM calculations on
    MD snapshots."
- Tools are numbered 1, 2, 3 … **within each stage**, so `1` means a different tool at each stage
  (Load: 1 = Structure Loader; Fix: 1 = PDB Filter; Simulate: 1 = Topology Generator).
- Input is not case-sensitive here (`N` worked as `n`, `X` as `x`). Invalid input:
  `Invalid choice: zz. Please try again.`
- `b` is not offered on the first stage (typing it: `Invalid choice: b`), and `n` is not offered on the last.
- Session log: `runs/_sessions/s1-home1/transcript.txt`

### Scenario 2: Moving between stages — RUN
- `n` (next) -> `→ Moved to Detect stage`. The stage you leave is marked ✓ **even if you did nothing in it**.
- `s` (jump to Simulate from Detect) -> `Jumped to Simulate stage`. The guide shows `Load:✓ → Detect:⊘ →
  Compare:⊘ → Fix:⊘ → Simulate:▶`, so stages jumped over are marked skipped.
- `s` again -> `Already on Simulate stage`.
- `l` -> `Jumped to Load stage`. Jumping back clears the marks.
- `b` (from Simulate) -> `← Moved back to Fix stage`. The guide then showed `Fix:✓` with no `▶` anywhere,
  because a stage already marked done keeps its ✓ when you return to it.

### Scenario 3: Loading 1UBQ and how stages unlock — RUN
- From Load stage: `1` (Structure Loader) -> `1` Load protein structure -> `1` RCSB -> `1` Enter PDB ID ->
  `1ubq` (lower case accepted) -> Enter (asymmetric unit) -> `✓ Downloaded 1UBQ` plus title/resolution/citation
  -> `6` Done -> back in the Structure Loader menu.
- Leaving the loader with **`b`** -> `✓ Structure loaded! Advanced to Detect stage`. The guide shows
  `Load:✓ → Detect:▶`, and Force Field Parameterizer turns from ○ to ✓.
- Leaving the loader with **`m`** (main menu) the first time did **not** auto-advance: still on Load:▶ with
  Biological Assembly Generator now ✓. Inconsistent (see notes).
- Walking on with `n`: Compare stage shows Homology Searcher ✓ and Structure Aligner ✓. Fix stage shows all ✓
  except `4. ○ Redox Site Preparer ⚠ Needs detected redox sites`. Simulate stage shows `1. ✓ Topology
  Generator`, `2. ○ Molecular Dynamics Manager ⚠ Needs a topology`, `3. ○ QM/MM Preparator ⚠ Needs detected
  redox sites`, `4. ✓ Molecular Docking`.
- Files: `runs/s1/wf/1UBQ.pdb`.

### Scenario 4: `v` Structure Viewer shortcut — RUN
- With nothing loaded, `v` opens the Structure Viewer menu with most items ○ (`⚠ Load a structure first`).
  "Save the scene shown in the open viewer" and "Load a saved scene" show ✓. The viewer itself is covered in
  another scope. `m` returns to the workflow screen.

### Scenario 5: `u` Utilities submenu — RUN
- Shows `w. Workspace Options`, `p. Preferences`, `h. Help & Workflows`, `f. Feedback & Support`, `c. Cancel`.
  Prompt `Select option (c):`, so Enter = cancel. `x` is rejected here (`Received 'x', which is not one of the
  available options`). w/p/f open the same screens as in Full Menu. h shows the Workflow-Mode help (below).

### Scenario 6: Switching between modes — RUN
- `a` (all tools) -> `Switching to Full Menu mode...` -> Full Menu. In this run the line
  `✓ Structure loaded! Advanced to Detect stage` was printed **twice** during the switch (spurious).
- `g` in Full Menu -> `Switching to Workflow Mode...`. The guide was back at `Load:✓ → Detect:▶`, so the progress
  made before switching (we had reached Fix) was lost.
- The switch lasts for this run only. The saved default is changed only in Preferences or with `--set-default`.

### Scenario 7: Reset workspace inside Workflow Mode — RUN
- `u` -> `w` -> `5` Reset workspace -> `y` -> `12` Back. ProPrep went to the **Full Menu**
  (`Navigation: Main Menu`) instead of Workflow Mode, because resetting also clears the remembered menu mode.
  (Session `s1-wf2`.)

### Options reference (Workflow Mode)
| Key | Action |
|---|---|
| `1`…`7` | Tools of the current stage (numbering restarts each stage) |
| `v` | Structure Viewer |
| `u` | Utilities submenu (w / p / h / f / c) |
| `n` | Mark this stage done and go to the next (not shown on Simulate) |
| `b` | Go back one stage (not shown on Load) |
| `l` `d` `c` `f` `s` | Jump to Load / Detect / Compare / Fix / Simulate (forward jumps mark skipped stages ⊘) |
| `a` | Switch to Full Menu |
| `x` | Exit (asks to confirm) |

---

## Feature: w. Workspace Options

**What it is for (plain English):** Tools for looking at or saving what ProPrep is currently holding in memory
(the "workspace": loaded structures, results), plus a debug switch and a command log. Most beginners will not
need it.
**Needs before you start:** Nothing (more is shown after loading a structure).
**Where to find it:** `w` in Full Menu, or `u` then `w` in Workflow Mode.

Menu: `1. Show workspace status`, `2. Show detailed workspace`, `3. Save workspace`, `4. Load workspace`,
`5. Reset workspace`, `6. Toggle debug mode (currently OFF)`, `7. Show command history`, `8. Save command history`,
`9. Load command history`, `10. Replay command history`, `11. Reset command history`, `12. Back to main menu`.
Prompt `Enter your choice (1):`. All were RUN in session `s1-first` (with 1UBQ loaded).

### Scenario 1: Show workspace status (1) — RUN
- A table "Workspace Contents" with Key / Type / Size. Empty session: `menu_mode`, `project_directory`. After
  loading 1UBQ: also `rcsb_metadata`, `rcsb_pdb_file`, `rcsb_structure`. These are internal names and mean little
  to a biologist.

### Scenario 2: Show detailed workspace (2) — RUN
- Output was **identical** to option 1. The "detailed" setting is not used by the code (`dump_rich(detailed=...)`
  ignores it).

### Scenario 3: Save workspace (3) — RUN
- Prompt `Enter filename to save workspace (workspace.json):` -> Enter -> `Workspace saved to workspace.json`.
- File `runs/s1/first/myproject/workspace.json` contained only `{"pdb_file": null}`. The loaded 1UBQ was **not**
  saved, because the code saves old key names (`pdb_file`, `filter_selections`, `completeness_results`). It
  does not restore your work. Use session files instead.

### Scenario 4: Load workspace (4) — RUN
- Enter (workspace.json) -> `Workspace loaded from workspace.json`. Only added `pdb_file: None` to the status table.
- `nofile.json` -> `File nofile.json not found`. Command history still listed this as ✓.

### Scenario 5: Reset workspace (5) — RUN
- `Are you sure you want to reset the workspace? [y/n] (n):`. `n` returns, but the breadcrumb then read
  `Navigation: Home`. `y` -> `Workspace has been reset`. The status table is completely empty, the loaded
  structure is forgotten (1b etc. go back to ○), and the files on disk stay. It also clears `menu_mode`
  and `project_directory` (see Workflow Scenario 7). The breadcrumb shows "Home" instead of "Main Menu" from then on.

### Scenario 6: Toggle debug mode (6) — RUN
- `Debug mode enabled`, and the label changes to `(currently ON)`. Pressing again: `Debug mode disabled`. Same
  effect as starting with `--debug`.

### Scenario 7: Show command history (7) — RUN
- `===== Recent Commands (last 10) =====` with times and ✓/○ marks, e.g. `5. [12:34:16] ✓ Save workspace`.
  Lists menu-level commands (show menu, save workspace...), not the answers typed inside tools.

### Scenario 8: Save command history (8) — RUN
- `Enter filename to save command history (history.json):` -> Enter -> `Command history saved to history.json`.
  File: `runs/s1/first/myproject/history.json` (JSON list of commands plus breadcrumb trail).

### Scenario 9: Load command history (9) — RUN
- Enter -> `This will replace current command history. Continue? [y/n] (n):` -> `y` ->
  `Command history loaded from history.json`, `Loaded 14 command records`, `History exported at: ...`.
- Afterwards option 7 showed only one command, so the 14 records were not actually put back. Only the
  breadcrumb trail is restored.

### Scenario 10: Replay command history (10) — RUN
- Enter (history.json) -> header (File, Commands: 14, Exported) and sub-options:
  `1. View all commands` (numbered list with ✓/✗ and timestamps, then `Press Enter to continue`),
  `2. Replay all functional commands` (lists 10 commands, asks `Execute 10 functional command(s)? [y/n] (n)`.
  Answered `n`), `3. Restore workspace state` (`No commands with captured workspace state found`),
  `4. Restore breadcrumb trail` (`Restore breadcrumb trail: Main Menu > Workspace Options? [y/n] (y)`. Answered `n`),
  `5. Back to workspace menu`.

### Scenario 11: Reset command history (11) — RUN
- Shows `Current command history: N commands`, `Current breadcrumb trail: M items`, then
  `Are you sure you want to clear all command history and breadcrumbs? [y/n] (n):`. `n` ->
  `Command history reset cancelled`. `y` -> `Command history and breadcrumbs cleared`.

### Scenario 12: Back (12) and invalid input — RUN
- `12` returns to the main menu. Invalid numbers are rejected with
  `Received '...', which is not one of the available options`.

---

## Feature: p. Preferences

**What it is for (plain English):** Settings that are remembered between runs, stored in
`~/.proprep/settings.json`.
**Needs before you start:** Nothing.
**Where to find it:** `p` in Full Menu, `u` then `p` in Workflow Mode.
All RUN with isolated HOME (session `s1-home1`, plus `s1-home2` for layout).

Menu as shown:
```
1. Change default menu mode (currently: Workflow Mode (guided))
2. Change full-menu layout (currently: List (single column))
3. Configure GitHub feedback (currently: not configured)
4. Check for new ProPrep releases at startup (currently: on)
5. Reset all preferences to defaults
6. Back to main menu
Enter your choice (1):
```

### Scenario P1: Change default menu mode (1) — RUN
- Shows `1. Workflow Mode (guided, step-by-step)`, `2. Full Menu Mode (all tools visible)`, `3. Ask each time`.
  Prompt `Select mode (1):` (the default is always 1, whatever the current setting).
- Picking the current mode -> `Mode unchanged`.
- `2` -> `✓ Default menu mode set to: full-menu` / `Note: Restart ProPrep for the change to take effect`.
  The current run keeps its mode.
- `3` -> `✓ Default menu mode set to: ask`. The "Choose Your Experience" question comes back on the next launch
  (confirmed in `s1-home2`).

### Scenario P2: Change full-menu layout (2) — RUN
- `1. List (single column, one tool per row)`, `2. Grid (2-column panels, good for screenshots)`, prompt
  `Select layout (1):`. `2` -> `✓ Full-menu layout set to: grid`. Choosing it again -> `Layout unchanged`.
- Unlike the menu mode, the layout change is used the next time the Full Menu is drawn, without restarting.

### Scenario P3: Configure GitHub feedback (3) — RUN (with a fake token, isolated HOME only)
- Shows a box explaining you need a GitHub Personal Access Token (classic, `repo` or `public_repo` scope, made at
  https://github.com/settings/tokens) and a repository `owner/repo`. Your token is stored locally in
  `~/.proprep/settings.json`.
- Prompt `GitHub repository (owner/repo) (Mag14011/ProPrep):`, then `Enter your GitHub Personal Access Token
  (Input will be hidden)` / `Token:` (typing is invisible).
  - Enter for repo + empty token -> `Token is required`. Nothing saved.
  - Repo `badrepo` + a token -> `Invalid repository format. Expected: owner/repo` (only after the token is typed).
  - Repo `someone/somerepo` + `faketoken123` -> `✓ GitHub feedback configuration saved`, `Repository:
    someone/somerepo`, `Token: (hidden)`, `You can now use the Feedback menu to submit issues!`. The menu then says
    `(currently: configured)`. The token is **not checked** with GitHub and is stored **as plain text** in
    `settings.json`.
  - Opening it again shows `Current repository: someone/somerepo`, `Current token: ******************** (hidden)`.
    Enter at both prompts -> `Keeping existing token`, saved again.
- The suggested default repository `Mag14011/ProPrep` differs from the one the update check uses
  (`mjgplab/proprep`). Ask the maintainers which is right before telling users.

### Scenario P4: Check for new releases at startup (4) — RUN
- Text: "At startup ProPrep can ask GitHub (github.com/mjgplab/proprep), at most once a day, whether a newer
  release exists, and print the update command if so. Nothing is downloaded or changed, and nothing about you
  is sent." / `The check is currently on.` / `1. On`, `2. Off`.
- `2` -> `✓ Startup check for new releases turned off`. Menu then reads `(currently: off)`.
- In a copy that came with AmberTools this option only prints an explanation that such copies never check
  (from source).

### Scenario P5: Reset all preferences (5) — RUN
- `This will reset all preferences to defaults. Continue? [y/n] (n):`. `n` -> `Reset cancelled`.
- `y` -> `✓ All preferences reset to defaults`. Layout back to List, GitHub not configured (token removed),
  update check on. The menu mode stayed "Ask each time", which is the default.

### Scenario P6: Back (6) and invalid input — RUN
- `9` -> `Received '9', which is not one of the available options`. `6` returns to the menu.

### Options reference (settings.json keys)
| Key | Values | Default | Changed by |
|---|---|---|---|
| `default_menu_mode` | workflow / full-menu / ask | ask | P1, `--set-default`, "save as default" at start-up |
| `menu_layout` | list / grid | list | P2 (per run: `--menu-grid` / `--menu-list`) |
| `github_repo`, `github_token` | text / null | null | P3 |
| `check_for_updates` | true / false | true | P4 |

---

## Feature: h. Help & Workflows

**What it is for (plain English):** A one-page reminder of typical tool orders.
**Where to find it:** `h` in Full Menu. In Workflow Mode `u` then `h` shows a different page.
There is only one page in each mode. Enter returns.

### Scenario 1: Full Menu help — RUN
- Title "Common Workflows", then:
  1. BASIC PROTEIN MD SIMULATION: `1a.Load → 3a.Filter → 3c.Fix → 3d.Protonate → 5b.tLEaP → 5c.MD`
  2. METALLOPROTEIN MD SIMULATION: `1a.Load → 3a.Filter → 3c.Fix → 4a.RedoxSite → 5a.FF-Params → 5b.tLEaP → 5c.MD`
  3. QM/MM CALCULATION: `1a.Load → 3a.Filter → 3c.Fix → 4a.RedoxSite → 5b.tLEaP → 5d.QM/MM`
  4. STRUCTURE COMPARISON: `1a.Load (multiple) → 2b.Aligner → ua.Viewer`
  5. CONTINUE EXISTING MD: `5c.MD (load existing topology/coordinates)`
  Then "Tool Categories" with one line per section. `Press Enter to return to main menu...`
- **The menu keys in this help are out of date** and do not match the current menu. The correct keys are:
  PDB Filter = 4a, Structure Fixer = 4c, Protonation = 4e, Redox Site Preparer = 4d, Force Field Parameterizer
  = 2c, Topology Generator (tLEaP) = 5a, MD Manager = 5b, QM/MM = 5c, Structure Aligner = 3b. Typing what the
  help says would open the wrong tool (e.g. `3a` opens Homology Searcher, `5b` opens MD Manager).

### Scenario 2: Workflow Mode help — RUN
- Title "Workflow Mode Help": how it works (`Use 'n' to advance, 'b' to go back, 's' to skip optional stages`,
  `Use 'a' to switch to Full Menu mode`), a list of stages Load / Detect / Compare / Fix / **Setup** / Simulate,
  and four workflows (Basic Protein MD, Metalloprotein MD, QM/MM, With Alignment).
- Mismatches: `s` does not "skip". It jumps to Simulate. There is no "Setup" stage in the guide bar (topology
  is in Simulate).

---

## Feature: f. Feedback & Support

**What it is for (plain English):** A form that turns your comment or bug report into a GitHub Issue for the
developers. Nothing was submitted in any test.
**Needs before you start:** A GitHub token set in Preferences option 3. Without one, submitting fails and the
feedback is saved to a local file instead (from source).
**Where to find it:** `f` in Full Menu, `u` then `f` in Workflow Mode. Developers can also be emailed at
mjgp@nmsu.edu (shown on the exit screen).

### Scenario 1: Look through the whole form, decline to submit — RUN
- Intro box: what you can send, "All feedback is submitted as a GitHub Issue (anonymously if you prefer)", and
  "First Time Setup" steps to create a token. Step 6 says `Main Menu → p. Preferences → 2. Configure GitHub`, but
  it is option **3** (option 2 is layout).
- Questions:
  1. `Select feedback type (3):` with 1 Bug Report, 2 Feature Request, 3 Suggestion, 4 Question, 5 Other.
  2. `Module number (or type module name) (0):` with 0 General, 1 PDB Loader, 2 Redox Site Detector, 3 Homology
     Searcher, 4 Structure Aligner, 5 PDB Filter, 6 Amino Acid Mutator, 7 Structure Fixer, 8 Protonation State
     Analyzer, 9 Redox Site Preparer, 10 MD Restraint Manager, 11 Force Field Parameterizer, 12 Topology
     Generator, 13 Molecular Dynamics Manager, 14 QM/MM Preparator, 15 Structure Viewer. Typing the name
     `Structure Viewer` also worked. The list uses the old name "PDB Loader" (= Structure Loader) and has no
     Docking, Membrane Builder, Orientator, Bio Assembly or Force Field Explorer.
  3. `Email address ():` blank = "Anonymous".
  4. `Summary:` one line (required).
  5. Detailed description: one line per prompt (` ():`), and an empty line finishes (required).
  6. "Include Session Context?" appears only when encryption is available and set up. It did **not** appear
     here.
- Then "Feedback Preview" (Type / Module / Email / Summary / Description) and
  `Submit this feedback to GitHub Issues? [y/n] (y):`. **The default is YES**, so pressing Enter would try to
  submit. Typed `n` -> `Feedback not submitted.`
- Session log: `runs/_sessions/s1-first/transcript.txt`

### Scenario 2: Cancel by leaving the summary empty — RUN
- Pressing Enter at `Summary:` -> `Summary is required.` / `Feedback cancelled.` Back to the menu. This is the
  easiest way out of the form. An empty description cancels the same way (from source).

### Scenario 3: Submitting — SOURCE-ONLY (not run on purpose)
- With `y`: if no token or repo is set -> `✗ GitHub integration not configured. Please run Preferences menu to set
  up GitHub token.` and the feedback is saved to `~/.proprep/feedback/feedback_<time>.json`
  ("You can manually submit this later or email it to the developer."). With a token, it posts an issue titled
  `[<Type>] <Summary>` to the configured repo, and reports success or an error (also saving locally on failure).

---

## Feature: Command-line flags for start-up

**Where to find it:** Typed after `proprep` in the terminal.

### Scenario 1: `proprep --help` — RUN
- Prints usage and grouped options (general; maintainer feedback tools; session options; menu mode options;
  workflow shortcuts; template and batch processing), then exits. No banner. Fast.

### Scenario 2: `proprep --version` — RUN
- Prints `ProPrep v1.23.0 - Proper Protein Preparation workflow manager` and exits (about 0.5 s).

### Scenario 3: `proprep --set-default workflow|full-menu|ask` — RUN (isolated HOME)
- Prints `✓ Default menu mode set to: workflow` (etc.) and exits without opening the menu. `settings.json` updated.
- `--set-default bogus` -> `proprep: error: argument --set-default: invalid choice: 'bogus' (choose from workflow,
  full-menu, ask)`.

### Scenario 4: `--workflow` / `--full-menu` — RUN
- Skip the menu-mode question for this run. Both together:
  `proprep: error: argument --full-menu: not allowed with argument --workflow`.

### Scenario 5: `proprep --debug` — RUN
- Starts normally. Workspace Options shows `Toggle debug mode (currently ON)`. Some tools print extra
  diagnostic lines, and if ProPrep crashes it shows the full technical error report (traceback) instead of the
  short `Error: ... Run with --debug flag for more details`. When loading the non-existent PDB `ZZZZ`, the
  user saw `PDB ZZZZ not found (404). May be a large structure only available in mmCIF format.` plus two warnings
  and `Failed to download PDB ZZZZ`, then back to source selection.
- Session `s1-debug`, files `runs/s1/dbg/dbgproj/1CRN.pdb`, `proprep_session_20261121_124604.json`.

### Scenario 6: `proprep --verbose` (or `-v`) — RUN
- Floods the screen with timestamped technical `DEBUG` lines (e.g. "Registered RedoxSite transformer: ...",
  "Registered structure type: ..."), both at start-up and around menus. Useful only when sending a problem
  report to developers.
- Session `s1-verbose`.

### Scenario 7: `--no-session` — RUN
- No "Session recording to" lines, no previous-session question, no session file, and `undo` does not work.

### Scenario 8: `--output-dir <folder>` — RUN
- Skips the Project Setup question and creates the folder (several levels allowed). Session file went to
  `runs/s1/cli/outdir_test/a/b/`. With a relative path, the log line and the stored project folder name are
  wrong: `Working directory set to: .../outdir_test/a/b/outdir_test/a/b`. This is a possible bug for tools that
  use the project folder. An absolute path avoids it.

---

## Notes for the manual writer

**Recommended beginner path:** `cd` to where you keep your work -> `proprep --workflow` -> type a short project
name without spaces or `/` -> Structure Loader (`1`) -> load a PDB ID -> leave the loader with `b` -> follow the
WORKFLOW GUIDE with `n`. To come back to the same project later, `cd` into the project folder first so ProPrep
offers to resume the recorded session.

**Jargon to explain:** workspace (what ProPrep holds in memory during a run), session file / session recording,
replay, rewind/`undo`, breadcrumb ("Navigation:" line), asymmetric unit vs biological assembly, topology, tLEaP,
redox site, PAT (GitHub token), debug/verbose.

**Common mistakes:**
- Pressing Enter at "Select menu mode": it has no default and repeats.
- Typing capital letters in the Full Menu (`W`, `X`): rejected. Workflow Mode accepts them.
- Typing workflow keys from the Help screen, which uses old keys (see Help Scenario 1).
- In Workflow Mode, number keys mean different tools at each stage.
- Pressing Enter at "Submit this feedback to GitHub Issues?" means YES.
- Expecting "Save workspace" to save your work. It does not. Session files are the real record.
- Launching from the parent folder and expecting the "Resume" offer.

**Things that looked like bugs (described, not fixed):**
1. Project name containing `/` (e.g. `sub/deeper`) crashes ProPrep at start-up with `[Errno 2] No such file or directory`.
2. Help & Workflows (Full Menu) lists outdated menu keys (3a Filter, 3c Fix, 3d Protonate, 5b tLEaP, 5c MD,
   5d QM/MM, 2b Aligner). Workflow help says `s` skips stages (it jumps to Simulate) and lists a non-existent
   "Setup" stage.
3. Feedback intro says Preferences → "2. Configure GitHub". It is option 3.
4. "Show detailed workspace" is identical to "Show workspace status".
5. "Save workspace" writes only `{"pdb_file": null}` and ignores the loaded structure. "Load command history"
   says "Loaded 14 command records" but restores none. A failed "Load workspace" is listed as ✓ in history.
6. "Reset workspace" also wipes the menu mode and project folder. In Workflow Mode it silently drops you into the
   Full Menu. Afterwards the breadcrumb says "Home".
7. Workflow Mode: leaving Structure Loader with `m` does not auto-advance to Detect, but `b` does. Switching with
   `a` printed "✓ Structure loaded! Advanced to Detect stage" twice. Returning with `g` lost the stage progress.
   After `b` the current stage shows ✓ instead of ▶. `n` marks a stage ✓ even if nothing was done.
8. "Session saved to: <file>" (intended at exit) never appeared in any run.
9. Cosmetic: the Structure Loader suggestion reads "or press  to return to the main menu" (key missing). Stray
   `[0m` in the undo table. An `INFO` log line is printed during rewind and on Ctrl-C.
10. `--output-dir` with a relative path stores a doubled, wrong project path.
11. GitHub token is stored unencrypted, is not validated, and the default repo (`Mag14011/ProPrep`) differs
    from the update-check repo (`mjgplab/proprep`).

**Surprises:**
- The MODELLER licence warning appears on every launch.
- The menu-mode, save-default, previous-session and project-name answers are not recorded, so `undo` can't
  reach them.
- A ○ tool can still be opened, so users can look at its options before its prerequisites are met.
- Exiting from a submenu and answering `n` lands you on the main menu, not where you were.
