# Sessions, replay, templates, batch processing, command-line shortcuts, and the web interface (scope s2)

All runs: ProPrep 1.23.0, macOS, env `/opt/homebrew/Caskroom/miniforge/base/envs/ProPrep`, 21 Nov 2026.
Working folders: `manual/runs/s2/...`. Driver transcripts: `manual/runs/_sessions/s2-*/transcript.txt`.
`P` below stands for `/opt/homebrew/Caskroom/miniforge/base/envs/ProPrep/bin/proprep`.

Note on the driver: a pending prompt line (e.g. `Select chains ... (all):`) only shows in the driver output
*after* the next answer is sent, so in transcripts the prompt text appears just above the `>>> SENT` line that answered it.

---

## Feature: Automatic session recording (every interactive run)

**What it is for (plain English):** ProPrep writes down every answer you type, in order, in a small file.
You can later play that file back to repeat the same preparation, pick up where you left off, fix an earlier
answer, or turn it into a template for many proteins.
**Needs before you start:** nothing; it is on by default.
**Where to find it:** automatic. Flags: `--no-session`, `--session-dir DIR`, `--session-description TEXT`.

### Scenario 1: Record a short real session (load 1UBQ, filter out water, exit) — RUN
- Starting point: empty folder `runs/s2/rec1`; command `P --full-menu --session-description "s2 first recording 1UBQ"`.
- Steps:
  1. Prompt "Project name (.):" -> typed `ubq` (folder `rec1/ubq` created).
  2. Screen prints `Session recording to: ./proprep_session_20261121_122825.json`, plus
     `Tip: Use --no-session to disable automatic session recording` and
     `Tip: Type 'undo' at any prompt to rewind to an earlier answer`.
  3. "Enter your choice (1a):" -> `1a` (Structure Loader); "Enter your choice (1):" -> `1` (Load protein structure);
     "Select source (1):" -> `1` (RCSB); "Select method (1):" -> `1` (Enter PDB ID directly);
     "Enter PDB ID (e.g., 1ABC, 7XYZ):" -> `1UBQ`; "Select structure type (1):" -> Enter (asymmetric unit).
     Download + metadata panel shown. "Select source (1):" -> `6` (Done); then `m` (main menu).
  4. "Enter your choice (1a):" -> `4a` (PDB Filter); `1` (Filter PDB structure);
     "Detect redox-active sites before filtering? [y/n] (y):" -> `n`;
     "Would you like to view the structure in the interactive 3D viewer? [y/n] (n):" -> Enter (no);
     "Select chains (comma-separated indices, or 'all') (all):" -> Enter;
     "Do you want to proceed with filtering? (y):" -> `y`; Chain A Protein "Choose option (r):" -> `r` (retain);
     Chain A Water "Choose option (r):" -> `d` (discard); "Do you want to apply these filters? (y):" -> `y`.
     Message: `Filtered structure automatically saved to: filtered_structure.pdb`.
  5. `x` then "Are you sure you want to exit? [y/n] (n):" -> `y`. Screen: `Goodbye!` / `[Session recording stopped]` / citation panel.
- Result / files produced (`runs/s2/rec1/ubq/`):
  - `1UBQ.pdb` — downloaded structure.
  - `filtered_structure.pdb` — protein only, waters removed.
  - `proprep_session_20261121_122825.json` — the session recording (19 answers).
- Time: ~5 min of typing; program steps are instant.
- Problems: the session file goes into the **project folder** (the folder named at "Project name"), not the folder you launched from. The project-name answer and the menu-mode choice are *not* recorded (they happen before recording starts).
- Session log: `runs/_sessions/s2-rec1/transcript.txt`

**Session file structure (JSON), as observed:**
```
{
  "session_format_version": "1.1",
  "start_time": "2026-11-21T12:28:25.706085",
  "end_time":   "2026-11-21T12:33:08.667112",      (null if the program was killed)
  "metadata": {
     "proprep_version": "1.23.0",
     "command_line": ".../proprep --full-menu --session-description s2 first recording 1UBQ",
     "working_directory": ".../runs/s2/rec1/ubq",
     "auto_recorded": true,
     "web_shell_mode": false,               (true when recorded inside proprep-web)
     "description": "s2 first recording 1UBQ",   (only with --session-description)
     "pdb_id": "1UBQ" / "pdb_file": "/abs/path.pdb"   (only with --pdbid / --pdbfile)
     -- later tools add: continued_from, continuation_time, edited_from, edit_time,
        changes_applied {edits, deletes}, truncated_at, truncation_time, removed_interactions
  },
  "interactions": [  one object per answer, in order:
     { "index": 4, "timestamp": "...", "type": "prompt" | "confirm",
       "prompt": "Enter PDB ID (e.g., 1ABC, 7XYZ)", "response": "1UBQ",
       "choices": null or ["r","s","d"],
       "context": { "default": ..., "module": "Structure Loader - RCSB",
                    "description": "Enter PDB ID to download",
                    "options_map": {"r": "Retain entire component", ...},   (menus only)
                    "option_label": "Discard entire component" } }      (label of the chosen option)
  ]
}
```
Yes/no answers are stored as `"yes"`/`"no"` with `type: "confirm"`. The file is written incrementally, so a run
killed by Ctrl-C or by closing proprep-web still leaves a usable file (with `end_time: null`).

### Scenario 2: `--session-description` — RUN
- Text is stored in `metadata.description` (above) and is shown in the "Notes" column of the session picker
  (Scenario 5): `│ 1 │ proprep_session_20261121_122825.json │ 2026-11-21 12:28 │ 19 │ s2 first recording 1UBQ │`.

### Scenario 3: `--no-session` — RUN
- Starting point: `runs/s2/nosess` containing an old session file; `P --full-menu --no-session`.
- What the user sees: no "Previous session files found!" question, no "Session recording to:" line. Project prompt as normal.
- Typing `undo` at the main menu: `Invalid choice: undo. Please try again.` (undo needs recording).
- Result: no new session file written.
- Transcript: `runs/_sessions/s2-nosess/transcript.txt`

### Scenario 4: `--session-dir DIR` — RUN (absolute path works; relative or missing folder crashes)
- 4a absolute, existing folder: `P --full-menu --session-dir .../runs/s2/sdir/sessions_abs`, project `projA`.
  Screen: `Session recording to: /Users/.../runs/s2/sdir/sessions_abs/proprep_session_20261121_124045.json`.
  File landed there; `projA/` stayed empty. Relaunching with the same `--session-dir` offers the
  "Previous session files found!" menu (discovery looks in `--session-dir`).
- 4b relative folder (`--session-dir sessions_abs`, launched from `runs/s2/sdir`), project `projB`:
  **FAILED** — `Error: [Errno 2] No such file or directory: 'sessions_abs/proprep_session_20261121_124110.json'`
  `Run with --debug flag for more details`, then (wrongly) `Session saved to: sessions_abs/...`. ProPrep exits.
  Cause: the relative path is resolved *after* ProPrep moves into the project folder.
- 4c absolute but non-existent folder: same `[Errno 2] No such file or directory` crash right after the project prompt.
- Workaround: always give an absolute path to a folder that already exists.
- Transcripts: `runs/_sessions/s2-sdir1`, `s2-sdir2`, `s2-sdir3`.

### Options reference
| Flag | Effect |
|---|---|
| (none) | Records to `proprep_session_YYYYMMDD_HHMMSS.json` in the project folder. |
| `--no-session` / `--no-auto-record` | No recording and no "previous sessions" question; `undo` unavailable. |
| `--session-dir DIR` | Where to look for old sessions and where to write the new one. Use an absolute, existing folder. |
| `--session-description TEXT` | Free-text note stored in metadata and shown in the session picker. |

---

## Feature: Starting ProPrep in a folder that already has a session file

**What it is for:** lets you pick up an earlier session instead of starting again.
**Needs before you start:** a `proprep_session_*.json` in the folder you launch from (or in `--session-dir`).
Launching from the *parent* folder of a project does **not** find it (verified: launching in `runs/s2/rec1`
went straight to the project prompt although `rec1/ubq/` had a session).
**Where to find it:** automatic at startup (skipped with `--no-session`, `--analysis`, `--pdbview`, `--force-replay`, `--resume-session`).

### Scenario 5: Startup menu and its choices — RUN
- Starting point: copies of `rec1/ubq` in `runs/s2/disc1..3`; `P --full-menu`.
- What the user sees:
  ```
  Previous session files found!
  What would you like to do?
    1. Start fresh session (record new interactions)
    2. Resume a previous session from any point
    3. Continue without session recording
  Select option  [1/2/3] (1):
  ```
- 5a `1` (fresh): project prompt, new `proprep_session_20261121_123418.json` alongside the old one. (`disc2`)
- 5b `3` (no recording): project prompt, nothing recorded, `undo` rejected as invalid. (`disc3`)
- 5c `2` (resume): table "Available Session Files" (#, Filename, Date/Time, Interactions, Notes; up to 10 newest),
  prompt `Select session to resume (0 to cancel)  [1/0] (0):`. Note the default is **0 = cancel**.
  - `0` -> treated as "start fresh" (project prompt). (`s2-disc4`)
  - `1` -> session panel (`Interactions: 19`) and
    ```
    Options:
      1. Continue from end (append new interactions)
      2. Open session editor (edit/delete multiple interactions)
      3. Cancel
    ```
    `1` = Scenario 6, `2` = Scenario 8, `3` = start fresh.
- Transcripts: `runs/_sessions/s2-disc1`, `s2-disc2`, `s2-disc3`, `s2-disc4`.

### Scenario 6: "Continue from end" on a session that ended with Exit — RUN (behaves badly)
- Steps: `2` -> `1` -> `1`.
- What the user sees: `[Session backup: proprep_session_20261121_122825.json.bak-20261121_123345]`,
  `Hybrid mode started: ... Mode: Replay existing → Record new`, then every recorded answer flashes by as
  `Enter your choice [REPLAY: 1a]`, `Enter PDB ID ... [REPLAY: 1UBQ]` ... including
  `Enter your choice [REPLAY: x]` / `Are you sure you want to exit? [REPLAY: Yes]` -> `Goodbye!`.
- Result: ProPrep **replays the recorded exit and quits**; the user never gets control. A `.bak-` copy is made; files are re-created (1UBQ re-downloaded, filter re-run).
- How to get past it: use the session editor (Scenario 8) and delete the last two answers (`d 18-end`) first, or
  only use "Continue" on sessions that were interrupted (Ctrl-C / closed window), which have no recorded exit (Scenario 7).
- No project prompt in this mode: ProPrep works in the folder you launched from.
- Transcript: `runs/_sessions/s2-disc1/transcript.txt`

### Scenario 7: Continue an interrupted session (Ctrl-C) with `--resume-session` and `--demo-delay` — RUN
- Recorded: `P --full-menu` in `runs/s2/ctrlc`, Enter at project prompt, `1a`, then Ctrl-C ->
  `Application interrupted by user` / `[Session recording stopped]`. File holds 1 answer (`1a`).
- Then `P --full-menu --resume-session proprep_session_20261121_130253.json --demo-delay 1`, option `1`:
  `Enter your choice [REPLAY: 1a]` (1 s pause), `[Replay complete - now recording new interactions]` (printed twice),
  then live prompt `Enter your choice (1):` inside the Structure Loader. New answers (`x`, `yes`) were appended to the same file: `['1a', 'x', 'yes']`.
- Transcript: `runs/_sessions/s2-ctrlc2`, `s2-ctrlc3`.

---

## Feature: `--resume-session FILE` and the interactive Session Editor

**What it is for:** reopen a recorded session to continue it, or correct/remove earlier answers and re-run.
**Needs before you start:** a session file. Missing file -> log line `ERROR ... Session file not found: nope.json` and exit.
**Where to find it:** `P --resume-session FILE` (or startup menu option 2 -> pick -> option 2).

### Scenario 8: Edit an earlier answer and delete the recorded exit — RUN
- Starting point: copy in `runs/s2/resume1`; `P --full-menu --resume-session proprep_session_20261121_122825.json`; option `2`.
- What the user sees: panel `Session Editor | File: ... | 19 interactions | no changes`, a warning that changing a
  menu choice may make later answers not match ("ProPrep will fall through to live input"), and a table:
  `# | Status | Module | Question | Answer` (Answer shown in words, e.g. "Query Protein Data Bank (RCSB)"),
  footer `[e]dit #  [d]elete #|#-#|#-end  [u]ndo #|#-#  [n]ext  [b]ack  [j]ump #  [p]review  save  quit  help`, prompt `> :`.
- Steps:
  1. `help` -> command list panel (see Options reference).
  2. `e 16` -> `Editing #16: Chain A - Water retention option`, lists `r/s/d` with `<-- current` on `d`,
     `New value (d):` -> `r` -> `Marked #16 for edit: Discard entire component -> Retain entire component`; Status column shows `E`, answer shows "(was: ...)".
  3. `d 18-end` -> `Deleted #18-#19 (2 interactions) (use 'u 18-19' to restore)`; Status `D`, answer `[deleted]`.
  4. `u 19` -> `Restored #19 to original`; `d 19` -> `Deleted #19 ...`.
  5. `p` -> "Pending Changes" table (EDIT #16, DELETE #18, DELETE #19).
  6. Invalid: `e 99` -> `Invalid index. Valid range: 1-19`; `foo` -> `Unknown command: foo. Type 'help' for commands.`; `n` -> `Already on the last page`.
  7. `save` -> reminder about path-changing edits, `Continue with save? [y/n] (y):` -> `y`.
- What happens next: `Session saved: proprep_session_20261121_122825_edited.json` / `1 edit(s), 2 deletion(s) applied`,
  backup made, hybrid replay of the edited file (water now `[REPLAY: r]`, filter summary shows `Water: 58 residues`),
  then `[Replay complete - now recording new interactions]` and a **live** PDB Filter menu. Option `2` (status) confirmed 58 waters kept.
- Files: original untouched; new `*_edited.json`; one `*.bak-<time>` per resume/rewind.
- Transcript: `runs/_sessions/s2-resume1/transcript.txt`

### Scenario 9: Resume -> option 3 Cancel — RUN
- `P --full-menu --resume-session s.json`, option `3` -> ProPrep starts a normal fresh session (project prompt). (`s2-resume2`)

### Options reference (Session Editor)
| Command | Meaning |
|---|---|
| `e <#>` | Change the answer of item # (shows the options for menus). |
| `d <#>`, `d <#>-<#>`, `d <#>-end` | Skip item(s) during replay. |
| `u <#>`, `u <#>-<#>` | Undo a pending edit/delete. |
| `v <#>` | (template mode only) mark as variable; in edit mode: "Variable marking is only available in template mode". |
| `n` / `b` / `j <#>` | Next page / previous page / jump to page holding #. |
| `p` | Preview pending changes. |
| `save` | Write `<name>_edited.json` (edit mode) and replay it, or `<name>_template.json` (template mode). With no changes: "No changes to save" + "Save unchanged copy anyway?". |
| `quit` | Leave; asks "Discard unsaved changes? (n)" if anything is pending. In resume mode, quitting = start fresh. |
| `help` | Command list. |
| Resume options | 1 Continue from end · 2 Open session editor · 3 Cancel (start fresh). |

---

## Feature: Typing `undo` at any prompt (rewind)

**What it is for:** go back to an earlier question without restarting ProPrep.
**Needs before you start:** session recording must be on (not with `--no-session`, not in `--force-replay`).
**Where to find it:** type `undo` at any recorded prompt.

### Scenario 10: `undo` -> pick an answer -> re-answer it live — RUN
- Starting point: live part of Scenario 8 (after viewing filter status).
- Steps: at `Enter your choice (1):` typed `undo` -> table "Recent answers in this session" (# / module / question / answer, 12 per page),
  prompt `Rewind to which answer? number to rewind to, [m]ore, or [x] cancel (17):` -> `15` (water option) ->
  `[r]e-answer it live, or [c]hange it and replay what followed? ([x] cancel) (r):` -> `r`.
- What the user sees: `[Rewinding: question 15 will be asked again; 3 answers from there on are set aside in a backup of the session log]`,
  `Rewinding session: ...`, backup line, fast replay up to the water question, then the question asked live `Choose option (r):` -> typed `d`, `y` to apply.
- Result: waters removed again; session file trimmed and continued.

### Scenario 11: `undo` -> [c]hange an answer and replay what followed — RUN (BUG: change ignored)
- Steps: `undo` -> `15` -> `c` -> `New answer for "Chain A - Water retention option":` -> `q` (rejected: `Answer must be one of: r, s, d`) -> `r`.
- What the user sees during replay: `Choose option [REPLAY: r]` immediately followed by
  `[REPLAY: 'r' now means a different option; using 'd' → Discard entire component]`.
- Result: the **old answer (discard) was used**; filter status shows no waters. In the saved file item 15 has
  `response: "r"` but `option_label: "Discard entire component"`, and ProPrep's label-matching on replay switches back to `d`.
  The Session Editor (`e #`) does update the label and works; `undo`+`c` does not (for menu-type questions).
- Workaround: use `undo` -> `r` (re-answer live), or the Session Editor.

### Scenario 12: `undo` edge cases — RUN
- `undo` typed at the **Project name** prompt (before recording starts) is taken literally: a project folder named `undo` was created.
- At the first menu after a typo: the typo itself (`u1`) had been recorded as answer 0.
- `x` -> `Rewind cancelled.` and the question is asked again; `m` -> shows the previous page (wraps around); `7` when only 0-1 exist -> `Enter a number between 0 and 1.`
- Rewind to `0` + `r` -> `[Session truncated to nothing - the first question will be asked live]`, `Edit at: interaction -1 (truncated)`.
- With `--no-session`: `Invalid choice: undo` (source message for other recorder-less cases: "Rewind needs an active session recording, and this run has none (replay-only or --no-session)").
- Transcript: `runs/_sessions/s2-undo/transcript.txt`

---

## Feature: `--force-replay FILE` and `--demo-delay SECONDS`

**What it is for:** play back a recorded session exactly, hands-off (e.g. to reproduce a preparation or demo it).
**Needs before you start:** session file; the structure file if the session was started with `--pdbfile`.
**Where to find it:** `P --force-replay FILE [--demo-delay 0.5]`.

### Scenario 13: Pure replay with a 0.5 s delay — RUN
- `P --force-replay proprep_session_20261121_122825.json --demo-delay 0.5` in `runs/s2/force1` (session file only).
- Without `--full-menu`/`--workflow` the menu-mode question still appears (`Select menu mode [1/2]:` -> `2`,
  `Save this as your default preference? (n)` -> `n`) — it is not part of the recording.
- Then: `[Session replay started: ...]`, `Replaying session from: ...`, `Demo mode: 0.5s delay between actions`;
  prompts show the recorded answer in brackets, e.g. `Enter your choice [1a]`, `Are you sure you want to exit? [Yes]`, then `Goodbye!` / `[Session replay stopped]`.
- Result: `1UBQ.pdb` and `filtered_structure.pdb` re-created in the launch folder; session file unchanged; no new session file. ~10 s.
- Transcript: `runs/_sessions/s2-force1/transcript.txt`

### Scenario 14: Replay a session that was started with `--pdbfile` — RUN
- File present: `Loading the structure this session was recorded with: my_protein.pdb` / `✓ Loaded structure from local file`, replay completes.
- File moved away: `The recorded session was started with --pdbfile my_protein.pdb, but that file is not here (looked for /.../my_protein.pdb and my_protein.pdb in the project directory). Copy it in, or replay with --pdbfile <path>.` Replay continues with no structure (all tools show "Needs a loaded structure").

### Scenario 15: Errors / orphan flags — RUN
- `--force-replay nope.json` -> log line `ERROR ... Session file not found: nope.json`, exits (citation panel still printed).
- `--demo-delay 1` without replay -> `Warning: --demo-delay has no effect without --force-replay or --resume-session`.

---

## Feature: Templates (`--create-template`, `--template-info`, `--validate-template`)

**What it is for:** turn one recorded session into a recipe with blanks ("variables"), e.g. the PDB ID, so the
same steps can be run on many proteins.
**Needs before you start:** a recorded session; for batch use, a list of values for every variable.
**Where to find it:** `P --create-template SESSION.json`, `P --template-info T.json`, `P --validate-template T.json [--batch-list LIST]`.

### Scenario 16: Create a two-variable template plus a batch CSV — RUN
- `P --create-template ubq_session.json` (copy of the 1UBQ session) in `runs/s2/tmpl`.
- What the user sees: "Template Creator" panel and the same table as the Session Editor, footer adds `[v]ariable #`.
  **Nothing is pre-marked** — the PDB ID is not detected automatically in 1.23.0 (auto-detection is switched off in the code); you must mark it yourself.
- Steps:
  1. `v 5` -> `Marking #5 as template variable / Interaction: Enter PDB ID to download / Original value: 1UBQ (1UBQ)`;
     `Variable name (input_protein):` -> Enter; `Description (Enter PDB ID to download):` -> Enter; `Required? [y/n] (y):` -> `y`. Row shows `V` and `{{ input_protein }}`.
  2. `v 16` -> shows options r/s/d (`<-- original` on d); `Variable name (var_16):` -> `water_option`;
     description -> `Keep (r) or discard (d) crystal waters`; required `y` -> `Marked as variable: {{ water_option }}`.
  3. `p` -> Pending Changes lists both as VARIABLE. `save` ->
     `Template saved: ubq_session_template.json` / `2 variable(s), 0 edit(s), 0 deletion(s)` /
     `*input_protein: Enter PDB ID to download (example: 1UBQ)` / `*water_option: ... (example: d)`.
  4. `Create a batch input file now? [y/n] (y):` -> `y`. For each run it asks each variable (choice variables list their options):
     Run 1 `input_protein (1UBQ):` Enter, `water_option (d):` `d`; Run 2 `1AKI`, `r`; Run 3 `1CRN`, `d`; Run 4 Enter on empty = finish.
     `Batch file saved: ubq_session_batch.csv (3 runs)` / `Run with: proprep --batch-replay ubq_session_template.json --batch-list ubq_session_batch.csv`.
- Files: `ubq_session_template.json` (session + `"template": true`, `"session_format_version": "1.2"`, `template_variables`
  {description, type text|choice, example, required, choices, options_map, original_label}, `template_metadata`;
  marked answers become `"{{ input_protein }}"`), `ubq_session_batch.csv` (`input_protein,water_option` header + 3 rows).
- Transcript: `runs/_sessions/s2-tmpl1/transcript.txt`

### Scenario 17: One-variable template, quit/save paths — RUN
- `save` with no variables -> `No variables marked. A template needs at least one variable.` / `Save as edited session instead? (y)` -> `n` -> program ends, nothing saved. (`s2-tmpl2`)
- `quit` with a marked variable -> `Discard unsaved changes? (n)` -> `n` returns to the editor; `save` -> `filter_session_template.json`;
  `Create a batch input file now?` -> `n` -> `You can create one later. Format: CSV with headers matching variable names.` (`s2-tmpl5`, `runs/s2/tmpl1var`)
- Pitfall seen: after `e 5` the next thing typed is the new value — typing `quit` there stored "quit" as the answer.

### Scenario 18: `--template-info` — RUN
- `P --template-info ubq_session_template.json`:
  ```
  Template Information
  File: ubq_session_template.json
  Version: unknown
  Created: 2026-11-21T12:43:15.939196
  Description: ProPrep batch processing template
  Template Variables:  input_protein | text | Unknown | 1UBQ ;  water_option | choice | Unknown | d
  Total Interactions: 19
  ```
  ("Version: unknown" and "Source: Unknown" always, for editor-made templates.)
- On a plain session file: same header with "Created: unknown", no variable table, exit 0 (no error).

### Scenario 19: `--validate-template` with and without `--batch-list` — RUN
| Command (`P --validate-template ubq_session_template.json ...`) | Output | Exit |
|---|---|---|
| (no list) | `✓ Template is valid: ubq_session_template.json` | 0 |
| `--batch-list ubq_session_batch.csv` | `✓ Template is valid` | 0 |
| `--batch-list ids.txt` (plain IDs: 1UBQ, 1CRN) | `✗ Template validation failed` / `• Variable 'water_option' (substituted into interaction 15): no run supplies a value` | 1 |
| `--batch-list missing_var.csv` (header `input_protein` only) | same "no run supplies a value" | 1 |
| `--batch-list partial.csv` (row 2 water blank) | `• Variable 'water_option' ...: no value in run(s) 2` | 1 |
| `--batch-list blank.csv` (only row has empty water) | "no run supplies a value" | 1 |
| `--batch-list extra.csv` (extra column) | `! Input column 'extra' matches no template variable and will be ignored (check the header for a typo)` then `✓ valid` | 0 |
| `--batch-list badchoice.csv` (water = `x`) | `✓ valid` — choice values are **not** checked | 0 |
| `--batch-list nonexist.csv` | `• Could not check nonexist.csv: [Errno 2] No such file or directory` | 1 |
| plain session file instead of template | `• Missing 'template' flag (should be true)` / `• Missing 'template_variables' section` | 1 |
| missing template file | `• Template file not found: nofile.json` | 1 |
- One-variable template + `ids.txt` -> valid (plain list lines fill `input_protein`).
- **BUG:** plain list of UniProt IDs (`P69905`, `P68871`) is mis-read as a CSV with header `P69905`:
  `! Input column 'P69905' matches no template variable ...` / `• Variable 'input_protein' ...: no run supplies a value`.
  Lists whose first line starts with a letter are taken as a header. Workaround: put `input_protein` as the first line.
- `--batch-list` without `--batch-replay`/`--validate-template` -> `Warning: --batch-list has no effect without --batch-replay or --validate-template`.

---

## Feature: Batch processing (`--batch-replay`, `--batch-list`, `--batch-continue`, `--web`/`--no-web`, `--output-dir`)

**What it is for:** run the same recorded recipe on a list of proteins without typing.
**Needs before you start:** a template; a list file (CSV with a header row of variable names, or plain text one PDB ID per line for one-variable templates; `#` lines are comments); an **existing** output folder.
**Where to find it:** `P --batch-replay T.json --batch-list LIST [--output-dir DIR] [--batch-continue] [--web|--no-web]`.

### Scenario 20: Three proteins, stop at first failure — RUN
- `--output-dir batch_stop` that did **not** exist -> **FAILED**: every run `Setup error: [Errno 2] No such file or directory: 'batch_stop/run_01'`, then a Python traceback `FileNotFoundError ... batch_report_...json`. Create the folder first (`mkdir`).
- With the folder created: `P --batch-replay ubq_session_template.json --batch-list ubq_session_batch.csv --output-dir batch_stop`.
- What the user sees: `═══ Batch Processing Mode ═══`, template, `Runs to process: 3`, variables list, `Replay mode: plain terminal`,
  a progress bar `Processing 1UBQ (1/3)`, a leaked traceback for the failing run, `Stopping batch due to error`, then the summary:
  `Planned: 3 / Successful: 2 (66.7%) / Failed: 1 (33.3%)`, results table (1UBQ ✓, 1AKI ✓, 1CRN ✗ "Replay diverged: no recorded answer for prompt 'Enter your choice'. At this point the recording expects a prompt here: 'Choose option'."),
  `Batch report saved to: batch_stop/batch_report_20261121_125028.json`,
  `Runs still owed (1): batch_stop/batch_retry_20261121_125028.csv` and the exact re-run command.
- Why 1CRN failed: crambin has no crystal waters, so the "Water" question the recipe expects never comes. Batch only works for proteins that ask the *same questions* as the recorded one.
- Files: `run_01/` (1UBQ.pdb, filtered_structure.pdb — 0 waters, run.log), `run_02/` (1AKI, water_option=r -> 79 HOH kept), `run_03/` (1CRN.pdb, run.log, error.log with traceback), `batch_report_*.json` (batch_info + per-run variables/success/error), `batch_retry_*.csv` (failed + not-attempted rows, ready for `--batch-list`).
- Time: ~1 s for 3 proteins (download + filter).
- If the failing run is first and `--batch-continue` is absent, the rest are "Not attempted" and the summary says `Batch stopped after run 1 of 3. Pass --batch-continue to work through failures instead of stopping at the first one.`

### Scenario 21: `--batch-continue` — RUN
- `cont.csv` = 1CRN, 1UBQ, 1AKI. Result: 1CRN ✗, 1UBQ ✓, 1AKI ✓; exit code 1 (any failure). Console: `runs/s2/tmpl/batch_cont_console.txt`.
- `--batch-continue` without `--batch-replay` -> `Warning: --batch-continue has no effect without --batch-replay`.

### Scenario 22: Plain-text ID list, all succeed — RUN
- One-variable template, `ids_ok.txt` (`# comment`, 1UBQ, 1AKI), `--output-dir b_ok` -> `Successful: 2 (100.0%)`, exit 0, no retry file.
- `run.log` in each run folder holds the full (colourless) ProPrep screen output for that protein.

### Scenario 23: `--web` / `--no-web` — RUN
- Terminal-recorded template with `--web`: `Replay mode: proprep-web (browser)`; 1UBQ fails:
  `Replay diverged: no recorded answer for prompt "Select chains ..."; the recording expects a confirm here: 'Would you like to view the structure in the interactive 3D viewer?'`
  (in browser mode that viewer question is skipped, so the recording gets out of step).
- `--no-web`: `Replay mode: plain terminal`, both succeed.
- Default: mode is taken from `metadata.web_shell_mode` of the recording (a session recorded in proprep-web had `"web_shell_mode": true`).
  The intended `[override]` note after "Replay mode:" is never shown (swallowed as formatting).
- `--web` without `--batch-replay` -> `Warning: --web/--no-web has no effect without --batch-replay`.

### Scenario 24: Batch input errors — RUN
| Situation | Message | Exit |
|---|---|---|
| no `--batch-list` | `Error: --batch-replay requires --batch-list` + usage line | 1 |
| plain session, not a template | `File is not a template: filter_session.json` | 1 |
| list missing a variable (`miss.csv`) | `Input list does not satisfy ubq_session_template.json:` / `• Variable 'water_option' ...: no run supplies a value` (nothing runs) | — |
| UniProt plain list | as Scenario 19 bug | — |
| empty list (only comments) | `No proteins found in input list: empty.txt` | 1 |
| invalid choice value `x` (1AKI) | run "succeeds": replay types `Choose option [x]`, ProPrep silently keeps the waters (treated like retain). No warning. | 0 |

### Options reference
| Flag | Meaning |
|---|---|
| `--batch-replay T.json` | Template to run. |
| `--batch-list FILE` | CSV (header = variable names) or plain list (one ID per line, binds `input_protein`). |
| `--output-dir DIR` | Parent folder for `run_01`, `run_02`, ... and reports (default: current folder). Must already exist. |
| `--batch-continue` | Keep going after a failed protein. |
| `--web` / `--no-web` | Force browser-mode / terminal-mode replay (default: as recorded). |

---

## Feature: Command-line shortcuts (`--pdbid`, `--pdbfile`, `--output-dir`, `--analysis`, `--pdbview`)

**What it is for:** skip the menus for common starts.
**Where to find it:** `P --help`, "workflow shortcuts" section.

### Scenario 25: `--pdbid 1UBQ` — RUN
- `P --full-menu --pdbid 1UBQ` in `runs/s2/short_pdbid`: **no project-name prompt** (current folder is used),
  `Session recording to: ./proprep_session_20261121_125144.json`, `Downloading 1UBQ...`, metadata panel, main menu with the structure loaded (most tools ticked ✓).
- Session metadata gets `"pdb_id": "1UBQ"`; the load itself is not a recorded answer, replay redoes it from metadata.
- Bad ID `--pdbid ZZZZ`: `PDB ZZZZ not found (404). May be a large structure only available in mmCIF format.` / `Failed to download PDB ZZZZ`, then the main menu with no structure (log line still says "Downloaded and loaded PDB: ZZZZ").
- `--pdbid` and `--pdbfile` together: `--pdbid` wins silently.

### Scenario 26: `--pdbid` + `--output-dir` — RUN
- `P --full-menu --pdbid 1CRN --output-dir crambin_project/sub`: nested folders created, `1CRN.pdb` and the session file saved there.
- Log line shows a doubled path: `Working directory set to: .../crambin_project/sub/crambin_project/sub` (the path stored as "project directory" inside ProPrep is wrong; files themselves are in the right place).

### Scenario 27: `--pdbfile my_protein.pdb` — RUN
- `✓ Loaded structure from local file / File: my_protein.pdb`, no project prompt, metadata `pdb_file` (absolute path).
- Missing file: log `ERROR ... PDB file not found: /.../nothere.pdb` and ProPrep exits.

### Scenario 28: `--analysis` — RUN
- `P --full-menu --analysis`: no project prompt, recording on, `Jumping to Simulation Analysis...`,
  `No simulations directory found` / `No completed simulations found in simulations/ directory`, then the "MD ANALYSIS FILE BROWSER"
  (help on what to select: .mdout/.out for energetics; .nc + .prmtop/.parm7 for trajectories; commands `select N`, `add N`, `done`, `clear`, `cd`, `find`, `q`).
  `q` -> main menu. (The analysis itself belongs to the MD scope.)
- `--analysis --pdbid 1UBQ` -> `Error: --analysis cannot be combined with --pdbid or --pdbfile`.

### Scenario 29: `--pdbview 1UBQ` / `--pdbview file.pdb` — RUN
- `P --pdbview 1UBQ`: still asks the menu-mode question unless `--full-menu`/`--workflow` is given. Then downloads 1UBQ and prints
  `✓ Structure viewer launched at http://localhost:8766/viewer (1 structure)` (tries to open the browser), and **drops into the full main menu** with the structure loaded.
  `curl http://localhost:8766/viewer` -> 200, page title "ProPrep Interactive Structure Viewer". No session recorded (view-only).
- Exiting ProPrep (`x`, `y`) stops the viewer (port no longer listening).
- `--pdbview 1UBQ.pdb` (local file) -> viewer on port 8767 (port number varies). Missing file -> `File not found: /.../nothere.pdb`, then main menu.

### Options reference
| Flag | Behaviour |
|---|---|
| `--pdbid ID` | Download from RCSB at start; current folder is the project; recorded in metadata. |
| `--pdbfile PATH` | Load a local PDB at start; must exist. |
| `--output-dir DIR` | Use/create DIR as project folder (no project prompt). |
| `--analysis` | Open the MD analysis file browser; no session discovery. |
| `--pdbview ID|FILE` | Load + open web viewer; no recording; remains in the menu. |

---

## Feature: `proprep-web` (browser-based ProPrep)

**What it is for:** the same ProPrep, but inside a browser tab, with the 3-D structure viewer docked beside the terminal
and a "Download project" button. Also a hosted "workshop seats" mode for teaching.
**Needs before you start:** ProPrep env; internet access in the browser (the page loads the xterm terminal library from cdn.jsdelivr.net).
**Where to find it:** command `proprep-web` (same bin folder).

### Scenario 30: `proprep-web --help` — RUN
Options: `THEME` (positional: dark/light/black/white), `--host` (default 127.0.0.1), `--port` (default 8000; scans upward if busy),
`--strict-port`, `--no-browser`, `--theme`, `--font-size PX` (8-32, default 13), `--high-contrast`, `--no-transcript`
(otherwise each session is logged to `~/.proprep/web_sessions/session_<timestamp>.log`), `--log-level`;
hosted mode: `--seats N`, `--seats-dir DIR` (default `~/.proprep/web_seats`), `--public-url URL`.

### Scenario 31: Start, open, use, and stop (local mode) — RUN
- Started in `runs/s2/web`: `proprep-web --no-browser --port 8420`:
  ```
  ProPrep Web Shell listening on http://127.0.0.1:8420/
    - Terminal pane runs the unmodified ProPrep CLI
    - Structure Viewer will dock on the right when launched from the CLI
    - Session transcript: ~/.proprep/web_sessions/session_<timestamp>.log
  INFO:     Uvicorn running on http://127.0.0.1:8420 (Press CTRL+C to quit)
  ```
  (Without `--no-browser` it opens the default browser at that URL.)
- `curl /` -> 200, 2.7 KB page "ProPrep Web Shell": top bar with "Download project", font A−/A+, "Detach viewer", status; left terminal pane; right pane placeholder "Structure Viewer — Launch the Interactive Structure Viewer from the CLI. It will appear here automatically."
  `/healthz` -> `ok`; `/seat-info` -> `{"name":"local","hosted":false,...,"cwd":".../runs/s2/web"}`; `/version` before a viewer -> 503 "viewer not running".
- Driving the terminal (via its websocket, standing in for a browser tab): the normal ProPrep start screens appear
  (menu-mode question, project name). I loaded 1UBQ and ran the PDB Filter in project `webproj`. Differences from the plain terminal:
  - on loading, `✓ Structure viewer launched at http://localhost:8765/viewer` — the viewer opens automatically and docks in the right pane (server log: `viewer server announced on port 8765`);
  - the PDB Filter does **not** ask "Would you like to view the structure...?"; it prints `Showing .../1UBQ.pdb in the docked viewer...`;
  - the session file is marked `"web_shell_mode": true`.
- With a tab open: `/viewer` through the web shell -> 200 (the docked viewer page), `/download` -> `proprep_local_<time>.zip` containing the project folder (e.g. `w3proj/1CRN.pdb`, session json).
- Stopping: in local mode the server **shuts itself down ~2 s after the last browser tab closes**
  (`No browser reconnected within 2.0s; shutting down.`), and the ProPrep inside it is ended (session file left with `end_time: null`). Ctrl-C in the launching terminal also stops it.
- Port 8420 still busy right after a shutdown -> `(port 8420 was busy; using 8421)`.
- Files: `runs/s2/web/webproj/` (1UBQ.pdb, filtered_structure.pdb, proprep_session_20261121_125705.json), web transcript copy `runs/s2/web/session_20261121_125649_108527.log`, client transcript `runs/s2/web/ws_client_transcript.txt`.
- Driver transcripts: `runs/_sessions/s2-web`, `s2-web2`, `s2-web3`.

### Scenario 32: `--strict-port` on a busy port — RUN
- Exit code 3, but it first prints `ProPrep Web Shell listening on http://127.0.0.1:8440/` and only then
  `ERROR: [Errno 48] error while attempting to bind on address ('127.0.0.1', 8440): address already in use` (misleading first line).

### Scenario 33: Hosted workshop seats — RUN (start/stop only)
- `proprep-web --seats 2 --seats-dir .../seatdir --port 8450 --strict-port --public-url https://workshop.example.org`:
  ```
  - Hosted mode: 2 seat(s); seat directories under .../seatdir
  - Sessions persist across browser disconnects; no auto-shutdown
  - Seat links (one per attendee; the token is the only access control):
      seat 01: https://workshop.example.org/?seat=<token>
      seat 02: https://workshop.example.org/?seat=<token>
  ```
  Created `seat_01/`, `seat_02/`, `seats.json`. `/` without token -> 403 "This ProPrep link is not valid. Ask the workshop host for your seat link."; with token -> 303 redirect into the seat. Stopped with Ctrl-C.

---

## Feature: `proprep-library` (user parameter library snapshots)

**What it is for:** keep, clear and bring back the personal library ProPrep builds up in `~/.proprep`
(`forcefield_params/` = saved force-field parameter sets, `transformers/`). Nothing is ever deleted: reset/restore first take a snapshot.
Settings, licence keys, templates, web sessions are not touched.
**Where to find it:** `proprep-library {status,snapshot,list,reset,restore}`, option `--root DIR` (default `$PROPREP_USER_DIR` or `~/.proprep`).

### Scenario 34: The `proprep-library` command itself — FAILED
- `proprep-library --help` -> `bad interpreter: .../envs/ProPrep/bin/python3.13: no such file or directory`.
  The launcher script points at Python 3.13 but the environment has 3.12 (`proprep-web` correctly points to 3.12).
- Workaround used for everything below: `.../envs/ProPrep/bin/python -m proprep.utils.user_library_snapshots ...` (identical behaviour).

### Scenario 35: Read-only commands on the real library — RUN
- `status`: `Library root: /Users/.../.proprep` / `forcefield_params 0 files 0.0 KB` / `transformers 0 files 0.0 KB` / `snapshots: 0 in .../library_snapshots`.
- `list`: `No snapshots in /Users/.../.proprep/library_snapshots`.

### Scenario 36: snapshot / reset / restore on a throw-away library (`--root runs/s2/libroot`) — RUN
- `snapshot before_test -m "sandbox test"` -> `Snapshot written: .../library_snapshots/before_test`.
- `list` -> `before_test  2026-11-21T12:54:58  forcefield_params: 1 file, transformers: 1 file  — sandbox test`.
- `reset -y -m clean` -> `Previous library saved to: .../before_reset_20261121_125458` / `Library is now empty.` / `Restart ProPrep so open sessions see the change.`
- `restore before_test -y` -> `Library restored from 'before_test'.` (an empty library is not snapshotted first).
- `restore nosuch -y` -> `error: no snapshot named 'nosuch' in ...`, exit 1.
- Without `-y`: `About to empty the library (its current contents are snapshotted first). Continue? [y/N]` -> `n` -> `Nothing changed.` (same for restore).
- Each snapshot folder: `forcefield_params/`, `transformers/`, `manifest.json`.

| Sub-command | Options | Effect |
|---|---|---|
| `status` | — | Show file counts/sizes and number of snapshots. |
| `list` | — | List snapshots with date, contents, note. |
| `snapshot [name]` | `-m NOTE` | Copy current library (default name = timestamp). |
| `reset` | `-m NOTE`, `-y` | Snapshot, then empty the library. |
| `restore name` | `-y` | Snapshot current (if not empty), then replace with snapshot. |

---

## Notes for the manual writer

**Jargon to explain:** session file / recording; replay; "hybrid mode" (replay old answers, then record new ones live);
template and "variable" (a blank filled per protein); batch list / CSV header; "replay diverged" (the program asked a question the
recording did not expect — usually because the new protein differs, e.g. no waters, extra chains or ligands); asymmetric unit; port/localhost; seat token.

**Recommended beginner path:**
1. Do one protein by hand (session is recorded automatically). Note the `Session recording to:` line.
2. `proprep --create-template that_session.json`, type `v <#>` on the "Enter PDB ID" row, `save`, and answer `y` to create the batch file.
3. `mkdir batch_out`, then `proprep --validate-template T.json --batch-list T_batch.csv`, then
   `proprep --batch-replay T.json --batch-list T_batch.csv --output-dir batch_out --batch-continue`.
4. Read `batch_out/run_XX/run.log` / `error.log` for failures; re-run with the printed `batch_retry_*.csv`.
Choose proteins similar to the recorded one (same number of chains, same kinds of components).

**Common mistakes / surprises:**
- Session files live in the project subfolder; launch ProPrep from *inside* that folder to be offered "Resume".
- "Continue from end" on a session you exited normally just replays your exit and quits. Delete the last two items in the editor first.
- In the resume picker the default answer is `0` (cancel) — pressing Enter starts a fresh session.
- `undo` does nothing before the project prompt (it becomes the project name) or with `--no-session`.
- The menu-mode question is not recorded; add `--full-menu` (or `--workflow`) to replay/`--pdbview` commands to avoid it.
- Batch output folder must exist beforehand.
- Templates recorded in proprep-web must be replayed in web mode (automatic) — mixing modes causes a divergence at the PDB Filter viewer question.
- PDB IDs in a plain list are fine; UniProt/AlphaFold IDs in a plain list are mistaken for a header.

**Looked like bugs (describe, not fixed):**
1. `undo` -> `[c]hange` on a menu question is reverted during replay (`'r' now means a different option; using 'd'`) because `option_label` is not updated. Session Editor edits work.
2. "Continue from end" does not strip the recorded Exit/confirm, so it exits immediately.
3. Relative or non-existent `--session-dir` crashes after the project prompt and then claims "Session saved to:".
4. `--batch-replay --output-dir X` crashes if X does not exist (no `mkdir -p`), with a raw Python traceback.
5. Plain-text batch list beginning with a letter (UniProt IDs) is parsed as CSV.
6. `proprep-library` launcher has a wrong interpreter (python3.13) and does not run.
7. `--output-dir a/b` stores a doubled path as the project directory (log shows `.../a/b/a/b`).
8. Batch replay does not validate choice values (`x` accepted silently); `--validate-template` does not check them either.
9. Batch failures leak full tracebacks/log lines onto the progress display; `[override]` note in "Replay mode" never shows.
10. `--strict-port` prints "listening on ..." before failing; `[Replay complete ...]` is printed twice; on a normal exit the "Session saved to:" line never appears (only "[Session recording stopped]").
11. An invalid menu entry (typo) is recorded as an answer.
12. One driver run of the template creator (`s2-tmpl4`) died silently mid-editor; could not reproduce on retry (`s2-tmpl5` worked) — probably a driver glitch, mentioned for completeness.
