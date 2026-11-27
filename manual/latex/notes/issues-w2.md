# User-facing problems mentioned by w2 (Ch. 4, 5, 6, 23, 26, App. B)

Format: Symptom | Why | What to do | Chapter label

## Starting, exiting, resuming
- **Pressing Enter at "Select menu mode [1/2]" just repeats the question** | The question has no default | Type 1 or 2 (2 = Full Menu, used in the manual), or start with `--full-menu` / `--workflow` | ch:starting
- **ProPrep stops at start-up with `Error: [Errno 2] No such file or directory: 'sub/deeper'`** | A project name containing `/` asks for two folder levels; only one can be created | Use a simple name without `/` (letters, digits, `-`, `_`); avoid spaces too | ch:starting
- **A project folder called `undo` appeared** | `undo` typed at "Project name" is taken literally (recording has not started yet) | Delete the folder; undo only works after the project question | ch:starting
- **Closing the browser tab ended ProPrep without asking** | proprep-web shuts down ~2 s after the last tab closes; the session file is left with `end_time: null` | Exit with `x` first, download the project, then close the tab | ch:starting
- **"Download project" zip is huge / contains unrelated files** | The zip holds everything in the folder `proprep-web` was started from (source: server.py `_zip_directory`) | Start proprep-web from a dedicated work folder, not the home folder | ch:starting
- **ProPrep did not offer to resume my project** | Session files are looked for only in the folder ProPrep is started from; launching from the parent folder does not find them | `cd` into the project folder first, then start ProPrep | ch:starting
- **Pressing Enter in the session picker started a fresh session** | Default of "Select session to resume (0 to cancel) [1/0] (0)" is 0 = cancel | Type the session number | ch:starting, ch:sessions
- **"Continue from end" replays my session and ProPrep quits immediately** | The recorded `x` / `y` exit is replayed too | Open the session editor, delete the last two answers (`d <n>-end`), save; or use Continue only for interrupted sessions | ch:starting, ch:sessions
- **"Session saved to: <file>" never appears at exit** | Message not printed in 1.23.0 (only "[Session recording stopped]") | Nothing to do; the session file is saved incrementally | ch:starting
- **MODELLER warning on every start** | Printed whenever MODELLER has no valid licence key | Ignore unless loop/residue repair is needed; add key (Ch. 3) | ch:starting

## Menus
- **`W`, `X`, `4E` rejected as "Invalid choice"** | Full Menu is case-sensitive (Workflow Mode is not) | Type lowercase | ch:menus
- **Typing `1a` inside a tool menu gives "Received '1a', which is not one of the available options"** | Main-menu codes do not work inside a tool menu | Press `m` first | ch:menus
- **Answering `n` to "Are you sure you want to exit?" inside a tool lands on the main menu** | Exit-cancel always returns to the main menu | Re-open the tool | ch:menus
- **Main menu still shows plain ✓ after a step was done** | ✓ means "available", not "done"; only tool menus show ● for done | Open the tool and look for ● | ch:menus
- **Help page (h) codes open the wrong tool (e.g. 3a opens Homology Searcher)** | Help & Workflows lists pre-1.23 codes (3a Filter, 3c Fix, 3d Protonate, 4a RedoxSite, 5a FF-Params, 5b tLEaP, 5c MD, 5d QM/MM, 2b Aligner) | Use the corrected table in Ch. 5 / Appendix B | ch:menus, app:menumap
- **Workflow help says `s` skips a stage and lists a "Setup" stage** | `s` jumps to Simulate; there is no Setup stage in the guide bar | Ignore; topology is on Simulate | ch:menus
- **`undo` → `[c]hange` on a menu question replays the old answer** (`'r' now means a different option; using 'd'`) | option label is not updated | Use `undo` → `r` and answer live, or the session editor | ch:menus, ch:sessions
- **`undo` rejected as "Invalid choice: undo"** | No session recording (`--no-session` or option 3 at start-up) | Start with recording on | ch:menus
- **"Save workspace" file contains only `{"pdb_file": null}`; "Load workspace" restores nothing** | Saves obsolete keys | Rely on session files and project-folder files | ch:menus
- **"Show detailed workspace" identical to "Show workspace status"** | detailed flag unused | — | ch:menus
- **"Load command history" says "Loaded 14 command records" but history shows one** | only breadcrumbs restored | — | ch:menus
- **After "Reset workspace" breadcrumb shows "Navigation: Home"; Workflow Mode drops to Full Menu** | Reset also clears menu mode and project directory | Press `g` to return to Workflow Mode | ch:menus
- **Workflow Mode: stage marks wrong** (leaving loader with `m` does not advance; `n` marks ✓ even if nothing done; after `b` stage shows ✓ not ▶; `a`→`g` loses progress; "Advanced to Detect stage" printed twice) | Guide-bar bookkeeping | Cosmetic; files unaffected | ch:menus
- **Feedback sent by accident** | "Submit this feedback to GitHub Issues? [y/n] (y)" defaults to YES | Type `n`; or press Enter at "Summary:" to cancel | ch:menus, ch:help
- **Cosmetic:** "or press  to return to the main menu" (missing `[m]`), stray `[0m` in undo table, INFO log line during rewind | — | — | ch:menus

## Project folder / file pickup
- **A tool works on an old or wrong structure file** | Each tool picks its input by a fixed priority, not by recency; ProPrep-made files (oriented, `_noH`, assembly, aligned, filtered…) outrank a fresh download | Read the "Using … file:" line; one protein per project; reset workspace or new project to switch | ch:projectfolder
- **4f keeps using an outdated oriented file after re-filtering / re-protonating** | "Oriented" outranks everything in 4f (and in 5a per source) | Orient last; re-run 4f after any change | ch:projectfolder
- **Oriented file names grow (`..._oriented_oriented_...pdb`)** | each 4f run starts from previous oriented file | Run 4f once per change | ch:projectfolder
- **`filtered_structure.pdb` / `structure_with_md_names.pdb` overwritten without asking** | fixed output names | Copy/rename versions you want to keep | ch:projectfolder
- **Residue numbers after 4c start at 1** (2RH1 29–1161 → 1–430 in 4e tables) | repaired file is renumbered | Use the new numbers | ch:projectfolder, ch:cleaning
- **`1UBQ.prmtop` is 0 bytes** | a second, failed tLEaP build wrote to the same name | Give each build its own output name | ch:projectfolder, ch:topology

## Command line
- **`--output-dir a/b` stores a doubled project path (`.../a/b/a/b`)** | relative path resolved after changing folder | Use a full path (`~/...`) | ch:cmdline
- **`--session-dir` crashes with `[Errno 2] No such file or directory`** (relative or non-existent folder), then wrongly prints "Session saved to:" | resolved after moving into project folder; no mkdir | `mkdir` the folder, give a full path | ch:cmdline
- **`--batch-replay --output-dir X` fails every run with "Setup error: [Errno 2]" plus traceback** | folder not created | `mkdir X` first | ch:cmdline, ch:batch
- **`--pdbview` / `--force-replay` still ask the menu-mode question** | menu choice not recorded / not implied | add `--full-menu` | ch:cmdline
- **`--pdbid` + `--pdbfile` together: file silently ignored** | pdbid wins | give only one | ch:cmdline
- **`--pdbid ZZZZ` fails but log says "Downloaded and loaded PDB: ZZZZ"** | misleading log line | trust "Failed to download"; check code | ch:cmdline
- **`proprep-library` → `bad interpreter: .../python3.13: no such file or directory`** | launcher points to Python 3.13; env has 3.12 | `python -m proprep.utils.user_library_snapshots <command>` with the ProPrep env active | ch:cmdline
- **`proprep-web --strict-port` prints "listening on …" and then "address already in use"** | message order | trust the error; choose another port | ch:cmdline

## Help / feedback
- **Feedback intro says "p. Preferences → 2. Configure GitHub"** | it is option 3 | Use p → 3 | ch:help
- **GitHub token stored in plain text in `~/.proprep/settings.json`** although the intro says "stored securely"; token not validated on entry | — | Minimal-scope token, not on shared accounts; Reset preferences removes it | ch:help
- **Feedback default repository `Mag14011/ProPrep` differs from release repository `mjgplab/proprep`** | — | Ask maintainers / email instead | ch:help
- **Feedback module list uses old names ("PDB Loader") and lacks newer tools (Docking, Membrane Builder, Orientator, Bio Assembly, FF Explorer)** | outdated list | choose 0 General and name the tool | ch:help
- **Exit-screen citation differs from CITATION.cff preferred citation** (title "An Interactive and Instructional Interface…" by Walker & Guberman-Pfeffer vs "(Metallo)protein MD Preparation That Keeps Users Informed and in Control" by Walker, Nielson, Cairns, Guberman-Pfeffer) | two references in circulation | Check the project page before publishing; maintainers should align | ch:help
