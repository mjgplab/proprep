# Brief for agents recording ProPrep scenarios

## Goal
We are writing a **user manual for ProPrep** (protein preparation for Amber MD,
v1.23.0) aimed at **non-technical scientists**. Your job is NOT to write the
manual. Your job is to **actually drive the application** through every path in
your assigned scope and **record what happens**, so a writer can later turn your
notes into plain-English instructions. Record real screens, real prompts, real
files produced, and real errors — never invent behaviour you did not see.
Reading source code is allowed to *find* options and understand them, but every
scenario you report must say whether it was **RUN** (you drove it) or
**SOURCE-ONLY** (you could not run it, and why).

## Environment facts (already verified)
- ProPrep env: `/opt/homebrew/Caskroom/miniforge/base/envs/ProPrep` (Python 3.12,
  ProPrep 1.23.0, AmberTools: tleap, sander, antechamber, parmchk2, MCPB.py,
  cpptraj, pdb4amber, reduce, propka3, pdb2pqr30, vina, packmol, packmol-memgen).
- **Not available:** MODELLER license key (MODELLER is installed but import
  fails with "Invalid license key" -> loop/residue repair cannot run), Gaussian
  (g16), ORCA, BLAST+ (blastp), pmemd/pmemd.cuda, xtb, crest, openbabel.
  When a feature needs these, run it as far as it goes, record the message the
  user sees, and mark what is needed.
- Internet works (RCSB, AlphaFold DB downloads are fine).
- Source code: `<repo>/src/proprep/`
  (menus: `application/menu_commands.py`; each tool lives in its own package).

## How to drive the app — `manual/tools/drive.py`
ProPrep is interactive (typed answers at prompts). The driver keeps it running
in a pseudo-terminal between your shell calls. Always run from
`<repo>/manual`:

```
D="python3 tools/drive.py"
$D start <session> --cwd runs/<your-scope>/<something> -- /opt/homebrew/Caskroom/miniforge/base/envs/ProPrep/bin/proprep --full-menu
$D send <session> "1a"            # types 1a + Enter, prints new output
$D send <session> "" --wait 3     # just press Enter (accept default)
$D send <session> "y" --wait 10 --timeout 600   # slow step (downloads, tleap, sander)
$D read <session> --wait 5 --timeout 600        # keep waiting for output
$D screen <session> --lines 80    # re-show the tail of the transcript
$D stop <session>
```
- Output printed is everything since your last command; `--wait` = seconds of
  silence that count as "done". Raise `--wait` for slow steps if a screen looks cut off.
- `<UP> <DOWN> <ENTER> <CTRL-C> <ESC>` are understood inside the text; `--no-enter` to not add Enter.
- **Prefix every session name with your scope id** (e.g. `s3-fixer-1`) so you
  never collide with other agents running at the same time. Use a cwd under
  `runs/<scope-id>/`. Never touch other agents' sessions or folders.
- Startup asks for a **project name**: give a short name; outputs go in that folder.
- Prompts show the default in parentheses, e.g. `(1a)`; Enter accepts it.
  Typing `undo` at a prompt rewinds. Typing something invalid shows the valid choices.
- Do NOT answer "y" to "Save this as your default preference?" — keep `n`, and
  do not change global Preferences unless that is your scope (other agents share `~/.proprep`).
- Do not submit feedback/GitHub issues, do not post anything online, do not
  install packages, and do not edit ProPrep source code.
- Every session writes a replayable `proprep_session_*.json` in the project
  folder — note its path for scenarios you complete; they are useful evidence.
- Good test structures: `1HRC`/`1AKK` (cytochrome c, heme), `1UBQ` (ubiquitin,
  small, clean), `4HHB` (hemoglobin, multi-chain, hemes), `1CRN` (crambin, tiny),
  `2RH1` (GPCR, membrane), `1HSG` (HIV protease + ligand MK1, docking),
  `3NIR`/`1ZNI` (zinc/metal sites), AlphaFold `P69905`. Pick small ones for slow steps.
- Keep shell calls under ~10 minutes each; for long runs use `read` repeatedly.

## What to deliver
Write ONE Markdown file: `manual/scenarios/<scope-id>-<short-name>.md`.
Structure it like this, one block per feature, and one sub-block per scenario
(every distinct path/option you exercised: success paths, alternative choices,
"what if I pick the other option", invalid input, missing prerequisite, errors):

```
# <Scope title>

## Feature: <Menu key + name>  (e.g. "1a. Structure Loader")
**What it is for (plain English):** one or two sentences a biologist understands.
**Needs before you start:** prerequisites (other steps, files, software, licence).
**Where to find it:** menu key(s), workflow-mode location, command-line flag if any.

### Scenario <n>: <short title>   — RUN | SOURCE-ONLY | FAILED
- Starting point: ...
- Steps (exactly what was typed and the prompt it answered):
  1. Prompt: "Enter your choice (1a):" -> typed `1a`
  2. Prompt: "..." -> typed `...`  (explain what the options meant)
- What the user sees: short quotes of key screens/messages (keep quotes short; summarise long tables).
- Result / files produced: paths + one-line description of each file.
- Time taken (rough), if notable.
- Problems, warnings, confusing bits, and how to get past them.
- Session log / transcript: runs/_sessions/<session>/transcript.txt

### Options reference
Table of every option/choice at each prompt of this feature and what it does
(from the screens, plus source code when a branch could not be run).

## Notes for the manual writer
Jargon to explain, common mistakes, recommended default path for a beginner,
things that surprised you, anything that looked like a bug (describe, don't fix).
```

Finish with a short final reply: file path written, count of scenarios RUN /
SOURCE-ONLY / FAILED, and the 3–5 most important findings (blockers, bugs,
missing software). Be thorough — "every possible scenario" is the user's ask —
but stay inside your scope.
