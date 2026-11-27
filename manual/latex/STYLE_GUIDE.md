# Writer's guide — ProPrep User Manual (LaTeX)

You are writing chapters of a printed user manual for **ProPrep 1.23.0**.
Readers are **scientists who are not programmers** (biochemists, chemists,
students). They are smart but have never used a terminal much and do not know
Amber jargon. The manual must be **accurate, plain, and practical**.

Folder: `<repo>/manual/latex/`
- `proprep-manual.sty` — the house style (do NOT edit)
- `main.tex` — includes every chapter (do NOT edit)
- `chapters/00-howto.tex` — **read this first**: it is the model chapter and
  shows every macro and box in use.
- `chapters/NN-*.tex`, `appendices/X-*.tex` — one file per chapter; edit ONLY
  the files you were assigned.
- `notes/` — write your side notes here (see "Side deliverables").

Source material (the truth):
- `../scenarios/s0…s8-*.md` — what was actually done in ProPrep, step by step,
  with real prompts, files and problems. Each scenario is marked RUN,
  SOURCE-ONLY or FAILED.
- `../runs/_sessions/<name>/transcript.txt` — full real screens (ANSI-stripped).
  Use these to copy exact screen text.
- `../runs/<scope>/...` — the files ProPrep produced.
- `../MANUAL_OUTLINE.md` — the agreed structure.
- ProPrep source (only to double-check an option): `../../src/proprep/`.

## 1. Accuracy rules (most important)
- Describe only what the recordings show. Copy prompts and menu wording
  **exactly** from the scenario files/transcripts.
- RUN scenarios → normal instructions.
- SOURCE-ONLY → you may mention the option briefly: "ProPrep also offers …
  (not tried for this manual)". Never invent its screens.
- FAILED / bug → a `watchout` box: what happens and the workaround, if one was
  found. Say "In version 1.23.0, …" — calm, factual, no blame, no word "bug"
  in body text unless needed.
- Missing software (MODELLER key, Gaussian, ORCA, pmemd, GPU, BLAST locally)
  → a `needsoftware` box saying what the step needs and how far ProPrep gets
  without it.
- Use the example proteins that were really used (1UBQ ubiquitin, 1CRN crambin,
  1HRC cytochrome c, 4HHB hemoglobin, 1HSG HIV protease + MK1, 1BL8 KcsA,
  2RH1, 1ZNI, P69905 AlphaFold, …) and the numbers really observed (times,
  atom counts, RMSDs).

## 2. Voice
- Second person ("you"), present tense, active voice. Short sentences
  (aim ≤ 20 words). One idea per paragraph. No jokes, no hype.
- Explain *why* a step matters in one sentence before *how*, especially for
  scientific choices (pH, water model, box size, charges).
- First use of every technical word: `\term{word}` plus a one-line plain
  explanation in the same sentence or the next. Also list it in your glossary
  note. Prefer everyday words: "file" not "artifact", "error message" not
  "traceback", "folder" not "directory" (but say "directory" once if ProPrep
  uses it on screen).
- Never assume the reader knows Linux, Python, conda, Amber file formats.
- Don't use "simply", "just", "obviously", "easy".
- UK or US spelling: **US** (color, behavior), matching ProPrep's screens.

## 3. Markup (all defined in proprep-manual.sty)
| Markup | Use for |
|---|---|
| `\key{4e}` `\key{y}` `\key{Enter}` | keys/menu codes typed in ProPrep |
| `\typed{1UBQ}` | free text the user types (escape `_` as `\_` inside) |
| `\prompt{Enter your choice}` | ProPrep's own words, inline (short) |
| `\file{filtered_structure.pdb}` | file/folder names (underscores OK) |
| `\menupath{4e}{Protonation State Analyzer}` | a main-menu tool |
| `\pdb{1UBQ}` | PDB/UniProt codes |
| `\term{topology}` | first use of a technical term (also indexes it) |
| `\begin{steps} \item … \end{steps}` | numbered steps (nested level allowed) |
| `\begin{screen} … \end{screen}` | real ProPrep screen text (verbatim) |
| `\begin{terminal} … \end{terminal}` | commands typed in the computer's Terminal |
| `tip`, `watchout`, `needsoftware`, `behindscenes`, `example` | boxes |
| `\begin{glance} \glanceitem{…}{…} \end{glance}` | chapter summary box |
| `tabularx` with `L` columns, `booktabs` rules | tables (`longtable` if long) |

Screens (`screen` environment):
- Verbatim text. **Max 86 characters per line** — trim wide tables (drop or
  shorten columns, cut long descriptions with `...`), rewrap if needed.
- Write `\begin{screen}` with nothing after it on the line. Mark what the user typed as `§ans«1a»` — `§ « »` are the only special
  characters inside a screen. Never use `\` inside a screen.
- Remove ANSI leftovers (`[0m`, `[37m`), spinner lines, repeated progress bars.
- Keep screens short (typically 3–20 lines); cut with `...`. Show the screens a
  user needs to recognise where they are and what to answer.
- Box characters (╭─│═ ✓ ○ ⚠ ● ━) are fine.

Escaping in normal text: `\_ \% \& \# \$ \{ \}`; `~` → `\textasciitilde{}`;
use `\file{}` for paths, `\url{}` for web addresses. Avoid raw `<` `>` in text
mode (use `\textless{}`/`\textgreater{}` or words).

## 4. Chapter skeleton (Part IV task chapters)
```
\chapter{Getting a structure}\label{ch:structure}
<2–4 sentence intro: what this is for, in plain words>
\begin{glance}
\glanceitem{What it does}{…}
\glanceitem{You need it when}{…}
\glanceitem{Before you start}{…}
\glanceitem{Where to find it}{\key{1a} Structure Loader …}
\glanceitem{You get}{\file{…}}
\end{glance}
\section{Step by step: <the recommended route>}   % the path that worked
\section{Other choices}                          % each option, briefly, with when to use it
\section{If something goes wrong}                % symptoms → what to do (refer to Ch. 24/25)
```
Optional: `\section{Understanding the results}` for science-heavy steps, and a
short `\section{Options at a glance}` table listing every choice at a menu.
Chapters in Parts I, II, V, VI and appendices don't need the glance box but
should follow the same voice.

Length guide: Part IV chapter 5–12 pages; others 2–8 pages. Better complete and
clear than long. Don't pad.

## 5. Labels and cross-references (use exactly these)
Chapters: `ch:howto ch:quickstart ch:whatis ch:requirements ch:installation
ch:starting ch:menus ch:projectfolder ch:bigpicture ch:tutorial ch:structure
ch:cleaning ch:protonation ch:redox ch:forcefields ch:topology ch:runningmd
ch:analysis ch:membrane ch:docking ch:qmmm ch:compare ch:sessions ch:batch
ch:cmdline ch:troubleshooting ch:knownissues ch:help`
Appendices: `app:glossary app:menumap app:software app:filetypes app:residues`
Sections: `sec:<chapterkey>-<short>` e.g. `sec:topology-water`.
Refer as `Chapter~\ref{ch:protonation}`, `Section~\ref{sec:…}`,
`Appendix~\ref{app:glossary}`. Refer freely to other chapters by these labels
even if not written yet. Only reference *section* labels you created yourself.

## 6. Menu-mode convention
Instructions use **Full Menu** keys (`\key{4e}`), which never change. Tell the
reader to start with `proprep --full-menu` (or choose 2 at the first question)
when following the manual. Mention Workflow Mode only where it matters
(Chapters 4–5 explain it).
How to start ProPrep (Chapter 4 explains in full): browser `proprep-web`
(recommended) or Terminal `proprep`. The screens are the same in both.

## 7. Side deliverables (write these too, in `notes/`)
- `notes/glossary-<your-id>.md` — every term you `\term{}`ed:
  `- **term** — plain definition, 1–2 sentences.`
- `notes/issues-<your-id>.md` — every user-facing problem you mention:
  `- **Symptom (what the user sees)** | Why | What to do | Chapter label`
  These feed Chapter 24 (Troubleshooting) and 25 (Known limitations).

## 8. Check your work
From `manual/latex`: `bash tools/check.sh chapters/09-getting-structure.tex`
(any number of files). Fix every error. Aim for 0 overfull hboxes (shorten
long `\file{}` paths, break long table cells, shorten screen lines). The
check prints the PDF path; render a page to look at it if you like:
`pdftoppm -r 60 -png -f 1 -l 3 <pdf> /tmp/<yourid>-pg` then view the PNGs.
Undefined-reference warnings for other chapters' labels are expected.

You may run ProPrep to capture a screen that is missing from the recordings
(see `../tools/AGENT_BRIEF.md` for the driver; session names must start with
your writer id, e.g. `w3-…`), but keep it short; the recordings should be
enough. Do not edit ProPrep, other writers' files, the .sty or main.tex.

## 9. Final reply
List the files you wrote, page counts from the check, anything you could not
document accurately, and any macro you wished existed.
