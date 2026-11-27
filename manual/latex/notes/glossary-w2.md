# Glossary terms introduced by w2 (Ch. 4, 5, 6, 23, 26, App. B)

- **Terminal** — The program where you type commands to your computer (Terminal on macOS, the Ubuntu window under WSL on Windows). ProPrep is started from it.
- **Full Menu** — ProPrep's main menu that shows every tool at once, each with a fixed code such as 1a or 4e. This manual uses Full Menu codes. Start it with `proprep --full-menu`.
- **Workflow Mode** — ProPrep's guided main menu, which shows only the tools of one stage (Load, Detect, Compare, Fix, Simulate) at a time; tools are numbered 1, 2, 3 … within each stage.
- **option (command-line option, flag)** — An extra instruction typed after a command, starting with two hyphens, such as `--full-menu`. It applies to that run only.
- **session file** — The file `proprep_session_<date>_<time>.json` in the project folder in which ProPrep records every answer you type, in order. It is used to undo, resume, replay and make templates.
- **breadcrumb** — The `Navigation:` line at the top of each ProPrep screen, showing where you are in the menus (e.g. `Main Menu > Structure Loader`).
- **default** — The answer shown in parentheses at the end of a question, e.g. `(y)`. Pressing Enter without typing anything accepts it.
- **workspace** — What ProPrep holds in its memory during a run: the structures loaded, the files it has made and the results of each tool. It starts empty every time ProPrep starts; files on disk stay.
- **project folder** — The folder in which ProPrep saves every file for one project; named at the "Project name" question.
- **PDB file** — A text file (`.pdb`) listing every atom of a structure with its position; the format used by the Protein Data Bank and by most ProPrep steps.
- **topology** — The file (`.prmtop`) describing every atom, bond and force-field parameter of the simulation system; made by the Topology Generator (5a) with Amber's tLEaP.
- **snapshot (parameter library)** — A saved copy of your personal force-field parameter library, made with the `proprep-library` command, which can be restored later.
- **traceback** — The full technical error report a program prints when it stops with an error: the list of program steps that led to it. ProPrep shows it with `--debug`.
- **issue (GitHub issue)** — A report (bug, question, suggestion) posted on a software project's public GitHub page. ProPrep's Feedback & Support form creates one.
- **personal access token** — A long code, like a password, made on the GitHub website, that lets another program (here ProPrep's feedback form) post on your behalf.
- **port** — The number after the colon in a local web address such as `http://127.0.0.1:8000/`; it tells the browser which program on your computer to talk to. `proprep-web` uses 8000 or the next free one.
- **127.0.0.1 / localhost** — The address meaning "this computer"; a page at this address is served by a program on your own machine, not from the internet.
