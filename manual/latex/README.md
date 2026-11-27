# ProPrep User Manual — LaTeX source

Build (from this folder):

```
latexmk            # → main.pdf  (XeLaTeX, runs makeindex, repeats until references settle)
latexmk -c         # tidy auxiliary files
```

Layout
- `main.tex` — book skeleton (parts, chapter order)
- `proprep-manual.sty` — house style: fonts, colours, boxes, `screen`/`terminal` environments, macros
- `chapters/` — one file per chapter (`00-*` front matter, `01`–`26`)
- `appendices/` — A glossary, B menu map, C software, D file types, E residue names
- `STYLE_GUIDE.md` — conventions for anyone editing the text
- `notes/` — writers' raw notes (glossary terms, user-facing issues) used for Ch. 24/25 and Appendix A
- `tools/check.sh` — compile one or more chapter files on their own

Requirements
- XeLaTeX with tcolorbox, tikz, fancyvrb, titlesec, fancyhdr, hyperref, url, makeidx, TeX Gyre fonts.
- The screen font is Menlo (macOS). On other systems it falls back to Latin Modern Mono for
  normal mono text, but `\screenfont` needs a font with box-drawing characters: change
  `\newfontfamily\screenfont{Menlo}` in the .sty to e.g. `DejaVu Sans Mono`.
- If XeLaTeX stops with "Mismatched LaTeX support files detected", the LaTeX files in your
  personal texmf folder are newer than your TeX installation's formats. Rebuild your own
  formats once with `fmtutil-user --byfmt xelatex` (and `pdflatex`, `lualatex` if you use them).

Writing rules in short (full version in STYLE_GUIDE.md)
- `\begin{screen}` … `\end{screen}`: real ProPrep output, ≤ 86 characters per line (≤ 80 inside a
  numbered step), nothing after `\begin{screen}` on its line; typed answers as `§ans«y»`.
- `\key{4e}`, `\typed{text}` (escape `_`), `\prompt{ProPrep's words}`, `\file{path}` (breaks
  automatically; don't escape `_`), `\term{word}` (indexed; define it in Appendix A).
- Boxes: `tip`, `watchout`, `needsoftware`, `behindscenes`, `example`, `glance`.
- Source of truth for behaviour: `../scenarios/*.md` and `../runs/_sessions/*/transcript.txt`.
