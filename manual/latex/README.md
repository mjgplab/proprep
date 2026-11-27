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
- This Mac's personal texmf tree has newer LaTeX files than its TeX Live 2025 base, which
  breaks the stock `xelatex` format. `latexmkrc` therefore uses a private format in
  `.texfmt/ppxelatex.fmt`. If it goes missing or stale, rebuild it with
  `(cd .texfmt && xelatex -ini -etex -jobname=ppxelatex xelatex.ini)`.
  The private format is used only when that file exists, so on a normal TeX installation
  (or Overleaf with XeLaTeX) the standard format is used automatically.

Writing rules in short (full version in STYLE_GUIDE.md)
- `\begin{screen}` … `\end{screen}`: real ProPrep output, ≤ 86 characters per line (≤ 80 inside a
  numbered step), nothing after `\begin{screen}` on its line; typed answers as `§ans«y»`.
- `\key{4e}`, `\typed{text}` (escape `_`), `\prompt{ProPrep's words}`, `\file{path}` (breaks
  automatically; don't escape `_`), `\term{word}` (indexed; define it in Appendix A).
- Boxes: `tip`, `watchout`, `needsoftware`, `behindscenes`, `example`, `glance`.
- Source of truth for behaviour: `../scenarios/*.md` and `../runs/_sessions/*/transcript.txt`.
