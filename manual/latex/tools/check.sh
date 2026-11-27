#!/bin/bash
# Compile one or more chapter files on their own, to check they build.
# Usage (from manual/latex):  bash tools/check.sh chapters/09-getting-structure.tex [more files]
# Prints errors; leaves the PDF in /tmp/ppcheck-XXXX/check.pdf for a look.
set -u
HERE="$(cd "$(dirname "$0")/.." && pwd)"
W="$(mktemp -d /tmp/ppcheck-XXXXXX)"
{
  echo '\documentclass[11pt,openany]{book}'
  echo "\\usepackage{$HERE/proprep-manual}"
  echo '\begin{document}'
  for f in "$@"; do echo "\\input{$HERE/${f%.tex}}"; done
  echo '\end{document}'
} > "$W/check.tex"
cd "$W"
TEXINPUTS="$HERE:$HERE/figures:" xelatex -interaction=nonstopmode -file-line-error check.tex >/dev/null 2>&1
TEXINPUTS="$HERE:$HERE/figures:" xelatex -interaction=nonstopmode -file-line-error check.tex > build.out 2>&1
grep -nE "^(/|\./).*:[0-9]+:|^! |Undefined control|Missing character" check.log | head -40
grep -c "Overfull \\\\hbox" check.log | sed 's/^/overfull hboxes: /'
grep -A1 "Overfull \\\\hbox" check.log | grep -oE "lines? [0-9]+(--[0-9]+)?" | head -15 | tr '\n' ' '; echo
echo "pages: $(grep -oE 'Output written on .* \(([0-9]+) page' check.log | grep -oE '[0-9]+ page')  pdf: $W/check.pdf"
