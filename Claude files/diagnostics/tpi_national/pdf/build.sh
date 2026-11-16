#!/bin/sh
# RENDER Claude files/docs/tpi_fix_explained.md TO PDF (pandoc 2.12 + LuaLaTeX + STIX fonts, as on this Mac).
#   sh "Claude files/diagnostics/tpi_national/pdf/build.sh"      (from the repo root)
set -e
K="Claude files/diagnostics/tpi_national/pdf"
B=$(mktemp -d)
python3 "$K/md_for_pandoc.py" "Claude files/docs/tpi_fix_explained.md" "$B/report.md"
cp "$K/header.tex" "$B/"
cd "$B"
pandoc -f markdown-superscript-subscript report.md -s -o report.tex -H header.tex -V mathfont="STIX Two Math" \
  -V geometry:margin=2.1cm -V fontsize=10pt -V colorlinks=true -V linkcolor=blue -V urlcolor=blue \
  -V toccolor=black --toc --toc-depth=2 --shift-heading-level-by=-1
lualatex -interaction=nonstopmode -halt-on-error report.tex > /dev/null
lualatex -interaction=nonstopmode -halt-on-error report.tex > /dev/null
cd - > /dev/null
cp "$B/report.pdf" "Claude files/docs/tpi_fix_explained.pdf"
echo "wrote Claude files/docs/tpi_fix_explained.pdf"
