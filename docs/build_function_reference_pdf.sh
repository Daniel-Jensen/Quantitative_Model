#!/usr/bin/env bash
# Build docs/function_reference.pdf from docs/function_reference.md.
#
# Requires: pandoc + xelatex (both checked below).
# The markdown's own H1 title and manual "Table of contents" section are
# stripped before conversion — pandoc generates the title block and a
# hyperlinked, page-numbered TOC instead.
set -euo pipefail
cd "$(dirname "$0")"

command -v pandoc >/dev/null || { echo "pandoc not found" >&2; exit 1; }
command -v xelatex >/dev/null || { echo "xelatex not found" >&2; exit 1; }

tmp="$(mktemp -t funcref_XXXX).md"
trap 'rm -f "$tmp"' EXIT

awk 'NR==1 && /^# /              {next}   # drop the H1 (title via metadata)
     /^## Table of contents/     {skip=1} # drop the manual TOC section …
     skip && /^---$/             {skip=0; next}  # … up to its closing rule
     !skip {gsub(/≪/, "<<"); print}       # ≪ has no glyph in Menlo
    ' function_reference.md > "$tmp"

pandoc "$tmp" -f gfm -o function_reference.pdf \
  --pdf-engine=xelatex \
  --shift-heading-level-by=-1 \
  --toc --toc-depth=3 \
  --metadata title="Runtime Function Reference" \
  --metadata subtitle="Two-country HANK–GK monetary union — code/global/" \
  --metadata author="Model documentation" \
  --metadata date="$(date +%Y-%m-%d)" \
  -V geometry:margin=2.2cm \
  -V fontsize=10pt \
  -V mainfont="STIXGeneral" \
  -V monofont="Menlo" \
  -V monofontoptions="Scale=0.78" \
  -V colorlinks=true -V linkcolor=blue -V urlcolor=blue -V toccolor=black

echo "Wrote $(pwd)/function_reference.pdf"
