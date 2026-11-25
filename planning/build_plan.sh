#!/bin/bash
# Rebuild planning/writeup_plan.pdf from writeup_plan.tex.
#
# Usage (from anywhere):
#   bash /home/sk999/Documents/ASME_2027/planning/build_plan.sh          # recompile the PDF only
#   bash /home/sk999/Documents/ASME_2027/planning/build_plan.sh --figs   # first redraw every figure, then recompile
#
# Without --figs the figure files are used as they are, so figures edited by hand are kept.
# With --figs the figures are redrawn from their scripts, which overwrites hand edits:
#   ASME_2027/figures/code/*.py      -> ASME_2027/figures/<same name>.pdf
set -euo pipefail

PLAN_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FIG_DIR="$PLAN_DIR/../figures"
TEX=writeup_plan

if [[ "${1:-}" == "--figs" ]]; then
    for script in "$FIG_DIR"/code/*.py; do
        echo "drawing $(basename "$script" .py)"
        (cd "$FIG_DIR/code" && python3 "$script")
    done
elif [[ -n "${1:-}" ]]; then
    echo "unknown option: $1 (use --figs or nothing)" >&2
    exit 2
fi

cd "$PLAN_DIR"
# two passes for cross-references; stop on the first LaTeX error
for pass in 1 2; do
    if ! pdflatex -interaction=nonstopmode -halt-on-error "$TEX.tex" > /dev/null; then
        echo "pdflatex failed on pass $pass; errors:" >&2
        grep -A3 '^!' "$TEX.log" >&2 || true
        exit 1
    fi
done

missing=$(grep -c "not found" "$TEX.log" || true)
if [[ "$missing" != "0" ]]; then
    echo "warning: $missing file(s) not found:" >&2
    grep "not found" "$TEX.log" >&2
fi
echo "built $PLAN_DIR/$TEX.pdf ($(pdfinfo "$TEX.pdf" 2>/dev/null | awk '/^Pages/ {print $2}') pages)"
