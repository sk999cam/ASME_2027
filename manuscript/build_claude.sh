#!/bin/bash
# Rebuild the draft manuscript manuscript/main_claude.tex into manuscript/main_claude.pdf on this machine.
#
# Usage (from anywhere):
#   bash /home/sk999/Documents/ASME_2027/manuscript/build_claude.sh           # local layout check (stand-in class)
#   bash /home/sk999/Documents/ASME_2027/manuscript/build_claude.sh --real    # try the real asmeconf class
#
# Why two modes: asmeconf v1.47 needs a recent TeX Live (Overleaf has it); the TeX Live 2018 on this machine
# does not compile it. The default mode therefore builds in manuscript/local_build/ with a stand-in class of
# the same page geometry, column separation, font and bibliography style
# (planning/local_tex/asmeconf_standin.cls, copied into local_build/ at each run). Text,
# figures, tables, equations and references are those of main_claude.tex; the page count is an estimate.
# The authoritative PDF is the Overleaf build (set the main document to manuscript/main_claude.tex there).
# --real compiles in manuscript/ with the real class and fails on TeX Live 2018.
#
# Nothing is redrawn: figures are read from ../figures as they are.
set -euo pipefail

MS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TEX=main_claude

if [[ "${1:-}" == "--real" ]]; then
    BUILD_DIR="$MS_DIR"
elif [[ -z "${1:-}" ]]; then
    BUILD_DIR="$MS_DIR/local_build"
    STANDIN="$MS_DIR/../planning/local_tex/asmeconf_standin.cls"
    mkdir -p "$BUILD_DIR"
    if [[ ! -f "$STANDIN" ]]; then
        echo "missing stand-in class $STANDIN" >&2
        exit 1
    fi
    # the stand-in must be named asmeconf.cls where the build runs; local_build/ is git-ignored
    cp "$STANDIN" "$BUILD_DIR/asmeconf.cls"
    # sources are linked, not copied, so the build always uses the current files
    for f in "$TEX.tex" references.bib references_claude.bib asmeconf.bst; do
        ln -sfn "$MS_DIR/$f" "$BUILD_DIR/$f"
    done
    ln -sfn "$MS_DIR/../figures" "$BUILD_DIR/figures"
else
    echo "unknown option: $1 (use --real or nothing)" >&2
    exit 2
fi

cd "$BUILD_DIR"
run_tex() {
    if ! pdflatex -interaction=nonstopmode -halt-on-error "$TEX.tex" > /dev/null; then
        echo "pdflatex failed ($1); errors:" >&2
        grep -A3 '^!' "$TEX.log" >&2 || true
        exit 1
    fi
}
run_tex "pass 1"
bibtex "$TEX" > /dev/null || { echo "bibtex reported errors; see $BUILD_DIR/$TEX.blg" >&2; }
run_tex "pass 2"
run_tex "pass 3"

if [[ "$BUILD_DIR" != "$MS_DIR" ]]; then
    cp "$BUILD_DIR/$TEX.pdf" "$MS_DIR/$TEX.pdf"
fi

undef=$(grep -c "undefined" "$TEX.log" || true)
missing=$(grep -c "not found" "$TEX.log" || true)
[[ "$undef" != "0" ]] && { echo "warning: $undef line(s) mention undefined references or citations:" >&2; grep "undefined" "$TEX.log" | head -10 >&2; }
[[ "$missing" != "0" ]] && { echo "warning: $missing file(s) not found:" >&2; grep "not found" "$TEX.log" >&2; }
pages=$(pdfinfo "$MS_DIR/$TEX.pdf" 2>/dev/null | awk '/^Pages/ {print $2}')
mode=$([[ "$BUILD_DIR" == "$MS_DIR" ]] && echo "real asmeconf class" || echo "stand-in class, page count is an estimate")
echo "built $MS_DIR/$TEX.pdf ($pages pages; $mode)"
