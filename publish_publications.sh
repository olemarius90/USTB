#!/bin/bash
#
# publish_publications.sh - Publish MATLAB publication scripts to HTML for the website
#
# Pre-built output is uploaded to the examples-v1 GitHub Release as
# publications-html.tar.gz; deploy-website.yml extracts it into website/examples/
# (same pattern as examples-html.tar.gz).
#
# Prerequisites: MATLAB with USTB on path, network for Zenodo dataset download
#
# Usage:
#   ./publish_publications.sh              # Build publications_html/ and tarball
#   ./publish_publications.sh --upload     # Also upload to GitHub Release examples-v1
#
# Run from the USTB repository root.

set -e
set -o pipefail 2>/dev/null || true

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
OUTPUT_DIR="${SCRIPT_DIR}/publications_html"
TARBALL="publications-html.tar.gz"

# Relative to repo root — keep in sync with publications/ tree
VRALSTAD_REL="publications/TUSON/Vralstad_et_al_2026_Retrospective_transmit_correction_of_blocked_arrays"
SRC_M="${VRALSTAD_REL}/Correction_of_simulated_blockage.m"

MATLAB_BATCH_EXTRA=()
_windows_env=0
if [ -n "${WINDIR:-}" ] || [ -n "${SYSTEMROOT:-}" ]; then
    _windows_env=1
fi
case "${OSTYPE:-}" in
    msys*|cygwin*|mingw*)
        _windows_env=1
        ;;
esac
if [ "$_windows_env" -eq 0 ]; then
    case "$(uname -s 2>/dev/null)" in
        Linux|Darwin)
            MATLAB_BATCH_EXTRA=(-nodisplay)
            ;;
    esac
fi

# Same as publish_examples.sh — Git Bash /c/... is not accepted by Windows MATLAB cd().
repo_path_for_matlab() {
    local p="$1"
    if [ "$_windows_env" -ne 1 ]; then
        printf '%s' "$p"
        return
    fi
    if command -v cygpath >/dev/null 2>&1; then
        if out=$(cygpath -m "$p" 2>/dev/null) && [ -n "$out" ]; then
            printf '%s' "$out"
            return
        fi
    fi
    case "$p" in
        /[a-zA-Z]/?*)
            local drive rest
            drive=$(printf '%s' "${p:1:1}" | tr '[:lower:]' '[:upper:]')
            rest="${p:3}"
            printf '%s:/%s' "$drive" "$rest"
            ;;
        *)
            printf '%s' "$p"
            ;;
    esac
}

SCRIPT_DIR_M=$(repo_path_for_matlab "$SCRIPT_DIR")
OUTPUT_DIR_M=$(repo_path_for_matlab "$OUTPUT_DIR")

echo "=== USTB publication HTML publisher ==="
echo "Source: ${SRC_M}"
echo "Output: ${OUTPUT_DIR}"
echo "(MATLAB paths: ${SCRIPT_DIR_M} -> ${OUTPUT_DIR_M})"
echo ""

rm -rf "${OUTPUT_DIR}"
mkdir -p "${OUTPUT_DIR}/${VRALSTAD_REL}"

if ! command -v matlab &> /dev/null; then
    echo "Error: MATLAB not found on PATH"
    exit 1
fi

echo "Publishing (evalCode)..."
# Single-line -batch — multiline parsing differs across platforms/shells.
matlab "${MATLAB_BATCH_EXTRA[@]}" -batch "addpath(genpath('${SCRIPT_DIR_M}')); src = fullfile('${SCRIPT_DIR_M}', '${SRC_M}'); out = fullfile('${OUTPUT_DIR_M}', '${VRALSTAD_REL}'); if ~isfile(src), error('Missing %s', src); end; opts = struct('outputDir', out, 'format', 'html', 'showCode', true, 'evalCode', true, 'catchError', true, 'createThumbnail', false, 'maxOutputLines', inf); publish(src, opts);" 2>&1 | tee publish_publications.log

echo ""
echo "=== Checking for errors in published HTML ==="
ERROR_COUNT=0
for f in $(find "${OUTPUT_DIR}" -name "*.html"); do
    if grep -q "Error using\|Error in " "$f" 2>/dev/null; then
        echo "  ERROR in $(basename "$f") — remove output or fix script"
        ERROR_COUNT=$((ERROR_COUNT + 1))
    fi
done
echo "Marked ${ERROR_COUNT} HTML file(s) with errors (see grep above)"

cd "${OUTPUT_DIR}" && tar -czf "${SCRIPT_DIR}/${TARBALL}" . && cd "${SCRIPT_DIR}"
echo ""
echo "Tarball: ${TARBALL} ($(du -h "${TARBALL}" | cut -f1))"

if [ "$1" = "--upload" ]; then
    echo ""
    echo "=== Uploading to GitHub Release examples-v1 ==="
    REPO="${2:-olemarius90/USTB}"
    gh release upload examples-v1 "${TARBALL}" --repo "${REPO}" --clobber
    echo "Uploaded ${TARBALL} to ${REPO} release examples-v1"
fi

echo "Done."
