#!/bin/bash
#
# publish_examples.sh - Generate published HTML examples for the USTB website
#
# This script runs MATLAB's publish() on all eligible examples, generates
# an index page, packages everything into a tarball, and optionally uploads
# it to a GitHub Release.
#
# Prerequisites:
#   - MATLAB R2024b or later with Signal_Processing_Toolbox
#   - Field II (optional, for field_II examples): install to /opt/field_ii
#   - Python 3 (for generate_examples_index.py)
#   - Recompiled MEX file (see +mex/build_mex.m or fix/recompile-mex branch)
#
# Usage:
#   ./publish_examples.sh              # Generate examples
#   ./publish_examples.sh --upload     # Generate and upload to GitHub Release
#
# The script must be run from the USTB repository root.

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
OUTPUT_DIR="${SCRIPT_DIR}/examples_html"
TARBALL="examples-html.tar.gz"
FIELD_II_PATH="/opt/field_ii"

echo "=== USTB Example Publisher ==="
echo "Output: ${OUTPUT_DIR}"
echo ""

# Clean previous output
rm -rf "${OUTPUT_DIR}"
mkdir -p "${OUTPUT_DIR}"

# Check MATLAB
if ! command -v matlab &> /dev/null; then
    echo "Error: MATLAB not found on PATH"
    exit 1
fi
echo "MATLAB: $(matlab -batch "disp(version)" 2>/dev/null | tail -1)"

# -nodisplay is Linux/macOS only; Windows MATLAB warns "Unrecognized command line option: nodisplay".
# Prefer env over uname: some shells report oddly, and detecting Windows avoids passing the flag even
# when MATLAB is the Windows binary (e.g. some WSL/Git Bash setups).
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

# Check MEX (platform-specific binary)
echo -n "MEX: "
if ls +mex/das_c.mexw64 >/dev/null 2>&1; then
    ls -la +mex/das_c.mexw64 2>/dev/null | awk '{print $6, $7, $8, $9}'
elif ls +mex/das_c.mexa64 >/dev/null 2>&1; then
    ls -la +mex/das_c.mexa64 2>/dev/null | awk '{print $6, $7, $8, $9}'
else
    echo "(no das_c mex found)"
fi
FIELD_II_CMD=""
if [ -d "${FIELD_II_PATH}" ]; then
    echo "Field II: ${FIELD_II_PATH}"
    FIELD_II_CMD="addpath('${FIELD_II_PATH}');"
else
    echo "Field II: not found (field_II examples will use fresnel fallback)"
fi

echo ""

# Publish examples. Single-line -batch: multiline breaks argument passing.
# Use MATLAB_BATCH_EXTRA (-nodisplay on Unix only; Windows MATLAB rejects it).
echo "Publishing examples..."
unset DISPLAY
matlab "${MATLAB_BATCH_EXTRA[@]}" -batch "cd('${SCRIPT_DIR}'); addpath(genpath(pwd)); ${FIELD_II_CMD}publish_all_examples('${OUTPUT_DIR}', true);" 2>&1 | tee publish_examples.log

# Check for errors in published HTML
echo ""
echo "=== Checking for errors in published HTML ==="
ERROR_COUNT=0
for f in $(find "${OUTPUT_DIR}" -name "*.html"); do
    if grep -q "Error using\|Error in " "$f" 2>/dev/null; then
        NAME=$(echo "$f" | sed "s|${OUTPUT_DIR}/||")
        echo "  ERROR: ${NAME}"
        # Remove error files
        DIR=$(dirname "$f")
        BASE=$(basename "$f" .html)
        rm -f "$f" "${DIR}/${BASE}"_*.png
        ERROR_COUNT=$((ERROR_COUNT + 1))
    fi
done
echo "Removed ${ERROR_COUNT} examples with errors"

# Generate index
echo ""
echo "=== Generating index ==="
python3 generate_examples_index.py "${OUTPUT_DIR}"

# Summary
HTML_COUNT=$(find "${OUTPUT_DIR}" -name "*.html" -not -name "index.html" | wc -l)
PNG_COUNT=$(find "${OUTPUT_DIR}" -name "*.png" | wc -l)
echo ""
echo "=== Summary ==="
echo "Examples: ${HTML_COUNT}"
echo "Figures:  ${PNG_COUNT}"
echo "Output:   ${OUTPUT_DIR}"

# Package
echo ""
echo "=== Packaging ==="
cd "${OUTPUT_DIR}" && tar -czf "${SCRIPT_DIR}/${TARBALL}" . && cd "${SCRIPT_DIR}"
echo "Tarball: ${TARBALL} ($(du -h ${TARBALL} | cut -f1))"

# Upload if requested
if [ "$1" = "--upload" ]; then
    echo ""
    echo "=== Uploading to GitHub Release ==="
    REPO="${2:-olemarius90/USTB}"
    gh release upload examples-v1 "${TARBALL}" --repo "${REPO}" --clobber
    echo "Uploaded to ${REPO} release examples-v1"
fi

echo ""
echo "Done!"
