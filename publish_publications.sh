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

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
OUTPUT_DIR="${SCRIPT_DIR}/publications_html"
TARBALL="publications-html.tar.gz"

# Relative to repo root — keep in sync with publications/ tree
VRALSTAD_REL="publications/TUSON/Vralstad_et_al_2026_Retrospective_transmit_correction_of_blocked_arrays"
SRC_M="${VRALSTAD_REL}/Correction_of_simulated_blockage.m"

echo "=== USTB publication HTML publisher ==="
echo "Source: ${SRC_M}"
echo "Output: ${OUTPUT_DIR}"
echo ""

rm -rf "${OUTPUT_DIR}"
mkdir -p "${OUTPUT_DIR}/${VRALSTAD_REL}"

if ! command -v matlab &> /dev/null; then
    echo "Error: MATLAB not found on PATH"
    exit 1
fi

echo "Publishing (evalCode)..."
matlab -nodisplay -batch "
    addpath(genpath('${SCRIPT_DIR}'));
    src = fullfile('${SCRIPT_DIR}', '${SRC_M}');
    out = fullfile('${OUTPUT_DIR}', '${VRALSTAD_REL}');
    if ~isfile(src), error('Missing %s', src); end
    opts = struct('outputDir', out, 'format', 'html', 'showCode', true, ...
        'evalCode', true, 'catchError', true, 'createThumbnail', false, ...
        'maxOutputLines', inf);
    publish(src, opts);
" 2>&1 | tee publish_publications.log

# Drop broken publishes (same idea as publish_examples.sh)
echo ""
echo "=== Checking for errors in published HTML ==="
for f in $(find "${OUTPUT_DIR}" -name "*.html"); do
    if grep -q "Error using\|Error in " "$f" 2>/dev/null; then
        echo "  ERROR in $(basename "$f") — remove output or fix script"
    fi
done

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
