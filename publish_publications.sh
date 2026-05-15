#!/usr/bin/env bash
#
# publish_publications.sh — MATLAB **publications-html.tar.gz** (TUSON / publish-m only)
#
# Example HTML gallery: ./publish_examples.sh
# Dataset previews: ./publish_datasets.sh
#
# Usage:
#   ./publish_publications.sh
#   ./publish_publications.sh --upload [OWNER/REPO]

set -e
set -o pipefail 2>/dev/null || true

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=publish_common.sh
. "${SCRIPT_DIR}/publish_common.sh"

publish_common_matlab_extra

OUTPUT_DIR="${SCRIPT_DIR}/publications_html"
TARBALL="publications-html.tar.gz"
VRALSTAD_REL="publications/TUSON/Vralstad_et_al_2026_Retrospective_transmit_correction_of_blocked_arrays"
SRC_M="${VRALSTAD_REL}/Correction_of_simulated_blockage.m"

SCRIPT_DIR_M=$(publish_repo_path_for_matlab "$SCRIPT_DIR")
OUTPUT_DIR_M=$(publish_repo_path_for_matlab "$OUTPUT_DIR")

echo "=== USTB publication HTML publisher ==="
echo "Source: ${SRC_M}"
echo "Output: ${OUTPUT_DIR}"
echo "(MATLAB paths: ${SCRIPT_DIR_M} -> ${OUTPUT_DIR_M})"
echo ""

rm -rf "${OUTPUT_DIR}"
mkdir -p "${OUTPUT_DIR}/${VRALSTAD_REL}"

publish_require_matlab

echo "Publishing (evalCode)..."
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
if [ "${ERROR_COUNT}" -gt 0 ]; then
    echo "Aborting: fix the script or published output before packaging (no tarball)." >&2
    exit 1
fi

cd "${OUTPUT_DIR}" && tar -czf "${SCRIPT_DIR}/${TARBALL}" . && cd "${SCRIPT_DIR}"
echo ""
echo "Tarball: ${TARBALL} ($(du -h "${SCRIPT_DIR}/${TARBALL}" | cut -f1))"

if [ "${1:-}" = "--upload" ]; then
    echo ""
    echo "=== Uploading to GitHub Release examples-v1 ==="
    REPO="${2:-olemarius90/USTB}"
    gh release upload examples-v1 "${SCRIPT_DIR}/${TARBALL}" --repo "${REPO}" --clobber
    echo "Uploaded ${TARBALL} to ${REPO} release examples-v1"
fi

echo "Done."
