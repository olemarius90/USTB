#!/usr/bin/env bash
#
# publish_publications.sh — MATLAB **publications-html.tar.gz**
#
# Publishes every script iframe-linked from website/publications.html via
# publish_all_publications.m → publications/<venue>/<slug>/*.html (+ PNGs).
# Example gallery: ./publish_examples.sh · Dataset previews: ./publish_datasets.sh
#
# Usage:
#   ./publish_publications.sh
#   ./publish_publications.sh --upload [OWNER/REPO]

set -e
set -o pipefail 2>/dev/null || true

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
# shellcheck source=publish_common.sh
. "${SCRIPT_DIR}/publish_common.sh"

publish_common_matlab_extra

OUTPUT_DIR="${REPO_ROOT}/publications_html"
TARBALL="publications-html.tar.gz"

REPO_ROOT_M=$(publish_repo_path_for_matlab "$REPO_ROOT")
OUTPUT_DIR_M=$(publish_repo_path_for_matlab "$OUTPUT_DIR")

echo "=== USTB publication HTML publisher ==="
echo "Scripts: publish_all_publications.m (website/publications.html)"
echo "Output: ${OUTPUT_DIR}"
echo "(MATLAB paths: repo ${REPO_ROOT_M} → ${OUTPUT_DIR_M})"
echo ""

rm -rf "${OUTPUT_DIR}"
mkdir -p "${OUTPUT_DIR}"

publish_require_matlab

echo "Publishing (evalCode, all publication iframes)..."
matlab "${MATLAB_BATCH_EXTRA[@]}" -batch "cd('${REPO_ROOT_M}'); addpath(genpath(pwd)); publish_all_publications('${OUTPUT_DIR_M}');" 2>&1 | tee publish_publications.log

echo ""
echo "=== Required HTML outputs ==="
MISSING=0
for rel in \
    publications/preprint/generalized_beamformer/CPWC_double_adaptive_redone.html \
    publications/TUSON/Vralstad_et_al_2026_Retrospective_transmit_correction_of_blocked_arrays/Correction_of_simulated_blockage.html \
    publications/TUFFC/Dynamic_range_2020/dynamic_range_test.html \
    publications/TUFFC/Prieur_fDMAS_2018/FI_UFF_FIeldII_simulations_Fig2_and_Fig3.html \
    publications/IUS/2018_virtual_source_model/Proceedings_FI_UFF_Verasonics_RTB_delay_models.html \
    publications/IUS/2017_dark_region_artifact/process_beamformed_experimental_data.html
do
    if [ ! -f "${OUTPUT_DIR}/${rel}" ]; then
        echo "  MISSING ${rel}"
        MISSING=$((MISSING + 1))
    fi
done
if [ "${MISSING}" -gt 0 ]; then
    echo "Aborting: ${MISSING} expected publication HTML file(s) missing." >&2
    exit 1
fi

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
