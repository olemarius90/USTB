#!/usr/bin/env bash
#
# publish_datasets.sh — Beamform PNG previews + datasets.html for the website
#
# Produces datasets-html.tar.gz (contents rooted at repo website/ ):
#   datasets.html
#   assets/images/datasets/*.png
#
# GitHub Actions: deploy-website.yml extracts this tarball with -C website
#
# Prerequisites: MATLAB + USTB, Python 3, network (Zenodo / ustb dataset hosts)
#
# Usage:
#   ./publish_datasets.sh
#   ./publish_datasets.sh --upload
#   ./publish_datasets.sh --upload unioslo/USTB

set -e
set -o pipefail 2>/dev/null || true

_THIS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=publish_common.sh
. "${_THIS_DIR}/publish_common.sh"

publish_common_matlab_extra
SCRIPT_DIR="${_THIS_DIR}"
TARBALL="datasets-html.tar.gz"
WEBSITE="${SCRIPT_DIR}/website"

publish_require_matlab
echo "MATLAB: $(matlab -batch "disp(version)" 2>/dev/null | tail -1)"

SCRIPT_DIR_M=$(publish_repo_path_for_matlab "$SCRIPT_DIR")
echo ""
echo "=== USTB dataset previews + datasets page ==="
echo "Website dir: ${WEBSITE}"
echo "(MATLAB repo path: ${SCRIPT_DIR_M})"
echo ""

mkdir -p "${WEBSITE}/assets/images/datasets"

echo "Exporting previews (MATLAB, may download many datasets)..."
unset DISPLAY
matlab "${MATLAB_BATCH_EXTRA[@]}" -batch "cd('${SCRIPT_DIR_M}'); addpath(genpath(pwd)); addpath(fullfile('${SCRIPT_DIR_M}','examples','dataset_smoke_tests')); export_dataset_previews_to_website();" \
    2>&1 | tee publish_datasets_matlab.log

echo ""
echo "Building datasets.html..."
python3 "${SCRIPT_DIR}/website/scripts/build_datasets_page.py"

if [ ! -f "${WEBSITE}/datasets.html" ]; then
    echo "Error: ${WEBSITE}/datasets.html missing after build" >&2
    exit 1
fi

echo ""
echo "=== Packaging ${TARBALL} ==="
(
    cd "${WEBSITE}"
    tar -czf "${SCRIPT_DIR}/${TARBALL}" datasets.html assets/images/datasets
)
echo "Tarball: ${TARBALL} ($(du -h "${SCRIPT_DIR}/${TARBALL}" | cut -f1))"

if [ "${1:-}" = "--upload" ]; then
    echo ""
    echo "=== Uploading to GitHub Release examples-v1 ==="
    REPO="${2:-olemarius90/USTB}"
    gh release upload examples-v1 "${SCRIPT_DIR}/${TARBALL}" --repo "${REPO}" --clobber
    echo "Uploaded ${TARBALL} to ${REPO} release examples-v1"
fi

echo "Done."
