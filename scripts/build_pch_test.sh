#!/usr/bin/env bash
#
# Build and run representative test TUs in build-pch/ (GBS_USE_PCH=ON).
# Run scripts/configure_pch.sh first.
#
# Usage: scripts/build_pch_test.sh [target ...]
set -euo pipefail
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="${REPO_ROOT}/build-pch"

DEFAULT_TARGETS=(tests_basisfunctions tests_bscurve tests_surfaces tests_bscapprox)
if [[ $# -gt 0 ]]; then
    TARGETS=("$@")
else
    TARGETS=("${DEFAULT_TARGETS[@]}")
fi

echo ">> build-pch: ${TARGETS[*]}"
time ninja -C "${BUILD_DIR}" "${TARGETS[@]}"

for t in "${TARGETS[@]}"; do
    bin="${BUILD_DIR}/tests/${t}"
    if [[ -x "${bin}" ]]; then
        echo ">> Running ${t}"
        "${bin}"
    fi
done
