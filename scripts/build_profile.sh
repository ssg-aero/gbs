#!/usr/bin/env bash
#
# Build representative test TUs in build-profile/ for compile-time analysis.
# Produces -ftime-trace JSON files alongside .o files in the build tree.
# Run scripts/configure_profile.sh first.
#
# Usage: scripts/build_profile.sh [target ...]
#   Defaults: tests_basisfunctions tests_bscurve tests_surfaces tests_bscapprox
set -euo pipefail
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="${REPO_ROOT}/build-profile"

DEFAULT_TARGETS=(tests_basisfunctions tests_bscurve tests_surfaces tests_bscapprox)
if [[ $# -gt 0 ]]; then
    TARGETS=("$@")
else
    TARGETS=("${DEFAULT_TARGETS[@]}")
fi

echo ">> build-profile: ${TARGETS[*]}"
time ninja -C "${BUILD_DIR}" "${TARGETS[@]}"
echo ""
echo ">> Trace files written to:"
find "${BUILD_DIR}" -name "*.json" -newer "${BUILD_DIR}/build.ninja" | sort
