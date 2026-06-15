#!/usr/bin/env bash
#
# Build and run gbs test targets.
#
# Usage:
#   scripts/build_test.sh                      # build + run the default eval/interp test set
#   scripts/build_test.sh tests_basisfunctions # build + run a single target
#   scripts/build_test.sh t1 t2 t3             # build + run several targets
#
# Honours the existing (already-configured) build/ directory and the ninja
# generator. Builds ONLY the requested test targets (no full re-link of the
# slow Python bindings). Set DOCTEST_FILTER to pass a doctest -tc=<filter>.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="${REPO_ROOT}/build"

# Default target set: the eval/interp coverage relevant to find_span (#42/#46).
DEFAULT_TARGETS=(tests_basisfunctions tests_bscurve tests_surface_derivatives tests_surfaces)

if [[ $# -gt 0 ]]; then
    TARGETS=("$@")
else
    TARGETS=("${DEFAULT_TARGETS[@]}")
fi

echo ">> Building targets: ${TARGETS[*]}"
ninja -C "${BUILD_DIR}" "${TARGETS[@]}"

for t in "${TARGETS[@]}"; do
    bin="${BUILD_DIR}/tests/${t}"
    if [[ ! -x "${bin}" ]]; then
        echo "!! No runnable binary at ${bin}, skipping run" >&2
        continue
    fi
    echo ">> Running ${t}"
    if [[ -n "${DOCTEST_FILTER:-}" ]]; then
        "${bin}" -tc="${DOCTEST_FILTER}"
    else
        "${bin}"
    fi
done
