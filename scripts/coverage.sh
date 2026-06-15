#!/usr/bin/env bash
#
# Reproduce the gbs test-coverage baseline (issue #50).
#
# Configures a DEDICATED build directory (build-cov/, never the main build/)
# with clang + source-based coverage instrumentation, builds the priority
# test targets, runs them and aggregates an llvm-cov report.
#
# Usage:
#   scripts/coverage.sh             # configure (if needed) + build + report
#   scripts/coverage.sh --report    # skip configure/build, just re-run + report
#
# Knobs (env):
#   GBS_COV_BUILD_DIR  build dir (default: build-cov)
#   CC / CXX           compiler (default: clang-21 / clang++-21)
#
# The matching `coverage` CMake target does the report step alone:
#   cmake --build build-cov --target coverage
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="${GBS_COV_BUILD_DIR:-${REPO_ROOT}/build-cov}"
CC_BIN="${CC:-clang-21}"
CXX_BIN="${CXX:-clang++-21}"

# Priority-perimeter test executables (eval/span/basis, interp, loft).
TARGETS=(
    tests_basisfunctions
    tests_bscurve
    tests_bssurf
    tests_curve_derivatives
    tests_surface_derivatives
    tests_curves
    tests_surfaces
    tests_bsinterp
    tests_knotsfunctions
    tests_bssbuild
    tests_bsctools
)

MODE="${1:-full}"

if [[ "${MODE}" != "--report" ]]; then
    if [[ ! -f "${BUILD_DIR}/CMakeCache.txt" ]]; then
        echo ">> configuring coverage build at ${BUILD_DIR}"
        cmake -S "${REPO_ROOT}" -B "${BUILD_DIR}" -G Ninja \
            -DCMAKE_BUILD_TYPE=RelWithDebInfo \
            -DCMAKE_C_COMPILER="${CC_BIN}" \
            -DCMAKE_CXX_COMPILER="${CXX_BIN}" \
            -DGBS_COVERAGE=ON \
            -DGBS_USE_PYTHON_BINDINGS=OFF \
            -DGBS_BUILD_STUBS=OFF
    fi
    echo ">> building priority test targets"
    ninja -C "${BUILD_DIR}" "${TARGETS[@]}"
fi

"${REPO_ROOT}/scripts/coverage_report.sh" "${BUILD_DIR}" "${TARGETS[@]}"
