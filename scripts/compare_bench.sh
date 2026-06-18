#!/usr/bin/env bash
#
# Configure + build + run the competitive comparison bench (#86), Front A (C++).
# Opt-in: enables GBS_BUILD_COMPARE_BENCH (OFF by default) and links OpenCASCADE.
# Runs inside the `gbs-occt8` micromamba env (mirrors the CI occt8 toolchain:
# clang + OCCT 8 + TBB). A dedicated build dir keeps it apart from the plain
# bench/build tree.
#
# Usage:
#   scripts/compare_bench.sh                 # build + run bench_compare_occt
#   scripts/compare_bench.sh ENV_NAME        # override env (default gbs-occt8)
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ENV_NAME="${1:-gbs-occt8}"
BUILD_DIR="${REPO_ROOT}/bench/build-compare"

micromamba run -n "${ENV_NAME}" cmake -S "${REPO_ROOT}/bench" -B "${BUILD_DIR}" -G Ninja \
    -DCMAKE_BUILD_TYPE=Release \
    -DGBS_BUILD_COMPARE_BENCH=ON \
    -DCMAKE_C_COMPILER=clang \
    -DCMAKE_CXX_COMPILER=clang++

micromamba run -n "${ENV_NAME}" ninja -C "${BUILD_DIR}" bench_compare_occt

echo ">> Running bench_compare_occt"
micromamba run -n "${ENV_NAME}" "${BUILD_DIR}/bench_compare_occt"
