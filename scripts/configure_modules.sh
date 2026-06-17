#!/usr/bin/env bash
#
# Configure build-modules/ with GBS_USE_MODULES=ON (C++20 modules POC).
# Measures the compile-time delta relative to the header-include baseline.
#
# Usage: scripts/configure_modules.sh
set -euo pipefail
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="${REPO_ROOT}/build-modules"

cmake -S "${REPO_ROOT}" -B "${BUILD_DIR}" \
    -G Ninja \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX=/home/sebastien/micromamba/envs/dev \
    -DCMAKE_CXX_COMPILER=clang++ \
    -DGBS_BUILD_TESTS=ON \
    -DGBS_USE_PYTHON_BINDINGS=OFF \
    -DGBS_USE_RENDER=ON \
    -DGBS_USE_MODULES=ON \
    -DGBS_USE_PCH=OFF
