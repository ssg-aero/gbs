#!/usr/bin/env bash
#
# Configure build-profile/ for compile-time profiling with -ftime-trace.
# Use this to capture clang JSON trace files for analysis.
# The -ftime-trace flag adds ~2-5 % overhead; profile data in the trace JSON files.
#
# Usage: scripts/configure_profile.sh
set -euo pipefail
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="${REPO_ROOT}/build-profile"

cmake -S "${REPO_ROOT}" -B "${BUILD_DIR}" \
    -G Ninja \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX=/home/sebastien/micromamba/envs/dev \
    -DCMAKE_CXX_COMPILER=clang++ \
    -DCMAKE_CXX_FLAGS="-ftime-trace" \
    -DGBS_BUILD_TESTS=ON \
    -DGBS_USE_PYTHON_BINDINGS=OFF \
    -DGBS_USE_RENDER=ON \
    -DGBS_USE_MODULES=OFF \
    -DGBS_USE_PCH=OFF
