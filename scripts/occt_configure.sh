#!/usr/bin/env bash
#
# Configure a dedicated OCCT build tree using a local micromamba env that
# mirrors the CI `occt` job (OCCT + clang + deps). Used to reproduce the
# occt8 / occt7 CI jobs locally.
#
# Usage:
#   scripts/occt_configure.sh [ENV_NAME] [BUILD_DIR]
# Defaults: ENV_NAME=gbs-occt8  BUILD_DIR=build-occt8
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ENV_NAME="${1:-gbs-occt8}"
BUILD_DIR="${2:-${REPO_ROOT}/build-occt8}"

micromamba run -n "${ENV_NAME}" cmake -B "${BUILD_DIR}" -S "${REPO_ROOT}" -G Ninja \
    -DCMAKE_BUILD_TYPE=Release \
    -DGBS_BUILD_TESTS=ON \
    -DGBS_USE_OCCT_UTILS=ON \
    -DGBS_USE_PYTHON_BINDINGS=OFF \
    -DGBS_BUILD_STUBS=OFF \
    -DGBS_USE_CUDA=OFF \
    -DGBS_USE_MODULES=OFF \
    -DGBS_BUILD_DOC=OFF \
    -DCMAKE_C_COMPILER=clang \
    -DCMAKE_CXX_COMPILER=clang++
