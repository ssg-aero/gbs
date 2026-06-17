#!/usr/bin/env bash
#
# Build & run the gbs-occt tests against a local OCCT env, mirroring the CI
# `occt` job: build target gbs-occt_tests, then run it.
#
# Usage:
#   scripts/occt_build_test.sh [ENV_NAME] [BUILD_DIR] [DOCTEST_FILTER]
# Defaults: ENV_NAME=gbs-occt8  BUILD_DIR=build-occt8  (no filter)
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ENV_NAME="${1:-gbs-occt8}"
BUILD_DIR="${2:-${REPO_ROOT}/build-occt8}"
FILTER="${3:-}"

micromamba run -n "${ENV_NAME}" cmake --build "${BUILD_DIR}" --target gbs-occt_tests --parallel

ENV_PREFIX="$(micromamba run -n "${ENV_NAME}" printenv CONDA_PREFIX)"
RUN_DIR="${BUILD_DIR}/gbs-occt/tests"
BIN="${RUN_DIR}/gbs-occt_tests"

mkdir -p "${RUN_DIR}/tests/out"   # IGES writers use relative tests/out/ paths
cd "${RUN_DIR}"
if [[ -n "${FILTER}" ]]; then
    LD_LIBRARY_PATH="${ENV_PREFIX}/lib:${LD_LIBRARY_PATH:-}" "${BIN}" -tc="${FILTER}"
else
    LD_LIBRARY_PATH="${ENV_PREFIX}/lib:${LD_LIBRARY_PATH:-}" "${BIN}"
fi
