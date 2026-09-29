#!/usr/bin/env bash
#
# Run the Python tests against the bindings of the local build/ tree, without
# installing pygbs into the active environment.
#
# Usage:
#   scripts/pytest_build_tree.sh                          # whole python/tests suite
#   scripts/pytest_build_tree.sh python/tests/test_x.py   # selected tests / pytest args
#
# Builds the pygbs target, assembles build/pytree/pygbs/ (compiled module +
# python/*.py, symlinked) and runs pytest with it first on PYTHONPATH.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="${REPO_ROOT}/build"
PKG_DIR="${BUILD_DIR}/pytree/pygbs"

echo ">> Building pygbs"
cmake --build "${BUILD_DIR}" --target pygbs

rm -rf "${PKG_DIR}"
mkdir -p "${PKG_DIR}"
ln -s "${BUILD_DIR}"/python/gbs*.so "${PKG_DIR}/"
ln -s "${REPO_ROOT}"/python/*.py "${PKG_DIR}/"

export PYTHONPATH="${BUILD_DIR}/pytree${PYTHONPATH:+:${PYTHONPATH}}"
export LD_LIBRARY_PATH="${BUILD_DIR}/gbs-render:${LD_LIBRARY_PATH:-}"
export PYVISTA_OFF_SCREEN=true

cd "${REPO_ROOT}"
if [[ $# -gt 0 ]]; then
    python -m pytest "$@"
else
    python -m pytest python/tests/
fi
