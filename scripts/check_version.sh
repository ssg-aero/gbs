#!/usr/bin/env bash
#
# Check that the project() VERSION propagates to the installed CMake package
# version file and to the Python bindings (pygbs.gbs.__version__).
#
# Usage:
#   scripts/check_version.sh
#
# Re-configures the existing build/ directory (keeps its cache), builds only
# the pygbs target, then prints both version strings.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="${REPO_ROOT}/build"

echo ">> Re-configuring ${BUILD_DIR}"
cmake -S "${REPO_ROOT}" -B "${BUILD_DIR}" > /dev/null

echo ">> CMake package version file"
grep 'set(PACKAGE_VERSION ' "${BUILD_DIR}/GBSConfigVersion.cmake"

echo ">> Building pygbs"
# build/ may have been configured from a shell whose compiler activation added
# -isystem $CONDA_PREFIX/include: CMake then treats it as implicit and drops it
# from the command line (pybind11 headers live there).
if [[ -n "${CONDA_PREFIX:-}" ]]; then
    export CPLUS_INCLUDE_PATH="${CONDA_PREFIX}/include${CPLUS_INCLUDE_PATH:+:${CPLUS_INCLUDE_PATH}}"
fi
cmake --build "${BUILD_DIR}" --target pygbs

echo ">> pygbs.gbs.__version__"
export LD_LIBRARY_PATH="${BUILD_DIR}/gbs-render:${LD_LIBRARY_PATH:-}"
PYTHONPATH="${BUILD_DIR}/python" python -c "import gbs; print(gbs.__version__)"
