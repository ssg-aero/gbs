#!/usr/bin/env bash
#
# Configure + build the standalone bench/ project (Release) and run a bench.
#
# Usage:
#   scripts/build_bench.sh                 # build all bench targets
#   scripts/build_bench.sh bench_loft      # build + run a single bench
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="${REPO_ROOT}/bench/build"

cmake -S "${REPO_ROOT}/bench" -B "${BUILD_DIR}" -G Ninja -DCMAKE_BUILD_TYPE=Release >/dev/null

if [[ $# -gt 0 ]]; then
    ninja -C "${BUILD_DIR}" "$@"
    for t in "$@"; do
        bin="${BUILD_DIR}/${t}"
        if [[ -x "${bin}" ]]; then
            echo ">> Running ${t}"
            "${bin}"
        fi
    done
else
    ninja -C "${BUILD_DIR}"
fi
