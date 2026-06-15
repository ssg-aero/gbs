#!/usr/bin/env bash
#
# Run already-built, coverage-instrumented test binaries and aggregate a
# source-based coverage report with llvm-cov.
#
# Usage:
#   scripts/coverage_report.sh <build_dir> <test_exe> [<test_exe> ...]
#
# Produces, under <build_dir>/coverage/:
#   - merged.profdata          (llvm-profdata merge of every run)
#   - report.txt               (llvm-cov report, per-file region/line/branch %)
#   - report_core.txt          (same, restricted to the priority perimeter)
#   - uncovered_core.txt       (llvm-cov show, only-regions, priority headers)
#
# The llvm tooling is auto-selected to match the clang that built the binaries
# (clang-21 -> llvm-profdata-21/llvm-cov-21); override with LLVM_PROFDATA /
# LLVM_COV. Binaries are run from the repo root so their relative input paths
# (tests/in/...) resolve.
set -euo pipefail

if [[ $# -lt 2 ]]; then
    echo "usage: $0 <build_dir> <test_exe> [<test_exe> ...]" >&2
    exit 2
fi

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="$1"; shift
TESTS=("$@")

# Skip micro-benchmarks: they only re-run already-covered eval paths millions of
# times (10M curve evals, 100M-element std::transform timing). Under coverage
# instrumentation they dominate wall-clock and add zero coverage. Override with
# DOCTEST_EXCLUDE="" to run everything.
# curve_reparametrized runs a full float arc-length integration (abs_curv over
# ~1080 samples + reparametrized re-evaluation); under instrumentation it alone
# takes minutes and only re-covers bsctools arc-length paths already hit cheaply
# elsewhere. Excluded from the baseline (documented in bench/COVERAGE_AUDIT.md).
DOCTEST_EXCLUDE="${DOCTEST_EXCLUDE-*perf*,cpp_algo.par_vs_seq,tests_bscurve.curve_reparametrized}"

COV_DIR="${BUILD_DIR}/coverage"
RAW_DIR="${COV_DIR}/raw"
mkdir -p "${RAW_DIR}"
rm -f "${RAW_DIR}"/*.profraw

# Pick the llvm tools that match the compiler version, falling back to unversioned.
detect_tool() {
    local base="$1" override="$2"
    if [[ -n "${override}" ]]; then echo "${override}"; return; fi
    for cand in "${base}-21" "${base}-20" "${base}"; do
        if command -v "${cand}" >/dev/null 2>&1; then echo "${cand}"; return; fi
    done
    echo "${base}"
}
PROFDATA="$(detect_tool llvm-profdata "${LLVM_PROFDATA:-}")"
COV="$(detect_tool llvm-cov "${LLVM_COV:-}")"

# Pin each run to a single core. The eval/interp paths use std::execution::par
# (TBB backend); with coverage instrumentation the per-region profile counters
# become cache-line-contended across all worker threads, so a parallel sweep
# that takes milliseconds balloons to minutes of 60-core thrash. Coverage is
# timing-independent, so serializing on one core gives identical coverage in a
# fraction of the wall time. Disable with GBS_COV_PIN=0.
PIN=()
if [[ "${GBS_COV_PIN:-1}" != "0" ]] && command -v taskset >/dev/null 2>&1; then
    PIN=(taskset -c 0)
fi

# Run each instrumented binary, one profraw per run. llvm-cov wants the FIRST
# binary positional and the rest behind -object, otherwise it treats a trailing
# positional source file as the main object ("not a valid object file").
OBJECTS=()
for t in "${TESTS[@]}"; do
    bin="${BUILD_DIR}/tests/${t}"
    if [[ ! -x "${bin}" ]]; then
        echo "!! missing binary ${bin}, skipping" >&2
        continue
    fi
    echo ">> running ${t}"
    DT_ARGS=()
    [[ -n "${DOCTEST_EXCLUDE}" ]] && DT_ARGS+=("-tce=${DOCTEST_EXCLUDE}")
    ( cd "${REPO_ROOT}" && LLVM_PROFILE_FILE="${RAW_DIR}/${t}.profraw" "${PIN[@]}" "${bin}" "${DT_ARGS[@]}" >/dev/null )
    if [[ ${#OBJECTS[@]} -eq 0 ]]; then
        OBJECTS+=("${bin}")          # main object (positional)
    else
        OBJECTS+=(-object "${bin}")  # additional objects
    fi
done

if [[ ${#OBJECTS[@]} -eq 0 ]]; then
    echo "!! no instrumented binaries ran" >&2
    exit 1
fi

echo ">> merging profiles with ${PROFDATA}"
"${PROFDATA}" merge -sparse "${RAW_DIR}"/*.profraw -o "${COV_DIR}/merged.profdata"

# Only attribute coverage to the library headers, never to test/3rd-party code.
# llvm-cov has no "only this dir" flag, so drop everything else by regex.
IGNORE_RE='(/micromamba/|/tests/|/testing/|/tests_3rdlibs/|/gbs-render/|/gbs-io/|/gbs-mesh/|/gbs-occt/|/gbs-cuda/|/Eigen/|/eigen|/usr/|/bench/)'
# Priority perimeter: the headers whose coverage this audit ranks.
CORE_FILES=(
    "${REPO_ROOT}/gbs/basisfunctions.ixx"
    "${REPO_ROOT}/gbs/bscurve.h"
    "${REPO_ROOT}/gbs/bssurf.h"
    "${REPO_ROOT}/gbs/bscinterp.h"
    "${REPO_ROOT}/gbs/bssinterp.h"
    "${REPO_ROOT}/gbs/knotsfunctions.ixx"
    "${REPO_ROOT}/gbs/bssbuild.h"
    "${REPO_ROOT}/gbs/bsctools.h"
)

echo ">> writing reports to ${COV_DIR}"
"${COV}" report "${OBJECTS[@]}" -instr-profile="${COV_DIR}/merged.profdata" \
    -show-branch-summary -ignore-filename-regex="${IGNORE_RE}" > "${COV_DIR}/report.txt"

# Focused per-file view of the priority perimeter.
"${COV}" report "${OBJECTS[@]}" -instr-profile="${COV_DIR}/merged.profdata" \
    -show-branch-summary "${CORE_FILES[@]}" > "${COV_DIR}/report_core.txt"

# Annotated source for the priority headers, regions + branch counts, so the
# audit can point at the exact uncovered branches.
"${COV}" show "${OBJECTS[@]}" -instr-profile="${COV_DIR}/merged.profdata" \
    -show-branches=count -show-regions=true -show-line-counts=true \
    "${CORE_FILES[@]}" > "${COV_DIR}/annotated_core.txt"

echo
echo "==== coverage (priority perimeter) ===="
cat "${COV_DIR}/report_core.txt"
echo
echo "full report:    ${COV_DIR}/report.txt"
echo "annotated core: ${COV_DIR}/annotated_core.txt"
