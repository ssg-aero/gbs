#!/usr/bin/env bash
# Diagnostic: time each priority coverage binary individually (pinned, perf
# excluded) with a per-binary timeout, to spot tests that blow up under
# instrumentation before committing to a full coverage run.
set -uo pipefail
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="${REPO_ROOT}/build-cov"
EXCLUDE='*perf*,cpp_algo.par_vs_seq,tests_bscurve.curve_reparametrized'
CAP="${1:-90}"
TARGETS=(tests_basisfunctions tests_bscurve tests_bssurf tests_curve_derivatives
         tests_surface_derivatives tests_curves tests_surfaces tests_bsinterp
         tests_knotsfunctions tests_bssbuild tests_bsctools)
for t in "${TARGETS[@]}"; do
    bin="${BUILD_DIR}/tests/${t}"
    start=$(date +%s)
    timeout "${CAP}" taskset -c 0 env LLVM_PROFILE_FILE=/tmp/_covt.profraw \
        "${bin}" -tce="${EXCLUDE}" >/dev/null 2>&1
    rc=$?
    end=$(date +%s)
    printf '%-32s %3ds  rc=%s%s\n' "${t}" "$((end-start))" "${rc}" \
        "$([[ ${rc} -eq 124 ]] && echo '  <-- TIMEOUT' || true)"
done
