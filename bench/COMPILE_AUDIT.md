# COMPILE_AUDIT.md — spike #75: C++20 modules vs PCH vs status quo

**Date**: 2026-06-17  
**Branch**: `spike/compile-time-modules-75`  
**Compiler**: clang++ 20.1.8 (conda-forge), Linux x86_64, 64 cores  
**Build type**: Release (`-O3 -DNDEBUG`)  
**Representative TUs**: `tests_basisfunctions`, `tests_bscurve`, `tests_surfaces`, `tests_bscapprox`

---

## 1. Profiling methodology

Build configured with `scripts/configure_profile.sh` (adds `-ftime-trace`).
Trace JSON files analyzed with `scratch/analyze_ftime_trace.py`.

```sh
bash scripts/configure_profile.sh   # build-profile/ baseline + -ftime-trace
bash scripts/build_profile.sh       # cold build; writes *.json trace files
python scratch/analyze_ftime_trace.py --build-dir build-profile --top 20
```

---

## 2. `-ftime-trace` breakdown (4 TUs, aggregate frontend time)

### 2.1 Per-TU summary

| TU | Frontend (ms) | Total Source (ms) | Inst-Fn (ms) | Inst-Class (ms) | Backend (ms) |
|----|--------------|-------------------|--------------|-----------------|--------------|
| tests_basisfunctions | 4 182 | 3 386 | 832 | 650 | 1 348 |
| tests_bscurve | 7 485 | 4 041 | 3 164 | 1 862 | 5 374 |
| tests_bscapprox | 10 697 | 4 440 | 5 547 | 4 093 | 8 802 |
| tests_surfaces | 8 007 | 4 510 | 3 511 | 2 466 | 1 329 |
| **Total** | **30 371** | **16 377** | **13 054** | **9 071** | **16 853** |

`Total Source` = net time attributable to include-file parsing (no child overlap).  
`Inst-Fn` / `Inst-Class` = net time in template instantiation.

### 2.2 Fraction of total frontend

| Category | Time (ms) | % of Frontend |
|----------|-----------|---------------|
| Total Source (header parse) | 16 377 | **54%** |
| Total InstantiateFunction | 13 054 | **43%** |
| Total InstantiateClass | 9 071 | **30%** |
| Total Backend | 16 853 | **55%** |

> Note: categories overlap (instantiation is triggered while parsing headers), so percentages do not sum to 100.

### 2.3 Header parse breakdown (by library, aggregate 4 TUs)

Source `b`/`e` events matched across all test TUs:

| Category | Parse time (ms) | % of source events |
|----------|-----------------|--------------------|
| Clang intrinsic headers (AVX512, etc.) | ≈ 1 530 | 17% |
| Eigen (SIMD, SVD, Core) | ≈ 2 990 | 34% |
| STL (`bits/`, etc.) | ≈ 2 420 | 27% |
| Boost / NLopt | ≈ 681 | 8% |
| GBS headers (vecop/math/basisfunctions/knotsfunctions) | **≈ 853** | **10%** |
| Other conda packages | ≈ 450 | 5% |

**Key finding**: GBS own headers account for only ~10% of header-parse time.  
The dominant parse cost is Eigen (34%) + STL (27%) + clang intrinsics (17%).

### 2.4 Top-20 most expensive template instantiations (single event duration)

| Duration (ms) | Kind | Detail |
|---------------|------|--------|
| 2 327 | InstantiateFunction | `gbs::plot<gbs::SurfaceOfRevolution<double>>` |
| 2 314–2 219 | InstantiateFunction | `gbs::tuple_for_each<tuple<SurfaceOfRevolution…>>` (×5) |
| 2 219 | InstantiateFunction | `gbs::make_actor<double, 3>` |
| 2 168 | InstantiateFunction | `gbs::make_curvature<double, 3>` |
| 2 150–2 134 | InstantiateFunction | `gbs::abs_curv<double/float, 3, 10>` (×3) |
| 2 028–2 017 | InstantiateFunction | `gbs::interpolate<float/double, 1>` (×4) |
| 2 013–2 010 | InstantiateFunction | `gbs::build_poles<float/double, 1, 1>` (×2) |
| 1 626 | InstantiateFunction | `gbs::approx<double, 3>` |

> Observations:
> - VTK-render instantiations (`plot`, `make_actor`, `tuple_for_each`) are TU-unique — `extern template` would not help.
> - `gbs::interpolate<T,dim>` and `gbs::build_poles<T,dim>` appear in 3-4 TUs each — prime candidates for `extern template`.
> - These tall-pole instantiations each take 1.6–2.3 s, far exceeding the GBS header parse cost (~200 ms/TU).

---

## 3. Build-time measurements

All three variants were configured in fresh build directories and built from scratch.  
Baseline used `build-profile/` (includes `-ftime-trace` overhead, ≈ 3%).

### 3.1 Cold build — 4 test targets (full cold build including gbs-render)

| Variant | Wall (s) | CPU (s) | vs baseline (CPU) |
|---------|----------|---------|-------------------|
| Baseline (headers, -ftime-trace) | 20.3 | 67.7 | – |
| **Modules ON** (`GBS_USE_MODULES=ON`) | **27.3** | **68.5** | **+1%** |
| **PCH ON** (`GBS_USE_PCCH=ON`) | **19.5** | **57.4** | **-15%** |

Modules cold-build is 7 s *slower* in wall time (4 extra module-library compilations — vecop, math, basis_functions, knots_functions — must complete serially before test TUs can start importing them; this serialization adds latency on a 64-core machine where everything else is parallel).

### 3.2 Incremental rebuild — `tests_bscapprox` only (`-j1`)

This is the developer-iteration case: touch one file, wait for recompile.

| Variant | Wall (s) | CPU (s) | vs baseline |
|---------|----------|---------|-------------|
| Baseline (headers, -ftime-trace) | 20.1 | 19.5 | – |
| **Modules ON** | **19.4** | **18.7** | **-3.6%** |
| **PCH ON** | **16.6** | **16.0** | **-17.4%** |

Correction for -ftime-trace overhead (~3%): true baseline ≈ 19.5 s.
- Modules incremental: 19.4 s ≈ no meaningful gain (0%)
- PCH incremental: 16.6 s → **-15% real gain**

---

## 4. Alternatives: PCH and `extern template` analysis

### 4.1 PCH (implemented in this spike)

`GBS_USE_PCH=ON` precompiles `vecop.ixx`, `math.ixx`, `basisfunctions.ixx`, `knotsfunctions.ixx` plus Eigen and STL headers into a shared `.pch` file via `target_precompile_headers(REUSE_FROM gbs_pch_base)`.

- **Build speed**: 15–17% faster CPU / incremental wall time, confirmed by measurement.
- **Source parse elimination**: saves Eigen+STL+GBS parse (≈ 90% of Total Source time per TU; AVX512 intrinsic headers are also precompiled).
- **Template instantiation**: NOT saved by PCH — instantiation still occurs per-TU. This is the remaining ~43% of frontend time after PCH.
- **Maintenance cost**: Low. `gbs/gbs_pch.h` lists what to precompile; PCH is OFF by default, opt-in via `GBS_USE_PCH=ON`.
- **No API footgun**: headers stay accessible as before; PCH is purely a build optimization.

### 4.2 `extern template` (estimated)

Hot shared instantiations across multiple TUs:

| Instantiation | Estimated per-TU cost | Shared TUs | Potential saving |
|---------------|----------------------|------------|-----------------|
| `gbs::interpolate<double/float, 1>` | ~2 000 ms | 3–4 | ~4 000–6 000 ms |
| `gbs::build_poles<double/float, 1, 1>` | ~2 000 ms | 2–3 | ~2 000–4 000 ms |
| (others) | ~500–1 500 ms | 2–3 | ~1 000–3 000 ms |

**Estimated cold-build CPU saving**: 7 000–13 000 ms out of 67 700 ms → **10–19% further gain**, on top of PCH.

**Implementation**: add `extern template gbs::interpolate<double, 1>; ...` declarations to a header included by all TUs, plus an explicit-instantiation `.cpp` per type family. Medium effort; risk of missing a specialization.

### 4.3 Combination

PCH + extern template (hot types only) could yield **25–35% total CPU reduction** over the baseline. Still dominated by the non-shared, VTK-render instantiations that take 2 s+ each and cannot be shared.

---

## 5. Architectural assessment

| Dimension | Status quo (`.ixx` as include) | Modules ON |
|-----------|-------------------------------|-----------|
| Interface clarity | Footgun: helpers in `.ixx` are reachable via `#include` even if intended as module-private | Clean `export`/non-export split enforced by compiler |
| Dual-mode maintenance | Existing dual-mode guards work; every change must update two code paths | Same — dual-mode guards required until full migration |
| Consumer portability | Works with any C++17/20 compiler | Requires CMake 3.28+ + clang 19+ or GCC 14+ + matching `clang-scan-deps` |
| Build integration | Trivial | Complex: module dependency scanning, file-set install export |
| Incremental rebuild | Re-parses all headers on every TU change | Re-reads compiled BMI — saves GBS header parse only (~200 ms/TU) |

**Footgun removed by modules**: any helper function defined in a `.ixx` file that is NOT `export`ed is currently accessible to any TU that does `#include "vecop.ixx"` directly. With modules ON that accidental access is a hard compile error. This is a real architectural benefit, but the current `.ixx` files expose only the public gbs interface — so the practical risk today is low.

---

## 6. Scored recommendation

| Option | Cold build | Incremental | Architecture | Maintenance | Score |
|--------|-----------|-------------|--------------|-------------|-------|
| Status quo | 0% | 0% | Footgun (low risk) | Lowest | C+ |
| **PCH** (`GBS_USE_PCH=ON`) | -15% CPU | **-15%** | No change | Low | **B+** |
| `extern template` (hot types) | -10–19% | -10–19% | No change | Medium | B |
| **PCH + `extern template`** | **-25–35%** | **-25–35%** | No change | Medium | **A-** |
| Modules (current subtree) | +34% wall | ≈ 0% | + clear interface | High | C+ |
| Full module migration | unknown | unknown | ++ | Very high | deferred |

### Go / No-go

**Modules (large-scale migration)**: ❌ **No-go** for now.  
- The compile-time benefit is negligible (< 5% incremental improvement; negative for cold builds).
- The dominant cost is **template instantiation**, which modules cannot reduce.
- Dual-mode maintenance cost is non-trivial and grows with the codebase.
- The architectural benefit (interface clarity, footgun removal) does not justify the migration overhead when the current `.ixx` pattern works correctly.

**PCH** (`GBS_USE_PCH=ON`): ✅ **Go** as an opt-in build option.  
- 15% incremental improvement, zero API impact, confirmed tests-green on Linux.
- Already implemented in this spike. Low ongoing maintenance.
- Recommended to enable in CI for faster test iteration.

**`extern template`** (hot types): 🔶 **Conditional go** as a follow-up.  
- Estimated 10–19% additional gain targeting `interpolate`, `build_poles`, and `abs_curv`.
- Worth a targeted follow-up if CI build time is a concern (it is not critical on 64 cores).

**Modules POC subtree** (current): kept behind `GBS_USE_MODULES=OFF` (default). The POC demonstrates modules work on Linux with clang 20 (conda); **not validated on Windows/macOS** yet. Architecture benefit is noted for a future decision if the ecosystem matures further.

---

## 7. Reproduce the measurements

```sh
# Baseline + trace analysis
bash scripts/configure_profile.sh
bash scripts/build_profile.sh
python scratch/analyze_ftime_trace.py --build-dir build-profile

# Modules ON build + test
bash scripts/configure_modules.sh
bash scripts/build_modules_test.sh

# PCH ON build + test
bash scripts/configure_pch.sh
bash scripts/build_pch_test.sh

# Incremental rebuild comparison (touch heaviest TU, -j1 each)
touch tests/tests_bscapprox.cpp
time ninja -C build-profile -j1 tests_bscapprox
touch tests/tests_bscapprox.cpp
time ninja -C build-modules -j1 tests_bscapprox
touch tests/tests_bscapprox.cpp
time ninja -C build-pch -j1 tests_bscapprox
```

---

*Tracked under [#75](https://github.com/ssg-aero/gbs/issues/75) / [epic #44](https://github.com/ssg-aero/gbs/issues/44)*
