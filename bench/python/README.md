# Front B — Python competitive benches (#86)

Compares the `pygbs` Python API against the Python B-spline ecosystem on the same
operations. See `bench/COMPARE_AUDIT.md` for the measured tables and analysis.

## Pinned environment

Measured on the `dev` micromamba env:

| package | version |
|---------|---------|
| python  | 3.11.13 |
| numpy   | 2.4.3   |
| scipy   | 1.16.1  |
| geomdl  | 5.4.0 (`pip install geomdl`; pure-Python reference, optional) |

`geomdl` is optional — the drivers skip its columns if it is absent.

## pygbs must carry the #88 parallel bulk evaluators

The bulk evaluator (`curve.value(u_list)`) was parallelized in **#88**. A pygbs
built before #88 measures the *serial* path and badly understates gbs. Build the
module from this branch and point `PYTHONPATH` at it:

```
ninja -C build pygbs
```

The drivers load the **build-tree** module preferentially (bare `import gbs` from
`build/python`) and fall back to an installed `pygbs.gbs` only if that fails — the
banner prints which `.so` was used, so always check it points at `build/python`.

## Run

```
PYTHONPATH=build/python:bench/python micromamba run -n dev python bench/python/compare_curve.py
PYTHONPATH=build/python:bench/python micromamba run -n dev python bench/python/compare_build.py
```

- `compare_curve.py` — curve evaluation (point d0, 1st/2nd deriv), batched vs
  per-call, on the *identical* spline mirrored into scipy/geomdl. Writes
  `results_curve.md`.
- `compare_build.py` — interpolation and least-squares approximation builds.
  Writes `results_build.md`.
- `_harness.py` — timing helpers (warm-up, median + IQR), env banner, module loader.

## Fairness notes

- Same knots/poles/degree feed every library for eval; same chord-length
  parametrization for interpolation (interpolation, not smoothing).
- Batched (`crv.value(U)` / `spl(U)`) and per-call (one point at a time) are
  reported separately — the pybind11 boundary cost differs from the vectorized
  path. A numpy array (not a Python list) is passed to the batched pygbs call.
- Median of 11 repeats, IQR tracked; warm-up call discarded.
