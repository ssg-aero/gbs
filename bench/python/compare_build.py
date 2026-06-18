"""Front B — curve construction: pygbs vs scipy.interpolate (vs geomdl).

Two build operations, full wall time of the build (the operation IS the build):

  * INTERPOLATION (deg 3, C2 through all points) — incl. the real 201-point
    case (#34). pygbs interpolate_cn vs scipy make_interp_spline, BOTH on a
    chord-length parametrization so it is like-for-like (interpolation, not
    smoothing). geomdl interpolate_curve if present.
  * APPROXIMATION (least-squares, FIXED ~N/4 poles) — pygbs approx vs scipy
    make_lsq_spline with a clamped knot vector giving the same coefficient
    count. This matches gbs's fixed-pole LSQ (not scipy splprep's smoothing-s
    fit, which solves a different problem; noted in the audit).

Run:  PYTHONPATH=build/python micromamba run -n dev python bench/python/compare_build.py
"""
import numpy as np
import _harness as H

gbs = H.gbs_module()
from scipy.interpolate import make_interp_spline, make_lsq_spline

try:
    from geomdl import fitting as gfit
    HAVE_GEOMDL = True
except Exception:
    HAVE_GEOMDL = False

SIZES = [50, 100, 201, 400, 800]


def sample_curve_pts(n):
    t = np.linspace(0.0, 1.0, n)
    return np.column_stack([np.cos(6 * t), np.sin(6 * t), 2 * t])


def chord_param(pts):
    d = np.linalg.norm(np.diff(pts, axis=0), axis=1)
    u = np.concatenate([[0.0], np.cumsum(d)])
    return u / u[-1]


def main():
    print("=== gbs vs scipy (vs geomdl) — Python front, CURVE build (#86) ===")
    for line in H.env_banner():
        print(" ", line)
    if not HAVE_GEOMDL:
        print("  (geomdl absent — its columns are skipped)")

    md = []
    k = 3

    # ---------------- interpolation ----------------
    tb = H.Table("CURVE interpolation deg 3 — full build (median ms, lower=better)",
                 (["N", "pygbs", "scipy", "geomdl"] if HAVE_GEOMDL else ["N", "pygbs", "scipy"]))
    for N in SIZES:
        pts = sample_curve_pts(N)
        ptl = [list(p) for p in pts]
        u = chord_param(pts)
        mg, _ = H.measure(lambda: gbs.interpolate_cn(ptl, k, gbs.KnotsCalcMode.CHORD_LENGTH))
        ms, _ = H.measure(lambda: make_interp_spline(u, pts, k=k))
        if HAVE_GEOMDL:
            mgd, _ = H.measure(lambda: gfit.interpolate_curve(ptl, k), repeats=3)
            tb.add(N, mg, ms, mgd)
        else:
            tb.add(N, mg, ms)
    md.append(tb.markdown())

    # ---------------- approximation (fixed poles, LSQ) ----------------
    ta = H.Table("CURVE approximation deg 3 (~N/4 poles) — LSQ fit (median ms)",
                 ["N", "pygbs", "scipy"])
    for N in [100, 201, 400, 800]:
        pts = sample_curve_pts(N)
        ptl = [list(p) for p in pts]
        u = chord_param(pts)
        n_poles = max(8, N // 4)
        n_interior = n_poles - k - 1
        interior = np.quantile(u, np.linspace(0, 1, n_interior + 2)[1:-1])
        t = np.r_[np.full(k + 1, u[0]), interior, np.full(k + 1, u[-1])]
        mg, _ = H.measure(lambda: gbs.approx(ptl, k, n_poles, gbs.KnotsCalcMode.CHORD_LENGTH))
        ms, _ = H.measure(lambda: make_lsq_spline(u, pts, t, k))
        ta.add(N, mg, ms)
    md.append(ta.markdown())

    with open("bench/python/results_build.md", "w") as f:
        f.write("\n\n".join(md) + "\n")
    print("\n[written] bench/python/results_build.md")


if __name__ == "__main__":
    main()
