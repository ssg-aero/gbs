"""Front B — curve evaluation: pygbs vs scipy.interpolate (vs geomdl).

We build ONE degree-3 curve in pygbs, then mirror the *identical* spline
(same knots/poles/degree) into scipy.BSpline (and geomdl if present), so the
eval comparison is on the same geometry — it isolates evaluator speed, not
construction. Two regimes, per the #86 fairness note:

  * BATCHED: pygbs crv.value(U)  (parallel bulk evaluator, #88)
             vs scipy spl(U)     (vectorized over the numpy array)
             vs geomdl evaluate_list (pure-python per-point loop)
  * PER-CALL: a Python loop calling value(u) one point at a time — this exposes
              the pybind11 / call-boundary overhead, where scipy's C inner loop
              and pygbs's per-call cost differ from the batched path.

Run:  PYTHONPATH=build/python micromamba run -n dev python bench/python/compare_curve.py
(or just `python compare_curve.py` if pygbs is installed in the env)
"""
import numpy as np
import _harness as H

gbs = H.gbs_module()
from scipy.interpolate import BSpline

try:
    from geomdl import BSpline as GBSpline
    HAVE_GEOMDL = True
except Exception:
    HAVE_GEOMDL = False

SIZES = [1_000, 10_000, 100_000, 1_000_000]


def sample_curve_pts(n):
    t = np.linspace(0.0, 1.0, n)
    return np.column_stack([np.cos(6 * t), np.sin(6 * t), 2 * t])


def build():
    pts = sample_curve_pts(60)
    crv = gbs.interpolate_cn([list(p) for p in pts], 3, gbs.KnotsCalcMode.CHORD_LENGTH)
    t = np.asarray(crv.knotsFlats(), dtype=float)
    c = np.asarray(crv.poles(), dtype=float)
    k = crv.degree()
    spl = BSpline(t, c, k)
    g = None
    if HAVE_GEOMDL:
        g = GBSpline.Curve()
        g.degree = k
        g.ctrlpts = c.tolist()
        # geomdl needs the knot vector normalized to [0,1]
        tn = (t - t[0]) / (t[-1] - t[0])
        g.knotvector = tn.tolist()
    return crv, spl, g, crv.bounds()


def main():
    print("=== gbs vs scipy (vs geomdl) — Python front, CURVE eval (#86) ===")
    for line in H.env_banner():
        print(" ", line)
    if not HAVE_GEOMDL:
        print("  (geomdl absent — its columns are skipped; install with `pip install geomdl`)")

    crv, spl, g, (u0, u1) = build()

    # sanity: identical geometry
    us = np.linspace(u0, u1, 97)
    a = np.asarray(crv.value(list(us)))
    b = spl(us)
    print(f"\n[sanity] max|pygbs - scipy| over 97 pts = {np.max(np.abs(a - b)):.2e}")

    md = []
    for d in (0, 1, 2):
        # ---- batched ----
        tb = H.Table(f"CURVE d={d} — BATCHED (median ms, lower=better)",
                     ["N", "pygbs", "scipy", "geomdl"] if HAVE_GEOMDL else ["N", "pygbs", "scipy"])
        for N in SIZES:
            U = np.linspace(u0, u1, N)
            Ul = U.tolist()
            mg, _ = H.measure(lambda: crv.value(U, d))
            ms, _ = H.measure(lambda: spl(U, nu=d))
            if HAVE_GEOMDL and d == 0 and N <= 100_000:
                mgd, _ = H.measure(lambda: g.evaluate_list(Ul), repeats=3)
                tb.add(N, mg, ms, mgd)
            elif HAVE_GEOMDL:
                tb.add(N, mg, ms, float("nan"))
            else:
                tb.add(N, mg, ms)
        md.append(tb.markdown())

        # ---- per-call (only d=0, smaller N: a python loop is O(N) interpreter) ----
        if d == 0:
            tc = H.Table("CURVE d=0 — PER-CALL python loop (median ms over N calls)",
                         ["N", "pygbs", "scipy"])
            for N in (1_000, 10_000):
                U = np.linspace(u0, u1, N).tolist()
                mg, _ = H.measure(lambda: [crv.value(u) for u in U], repeats=5)
                msc, _ = H.measure(lambda: [spl(u) for u in U], repeats=5)
                tc.add(N, mg, msc)
            md.append(tc.markdown())

    with open("bench/python/results_curve.md", "w") as f:
        f.write("\n\n".join(md) + "\n")
    print("\n[written] bench/python/results_curve.md")


if __name__ == "__main__":
    main()
