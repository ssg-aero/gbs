"""Shared timing harness for the gbs Python competitive benches (#86, Front B).

Fairness rules baked in here:
  * Warm up every callable once before timing (JIT/interpreter/page-fault warmup).
  * Repeat each measurement; report the MEDIAN and the inter-quartile range (IQR)
    rather than a single best, so dispersion is visible.
  * Time only the operation; build/setup is done outside the timed callable.
  * Pin the environment and record exact versions (see env_banner()).

This module has no third-party deps beyond the standard library + numpy, so it
imports even when scipy/geomdl/pygbs are missing (the drivers gate on those).
"""
import sys
import time
import statistics
import platform


def env_banner():
    """Return a list of 'name: version' strings for everything we can see."""
    lines = [f"python: {platform.python_version()} ({sys.executable})",
             f"platform: {platform.platform()}"]
    for mod in ("numpy", "scipy", "geomdl"):
        try:
            m = __import__(mod)
            lines.append(f"{mod}: {getattr(m, '__version__', '?')}")
        except Exception:
            lines.append(f"{mod}: NOT INSTALLED")
    lines.append(f"pygbs: {gbs_module().__file__}")
    return lines


def gbs_module():
    """Load the gbs Python module, preferring a fresh build-tree module (bare
    `gbs`, via PYTHONPATH=build/python) over the installed `pygbs.gbs` package.
    The build-tree one is what carries the post-#88 parallel bulk evaluators."""
    try:
        import gbs  # build/python/gbs.*.so
        return gbs
    except ImportError:
        import pygbs.gbs as gbs
        return gbs


def measure(fn, repeats=11, inner=1):
    """Time `fn` `repeats` times (each doing `inner` calls). Return (median_ms,
    iqr_ms) of the per-batch wall time in milliseconds. Warms up once first."""
    fn()  # warmup (not timed)
    samples = []
    for _ in range(repeats):
        t0 = time.perf_counter()
        for _ in range(inner):
            fn()
        t1 = time.perf_counter()
        samples.append((t1 - t0) * 1e3 / inner)
    samples.sort()
    med = statistics.median(samples)
    n = len(samples)
    q1 = samples[n // 4]
    q3 = samples[(3 * n) // 4]
    return med, (q3 - q1)


class Table:
    """A small fixed-width results table that also remembers its rows so a
    driver can dump a markdown version for the audit."""

    def __init__(self, title, cols):
        self.title = title
        self.cols = cols
        self.rows = []
        print(f"\n{title}")
        print("  ".join(f"{c:>14}" for c in cols))
        print("-" * (16 * len(cols)))

    def add(self, *vals):
        self.rows.append(vals)
        cells = []
        for v in vals:
            cells.append(f"{v:>14.4f}" if isinstance(v, float) else f"{v:>14}")
        print("  ".join(cells))

    def markdown(self):
        out = [f"#### {self.title}", "",
               "| " + " | ".join(self.cols) + " |",
               "|" + "|".join("---" for _ in self.cols) + "|"]
        for r in self.rows:
            cells = [f"{v:.4f}" if isinstance(v, float) else str(v) for v in r]
            out.append("| " + " | ".join(cells) + " |")
        return "\n".join(out)
