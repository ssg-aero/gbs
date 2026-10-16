"""Compares two runs of bench_extrema query by query.

    bench_extrema before.csv   # built against the headers of the reference version
    bench_extrema after.csv    # built against the headers of the new version
    python bench/python/compare_extrema.py before.csv after.csv

For each case: the queries where each version reaches a smaller true distance (by more
than 1e-9), the equal ones, the largest gains and losses, and the queries that threw.
"""
import collections
import csv
import math
import sys


def load(path):
    return {(r[0], int(r[1])): float(r[2]) for r in csv.reader(open(path))}


def main(before_path, after_path):
    before, after = load(before_path), load(after_path)
    stats = collections.defaultdict(lambda: {"after": 0, "before": 0, "equal": 0, "gain": 0.0, "loss": 0.0, "threw": 0})
    for key, d_after in after.items():
        s = stats[key[0]]
        d_before = before.get(key, math.nan)
        if math.isnan(d_before) or math.isnan(d_after):
            s["threw"] += 1
        elif d_after < d_before - 1e-9:
            s["after"] += 1
            s["gain"] = max(s["gain"], d_before - d_after)
        elif d_before < d_after - 1e-9:
            s["before"] += 1
            s["loss"] = max(s["loss"], d_after - d_before)
        else:
            s["equal"] += 1
    print(f"{'case':28}{'after better':>13}{'before better':>14}{'equal':>7}{'max gain':>10}{'max loss':>10}{'threw':>7}")
    for name, s in stats.items():
        print(f"{name:28}{s['after']:>13}{s['before']:>14}{s['equal']:>7}{s['gain']:>10.2e}{s['loss']:>10.2e}{s['threw']:>7}")


if __name__ == "__main__":
    main(*sys.argv[1:3])
