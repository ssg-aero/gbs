"""
Analyze clang -ftime-trace JSON files from build-profile/.

Usage:
    python scratch/analyze_ftime_trace.py [--build-dir build-profile] [--top N]

Reads all *.json trace files under <build-dir>/CMakeFiles/,
summarizes time by event category, and lists the top-N most expensive
template instantiations.

Output columns (all in ms):
  total_s  : sum of event durations for that category (may overlap)
  count    : number of events
  top_pct  : fraction of the wall-clock Frontend time this category takes

Categories that matter:
  Frontend             - total clang frontend time (parse + sema + codegen)
  Backend              - LLVM IR lowering + optimization + MC emit
  ParseTemplate        - parsing template syntax
  InstantiateTemplate  - template instantiation
  Source               - include-file parse time (also reports source file)
  IncludeFile          - recursive header parsing
"""

import argparse
import json
import pathlib
import sys
from collections import defaultdict

def load_traces(build_dir: pathlib.Path):
    events = []
    files = sorted(build_dir.rglob("*.json"))
    if not files:
        print(f"No JSON trace files found under {build_dir}", file=sys.stderr)
        sys.exit(1)
    print(f"Found {len(files)} trace file(s):")
    for f in files:
        rel = f.relative_to(build_dir.parent) if build_dir.parent in f.parents else f
        print(f"  {rel}")
        with open(f) as fh:
            data = json.load(fh)
        for ev in data.get("traceEvents", []):
            if ev.get("ph") == "X":
                ev["_file"] = str(f)
                events.append(ev)
    return events


def summarize(events, top_n: int):
    by_category = defaultdict(list)
    for ev in events:
        by_category[ev["name"]].append(ev["dur"])

    # Frontend wall time per file
    frontend_by_file = defaultdict(int)
    for ev in events:
        if ev["name"] == "Frontend":
            frontend_by_file[ev["_file"]] += ev["dur"]

    total_frontend_us = sum(frontend_by_file.values())
    total_backend_us = sum(d for ev in events if ev["name"] == "Backend" for d in [ev["dur"]])

    print("\n=== Per-file Frontend time (ms) ===")
    for f, us in sorted(frontend_by_file.items(), key=lambda x: -x[1]):
        fname = pathlib.Path(f).stem.replace(".cpp", "")
        print(f"  {fname:50s} {us/1000:8.1f} ms")

    print(f"\n=== Category summary (total across all TUs, ms) ===")
    print(f"  {'Category':<40} {'Sum(ms)':>10} {'Count':>8} {'%Frontend':>10}")
    print(f"  {'-'*40} {'-'*10} {'-'*8} {'-'*10}")
    ordered = sorted(by_category.items(), key=lambda x: -sum(x[1]))
    for name, durs in ordered:
        total_ms = sum(durs) / 1000
        pct = (sum(durs) / total_frontend_us * 100) if total_frontend_us else 0
        print(f"  {name:<40} {total_ms:>10.1f} {len(durs):>8}   {pct:>8.1f}%")

    print(f"\n  Frontend total : {total_frontend_us/1000:,.1f} ms")
    print(f"  Backend total  : {total_backend_us/1000:,.1f} ms")

    # Top instantiations
    inst_events = [(ev["args"].get("detail", "?"), ev["dur"])
                   for ev in events
                   if ev["name"] == "InstantiateTemplate"]
    if inst_events:
        print(f"\n=== Top {top_n} most expensive template instantiations (ms) ===")
        top = sorted(inst_events, key=lambda x: -x[1])[:top_n]
        for detail, dur_us in top:
            print(f"  {dur_us/1000:8.2f} ms  {detail[:120]}")

    # Top include files
    inc_events = [(ev["args"].get("detail", "?"), ev["dur"])
                  for ev in events
                  if ev["name"] == "Source"]
    if inc_events:
        # Aggregate by file
        inc_agg = defaultdict(int)
        inc_cnt = defaultdict(int)
        for detail, dur in inc_events:
            inc_agg[detail] += dur
            inc_cnt[detail] += 1
        print(f"\n=== Top {top_n} most expensive source/header files (aggregate ms) ===")
        top_inc = sorted(inc_agg.items(), key=lambda x: -x[1])[:top_n]
        for detail, total_us in top_inc:
            print(f"  {total_us/1000:8.2f} ms  (x{inc_cnt[detail]})  {detail[-100:]}")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--build-dir", default="build-profile",
                   help="Build directory containing trace JSON files (default: build-profile)")
    p.add_argument("--top", type=int, default=20,
                   help="Number of top entries to show (default: 20)")
    args = p.parse_args()

    build_dir = pathlib.Path(args.build_dir)
    if not build_dir.exists():
        print(f"ERROR: build directory '{build_dir}' does not exist. Run scripts/configure_profile.sh first.")
        sys.exit(1)

    events = load_traces(build_dir)
    summarize(events, args.top)


if __name__ == "__main__":
    main()
