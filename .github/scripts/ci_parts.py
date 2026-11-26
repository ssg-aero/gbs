#!/usr/bin/env python3
"""Split the build and the C++ tests of a CI job into parts, and check that the
parts together build and test everything, once.

The reference is what CMake configured: every target of the build tree with an
artifact (executables, libraries, the Python module), read from the CMake File
API. It does not depend on the split.

  plan      after the configuration: the targets of this part (deterministic);
            every tests/tests_*.cpp of the repository must be a target
  select    after the build: the ctest cases of the targets of this part;
            a test executable of the part must register some
  manifest  at the end of the part: what it built, ran, how it ended, and
            its configuration (compilers, CMake options, conda packages)
  check     in the final job of a platform: the manifests of all its parts,
            every case run and passed, the same configuration in each

The split: the targets whose sources include the BREP kernel (the same rule as
the gbs_heavy pool of the CMake files) are dealt to the parts by weight (lines
of their sources, heaviest first, to the lightest part); the Python module goes
to part 1 and weighs as much as BINDINGS_WEIGHT of them; every other target
goes to the last part. A new tests_*.cpp needs no change here.

Standard library only.
"""

import argparse
import glob
import hashlib
import json
import os
import re
import subprocess
import sys
import xml.etree.ElementTree as ET

HEAVY_INCLUDE = re.compile(r"#include <gbs-brep/")
NOT_BUILT = re.compile(r"^(.+)_NOT_BUILT-[0-9a-f]+$")
BINDINGS_TARGET = "pygbs"
BINDINGS_WEIGHT = 6  # in average BREP tests: the bindings and their link
# Test executables that register no case on purpose (any other one that
# registers none fails the part: a test file emptied by mistake)
NO_CASES = {"tests_main"}  # tests/tests_main.cpp: the former gtest runner, empty
# Targets the CI configuration must have, whatever the File API returns
REQUIRED_TARGETS = {"gbs-render", "gbs_testing_main", BINDINGS_TARGET}
# The configuration the parts of a platform must share (CMake cache)
CACHE_ENTRIES = re.compile(r"^(GBS_.*|CMAKE_BUILD_TYPE|CMAKE_(C|CXX)_FLAGS(_RELEASE)?|CMAKE_(C|CXX)_COMPILER(_LAUNCHER)?"
                           r"|CMAKE_(EXE|SHARED|MODULE)_LINKER_FLAGS(_RELEASE)?|CMAKE_CXX_STANDARD)$")
ARTIFACT_TYPES = {"EXECUTABLE", "MODULE_LIBRARY", "SHARED_LIBRARY", "STATIC_LIBRARY"}


def fail(message):
    print(f"::error::{message}")
    sys.exit(1)


def norm(path):
    return os.path.normcase(os.path.normpath(os.path.abspath(path)))


def api_reply(build, kind):
    reply = os.path.join(build, ".cmake", "api", "v1", "reply")
    indexes = sorted(glob.glob(os.path.join(reply, "index-*.json")))
    if not indexes:
        fail(f"no CMake File API reply in {reply}: query {kind} before the configuration")
    index = json.load(open(indexes[-1], encoding="utf-8"))
    if kind not in index["reply"]:
        fail(f"no {kind} in the CMake File API reply: query it before the configuration")
    return reply, json.load(open(os.path.join(reply, index["reply"][kind]["jsonFile"]), encoding="utf-8"))


def configuration(build):
    """What the parts of a platform must share: compilers, CMake options, conda packages."""
    _, toolchains = api_reply(build, "toolchains-v1")
    _, cache = api_reply(build, "cache-v2")
    compilers = {t["language"]: f"{t['compiler'].get('id')} {t['compiler'].get('version')} "
                 f"{os.path.basename(t['compiler'].get('path', ''))}" for t in toolchains["toolchains"]}
    entries = {e["name"]: e["value"] for e in cache["entries"] if CACHE_ENTRIES.match(e["name"])}
    prefix = os.environ.get("CONDA_PREFIX", "")
    packages = sorted(os.path.basename(f)[:-5] for f in glob.glob(os.path.join(prefix, "conda-meta", "*.json")))
    if not packages:
        fail(f"no conda package found in {prefix}/conda-meta (CONDA_PREFIX)")
    digest = hashlib.sha256("\n".join(packages).encode()).hexdigest()[:16]
    return {"compilers": compilers, "cmake": entries, "conda": {"packages": len(packages), "digest": digest}}


def code_model(build):
    reply, model = api_reply(build, "codemodel-v2")
    source_dir = model["paths"]["source"]
    targets = {}
    for entry in model["configurations"][0]["targets"]:
        target = json.load(open(os.path.join(reply, entry["jsonFile"]), encoding="utf-8"))
        if target["type"] not in ARTIFACT_TYPES:
            continue
        sources = [s["path"] for s in target.get("sources", []) if not s.get("isGenerated")]
        sources = [p if os.path.isabs(p) else os.path.join(source_dir, p) for p in sources]
        artifacts = [a["path"] if os.path.isabs(a["path"]) else os.path.join(build, a["path"])
                     for a in target.get("artifacts", [])]
        targets[target["name"]] = {"type": target["type"], "sources": sources, "artifacts": artifacts}
    return targets


def read(path):
    try:
        with open(path, encoding="utf-8", errors="replace") as f:
            return f.read()
    except OSError:
        return ""


def split(targets, parts):
    """{target: part} for parts 1..n"""
    heavy, weight = [], {}
    for name, t in targets.items():
        texts = [read(p) for p in t["sources"] if p.endswith((".cpp", ".cxx", ".cc", ".cu"))]
        if name != BINDINGS_TARGET and any(HEAVY_INCLUDE.search(x) for x in texts):
            heavy.append(name)
            weight[name] = max(1, sum(x.count("\n") for x in texts))
    load = [0] * parts
    assign = {}
    if BINDINGS_TARGET in targets:
        assign[BINDINGS_TARGET] = 1
        mean = sum(weight.values()) / len(weight) if weight else 1
        weight[BINDINGS_TARGET] = round(BINDINGS_WEIGHT * mean)
        load[0] += weight[BINDINGS_TARGET]
    for name in sorted(heavy, key=lambda n: (-weight[n], n)):
        i = min(range(parts), key=lambda i: (load[i], i))
        assign[name] = i + 1
        load[i] += weight[name]
    for name in targets:
        assign.setdefault(name, parts)
    return assign, weight, load


def cmd_plan(a):
    targets = code_model(a.build)
    missing = sorted(REQUIRED_TARGETS - set(targets))
    source_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    files = sorted(os.path.basename(f)[:-4] for f in glob.glob(os.path.join(source_dir, "tests", "tests_*.cpp")))
    missing += [f"{f} (tests/{f}.cpp)" for f in files if f not in targets]
    if missing:
        fail(f"targets missing from the configuration: {missing}")
    config = configuration(a.build)
    assign, weight, load = split(targets, a.parts)
    mine = sorted(n for n, p in assign.items() if p == a.part)
    os.makedirs(a.out, exist_ok=True)
    with open(os.path.join(a.out, "targets.txt"), "w", encoding="utf-8") as f:
        f.write("\n".join(mine) + "\n")
    json.dump({"part": a.part, "parts": a.parts, "reference": sorted(targets), "assigned": mine,
               "executables": sorted(n for n in targets if targets[n]["type"] == "EXECUTABLE"),
               "configuration": config, "artifacts": {n: targets[n]["artifacts"] for n in targets}},
              open(os.path.join(a.out, "plan.json"), "w", encoding="utf-8"), indent=1)
    print(f"Split of {len(targets)} targets in {a.parts} part(s); weights in lines of source:")
    for p in range(1, a.parts + 1):
        names = sorted((n for n, q in assign.items() if q == p), key=lambda n: (-weight.get(n, 0), n))
        heavy = [f"{n} ({weight[n]})" for n in names if n in weight]
        print(f"  part {p}: {len(names)} targets, heavy weight {load[p - 1]}: {', '.join(heavy) or '-'}"
              f"; {len(names) - len(heavy)} light{' (this part)' if p == a.part else ''}")
    print(f"This part builds: {' '.join(mine)}")
    print(f"Configuration: {json.dumps(config['compilers'])}, {config['conda']['packages']} conda packages "
          f"({config['conda']['digest']}), {len(config['cmake'])} CMake entries compared between the parts")
    want_bindings = "true" if BINDINGS_TARGET in mine else "false"
    if os.environ.get("GITHUB_OUTPUT"):
        with open(os.environ["GITHUB_OUTPUT"], "a", encoding="utf-8") as f:
            f.write(f"bindings={want_bindings}\n")


def ctest_list(build):
    out = subprocess.run(["ctest", "--test-dir", build, "--show-only=json-v1"],
                         check=True, capture_output=True, text=True).stdout
    return json.loads(out)["tests"]


def cmd_select(a):
    plan = json.load(open(os.path.join(a.plan, "plan.json"), encoding="utf-8"))
    mine = set(plan["assigned"])
    by_artifact = {norm(p): n for n, ps in plan["artifacts"].items() for p in ps}
    exe_of = {}
    for n, ps in plan["artifacts"].items():
        for p in ps:
            exe_of.setdefault(n, p)
    selected, problems = [], []
    cases_of = set()
    seen_targets = set()
    for test in ctest_list(a.build):
        name = test["name"]
        m = NOT_BUILT.match(name)
        if m:
            target = m.group(1)
            if target not in plan["reference"]:
                problems.append(f"{name}: placeholder of an unknown target")
            elif target in mine:
                problems.append(f"{target} is in this part but was not built ({name})")
            seen_targets.add(target)
            continue
        command = test.get("command") or []
        if not command and name in plan["reference"]:
            # add_test(NAME <target> COMMAND <target>) of a target not built here:
            # ctest cannot resolve its executable
            if name in mine:
                problems.append(f"{name} is in this part but was not built")
            seen_targets.add(name)
            continue
        target = by_artifact.get(norm(command[0])) if command else None
        if target is None:
            problems.append(f"{name}: its command {command[:1]} is no artifact of the build tree")
            continue
        seen_targets.add(target)
        if target in mine:
            if not os.path.exists(command[0]):
                problems.append(f"{name}: {command[0]} was not built")
            selected.append(name)
            cases_of.add(target)
    problems += [f"{t} registers no ctest case" for t in sorted(mine & set(plan["executables"]))
                 if t not in cases_of and t not in NO_CASES]
    if problems:
        for p in problems:
            print(f"::error::{p}")
        sys.exit(1)
    if len(set(selected)) != len(selected):
        fail("a ctest case is listed twice")
    with open(os.path.join(a.plan, "tests.txt"), "w", encoding="utf-8") as f:
        f.write("\n".join(selected) + "\n")
    print(f"{len(selected)} ctest cases in this part, from "
          f"{len(seen_targets & mine)} of its {len(mine)} targets (the others register none)")
    if not selected:
        fail("no ctest case in this part")


def junit_cases(path):
    cases = {}
    if not path or not os.path.exists(path):
        return cases
    for case in ET.parse(path).getroot().iter("testcase"):
        cases[case.get("name")] = case.get("status") or "run"
    return cases


def cmd_manifest(a):
    plan = json.load(open(os.path.join(a.plan, "plan.json"), encoding="utf-8"))
    tests_file = os.path.join(a.plan, "tests.txt")
    selected = [t for t in read(tests_file).splitlines() if t]
    cases = junit_cases(a.junit)
    manifest = {
        "platform": a.platform, "part": plan["part"], "parts": plan["parts"], "status": a.status,
        "reference": plan["reference"], "assigned": plan["assigned"], "selected": selected,
        "ran": sorted(cases), "failed": sorted(n for n, s in cases.items() if s not in ("run", "passed")),
        "python_tests": a.python_tests, "configuration": plan["configuration"],
    }
    json.dump(manifest, open(a.out, "w", encoding="utf-8"), indent=1)
    print(f"part {plan['part']}/{plan['parts']} ({a.status}): {len(plan['assigned'])} targets, "
          f"{len(selected)} cases selected, {len(cases)} run, {len(manifest['failed'])} failed, "
          f"Python tests: {a.python_tests}")


def cmd_check(a):
    errors = []
    if a.needs_result != "success":
        errors.append(f"the parts ended with '{a.needs_result}', not 'success'")
    files = sorted(glob.glob(os.path.join(a.manifests, "**", "*.json"), recursive=True))
    manifests = [json.load(open(f, encoding="utf-8")) for f in files]
    manifests = [m for m in manifests if m.get("platform") == a.platform]
    parts = sorted(m["part"] for m in manifests)
    if parts != list(range(1, a.parts + 1)):
        errors.append(f"manifests of parts {parts}, expected 1 to {a.parts}")
    lines = []
    for m in sorted(manifests, key=lambda m: m["part"]):
        lines.append(f"part {m['part']}/{m['parts']}: {m['status']}, {len(m['assigned'])} targets, "
                     f"{len(m['ran'])} ctest cases, Python tests {m['python_tests']}")
        for name in m["failed"]:
            lines.append(f"    failed: {name}")
        if m["failed"]:
            errors.append(f"part {m['part']}: {len(m['failed'])} cases failed, not run or disabled {m['failed'][:5]}")
        if m["status"] != "success":
            errors.append(f"part {m['part']} ended '{m['status']}'")
        if sorted(m["ran"]) != sorted(m["selected"]):
            missing = sorted(set(m["selected"]) - set(m["ran"]))
            extra = sorted(set(m["ran"]) - set(m["selected"]))
            errors.append(f"part {m['part']}: {len(missing)} selected cases not run {missing[:5]}, "
                          f"{len(extra)} run but not selected {extra[:5]}")
    if manifests:
        references = {tuple(m["reference"]) for m in manifests}
        if len(references) != 1:
            errors.append("the parts did not configure the same targets")
        configurations = {json.dumps(m.get("configuration"), sort_keys=True) for m in manifests}
        if len(configurations) != 1:
            for m in manifests:
                lines.append(f"part {m['part']} configuration: {json.dumps(m.get('configuration'), sort_keys=True)[:400]}")
            errors.append("the parts did not build with the same compilers, CMake options and conda packages "
                          "(a package published between their environment solves: re-run the jobs)")
        reference = set(manifests[0]["reference"])
        owners = {}
        for m in manifests:
            for t in m["assigned"]:
                owners.setdefault(t, []).append(m["part"])
        twice = {t: p for t, p in owners.items() if len(p) > 1}
        if twice:
            errors.append(f"targets built by several parts: {twice}")
        missing = sorted(reference - set(owners))
        unknown = sorted(set(owners) - reference)
        if missing:
            errors.append(f"targets built by no part: {missing}")
        if unknown:
            errors.append(f"targets outside the configuration: {unknown}")
        ran = [c for m in manifests for c in m["ran"]]
        if len(ran) != len(set(ran)):
            errors.append("ctest cases run by several parts")
        python = [m["part"] for m in manifests if m["python_tests"] == "success"]
        if BINDINGS_TARGET in reference and len(python) != 1:
            errors.append(f"Python tests passed in parts {python}, expected exactly one")
        lines.append(f"total: {len(reference)} targets, each built once; {len(ran)} ctest cases, each run once")
    print("\n".join(lines))
    if os.environ.get("GITHUB_STEP_SUMMARY"):
        with open(os.environ["GITHUB_STEP_SUMMARY"], "a", encoding="utf-8") as f:
            f.write(f"### {a.platform}\n```\n" + "\n".join(lines + [f"error: {e}" for e in errors]) + "\n```\n")
    for e in errors:
        print(f"::error::{e}")
    sys.exit(1 if errors else 0)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("plan")
    p.add_argument("--build", required=True)
    p.add_argument("--part", type=int, required=True)
    p.add_argument("--parts", type=int, required=True)
    p.add_argument("--out", required=True)
    p = sub.add_parser("select")
    p.add_argument("--build", required=True)
    p.add_argument("--plan", required=True)
    p = sub.add_parser("manifest")
    p.add_argument("--plan", required=True)
    p.add_argument("--junit")
    p.add_argument("--status", required=True)
    p.add_argument("--python-tests", required=True)
    p.add_argument("--platform", required=True)
    p.add_argument("--out", required=True)
    p = sub.add_parser("check")
    p.add_argument("--manifests", required=True)
    p.add_argument("--parts", type=int, required=True)
    p.add_argument("--platform", required=True)
    p.add_argument("--needs-result", required=True)
    a = parser.parse_args()
    if a.command == "plan":
        if not 1 <= a.part <= a.parts:
            fail(f"part {a.part} of {a.parts}")
        cmd_plan(a)
    else:
        {"select": cmd_select, "manifest": cmd_manifest, "check": cmd_check}[a.command](a)


if __name__ == "__main__":
    main()
