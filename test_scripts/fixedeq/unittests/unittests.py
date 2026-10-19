#!/usr/bin/env python3
"""Build and run the S1 unit tests (unittest/) in one variant and write the residual table
(docs/agent/CODE_PLAN.md, sections 3.2 and 3.6).

  python test_scripts/fixedeq/unittests/unittests.py run --variant release|asan-ubsan \
      --build-dir DIR --out DIR [--git git.exe] [--allow-dirty]

--out must be a new directory. It receives the cmake and test logs, the googletest XML,
residuals.tsv (one row per fixture) and provenance.json. Exit status 1 if the configure, the build
or any test fails, if the XML is missing or unreadable, or if a fixture has no row.
"""
import argparse
import datetime
import json
import os
import platform
import re
import subprocess
import sys
import time
import xml.etree.ElementTree as ET
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "regression"))
sys.path.insert(0, str(HERE.parent / "sanitizer"))
from regress import REPO, capture, read, write  # noqa: E402
from sanitize import RUN_ENV, VARIANTS, build_env  # noqa: E402

SCRIPT = Path(__file__).resolve().relative_to(REPO).as_posix()
SOURCE = REPO / "unittest"
FIXTURES = SOURCE / "fixtures"
BINARY = "fixedeqchart_tests"
COLUMNS = ["fixture", "status", "truth", "n", "kappa", "q_limit", "build_error", "roundtrip_error", "seed_error",
           "stationarity", "row_sums", "mean_rate_error", "nu_residual", "exit_rate_mismatch", "min_rate"]
STANDARD = {"name", "status", "result", "time", "timestamp", "classname", "value_param", "type_param", "file", "line"}
MAXIMA = ["build_error", "roundtrip_error", "seed_error", "stationarity", "row_sums", "mean_rate_error", "nu_residual",
          "exit_rate_mismatch"]


def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")


def configure_args(variant, build):
    """Release takes IQ-TREE's Release flags from unittest/CMakeLists.txt. A sanitizer variant takes
    decision 022's flags at -g -O1, the flags of IQ-TREE's Mem build type, which this project lacks."""
    base = ["cmake", "-S", str(SOURCE), "-B", str(build), "-DCMAKE_CXX_COMPILER=clang++"]
    if variant == "release":
        return base + ["-DCMAKE_BUILD_TYPE=Release"]
    san = VARIANTS[variant][0]
    return base + ["-DCMAKE_BUILD_TYPE=", f"-DCMAKE_CXX_FLAGS=-g -O1 {san} -fno-omit-frame-pointer",
                   f"-DCMAKE_EXE_LINKER_FLAGS={san}"]


def fixture_names():
    return sorted(p.stem for p in FIXTURES.glob("*.txt") if p.name != "targets.txt")


def results(xml_path, expected):
    """(problems, rows) from googletest's XML: one row per test case that recorded a fixture.
    Properties are read both as <property> children and as test-case attributes."""
    try:
        root = ET.parse(xml_path).getroot()
    except (OSError, ET.ParseError) as e:
        return [f"no readable googletest XML ({e})"], []
    problems = []
    failed = int(root.get("failures", 0)) + int(root.get("errors", 0))
    if failed:
        problems.append(f"{failed} of {root.get('tests')} tests failed")
    rows = []
    for tc in root.iter("testcase"):
        props = {k: v for k, v in tc.attrib.items() if k not in STANDARD}
        props.update({p.get("name"): p.get("value") for p in tc.iter("property")})
        if "fixture" in props:
            props["status"] = "failed" if tc.find("failure") is not None else "passed"
            rows.append(props)
    missing = sorted(set(expected) - {r["fixture"] for r in rows})
    if missing:
        problems.append(f"no row for {len(missing)} fixtures: {', '.join(missing)}")
    return problems, sorted(rows, key=lambda r: r["fixture"])


def maxima(rows):
    out = {}
    for k in MAXIMA:
        vals = [float(r[k]) for r in rows if r.get(k)]
        out[k] = max(vals) if vals else None
    ratios = [float(r["build_error"]) / float(r["q_limit"]) for r in rows if r.get("build_error")]
    out["build_error_over_limit"] = max(ratios) if ratios else None
    return out


def eigen_version():
    for path in ("/usr/include/eigen3/Eigen/src/Core/util/Macros.h",):
        try:
            text = read(path)
        except OSError:
            continue
        parts = [re.search(rf"#define EIGEN_{k}_VERSION (\d+)", text) for k in ("WORLD", "MAJOR", "MINOR")]
        if all(parts):
            return ".".join(m.group(1) for m in parts)
    return "unavailable"


def cmd_run(a):
    prov = {"script": SCRIPT, "command": sys.argv, "variant": a.variant, "started_utc": now(),
            "git_commit": capture(a.git, "rev-parse", "HEAD", cwd=REPO),
            "git_clean": capture(a.git, "status", "--porcelain", cwd=REPO) == ""}
    if not prov["git_clean"] and not a.allow_dirty:
        sys.exit(f"the repository is not clean ({prov['git_commit']}); commit first or pass --allow-dirty")
    out, build = Path(a.out).expanduser(), Path(a.build_dir).expanduser()
    if out.exists():
        sys.exit(f"{out} exists; results go into a new directory")
    out.mkdir(parents=True)
    cmake_text = read(SOURCE / "CMakeLists.txt")
    prov.update({
        "build_dir": str(build), "build_dir_was_new": not build.exists(),
        "platform": platform.platform(), "cmake": capture("cmake", "--version").splitlines()[0],
        "compiler": capture("clang++", "--version").splitlines()[0], "eigen": eigen_version(),
        "googletest": {"url": re.search(r"URL (\S+)", cmake_text).group(1),
                       "sha256": re.search(r"URL_HASH SHA256=(\S+)", cmake_text).group(1)},
        "run_env": RUN_ENV if a.variant != "release" else {}, "steps": []})
    env = build_env()
    run_env = {**os.environ, **prov["run_env"]}
    xml = out / "gtest.xml"
    steps = [("configure", configure_args(a.variant, build), env),
             ("build", ["cmake", "--build", str(build), "-j", str(a.jobs)], env),
             ("test", [str(build / BINARY), f"--gtest_output=xml:{xml}"], run_env)]
    for name, argv, step_env in steps:
        t0 = time.perf_counter()
        with open(out / f"{name}.log", "w", encoding="utf-8") as fh:
            fh.write(f"$ {' '.join(argv)}\n")
            fh.flush()
            code = subprocess.run(argv, env=step_env, stdout=fh, stderr=subprocess.STDOUT).returncode
        prov["steps"].append({"step": name, "argv": argv, "exit_code": code, "wall_s": round(time.perf_counter() - t0, 1)})
        print(f"{name}: exit {code}", flush=True)
        if code != 0 and name != "test":
            break
    problems, rows = results(xml, fixture_names())
    problems = [f"{s['step']} exited {s['exit_code']}" for s in prov["steps"] if s["exit_code"] != 0] + problems
    if len(prov["steps"]) < len(steps):
        problems.append("the tests did not run")
    write(out / "residuals.tsv", "\n".join(["\t".join(COLUMNS)] + ["\t".join(r.get(c, "") for c in COLUMNS)
                                                                    for r in rows]) + "\n")
    prov.update({"rows": len(rows), "maxima": maxima(rows), "problems": problems, "finished_utc": now()})
    write(out / "provenance.json", json.dumps(prov, indent=1))
    print(f"{len(rows)} fixture rows; maxima: " + ", ".join(f"{k} {v:.2e}" for k, v in prov["maxima"].items()
                                                            if v is not None))
    for p in problems:
        print(f"PROBLEM: {p}")
    sys.exit(1 if problems else 0)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run", help="configure, build and run one variant")
    r.add_argument("--variant", required=True, choices=["release", "asan-ubsan", "ubsan"])
    r.add_argument("--build-dir", required=True)
    r.add_argument("--out", required=True)
    r.add_argument("--jobs", type=int, default=2)
    r.add_argument("--git", default="git", help="git executable (git.exe under WSL)")
    r.add_argument("--allow-dirty", action="store_true")
    a = ap.parse_args()
    cmd_run(a)


if __name__ == "__main__":
    main()
