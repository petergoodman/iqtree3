#!/usr/bin/env python3
"""Sanitizer build and run of IQ-TREE (decision 022; docs/agent/CODE_PLAN.md section 3.6).

  build  configure and build into a new directory with decision 022's flags; write sanitizer-build.json
  run    run the regression list, the differential cases and the probes with a sanitizer binary, then scan
  scan   group the sanitizer reports found in a run directory; write findings.json and findings.md

Findings in upstream code are recorded, not fixed.
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
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "regression"))
from regress import PACKAGES, REPO, capture, read, sha256, write  # noqa: E402

SCRIPT = Path(__file__).resolve().relative_to(REPO).as_posix()
DRIVERS = HERE.parent
# decision 022, in fallback order
VARIANTS = {
    "asan-ubsan": ("-fsanitize=address,undefined", []),
    "asan-ubsan-nocmaple": ("-fsanitize=address,undefined", ["-DUSE_CMAPLE=OFF"]),
    "ubsan": ("-fsanitize=undefined", []),
}
RUN_ENV = {"ASAN_OPTIONS": "detect_leaks=0:allow_user_segv_handler=0", "UBSAN_OPTIONS": "print_stacktrace=1"}

ASAN = re.compile(r"^==\d+==ERROR: AddressSanitizer: (\S+)")
FRAME = re.compile(r"^\s+#\d+ 0x[0-9a-f]+ in (.+?) (/\S+?):(\d+)(?::\d+)?$")
UBSAN = re.compile(r"^(\S+?):(\d+):(\d+): runtime error: (.*)$")
CRASH = re.compile(r"CRASHES WITH SIGNAL (\S+)")


def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")


def build_env():
    """The environment without conda, whose Boost, OpenMP, zlib and libstdc++ CMake could pick up (CLAUDE.md)."""
    env = {k: v for k, v in os.environ.items() if not k.startswith("CONDA") and k != "CMAKE_PREFIX_PATH"}
    env["PATH"] = ":".join(p for p in env.get("PATH", "").split(":") if "conda" not in p)
    return env


def cmd_build(a):
    build = Path(a.build_dir).expanduser()
    if build.exists():
        sys.exit(f"{build} exists; a sanitizer build goes into a new directory")
    san, extra = VARIANTS[a.variant]
    flags = f"{san} -fno-omit-frame-pointer"
    configure = ["cmake", "-S", str(REPO), "-B", str(build), "-DCMAKE_BUILD_TYPE=Mem",
                 "-DCMAKE_C_COMPILER=clang", "-DCMAKE_CXX_COMPILER=clang++",
                 f"-DCMAKE_C_FLAGS={flags}", f"-DCMAKE_CXX_FLAGS={flags}", f"-DCMAKE_EXE_LINKER_FLAGS={san}", *extra]
    compile_ = ["cmake", "--build", str(build), "-j", str(a.jobs)]
    env = build_env()
    rec = {"script": SCRIPT, "command": sys.argv, "variant": a.variant, "started_utc": now(),
           "git_commit": capture(a.git, "rev-parse", "HEAD", cwd=REPO),
           "git_clean": capture(a.git, "status", "--porcelain", cwd=REPO) == "",
           "platform": platform.platform(), "cmake": capture("cmake", "--version").splitlines()[0],
           "compiler": capture("clang++", "--version").splitlines()[0],
           "packages": {p: capture("dpkg-query", "-W", "-f", "${Version}", p) for p in PACKAGES},
           "path": env["PATH"], "steps": []}
    build.mkdir(parents=True)
    log = build / "sanitizer-build.log"
    with open(log, "w", encoding="utf-8") as fh:
        for argv in (configure, compile_):
            t0 = time.perf_counter()
            fh.write(f"$ {' '.join(argv)}\n")
            fh.flush()
            code = subprocess.run(argv, env=env, stdout=fh, stderr=subprocess.STDOUT).returncode
            rec["steps"].append({"argv": argv, "exit_code": code, "wall_s": round(time.perf_counter() - t0, 1)})
            print(f"{argv[1]}: exit {code}", flush=True)
            if code != 0:
                break
    binary = build / "iqtree3"
    rec["binary"] = {"path": str(binary), "sha256": sha256(binary)} if binary.exists() else None
    rec["finished_utc"] = now()
    write(build / "sanitizer-build.json", json.dumps(rec, indent=1))
    tail = read(log).splitlines()[-15:]
    print("\n".join(tail))
    sys.exit(0 if rec["binary"] and all(s["exit_code"] == 0 for s in rec["steps"]) else 1)


def scan_text(text):
    """Sanitizer and crash reports in one output file: [(kind, location, message)]."""
    found, lines = [], text.splitlines()
    for i, line in enumerate(lines):
        m = ASAN.match(line)
        if m:
            loc = None
            for frame in lines[i + 1:i + 60]:
                f = FRAME.match(frame)
                if f and not f.group(2).startswith(("/usr/", "../")) and "compiler-rt" not in f.group(2):
                    loc = f"{f.group(2)}:{f.group(3)} ({f.group(1)})"
                    break
            found.append(("AddressSanitizer", loc, m.group(1)))
            continue
        m = UBSAN.match(line)
        if m:
            found.append(("UndefinedBehaviorSanitizer", f"{m.group(1)}:{m.group(2)}", re.sub(r"-?\d[\d.e+-]*", "#", m.group(4))))
            continue
        m = CRASH.search(line)
        if m:
            found.append(("IQ-TREE crash", None, f"signal {m.group(1)}"))
    return found


def cmd_scan(a):
    root = Path(a.dir).expanduser()
    groups, scanned = {}, 0
    for p in sorted(root.rglob("*")):
        if not p.is_file() or not (p.name == "stdout.txt" or p.suffix == ".log"):
            continue
        scanned += 1
        for kind, loc, msg in scan_text(read(p)):
            g = groups.setdefault((kind, loc, msg), {"kind": kind, "location": loc, "message": msg, "count": 0, "files": []})
            g["count"] += 1
            rel = str(p.relative_to(root))
            if rel not in g["files"]:
                g["files"].append(rel)
    meta = root / "sanitizer.json"
    result = {"script": SCRIPT, "command": sys.argv, "scanned_utc": now(), "files_scanned": scanned,
              "run": json.loads(read(meta)) if meta.exists() else None, "groups": list(groups.values())}
    write(root / "findings.json", json.dumps(result, indent=1))
    lines = ["# Sanitizer findings", "", f"Scanned {scanned} output files under `{root}` on {result['scanned_utc']}.", ""]
    if result["run"]:
        r = result["run"]
        lines += [f"Binary SHA-256 `{r['binary']['sha256']}`; options {r['options']}.", "",
                  "| Step | Exit | Wall time (s) |", "|---|---|---|"]
        lines += [f"| {s['name']} | {s['exit_code']} | {s['wall_s']} |" for s in r["steps"]] + [""]
    if groups:
        lines += ["| Kind | Location | Message | Count | Files |", "|---|---|---|---|---|"]
        lines += [f"| {g['kind']} | {g['location']} | {g['message']} | {g['count']} | {', '.join(g['files'][:5])}"
                  f"{' ...' if len(g['files']) > 5 else ''} |" for g in groups.values()]
    else:
        lines.append("No sanitizer reports and no IQ-TREE crash messages were found.")
    write(root / "findings.md", "\n".join(lines) + "\n")
    print(f"{len(groups)} distinct findings in {scanned} files; wrote {root / 'findings.md'}")


def cmd_run(a):
    out = Path(a.out).expanduser()
    if out.exists():
        sys.exit(f"{out} exists; every sanitizer run goes into a new directory")
    binary = Path(a.binary).expanduser().resolve()
    out.mkdir(parents=True)
    env = {**os.environ, **RUN_ENV}
    common = ["--git", a.git, "--timeout", str(a.timeout)] + (["--allow-dirty"] if a.allow_dirty else [])
    py = sys.executable
    steps = [
        ("regression", [py, str(DRIVERS / "regression" / "regress.py"), "run", "--binary", str(binary),
                        "--out", str(out / "regression"), "--build-dir", a.build_dir, *common]),
        ("differential", [py, str(DRIVERS / "differential" / "differential.py"), "--binary", str(binary),
                          "--out", str(out / "differential"), *common]),
        ("probes", [py, str(DRIVERS / "probes" / "probes.py"), "--binary", str(binary),
                    "--out", str(out / "probes"), *common]),
    ]
    meta = {"script": SCRIPT, "command": sys.argv, "started_utc": now(),
            "git_commit": capture(a.git, "rev-parse", "HEAD", cwd=REPO),
            "git_clean": capture(a.git, "status", "--porcelain", cwd=REPO) == "",
            "binary": {"path": str(binary), "sha256": sha256(binary)}, "build_dir": a.build_dir,
            "options": RUN_ENV, "steps": []}
    for name, argv in steps:
        t0 = time.perf_counter()
        print(f"{name}: started {datetime.datetime.now():%H:%M:%S}", flush=True)
        with open(out / f"{name}.log", "w", encoding="utf-8") as fh:
            code = subprocess.run(argv, env=env, cwd=REPO, stdout=fh, stderr=subprocess.STDOUT).returncode
        meta["steps"].append({"name": name, "argv": argv, "exit_code": code, "wall_s": round(time.perf_counter() - t0, 1)})
        print(f"{name}: exit {code}", flush=True)
    meta["finished_utc"] = now()
    write(out / "sanitizer.json", json.dumps(meta, indent=1))
    cmd_scan(argparse.Namespace(dir=str(out)))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    b = sub.add_parser("build")
    b.add_argument("--build-dir", required=True)
    b.add_argument("--variant", choices=list(VARIANTS), default="asan-ubsan")
    b.add_argument("--jobs", type=int, default=2)
    b.add_argument("--git", default="git")
    r = sub.add_parser("run")
    r.add_argument("--binary", required=True)
    r.add_argument("--build-dir", required=True)
    r.add_argument("--out", required=True)
    r.add_argument("--timeout", type=float, default=10800, help="seconds per IQ-TREE run")
    r.add_argument("--git", default="git")
    r.add_argument("--allow-dirty", action="store_true")
    s = sub.add_parser("scan")
    s.add_argument("dir")
    a = parser.parse_args()
    {"build": cmd_build, "run": cmd_run, "scan": cmd_scan}[a.command](a)


if __name__ == "__main__":
    main()
