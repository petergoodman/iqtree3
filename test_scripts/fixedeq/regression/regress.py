#!/usr/bin/env python3
"""Regression driver for the fixedeq work (docs/agent/CODE_PLAN.md, section 3.4).

  freeze    copy an IQ-TREE binary into a directory and record its SHA-256
  run       run the baseline list with a binary into a new directory; write extract.json
  baseline  merge repeated extracts into baseline.json and baseline.md
  compare   check an extract.json against baseline.json; exit status 1 on any difference
"""
import argparse
import datetime
import hashlib
import json
import os
import platform
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
SCRIPT = Path(__file__).resolve().relative_to(REPO).as_posix()

AA = "@example/aa_example.phy"
DNA = "@example/example.phy"
TURTLE = ["-s", "@test_scripts/test_data/turtle_aa.fasta"]
NEX = "@test_scripts/test_data/turtle_aa.nex"
FREQ = ("0.08,0.06,0.04,0.05,0.02,0.04,0.07,0.07,0.02,0.05,"
        "0.10,0.06,0.02,0.04,0.05,0.07,0.05,0.01,0.03,0.07")

# name: (arguments, threads); "@" marks a path relative to the repository root
RUNS = {
    "run01": (["-s", AA, "-m", "LG+G4"], 1),
    "run02": (["-s", AA, "-m", "NONREV"], 1),
    "run03": (["-s", AA, "-m", "NONREV+F{" + FREQ + "}"], 1),
    "run04": (["-s", AA, "-m", "NQ.pfam"], 1),
    "run05": (["-s", AA, "-m", "GTR20"], 1),
    "run06": (["-s", DNA, "-m", "UNREST"], 1),
    "run07": (["-s", DNA, "-m", "12.12"], 1),
    "run08": (TURTLE + ["-p", NEX, "--model-joint", "NONREV"], 1),
    "run09": (TURTLE + ["-S", NEX, "--model-joint", "NONREV"], 1),
    # relative names, because IQ-TREE prints the Q file's path inside the compared model section
    "run10": (TURTLE + ["-p", NEX, "-m", "run08.Q.txt", "-te", "run08.treefile"], 1),
    "run11": (TURTLE + ["-p", NEX, "--model-joint", "NONREV"], 4),
}

# report parts that hold paths, the build date or clock times
SKIPPED = {"(header)", "ALISIM COMMAND", "TIME STAMP"}
Q_HEADER = "Full Q matrix and state frequencies (can be used as input for IQ-TREE):"
NUMBER = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?")
SCALARS = {
    "log_likelihood": r"^Log-likelihood of the tree: (\S+)",
    "free_parameters": r"^Number of free parameters \(#branches \+ #model parameters\): (\S+)",
    "tree_length": r"^Total tree length \(sum of branch lengths\): (\S+)",
    "wall_clock_s": r"^Total wall-clock time used: (\S+) seconds",
}
PACKAGES = ["libc6", "libstdc++6", "libgcc-s1", "zlib1g", "libomp-14-dev", "clang-14", "lld-14"]


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read(path):
    return Path(path).read_text(encoding="utf-8", errors="replace")


def write(path, text):
    Path(path).write_text(text, encoding="utf-8", newline="\n")


def resolve(args):
    return [str(REPO / a[1:]) if a.startswith("@") else a for a in args]


def split_sections(text):
    """Map each report section title to its text; the part before the first title is '(header)'."""
    lines = text.splitlines()
    heads = [i for i in range(len(lines) - 1)
             if lines[i].strip() and set(lines[i + 1]) == {"-"} and len(lines[i + 1]) == len(lines[i])]
    sections = {"(header)": "\n".join(lines[:heads[0] if heads else len(lines)])}
    for k, i in enumerate(heads):
        end = heads[k + 1] if k + 1 < len(heads) else len(lines)
        key, n = lines[i], 2
        while key in sections:
            key, n = f"{lines[i]} #{n}", n + 1
        sections[key] = "\n".join(lines[i:end])
    return sections


def q_block(report_text):
    """The 20 Q rows and the frequency row that IQ-TREE prints for re-import with -m FILE."""
    lines = report_text.splitlines()
    start = next((i for i, line in enumerate(lines) if line.strip() == Q_HEADER), None)
    if start is None:
        raise ValueError("no full Q matrix block in the report")
    rows = []
    for line in lines[start + 1:]:
        if line.strip():
            rows.append(line)
        elif rows:
            break
    if len(rows) != 21 or any(len(r.split()) != 20 for r in rows):
        raise ValueError(f"Q block has {len(rows)} rows; expected 21 rows of 20 numbers")
    for t in " ".join(rows).split():
        float(t)
    return "\n".join(rows) + "\n"


def split_numbers(text):
    return NUMBER.sub("#", text), NUMBER.findall(text)


def capture(*argv, cwd=None):
    try:
        return subprocess.run(argv, cwd=cwd, capture_output=True, text=True, check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError) as e:
        return f"unavailable ({e})"


def provenance(git, binary=None, build_dir=None):
    debian = Path("/etc/debian_version")
    prov = {
        "script": SCRIPT,
        "command": sys.argv,
        "started_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
        "git_commit": capture(git, "rev-parse", "HEAD", cwd=REPO),
        "git_clean": capture(git, "status", "--porcelain", cwd=REPO) == "",
        "python": sys.version,
        "platform": platform.platform(),
        "debian": debian.read_text().strip() if debian.exists() else None,
        "cpu_count": os.cpu_count(),
        "packages": {p: capture("dpkg-query", "-W", "-f", "${Version}", p) for p in PACKAGES},
        "environment": {k: v for k, v in os.environ.items()
                        if k == "LANG" or k.startswith(("LC_", "OMP_", "KMP_"))},
        "inputs": {a[1:]: sha256(REPO / a[1:])
                   for args, _ in RUNS.values() for a in args if a.startswith("@")},
    }
    if binary:
        prov["binary"] = {"path": str(binary), "sha256": sha256(binary)}
    if build_dir:
        b = Path(build_dir).expanduser()
        cache = dict(line.split("=", 1) for line in read(b / "CMakeCache.txt").splitlines()
                     if "=" in line and not line.startswith(("#", "//")))
        compiler = cache.get("CMAKE_CXX_COMPILER:STRING") or cache.get("CMAKE_CXX_COMPILER:FILEPATH")
        target = b / "CMakeFiles" / "iqtree3.dir"
        prov["build"] = {
            "dir": str(b),
            "build_type": cache.get("CMAKE_BUILD_TYPE:STRING"),
            "compiler": compiler,
            "compiler_version": capture(compiler, "--version").splitlines()[0] if compiler else None,
            "cxx_flags": [l for l in read(target / "flags.make").splitlines() if l.startswith("CXX_FLAGS")],
            "link_flags": sorted({t for t in read(target / "link.txt").split()
                                  if t.startswith("-") and t != "-o"}),
        }
    return prov


def extract_run(rundir, name, exit_code):
    """Everything a later run is compared on, plus summary scalars and diagnostics."""
    rec = {"exit_code": exit_code, "items": {}, "scalars": {}}
    report = rundir / f"{name}.iqtree"
    if report.exists():
        text = read(report)
        rec["items"] = {k: v for k, v in split_sections(text).items() if k not in SKIPPED}
        rec["scalars"] = {k: (m.group(1) if (m := re.search(p, text, re.M)) else None)
                          for k, p in SCALARS.items()}
    tree = rundir / f"{name}.treefile"
    if tree.exists():
        rec["items"]["(treefile)"] = read(tree)
    rec["items"]["(exit code)"] = str(exit_code)
    log = rundir / f"{name}.log"
    rec["log_header"] = ([l for l in read(log).splitlines()
                          if l.startswith(("IQ-TREE version", "Host:", "Kernel:"))]
                         if log.exists() else [])
    if exit_code != 0:
        out = rundir / "stdout.txt"
        rec["stdout_tail"] = read(out).splitlines()[-30:] if out.exists() else []
    return rec


def make_baseline(extracts):
    """Merge repeated extracts. An item identical in every repeat is compared exactly; one whose
    text differs only in its numbers gets a per-number spread (max minus min); anything else is
    'unstable' and cannot be compared."""
    runs = {}
    for name, (_, threads) in RUNS.items():
        reps = [(label, ex["runs"][name]) for label, ex in extracts if name in ex["runs"]]
        if not reps:
            continue
        items = {}
        for key in dict.fromkeys(k for _, r in reps for k in r["items"]):
            texts = [r["items"].get(key) for _, r in reps]
            if all(t == texts[0] for t in texts):
                items[key] = {"mode": "exact", "text": texts[0]}
                continue
            if None not in texts:
                parts = [split_numbers(t) for t in texts]
                if all(p[0] == parts[0][0] and len(p[1]) == len(parts[0][1]) for p in parts):
                    columns = zip(*([float(x) for x in p[1]] for p in parts))
                    items[key] = {"mode": "spread", "skeleton": parts[0][0], "numbers": parts[0][1],
                                  "spread": [max(c) - min(c) for c in columns]}
                    continue
            items[key] = {"mode": "unstable", "texts": texts}
        runs[name] = {
            "arguments": reps[0][1].get("arguments"),
            "threads": threads,
            "repeats": [label for label, _ in reps],
            "exit_codes": [r.get("exit_code") for _, r in reps],
            "scalars": [r.get("scalars", {}) for _, r in reps],
            "driver_wall_s": [r.get("driver_wall_s") for _, r in reps],
            "log_header": reps[0][1].get("log_header", []),
            "items": items,
        }
    return runs


def first_difference(a, b):
    for i, (x, y) in enumerate(zip(a.splitlines(), b.splitlines()), 1):
        if x != y:
            return f"line {i}: {x.strip()[:70]!r} became {y.strip()[:70]!r}"
    return "one text is longer"


def compare_run(base, new):
    """The differences between one run's baseline entry and a new extract record."""
    problems = []
    new_items = new.get("items", {})
    for key in new_items:
        if key not in base["items"]:
            problems.append(f"{key}: not in the baseline")
    for key, b in base["items"].items():
        text = new_items.get(key)
        if text is None:
            problems.append(f"{key}: missing")
        elif b["mode"] == "exact":
            if text != b["text"]:
                problems.append(f"{key}: {first_difference(b['text'], text)}")
        elif b["mode"] == "spread":
            skeleton, numbers = split_numbers(text)
            if skeleton != b["skeleton"] or len(numbers) != len(b["numbers"]):
                problems.append(f"{key}: text differs beyond its numbers")
                continue
            outside = [(x, y, s) for x, y, s in zip(numbers, b["numbers"], b["spread"])
                       if abs(float(x) - float(y)) > s]
            if outside:
                x, y, s = outside[0]
                problems.append(f"{key}: {len(outside)} numbers outside their spread, "
                                f"first {x} against {y} (spread {s:.3g})")
        else:
            problems.append(f"{key}: unstable in the baseline, cannot be compared")
    return problems


def cmd_freeze(a):
    src = Path(a.binary).expanduser().resolve()
    dest = Path(a.dest).expanduser()
    target = dest / src.name
    if target.exists():
        sys.exit(f"{target} exists; a frozen binary is never overwritten")
    dest.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, target)
    digest = sha256(target)
    if digest != sha256(src):
        sys.exit("the copy differs from its source")
    write(dest / f"{src.name}.sha256", f"{digest}  {src.name}\n")
    print(f"{target}\n{digest}")


def cmd_run(a):
    out = Path(a.out).expanduser()
    if out.exists():
        sys.exit(f"{out} exists; every run set goes into a new directory")
    names = a.runs.split(",") if a.runs else list(RUNS)
    unknown = [n for n in names if n not in RUNS]
    if unknown:
        sys.exit(f"unknown runs: {unknown}")
    binary = Path(a.binary).expanduser().resolve()
    prov = provenance(a.git, binary, a.build_dir)
    if not prov["git_clean"] and not a.allow_dirty:
        sys.exit(f"the repository is not clean ({prov['git_commit']}); commit first or pass --allow-dirty")
    out.mkdir(parents=True)
    runs = {}
    for name in names:
        args, threads = RUNS[name]
        rundir = out / name
        rundir.mkdir()
        if name == "run10":
            src = out / "run08"
            try:
                write(rundir / "run08.Q.txt", q_block(read(src / "run08.iqtree")))
                shutil.copyfile(src / "run08.treefile", rundir / "run08.treefile")
            except (OSError, ValueError) as e:
                runs[name] = {"exit_code": "skipped", "scalars": {},
                              "items": {"(skipped)": f"run08 output unusable: {e}"}}
                print(f"{name}: skipped, run08 output unusable: {e}", flush=True)
                continue
        common = ["-seed", "1", "-T", str(threads), "--prefix", name]
        argv = [str(binary)] + resolve(args) + common
        print(f"{name}: started {datetime.datetime.now():%H:%M:%S}", flush=True)
        t0 = time.perf_counter()
        with open(rundir / "stdout.txt", "w", encoding="utf-8") as fh:
            try:
                code = subprocess.run(argv, cwd=rundir, stdout=fh, stderr=subprocess.STDOUT,
                                      timeout=a.timeout).returncode
            except subprocess.TimeoutExpired:
                code = "timeout"
        rec = extract_run(rundir, name, code)
        rec.update(arguments=args + common, driver_wall_s=round(time.perf_counter() - t0, 1))
        runs[name] = rec
        print(f"{name}: exit {code}, log-likelihood {rec['scalars'].get('log_likelihood')}, "
              f"{rec['driver_wall_s']} s", flush=True)
    prov["finished_utc"] = datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")
    write(out / "extract.json", json.dumps({"provenance": prov, "runs": runs}, indent=1, ensure_ascii=False))
    print(f"wrote {out / 'extract.json'}")


def markdown(baseline):
    p = baseline["provenance"]
    first = p["extracts"][0]
    build = first.get("build", {})
    lines = [
        "# S0 regression baseline",
        "",
        f"Generated by `{SCRIPT} baseline` on {p['baseline']['started_utc']}. Do not edit by hand;",
        "`baseline.json` holds the compared texts and the full provenance.",
        "",
        f"- Driver commit: `{first['git_commit']}` (clean working tree)",
        f"- Binary: `{first['binary']['path']}`, SHA-256 `{first['binary']['sha256']}`",
        f"- Built with: {build.get('compiler_version')}, {build.get('build_type')}",
        f"- System: {first['platform']}, Debian {first['debian']}",
        f"- IQ-TREE log header: {'; '.join(next(iter(baseline['runs'].values()))['log_header'])}",
        f"- Repeats: {', '.join(e['label'] + ' ' + e['started_utc'] for e in p['extracts'])}",
        "",
        "| Run | Arguments | Threads | Repeats | Log-likelihood | Free parameters | Tree length "
        "| Items varying across repeats | Wall time per repeat (s) |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for name, r in baseline["runs"].items():
        def values(key):
            vals = list(dict.fromkeys(str(s.get(key)) for s in r["scalars"]))
            return " / ".join(vals)
        shown = " ".join(Path(a[1:]).name if a.startswith("@") else a
                         for a in (r["arguments"] or [])[:-6])
        varying = [k for k, v in r["items"].items() if v["mode"] != "exact"]
        lines.append(f"| {name} | `{shown}` | {r['threads']} | {len(r['repeats'])} "
                     f"| {values('log_likelihood')} | {values('free_parameters')} "
                     f"| {values('tree_length')} | {', '.join(varying) or 'none'} "
                     f"| {', '.join(str(t) for t in r['driver_wall_s'])} |")
    return "\n".join(lines) + "\n"


def cmd_baseline(a):
    out = Path(a.out)
    if (out / "baseline.json").exists():
        sys.exit(f"{out / 'baseline.json'} exists; a baseline is never overwritten")
    extracts = [(Path(p).expanduser().parent.name, json.loads(read(Path(p).expanduser())))
                for p in a.extracts]
    provs = [ex["provenance"] for _, ex in extracts]
    for field in ("git_commit", "inputs"):
        if any(pr[field] != provs[0][field] for pr in provs):
            sys.exit(f"the extracts differ in {field}")
    if any(pr["binary"]["sha256"] != provs[0]["binary"]["sha256"] for pr in provs):
        sys.exit("the extracts come from different binaries")
    if not all(pr["git_clean"] for pr in provs):
        sys.exit("an extract was recorded from an unclean working tree")
    baseline = {
        "provenance": {
            "baseline": {"script": SCRIPT, "command": sys.argv,
                         "started_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
                         "git_commit": capture(a.git, "rev-parse", "HEAD", cwd=REPO),
                         "git_clean": capture(a.git, "status", "--porcelain", cwd=REPO) == ""},
            "extracts": [dict(pr, label=label) for (label, _), pr in zip(extracts, provs)],
        },
        "runs": make_baseline(extracts),
    }
    out.mkdir(parents=True, exist_ok=True)
    write(out / "baseline.json", json.dumps(baseline, indent=1, ensure_ascii=False))
    write(out / "baseline.md", markdown(baseline))
    print(f"wrote {out / 'baseline.json'} and baseline.md")


def cmd_compare(a):
    baseline = json.loads(read(Path(a.baseline).expanduser()))
    new = json.loads(read(Path(a.extract).expanduser()))
    failed = False
    for name, rec in new["runs"].items():
        if name not in baseline["runs"]:
            print(f"{name}: FAIL, not in the baseline")
            failed = True
            continue
        problems = compare_run(baseline["runs"][name], rec)
        print(f"{name}: {'FAIL' if problems else 'pass'}")
        for problem in problems:
            print(f"    {problem}")
        failed = failed or bool(problems)
    sys.exit(1 if failed else 0)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    f = sub.add_parser("freeze", help="copy a binary and record its SHA-256")
    f.add_argument("--binary", required=True)
    f.add_argument("--dest", required=True)
    r = sub.add_parser("run", help="run the baseline list into a new directory")
    r.add_argument("--binary", required=True)
    r.add_argument("--out", required=True)
    r.add_argument("--runs", help="comma-separated run names (default: all)")
    r.add_argument("--build-dir", help="CMake build directory, recorded in the provenance")
    r.add_argument("--git", default="git", help="git executable (git.exe under WSL)")
    r.add_argument("--allow-dirty", action="store_true")
    r.add_argument("--timeout", type=float, default=7200, help="seconds per run")
    b = sub.add_parser("baseline", help="merge extract.json files into a baseline")
    b.add_argument("extracts", nargs="+")
    b.add_argument("--out", required=True)
    b.add_argument("--git", default="git", help="git executable (git.exe under WSL)")
    c = sub.add_parser("compare", help="compare an extract.json with a baseline")
    c.add_argument("--baseline", required=True)
    c.add_argument("extract")
    a = parser.parse_args()
    {"freeze": cmd_freeze, "run": cmd_run, "baseline": cmd_baseline, "compare": cmd_compare}[a.command](a)


if __name__ == "__main__":
    main()
