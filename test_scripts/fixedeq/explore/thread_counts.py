"""Exploratory: how upstream's joint NONREV fit and fixed-point evaluations behave across thread
counts, the measurements behind decision 038 (PLAN.md risk 19, retired). Not a test; nothing here is
a pass line.

  python test_scripts/fixedeq/explore/thread_counts.py run --binary ~/iqtree3-baseline/iqtree3 \
      --runs ~/iqtree3-runs --out ~/iqtree3-runs/thread-counts/NEW --git git.exe
  python test_scripts/fixedeq/explore/thread_counts.py summarize DIR [DIR ...] --git git.exe > OUT.txt

`run` checks the binary against its SHA-256 file and runs, on turtle with baseline run 8's tree (-p)
or run 9's trees (-S) and run 10's copy of run 8's Q block, each with -seed 1 -pre run:
  pL2   -p --model-joint NONREV -te                    -T 1 x2, -T 4 x5
  SL2   -S --model-joint NONREV -te                    -T 1 x2, -T 4 x5
  SL1   -S --model-joint NONREV -te -blfix             -T 1 x2, -T 4 x5
  evS   -S -m Q -te -blfix --show-lh -vvv              -T 1 x3, -T 4 x10
  evP   -p -m Q -te -blfix --show-lh -vvv              -T 1 x3, -T 4 x10
  fitS1 -S --model-joint NONREV -te -blfix -vvv        -T 1 x2, -T 4 x3
`summarize` reads such directories (the 2026-10-22 records, made by one-off shell scripts running
these same commands, split them over two) and prints, per
setup: the printed log-likelihood of each run, the checkpoint values that differ between runs, and
the 17-digit initial log-likelihoods that -vvv prints (model/modelfactory.cpp, one per partition per
round under -S).
"""
import argparse
import datetime
import hashlib
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / "test_scripts/fixedeq"))
sys.path.insert(0, str(REPO / "test_scripts/fixedeq/probes"))
from probes import read_checkpoint  # noqa: E402

DATA = REPO / "test_scripts/test_data"
SETUPS = {  # name: (partition option, tree file, extra arguments, repeats at -T 1, repeats at -T 4)
    "pL2": ("-p", "p.tree", ["--model-joint", "NONREV"], 2, 5),
    "SL2": ("-S", "S.tree", ["--model-joint", "NONREV"], 2, 5),
    "SL1": ("-S", "S.tree", ["--model-joint", "NONREV", "-blfix"], 2, 5),
    "evS": ("-S", "S.tree", ["-m", "Q.txt", "-blfix", "--show-lh", "-vvv"], 3, 10),
    "evP": ("-p", "p.tree", ["-m", "Q.txt", "-blfix", "--show-lh", "-vvv"], 3, 10),
    "fitS1": ("-S", "S.tree", ["--model-joint", "NONREV", "-blfix", "-vvv"], 2, 3),
}
SKIP = re.compile(r"time|command|StopRule|seed")
NUM = re.compile(r"-?\d+\.?\d*(?:e[-+]?\d+)?")
INIT = re.compile(r"Initial log-likelihood: (-?\d+\.\d+)")
FINAL = re.compile(r"^Log-likelihood of the tree: (\S+)", re.M)


def git(exe, *args):
    try:
        return subprocess.run([exe, *args], cwd=REPO, capture_output=True, text=True, check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError) as e:
        return f"unavailable ({e})"


def provenance(a, extra):
    print(f"script test_scripts/fixedeq/explore/thread_counts.py; commit {git(a.git, 'rev-parse', 'HEAD')}; "
          f"clean {git(a.git, 'status', '--porcelain') == ''}; "
          f"{datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds')}; {extra}")


def run(a):
    binary, runs, out = Path(a.binary).expanduser(), Path(a.runs).expanduser(), Path(a.out).expanduser()
    want = (binary.parent / (binary.name + ".sha256")).read_text().split()[0]
    have = hashlib.sha256(binary.read_bytes()).hexdigest()
    if have != want:
        sys.exit(f"SHA-256 mismatch for {binary}: {have}")
    out.mkdir(parents=True)
    provenance(a, f"binary {binary} sha256 {have}; out {out}")
    base = runs / "baseline/rep1"
    shutil.copy(base / "run08/run08.treefile", out / "p.tree")
    shutil.copy(base / "run09/run09.treefile", out / "S.tree")
    shutil.copy(base / "run10/run08.Q.txt", out / "Q.txt")
    for name, (part, tree, extra, n1, n4) in SETUPS.items():
        for threads, n in ((1, n1), (4, n4)):
            for rep in range(1, n + 1):
                d = out / name / f"T{threads}-rep{rep}"
                d.mkdir(parents=True)
                args = [str(binary), "-s", str(DATA / "turtle_aa.fasta"), part, str(DATA / "turtle_aa.nex"),
                        "-te", str(out / tree)] + [str(out / x) if x == "Q.txt" else x for x in extra] + \
                       ["-seed", "1", "-T", str(threads), "-pre", "run"]
                t0 = time.time()
                with open(d / "stdout.txt", "w") as fh:
                    code = subprocess.run(args, cwd=d, stdout=fh, stderr=subprocess.STDOUT).returncode
                print(f"{name} T{threads} rep{rep}: exit {code}, wall {time.time() - t0:.1f} s", flush=True)


def numbers(text):
    return [float(x) for x in NUM.findall(text or "")]


def summarize_setup(name, dirs):
    print(f"\n== {name}")
    final = {d.name: (FINAL.search((d / "run.iqtree").read_text()) or [None, "missing"])[1] for d in dirs}
    print("  printed log-likelihood: " + ", ".join(f"{k} {v}" for k, v in final.items()))
    t1 = {v for k, v in final.items() if k.startswith("T1")}
    t4 = [v for k, v in final.items() if k.startswith("T4")]
    print(f"  -T 1 distinct {len(t1)}; -T 4 distinct {len(set(t4))} of {len(t4)}, equal to -T 1: "
          f"{sum(v in t1 for v in t4)}; -T 4 range {max(map(float, t4)) - min(map(float, t4)):.4f}")

    ckp = {d.name: {k: v for k, v in read_checkpoint(d / "run.ckp.gz").items() if not SKIP.search(k)}
           for d in dirs if (d / "run.ckp.gz").exists()}
    ref = next(k for k in ckp if k.startswith("T1"))
    keys = sorted(set().union(*ckp.values()))
    differing = [k for k in keys if len({c.get(k) for c in ckp.values()}) > 1]
    print(f"  checkpoint values compared {len(keys)}, differing between runs {len(differing)}")
    for k in differing:
        r = numbers(ckp[ref].get(k))
        diffs = [max(abs(p - q) for p, q in zip(r, x)) for x in (numbers(c.get(k)) for c in ckp.values())
                 if len(x) == len(r) and r]
        print(f"    {k}: max |difference| from {ref} {max(diffs) if diffs else float('nan'):.3g}")

    init = {d.name: INIT.findall((d / "stdout.txt").read_text()) for d in dirs}
    if not any(init.values()):
        return
    print(f"  initial log-likelihoods printed per run: {sorted({len(v) for v in init.values()})}")
    for t in ("T1", "T4"):
        group = [v for k, v in init.items() if k.startswith(t)]
        print(f"    -{t[0]} {t[1:]}: {len(group)} runs, {len({tuple(v) for v in group})} distinct sequences")
    one = next(v for k, v in init.items() if k.startswith("T1"))
    four = [v for k, v in init.items() if k.startswith("T4")]
    for i in range(min([len(one)] + [len(v) for v in four] + [6])):
        vals = sorted({v[i] for v in four})
        rel = max(abs(float(x) - float(one[i])) / abs(float(one[i])) for x in vals)
        print(f"    value {i}: -T 1 {one[i]}; -T 4 distinct {len(vals)}; largest relative difference "
              f"from -T 1 {rel:.2e}")


def summarize(a):
    provenance(a, "dirs " + " ".join(a.dirs))
    for top in (Path(x).expanduser() for x in a.dirs):
        for name in SETUPS:
            dirs = sorted((top / name).glob("T*-rep*"))
            if dirs:
                summarize_setup(f"{name} ({top.name})", dirs)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--binary", required=True, help="the frozen unmodified binary, with NAME.sha256 beside it")
    r.add_argument("--runs", required=True, help="the run directory holding baseline/rep1")
    r.add_argument("--out", required=True, help="a new directory for the outputs")
    s = sub.add_parser("summarize")
    s.add_argument("dirs", nargs="+", help="directories written by run")
    for p in (r, s):
        p.add_argument("--git", default="git", help="git executable (git.exe under WSL)")
    a = ap.parse_args()
    run(a) if a.cmd == "run" else summarize(a)


if __name__ == "__main__":
    main()
