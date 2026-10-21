"""Exploratory: the measurement behind decision 030's final-score line for the toy optimizer
comparison. Not a test; the line it prints is set before the compiled comparison runs.

  python test_scripts/fixedeq/explore/toy_calibration.py --alisim ~/iqtree3-baseline/iqtree3 \
      --out DIR --git git.exe

Writes the toy problem (nqc/toy.py) into the new directory DIR, then runs three oracle fits that
differ only by rounding: the plain fit; every start coordinate moved one unit in the last place;
and the LU route to nu of explore/thresholds.py in place of the pivoted QR. The spread is the
largest final-score difference from the plain fit, and the line is max(10 x spread, 1e-8 |lnL|).
Writes results.json in DIR, which nqc/nqcchecks.py reads, and prints a summary with provenance.
"""
import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse  # noqa: E402
import datetime  # noqa: E402
import json  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
from concurrent.futures import ProcessPoolExecutor  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
import scipy  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE.parent / "nqc"))
sys.path.insert(0, str(HERE.parent / "probes"))
from oracle import chart, nqc, readers, target  # noqa: E402
from probes import turtle_target  # noqa: E402
from regress import REPO, capture, sha256, write  # noqa: E402
from thresholds import build_alt  # noqa: E402
import toy  # noqa: E402

SCRIPT = Path(__file__).resolve().relative_to(REPO).as_posix()
FACTOR, FLOOR = 10.0, 1e-8  # decision 030


def lu_build(z, pi, ref):
    return build_alt(z, pi, ref)[:2]


def one_fit(args):
    label, out, pi, ref, z0, use_lu = args
    f = nqc.Objective(readers.read_phylip(Path(out) / "toy.phy"), readers.read_newick(Path(out) / "toy.nwk"),
                      pi, ref, build=lu_build if use_lu else chart.build)
    t0 = time.perf_counter()
    start = -f(np.array(z0))
    r = nqc.fit(f, np.array(z0))
    Q = f.build(r.x, pi, ref)[0]
    return label, {"start_lnl": start, "final_lnl": float(r.score), "iterations": r.result.iterations,
                   "stop": r.result.reason, "failed_searches": r.result.failed_searches,
                   "evaluations": f.evaluations, "wall_s": round(time.perf_counter() - t0, 1),
                   "x": [float(v) for v in r.x], "Q": Q.tolist()}


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--alisim", required=True, help="the frozen unmodified binary")
    ap.add_argument("--out", required=True, help="new directory outside the repository")
    ap.add_argument("--git", default="git", help="git executable (git.exe under WSL)")
    a = ap.parse_args()
    out = Path(a.out).expanduser().resolve()
    if out.exists() or out.is_relative_to(REPO):
        sys.exit(f"{out}: needs a new directory outside the repository")
    out.mkdir(parents=True)
    alisim = Path(a.alisim).expanduser().resolve()
    prov = {"script": SCRIPT, "command": sys.argv,
            "started_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
            "git_commit": capture(a.git, "rev-parse", "HEAD", cwd=REPO),
            "git_clean": capture(a.git, "status", "--porcelain", cwd=REPO) == "",
            "alisim": {"path": str(alisim), "sha256": sha256(alisim)},
            "numpy": np.__version__, "scipy": scipy.__version__, "python": sys.version,
            "threads": {v: os.environ[v] for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")},
            "seed": toy.SEED, "factor": FACTOR, "floor": FLOOR}
    target_text = turtle_target()
    pi = target.parse_target(target_text)[0]
    prov["toy"] = toy.make(out, alisim, pi)
    lg, ref = nqc.lg_seed(pi)
    z0 = chart.inverse(lg, ref)
    jobs = [("plain", str(out), pi, ref, z0, False),
            ("start_ulp", str(out), pi, ref, np.nextafter(z0, np.inf), False),
            ("lu_nu", str(out), pi, ref, z0, True)]
    with ProcessPoolExecutor(max_workers=len(jobs)) as pool:
        fits = dict(pool.map(one_fit, jobs))
    plain = fits["plain"]
    diffs = {k: abs(v["final_lnl"] - plain["final_lnl"]) for k, v in fits.items() if k != "plain"}
    qdiff = {k: float(np.max(np.abs(np.array(v["Q"]) - plain["Q"]) / np.abs(plain["Q"])))
             for k, v in fits.items() if k != "plain"}
    spread = max(diffs.values())
    line = max(FACTOR * spread, FLOOR * abs(plain["final_lnl"]))
    prov["finished_utc"] = datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")
    result = {"provenance": prov, "target_text": target_text, "references": [int(r) for r in ref],
              "start": [float(v) for v in z0], "fits": fits, "final_lnl_differences": diffs,
              "q_relative_differences": qdiff, "spread": spread, "line": line}
    write(out / "results.json", json.dumps(result, indent=1))
    print(f"script {SCRIPT}; commit {prov['git_commit']}; clean {prov['git_clean']}; {prov['started_utc']}; "
          f"seed {toy.SEED}; numpy {np.__version__}; scipy {scipy.__version__}")
    for k, v in fits.items():
        print(f"  {k}: start {v['start_lnl']!r}, final {v['final_lnl']!r}, {v['iterations']} iterations "
              f"({v['stop']}), {v['failed_searches']} failed line searches, {v['evaluations']} evaluations, "
              f"{v['wall_s']} s")
    for k in diffs:
        print(f"  {k} against plain: final lnL difference {diffs[k]:.3e}; largest relative Q difference {qdiff[k]:.3e}")
    print(f"spread {spread:.3e}; line max({FACTOR:g} x spread, {FLOOR:g} |lnL|) = {line:.3e}")
    print(f"wrote {out / 'results.json'}")


if __name__ == "__main__":
    main()
