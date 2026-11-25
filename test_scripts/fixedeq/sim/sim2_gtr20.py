"""SIM-2's seeded GTR20 arm (decision 046; CODE_PLAN.md 3.8, "SIM-2's seeded GTR20 arm"), a study of its own beside
the five shards of sim2.py, on round 1's data, binaries and settings. Python only; IQ-TREE and AliSim run as the
frozen and gate binaries given, refused unless they are round 1's.

  sim2_gtr20.py design [--smoke] [--out DIR] --git git.exe [--allow-dirty]
  sim2_gtr20.py run --out RUNDIR --shard K --cores N --frozen BIN --gate BIN [--deadline-hours 46] [--curve FILE]
                    [--smoke] [--skip-checks] --git git.exe [--allow-dirty]
  sim2_gtr20.py collect --out RUNDIR

design: cell D's 176 GTR20 tasks of sim2.design (160 D fits and Dr's 16 relaxed refits), ids suffixed -seeded,
each starting from --init-model seed.paml (every exchangeability 10, Tyr-Val P, then pi*: pipeline.gtr20_seed), with
a 40 h fit cap and estimates from round 1's default-seed medians times the laptop calibration's seeded-to-default
wall ratio (explore/gtr20_window.py); three shards of one wave each (R0's D fits, NQ's D fits, Dr); the first-wave
check and the charged-cost projection, cumulative with the five shards' 19,614, against 40,000. run, phase and
collect are sim2.py's (checks first, then phases longest first, each in its own process; a seeded fit adds
decision 046's checkpoint rows to the shard's report), with this study's design; collect writes
collect-gtr20-<stamp>.tar.gz. --smoke: the main smoke design's 5 GTR20 tasks (3 genes, 2 rounds), one shard.
"""
import argparse
import hashlib
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import sim2  # noqa: E402  (BLAS threads at 1; puts the drivers on sys.path)
import common  # noqa: E402
import pipeline as pl  # noqa: E402
from nqcchecks import check  # noqa: E402
from regress import REPO, sha256, write  # noqa: E402

SCRIPT = Path(__file__).resolve().relative_to(REPO).as_posix()
RUNS = "~/iqtree3-runs/sim2-gtr20"
STUDY = "gtr20-seeded"
# decision 046: p = 0.05, chosen by Peter on 2026-11-24 against the laptop rule's 0.02, whose headroom test both
# fits failed alike (each with I-V at the ceiling); 0.05 met the rule's lnL test and stopped by its gain test in 17
# rounds, 0.02 at the 20-round cap
P = 0.05
# the calibration's seeded-to-default wall ratio at p = 0.05 on R0 at 30 genes, 10,852.0 s over 10,701.5 s
# (~/iqtree3-runs/sim2-gtr20/20261124T192220Z-182a0dd7-window; 20261124-smoke-diagnosis/window_bounds.txt, section 3)
RATIO = 1.0141
FIT_CAP_H = 40.0    # a resource limit, not the 20-round setting: round 1's longest default-seed fit took 27.50 h
# round 1 and 2's default-seed GTR20 wall-time medians, hours, by (method, relaxed), 1 thread (relaxed 2):
# ~/iqtree3-runs/sim2/20261124-round2-read/window_checks.txt, (b)
EST_H = {("gtr20pi", False): 13.26, ("gtr20f", False): 16.01, ("gtr20pi", True): 13.67, ("gtr20f", True): 15.13}
LONGEST_H = 27.50
SHARD_CORES = {1: 84, 2: 84, 3: 34}   # R0's 80 D fits, NQ's 80 D fits, Dr's 16 at 2 threads, and a few spare
SMOKE_CORES = 12
DEADLINE_H = 46.0
PRIOR_CHARGED = 19614   # the five shards' charged core-hours, cores x wall (CHANGELOG.md, 2026-11-24)
CEILING = 40000         # Peter's ceiling (PLAN.md risk 26)
FROZEN_SHA256 = "9a72950bf61187f9923372f49badde3aa527864f57de795c5ad751d7c04fb501"
GATE_SHA256 = "25d78d6ee2391bd36f684a21bc6000959f5cb100e0e8d58b62020535956c7294"


def design(smoke=False):
    """The seeded arm's datasets, tasks and shards, deterministic: sim2.design(smoke)'s cell D and Dr GTR20
    tasks with ids suffixed -seeded and gtr20_start {"p": P}, otherwise unchanged except, outside the smoke
    design, the fit phase's cap and estimate."""
    main = sim2.design(smoke)
    tasks = []
    for t in main["tasks"]:
        if t["cell"] not in ("D", "Dr") or t["method"] == "nqc":
            continue
        t = {**t, "id": f"{t['id']}-seeded", "gtr20_start": {"p": P}}
        if not smoke:
            fit, rest = t["phases"][0], t["phases"][1:]
            t["phases"] = [{**fit, "cap_h": FIT_CAP_H, "est_h": EST_H[(t["method"], t["relaxed"])] * RATIO}] + rest
        tasks.append(t)
    need, stack = set(), [d for t in tasks for d in (t["dataset"], t["heldout"])]
    while stack:
        d = stack.pop()
        if d not in need:
            need.add(d)
            if main["datasets"][d]["source"]:
                stack.append(main["datasets"][d]["source"])

    def shard(t):
        if smoke:
            return 1
        return 3 if t["cell"] == "Dr" else 1 if main["datasets"][t["dataset"]]["truth"] == sim2.TRUTHS["R0"] else 2

    seed = pl.gtr20_seed(P, common.target_text())
    return {"study": STUDY, "smoke": smoke, "base": main["base"], "shards": 1 if smoke else len(SHARD_CORES), "p": P,
            "ratio": RATIO, "seed_sha256": hashlib.sha256(seed.encode()).hexdigest(),
            "datasets": {d: main["datasets"][d] for d in sorted(need)}, "tasks": tasks,
            "assignment": {t["id"]: shard(t) for t in tasks}}


def data_hours(phases):
    """The longest estimated chain of a shard's data phases (a simulation, then its masks)."""
    data = {pid: p for pid, p in phases.items() if pid.startswith("data/")}
    return max((p["est_h"] + sum(data[d]["est_h"] for d in p["deps"]) for p in data.values()), default=0.0)


def first_wave(des, shard, cores, deadline_h=DEADLINE_H):
    """Whether every fit of the shard starts in its first wave: once its data are made, one sim2.startable call
    with every core free (DEADLINE_FACTOR times each estimate within the deadline) starts them all."""
    phases = sim2.shard_phases(des, shard)
    data = [pid for pid in phases if pid.startswith("data/")]
    started = {p["id"] for p in sim2.startable(phases, dict.fromkeys(data, "done"), {}, cores, data_hours(phases),
                                               deadline_h, cores)}
    fits = [pid for pid, p in phases.items() if p["kind"] == "fit"]
    return all(f in started for f in fits), len(fits)


def projection(des, cores):
    """sim2.project on each shard at 1 and 1.5 times the estimates, with the cumulative total; and, reported, the
    charge if every fit took round 1's longest times RATIO, and if every job ran to the fit cap."""
    out = {}
    for name, mult in (("x1", 1.0), ("x1.5", 1.5)):
        shards = {k: sim2.project(sim2.shard_phases(des, k), c, DEADLINE_H, mult) for k, c in cores.items()}
        charged = sum(s["charged_core_h"] for s in shards.values())
        left = sum(s["left_core_h"] for s in shards.values())
        out[name] = {"shards": shards, "charged_core_h": charged, "round2_core_h": left,
                     "cumulative_core_h": PRIOR_CHARGED + charged + left}
    for name, fit_h in (("longest", LONGEST_H * RATIO), ("cap", FIT_CAP_H)):
        charged = sum(c * (data_hours(sim2.shard_phases(des, k)) + fit_h + sim2.METRICS_H) for k, c in cores.items())
        out[name] = {"fit_h": fit_h, "charged_core_h": charged, "cumulative_core_h": PRIOR_CHARGED + charged}
    return out


def binary_check(frozen, gate):
    """The binaries that differ from round 1's (decision 046: the seed is the only change), as (name, sha256)."""
    found = {"frozen": sha256(Path(frozen).expanduser().resolve()), "gate": sha256(Path(gate).expanduser().resolve())}
    want = {"frozen": FROZEN_SHA256, "gate": GATE_SHA256}
    return [(k, found[k]) for k in found if found[k] != want[k]]


def cmd_design(a):
    prov = common.provenance(a, SCRIPT, "design-smoke" if a.smoke else "design", "design", smoke=a.smoke)
    out = common.new_out(a.out, prov, runs=RUNS)
    des = design(a.smoke)
    write(out / "design.json", sim2.design_body(des))
    cores = {1: SMOKE_CORES} if a.smoke else SHARD_CORES
    proj = projection(des, cores)
    counts = {}
    for t in des["tasks"]:
        key = f"{t['cell']} {t['method']}"
        counts[key] = counts.get(key, 0) + 1
    waves = {k: first_wave(des, k, c) for k, c in cores.items()}
    rows = [check("tasks by cell and method", None, **counts),
            check("datasets", None, n=len(des["datasets"])),
            check("the seed (046)", None, p=des["p"], ratio=des["ratio"], seed_sha256=des["seed_sha256"]),
            check(f"every fit starts in its shard's first wave (sim2.startable: {sim2.DEADLINE_FACTOR:g} x the estimate "
                  f"within {DEADLINE_H:g} h)", all(w[0] for w in waves.values()),
                  **{f"shard{k}": {"cores": cores[k], "fits": w[1], "all_start": w[0]} for k, w in waves.items()}),
            check("projection: charged core-hours by shard at 1 and 1.5 times the estimates", None,
                  **{k: proj[k] for k in ("x1", "x1.5")}),
            check(f"cumulative charged core-hours within {CEILING:,}, the five shards' {PRIOR_CHARGED:,} included "
                  f"(PLAN.md risk 26)", all(proj[k]["cumulative_core_h"] <= CEILING for k in ("x1", "x1.5")),
                  **{k: proj[k]["cumulative_core_h"] for k in ("x1", "x1.5")}),
            check("charged core-hours if every fit took round 1's longest times the ratio, and if every job ran to the "
                  "fit cap", None, **{k: proj[k] for k in ("longest", "cap")}),
            check(f"round 1's longest default-seed fit times the ratio against the {FIT_CAP_H:g} h fit cap", None,
                  longest_h=LONGEST_H * RATIO, cap_h=FIT_CAP_H, within=LONGEST_H * RATIO <= FIT_CAP_H)]
    common.finish(out, prov, "SIM-2 seeded GTR20 design", rows, {"projection": proj})


def cmd_run(a):
    bad = binary_check(a.frozen, a.gate)
    if bad:
        sys.exit(f"not round 1's binaries (decision 046): {bad}")
    a.shards = design(a.smoke)["shards"]
    sim2.cmd_run(a, design, Path(__file__).resolve(), "SIM-2 seeded GTR20 run")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="step", required=True)
    d = sub.add_parser("design")
    d.add_argument("--smoke", action="store_true")
    d.add_argument("--out", help=f"new directory outside the repository (default under {RUNS}/)")
    r = sub.add_parser("run")
    r.add_argument("--shard", type=int, required=True)
    r.add_argument("--cores", type=int, required=True)
    r.add_argument("--deadline-hours", type=float, default=DEADLINE_H)
    r.add_argument("--smoke", action="store_true")
    r.add_argument("--skip-checks", action="store_true", help="development only (the smoke run)")
    k = sub.add_parser("collect")
    ph = sub.add_parser("phase", help="internal: one phase, started by run")
    ph.add_argument("--shard", type=int, required=True)
    ph.add_argument("--id", required=True)
    ph.add_argument("--smoke", action="store_true")
    ph.add_argument("--cores", type=int, required=True, help="the shard's cores, a cap on threads and workers")
    for p in (r, ph):
        p.add_argument("--frozen", required=True, help="the frozen binary (AliSim)")
        p.add_argument("--gate", required=True, help="the gate binary (IQ-TREE fits)")
        p.add_argument("--curve", default=str(sim2.CURVE), help="CLOAK's binned removal curve")
    for p in (d, r):
        p.add_argument("--git", default="git", help="git executable (git.exe under WSL)")
        p.add_argument("--allow-dirty", action="store_true")
    for p in (r, k, ph):
        p.add_argument("--out", required=True, help="the run directory")
    a = ap.parse_args()
    {"design": cmd_design, "run": cmd_run, "collect": lambda a: sim2.cmd_collect(a, "collect-gtr20"),
     "phase": lambda a: sim2.cmd_phase(a, design)}[a.step](a)


if __name__ == "__main__":
    main()
