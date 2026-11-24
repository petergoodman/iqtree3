"""SIM-2's driver (CODE_PLAN.md 3.8, "SIM-2's plan"; decisions 039 to 041 and 043 to 045). Python only; IQ-TREE
and AliSim run as the frozen and gate binaries given.

  sim2.py design [--smoke] [--out DIR] --git git.exe
  sim2.py checks --frozen BIN --gate BIN [--write-fixture] [--out DIR] --git git.exe [--allow-dirty]
  sim2.py run --out RUNDIR --shard K [--shards 5] --cores N --frozen BIN --gate BIN [--deadline-hours 46]
              [--curve FILE] [--smoke] [--skip-checks] --git git.exe [--allow-dirty]
  sim2.py collect --out RUNDIR

design: the datasets, the tasks (308 NQC fits, each with stage B, and 176 GTR20 fits), their phases with cores,
caps and estimates, the shards (cells D and F in the last, CLOAK's contingency) and the charged-cost projection
from the scheduler run on the estimates; deterministic. checks: decision 021's AliSim rows (alisim_check.py), the
pipeline and wiring rows for R0, D1, NQ and QB, and the cross-machine row against sim/fixtures/
sim2_crossmachine.json, which --write-fixture writes instead. run: one shard into RUNDIR/shard-K, checks first
(stopping on a failure), then phases longest first under --cores, each in its own process; a phase starts only if
1.5 times its estimate fits before the deadline, works in a .partial directory renamed on success, and on a later
run a leftover .partial becomes .interrupted-N and reruns. collect: one small tarball of RUNDIR's records.
--smoke: a miniature design (5 genes, replicate 1, 2 rounds, 60 s climbs, one shard) for the laptop.

Stage B, the oracle's climb from a capped IQ-TREE fit's endpoint (decision 039), is stage_b; the calibration
(explore/stage_b_calibration.py) times it on the pilot's endpoints.
"""
import argparse
import datetime
import heapq
import json
import os
import socket
import subprocess
import sys
import tarfile
import time
import traceback
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")   # each worker process one BLAS thread, before NumPy loads

import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
for _p in ("", "regression", "probes", "nqc", "differential", "explore", "profile", "reference"):
    sys.path.insert(0, str(HERE.parent / _p))
sys.path.insert(0, str(HERE))
from oracle import chart, flux, nqc, optimize, readers, simulate  # noqa: E402
from regress import REPO, read, sha256, write  # noqa: E402
from probes import PREFIX, read_checkpoint  # noqa: E402
from nqcchecks import FIT, check  # noqa: E402
from differential import write_matrix  # noqa: E402
import common  # noqa: E402
import pipeline as pl  # noqa: E402

SCRIPT = Path(__file__).resolve().relative_to(REPO).as_posix()
RUNS = "~/iqtree3-runs/sim2"
CERTIFICATE = 1e-9     # 021's certificate line, times max(1, |lnL|)
TOL_ACCEPTED = 1e-8    # 021's Level 1 row: the accepted score against the oracle's, times |lnL|

# ---------------------------------------------------------------- the design (CODE_PLAN 3.8)

BASE = 20261119                  # SIM-2's seed base, apart from the pilot's (20261117), whose SDs set n
TAXA26_SEED = pl.SEED            # cell E keeps the pilot's seeded 26-taxon set
ROUNDS = 20
HELDOUT_GENES, HELDOUT_REPLICATE = 10, 1000
TRUTHS = {"R0": "dial:0", "D1": "dial:0.1", "NQ": "NQ.bird", "QB": "Q.bird"}
SHIFTS = (("LG", "0.06"), ("Q.MAMMAL", "0.06"), ("Q.INSECT", "0.06"), ("WAG", "0.06"), ("Q.MAMMAL", "own"))
PROFILES = "profiles:C20"
STRENGTHS = (0.5, 1.0, 2.0)
THREADS = {30: 1, 100: 1, 300: 4, 1000: 16}
THREADS_26, THREADS_G, THREADS_H, THREADS_RELAXED = 2, 2, 16, 2
WORKERS = {30: 1, 100: 4, 300: 8, 1000: 32}
P_CAP_H = {30: 8, 100: 30, 300: 40, 1000: 40}
B_CAP_H = {30: 2, 100: 4, 300: 6, 1000: 6}      # provisional; set by the calibration (Peter's rule, CODE_PLAN 3.8)
CERTIFY = False                                  # the certified subset, decided by the calibration's probe
CERT_WORKERS, CERT_CAP_H = 8, 8.0
SHARDS, PINNED = 5, ("D", "Dr", "F")   # Dr: cell D's GTR20 convergence check
CURVE = HERE / "masking" / "cloak_curve.tsv"
FIXTURE = HERE / "fixtures" / "sim2_crossmachine.json"
SCAFFOLD = HERE / "scaffold"
# the smoke run: at 5 genes a GTR20 fit, or NQC's first linked round from LG, took over 15 min on two threads, so
# 3 genes; its relaxed check runs at the smoke's 2 rounds rather than 300 (the relaxed options are SIM-1's)
SMOKE = {"genes": 3, "rounds": 2, "p_cap_h": 1.0, "b_cap_h": 60 / 3600, "cert_cap_h": 60 / 3600}

# cost model (archive SIM2_plan, A2), with the calibration's stage B constants
PUMA = 1.3                                               # Puma's slowdown per core against the laptop, assumed
EVALUATIONS = {30: 29.5e3, 100: 39.7e3, 300: 45e3, 1000: 50e3}   # per capped fit; 300 and 1000 assumed
EFFICIENCY = {1: 1.0, 2: 0.85, 4: 0.65, 8: 0.55, 16: 0.45, 32: 0.35}
COLUMNS = 520
MULT = {"G": 2.5, "H": 3.3, "relaxed": 3.0}
SIMULATE_S_PER_GENE, DERIVE_S_PER_GENE, METRICS_H = 1.0, 0.3, 0.05
DEADLINE_FACTOR = 1.5


def ds_label(truth):
    return next((k for k, v in TRUTHS.items() if v == truth), truth.replace(":", "-"))


def fit_hours(genes, taxa, threads, mult=1.0):
    ev = EVALUATIONS.get(genes, EVALUATIONS[30])
    return (2.55e3 / 500 * (2 * taxa - 3) / 37 * 1e-6 * genes * COLUMNS * ev * mult * PUMA / 3600
            / EFFICIENCY.get(threads, 0.35) / threads)


def climb_hours(cap_h):
    """Stage B's estimate: its cap, since under the stall rule's low-flux condition no calibration climb stalled
    within its cap (CODE_PLAN 3.8, the calibration)."""
    return cap_h


def design(smoke=False):
    """Datasets and tasks, deterministic. A task's phases: NQC P, B, (cert), metrics; GTR20 fit, metrics."""
    g_ = (lambda g: SMOKE["genes"]) if smoke else (lambda g: g)
    reps = (lambda n: [1]) if smoke else (lambda n: list(range(1, n + 1)))
    hg = SMOKE["genes"] if smoke else HELDOUT_GENES
    datasets, tasks = {}, []

    def dataset(truth, genes, rep, kind="simulate", source=None, s=None):
        base = f"{ds_label(truth)}-g{genes}-r{rep}"
        did = base if kind == "simulate" else f"{source}-s{s:g}" if kind == "mask" else f"{source}-t26"
        datasets.setdefault(did, {"id": did, "kind": kind, "truth": truth, "genes": genes, "replicate": rep,
                                  "source": source, "s": s})
        return did

    def heldout(truth, pruned=False):
        did = dataset(truth, hg, HELDOUT_REPLICATE)
        return dataset(truth, hg, HELDOUT_REPLICATE, "prune", did) if pruned else did

    def task(cell, method, did, held, size, **kw):
        d = datasets[did]
        t = {"id": f"{cell}-{method}-{did}" + (f"-{kw['init']}" if kw.get("init") else ""), "cell": cell,
             "method": method, "dataset": did, "heldout": held, "size": size, "genes": d["genes"],
             "taxa": 26 if d["kind"] == "prune" else 52, "threads": THREADS[size], "workers": WORKERS[size],
             "rounds": SMOKE["rounds"] if smoke else ROUNDS, "relaxed": False, "level2": False, "init": None,
             "init_from": None, "certify": False}
        t.update(kw)
        if smoke:   # the laptop's smoke run: the same code paths on two cores a phase
            t["threads"], t["workers"] = min(t["threads"], 2), min(t["workers"], 2)
        tasks.append(t)
        return t

    for name, truth in TRUTHS.items():                              # A: the ladder
        if smoke and name not in ("R0", "NQ"):
            continue
        for size in ((30,) if smoke else (30, 100, 300)):
            for r in reps(10):
                task("A", "nqc", dataset(truth, g_(size), r), heldout(truth), size,
                     certify=(CERTIFY or smoke) and size == 30 and name in ("R0", "NQ") and r <= 3)
        if name != "QB" and not smoke:   # at the smoke's one gene count it would repeat the 30-gene task
            for r in reps(5):
                task("A", "nqc", dataset(truth, g_(1000), r), heldout(truth), 1000)
    for model, tv in (SHIFTS[:1] if smoke else SHIFTS):              # B: shifted targets
        truth = f"shift:{model}:{tv}"
        for r in reps(10):
            task("B", "nqc", dataset(truth, g_(300), r), heldout(truth), 300)
    for r in reps(10):                                               # C: site profiles
        task("C", "nqc", dataset(PROFILES, g_(300), r), heldout(PROFILES), 300)
    for name in (("R0",) if smoke else ("R0", "NQ")):
        truth = TRUTHS[name]
        for r in reps(10):                                           # D: masking
            base = dataset(truth, g_(100), r)
            for s in (STRENGTHS[1:2] if smoke else STRENGTHS):
                did = dataset(truth, g_(100), r, "mask", base, s)
                for method in ("nqc", "gtr20pi", "gtr20f"):
                    task("D", method, did, heldout(truth), 100)
            for method in ("gtr20pi", "gtr20f"):
                task("D", method, base, heldout(truth), 100)
        for r in reps(2):                                            # D: GTR20 convergence (039)
            base = dataset(truth, g_(100), r)
            for did in ((base,) if smoke else (base, f"{base}-s1")):
                for method in (("gtr20pi",) if smoke else ("gtr20pi", "gtr20f")):
                    task("Dr", method, did, heldout(truth), 100, relaxed=True, threads=THREADS_RELAXED)
        for r in reps(10):                                           # E: 26 taxa
            did = dataset(truth, g_(300), r, "prune", dataset(truth, g_(300), r))
            task("E", "nqc", did, heldout(truth, pruned=True), 300, threads=THREADS_26)
        for r in reps(5):                                            # F: starts
            base = dataset(truth, g_(100), r)
            task("F", "nqc", base, heldout(truth), 100, init="gtr20", init_from=f"D-gtr20pi-{base}")
            task("F", "nqc", base, heldout(truth), 100, init="transposed")
        for r in reps(5):                                            # G: Level 2
            task("G", "nqc", dataset(truth, g_(100), r), heldout(truth), 100, level2=True, threads=THREADS_G)
    for name in (("R0",) if smoke else ("R0", "D1", "NQ")):          # H: default rounds
        truth = TRUTHS[name]
        task("H", "nqc", dataset(truth, g_(300), 1), heldout(truth), 300, rounds=None, threads=THREADS_H)
    for t in tasks:
        t["phases"] = task_phases(t, smoke)
    shards = 1 if smoke else SHARDS
    return {"smoke": smoke, "base": BASE, "shards": shards, "datasets": datasets, "tasks": tasks,
            "assignment": assign(tasks, shards)}


def task_phases(t, smoke):
    """Each phase's name, cores, cap and estimate (hours); the D-cell relaxed fits are GTR20 only."""
    size, genes, taxa = t["size"], t["genes"], t["taxa"]
    p_cap = SMOKE["p_cap_h"] if smoke else P_CAP_H[size]
    if t["method"] != "nqc":
        mult = MULT["relaxed"] if t["relaxed"] else 1.0
        return [{"name": "fit", "cores": t["threads"], "cap_h": p_cap, "est_h": fit_hours(genes, taxa, t["threads"], mult)},
                {"name": "metrics", "cores": 1, "cap_h": None, "est_h": METRICS_H}]
    mult = MULT["G"] if t["level2"] else MULT["H"] if t["rounds"] is None else 1.0
    b_cap = SMOKE["b_cap_h"] if smoke else B_CAP_H[size]
    out = [{"name": "P", "cores": t["threads"], "cap_h": p_cap, "est_h": min(p_cap, fit_hours(genes, taxa, t["threads"], mult))},
           {"name": "B", "cores": t["workers"], "cap_h": b_cap, "est_h": climb_hours(b_cap)}]
    if t["certify"]:
        cap = SMOKE["cert_cap_h"] if smoke else CERT_CAP_H
        out.append({"name": "cert", "cores": min(CERT_WORKERS, 2) if smoke else CERT_WORKERS, "cap_h": cap,
                    "est_h": cap})
    return out + [{"name": "metrics", "cores": 1, "cap_h": None, "est_h": METRICS_H}]


def task_core_hours(t):
    return sum(p["cores"] * p["est_h"] for p in t["phases"])


def assign(tasks, shards):
    """Cells D and F (with the GTR20 fits F starts from) alone in the last shard, which waits for CLOAK's curve;
    the rest longest first to the least loaded of the others, ties to the lowest."""
    load, out = [0.0] * shards, {}
    for t in tasks:
        if t["cell"] in PINNED:
            out[t["id"]] = shards
            load[-1] += task_core_hours(t)
    free = range(shards - 1) if shards > 1 else range(1)
    for t in sorted((t for t in tasks if t["id"] not in out), key=lambda t: (-task_core_hours(t), t["id"])):
        k = min(free, key=lambda i: (load[i], i))
        out[t["id"]] = k + 1
        load[k] += task_core_hours(t)
    return out


def dataset_hours(d):
    per = SIMULATE_S_PER_GENE if d["kind"] == "simulate" else DERIVE_S_PER_GENE
    return d["genes"] * per / 3600


def shard_phases(des, shard):
    """The shard's phases by id: its datasets (data/ID, with their sources) and its tasks' phases (TASK/NAME),
    each with deps, cores, estimate and tail, the longest estimated path from its start to the end."""
    tasks = {t["id"]: t for t in des["tasks"] if des["assignment"][t["id"]] == shard}
    need, stack = set(), [d for t in tasks.values() for d in (t["dataset"], t["heldout"])]
    while stack:
        d = stack.pop()
        if d not in need:
            need.add(d)
            if des["datasets"][d]["source"]:
                stack.append(des["datasets"][d]["source"])
    phases = {}
    for d in sorted(need):
        ds = des["datasets"][d]
        phases[f"data/{d}"] = {"id": f"data/{d}", "kind": ds["kind"], "size": ds["genes"], "cores": 1,
                               "cap_h": None, "est_h": dataset_hours(ds), "dataset": d,
                               "deps": [f"data/{ds['source']}"] if ds["source"] else []}
    for t in tasks.values():
        first, prev, between = None, None, []
        for p in t["phases"]:
            pid = f"{t['id']}/{p['name']}"
            after = []
            if prev is None:   # P or the GTR20 fit
                deps = [f"data/{t['dataset']}"] + ([f"{t['init_from']}/fit"] if t["init_from"] else [])
                first = pid
            elif p["name"] == "metrics":   # needs the fit; waits for B (and cert), done or failed
                deps, after = [first, f"data/{t['heldout']}"], list(between)
            else:
                deps = [prev]
                between.append(pid)
            phases[pid] = {"id": pid, "kind": p["name"], "size": t["size"], "cores": p["cores"], "cap_h": p["cap_h"],
                           "est_h": p["est_h"], "task": t["id"], "deps": deps, "after": after}
            prev = pid
    children = {pid: [] for pid in phases}
    for p in phases.values():
        for d in p["deps"] + p.get("after", []):
            children[d].append(p["id"])
    tails = {}

    def tail(pid):
        if pid not in tails:
            tails[pid] = phases[pid]["est_h"] + max((tail(c) for c in children[pid]), default=0.0)
        return tails[pid]

    for pid in phases:
        phases[pid]["tail"] = tail(pid)
    return phases


# ---------------------------------------------------------------- the scheduler

def startable(phases, state, running, free, now_h, deadline_h, total, ratio=lambda p: 1.0):
    """Phases to start now, longest tail first: every dependency done, its cores free (a phase asking for more
    than the total takes the total), and DEADLINE_FACTOR times its (corrected) estimate before the deadline.
    running maps each running phase to its expected end (hours). EASY backfilling, so that wide phases are not
    starved by narrow ones: the first phase in order that fits the deadline but not the free cores reserves the
    earliest time enough cores free up, and a later phase starts only if it ends before then or uses cores the
    reservation leaves spare. A phase whose dependency failed or was blocked is marked blocked in state."""
    changed = True
    while changed:
        changed = False
        for p in phases.values():
            if p["id"] not in state and any(state.get(d) in ("failed", "blocked") for d in p["deps"]):
                state[p["id"]] = "blocked"
                changed = True
    ready = sorted((p for p in phases.values() if p["id"] not in state and p["id"] not in running
                    and all(state.get(d) == "done" for d in p["deps"])
                    and all(d in state for d in p.get("after", []))), key=lambda p: (-p["tail"], p["id"]))
    out, ends, shadow, spare = [], dict(running), None, 0
    for p in ready:
        need, dur = min(p["cores"], total), p["est_h"] * ratio(p)
        if now_h + DEADLINE_FACTOR * dur > deadline_h:
            continue
        if need > free:
            if shadow is None:
                avail, shadow = free, now_h
                for end, pid in sorted((e, i) for i, e in ends.items()):
                    if avail >= need:
                        break
                    avail, shadow = avail + min(phases[pid]["cores"], total), end
                spare = avail - need
            continue
        if shadow is not None and now_h + dur > shadow:
            if need > spare:
                continue
            spare -= need
        out.append(p)
        free -= need
        ends[p["id"]] = now_h + dur
    return out


def project(phases, cores, deadline_h, mult=1.0):
    """The scheduler on the estimates times mult: makespan, busy and charged core-hours (cores times the job's
    wall time), and the phases left for a later round."""
    state, running, now, busy, events = {}, {}, 0.0, 0.0, []
    while True:
        free = cores - sum(min(phases[i]["cores"], cores) for i in running)
        for p in startable(phases, state, running, free, now, deadline_h, cores):
            dur = p["est_h"] * mult
            running[p["id"]] = now + dur
            busy += min(p["cores"], cores) * dur
            heapq.heappush(events, (now + dur, p["id"]))
        if not running:
            break
        now, pid = heapq.heappop(events)
        del running[pid]
        state[pid] = "done"
    left = [p for p in phases if p not in state]
    return {"makespan_h": now, "busy_core_h": busy, "charged_core_h": cores * now, "phases_left": len(left),
            "left_core_h": sum(phases[p]["cores"] * phases[p]["est_h"] * mult for p in left)}


def projection(des, cores=94, deadline_h=46.0):
    out = {}
    for mult in (1.0, 1.5):
        shards = {k: project(shard_phases(des, k), cores, deadline_h, mult) for k in range(1, des["shards"] + 1)}
        out[f"x{mult:g}"] = {"shards": shards, "charged_core_h": sum(s["charged_core_h"] for s in shards.values()),
                             "round2_core_h": sum(s["left_core_h"] for s in shards.values())}
    return out


# ---------------------------------------------------------------- stage B

def nqc_state(fitdir):
    """A fit's linked NQC state from its checkpoint: coordinates (17 digits), accepted score, references and
    domain."""
    ck = read_checkpoint(Path(fitdir) / f"{PREFIX}.ckp.gz")

    def get(key):
        hits = [v for k, v in ck.items() if k.endswith(f"!NQC!ModelNonrevFixedEq!{key}")]
        if len(hits) != 1:
            raise ValueError(f"{fitdir}: {len(hits)} linked NQC values of {key} in the checkpoint")
        return [float(v) for v in hits[0].replace(" ", "").split(",")]

    return {"z": np.array(get("coordinates")), "accepted_score": get("accepted_score")[0],
            "references": [int(v) for v in get("references")], "domain": get("domain")}


def objective_for(data, tree_text, pi, level2=None):
    """The oracle's linked objective of a replicate's data: Level 1 as the wiring row builds it (gene k on the
    tree at partition rate 1 scaled by rho_k, its true shape); with level2, every gene on that fit's tree at its
    partition rate and shape (pl.level2_state)."""
    man = json.loads(read(Path(data) / "data.json"))
    alns = pl.gene_alignments_of(data, len(man["genes"]))
    if level2 is None:
        return pl.wiring_objective(alns, tree_text, man["genes"], pi), man
    ref = nqc.lg_seed(pi)[1]
    tree = readers.parse_newick(level2["tree"])
    return nqc.LinkedObjective([nqc.Objective(a, tree, pi, ref, shape=s, part_rate=r)
                                for a, s, r in zip(alns, level2["shapes"], level2["rates"])]), man


def stage_b(data, fitdir, out, workers, cap, tree_text, level2=None, **climb):
    """Decision 039's stage B from the fit in fitdir: the oracle re-scores P's checkpointed coordinates,
    optimize.stage_b_climb climbs from them over workers gene-parallel processes, and flux.certificate labels
    the endpoint "certified peak" or "highest found". Writes B.qmat into out; returns the record and its rows.
    The truth is scored on the training data when it lies at pi* (data.json's at_target, true when absent)."""
    data, out = Path(data), Path(out)
    man = json.loads(read(data / "data.json"))
    pi = np.array([float(v) for v in man["target"].split(",")])
    p = nqc_state(fitdir)
    obj, man = objective_for(data, tree_text, pi, level2)
    if [int(v) for v in obj.ref] != p["references"]:
        raise ValueError("P's references differ from LG's at pi*'s")
    lo, hi = p["domain"]
    at_target = man.get("at_target", True)
    Qt = readers.read_qmat(read(data / "truth.qmat"))[0]
    with nqc.GeneParallel(obj.parts, workers) as par:
        r = optimize.stage_b_climb(par.gradient, p["z"], lo, hi, nqc.climb_watch(pi, obj.ref), cap, **climb)
        Q = chart.build(r.x, pi, obj.ref)[0]
        value, gf = flux.LikelihoodFlux(par, pi)(flux.flux_of(Q, pi))
        truth_lnl = par.log_likelihood_gradient(Qt)[0] if at_target else None
    cert = flux.certificate(flux.flux_of(Q, pi), gf, flux.matched(obj.ref, lo, hi))
    line = CERTIFICATE * max(1.0, abs(value))
    lnl_p, lnl_b = -r.f_start, -r.f
    write_matrix(out / "B.qmat", Q, pi)
    rec = {"workers": workers, "cap_s": cap, "lnl_p_oracle": lnl_p, "accepted_score": p["accepted_score"],
           "lnl_b": lnl_b, "gain": r.gain, "iterations": r.iterations, "gradients": r.gradients,
           "seconds": round(r.seconds, 1), "stop": r.stop, "scipy_message": r.message, "pg_max": r.pg_max,
           "at_lower": len(r.at_lower), "at_upper": len(r.at_upper), "window": r.window,
           "C_p": simulate.circulation(chart.build(p["z"], pi, obj.ref)[0], pi), "C_b": simulate.circulation(Q, pi),
           "certificate": cert, "certificate_line": line,
           "label": "certified peak" if cert["residual"] <= line else "highest found",
           "truth_lnl": truth_lnl, "level2": level2 is not None, "z_b": r.x.tolist(), "trace": r.trace}
    rows = [check("B not below P (021)", lnl_b >= lnl_p, lnl_b=lnl_b, lnl_p=lnl_p)]
    if level2 is None:
        gap = abs(lnl_p - p["accepted_score"])
        rows.append(check("Level 1 accepted score against the oracle's, <= 1e-8 |lnL| (021)",
                          gap <= TOL_ACCEPTED * abs(lnl_p), difference=lnl_p - p["accepted_score"]))
    else:
        rec["rescore_gap"] = lnl_p - p["accepted_score"]   # G: reported, outside 021's row (CODE_PLAN 3.5)
    return rec, rows


def certify(data, fitdir, bdir, workers, cap, tree_text, level2=None):
    """The certified subset: the flux reference (decision 042) from B's endpoint, as exact.fit_and_check runs it,
    and its certificate."""
    from reference import Pool, fit_record

    man = json.loads(read(Path(data) / "data.json"))
    pi = np.array([float(v) for v in man["target"].split(",")])
    p = nqc_state(fitdir)
    obj = objective_for(data, tree_text, pi, level2)[0]
    dom = flux.matched(obj.ref, *p["domain"])
    z_b = json.loads(read(Path(bdir) / "record.json"))["record"]["z_b"]
    f0, moved = flux.inward(chart.build(np.asarray(z_b), pi, obj.ref)[0], pi, dom)
    pool = Pool(obj, workers)
    try:
        r = flux.fit(flux.LikelihoodFlux(obj, pi), f0, dom, objective_map=pool.fluxes, max_seconds=cap)
    finally:
        pool.close()
        if pool.pool is not None:
            pool.pool.join()
    rec = {"fit": fit_record(r), "start_moved": moved}
    if r.status != "failed":
        value, gf = flux.LikelihoodFlux(obj, pi)(r.f)
        cert = flux.certificate(r.f, gf, dom)
        rec.update({"lnl": value, "certificate": cert, "certified": cert["residual"] <= CERTIFICATE * max(1, abs(value))})
    return rec


# ---------------------------------------------------------------- phases

def phase_dir(shard_dir, p):
    if p["id"].startswith("data/"):
        return shard_dir / "data" / p["dataset"]
    return shard_dir / "tasks" / p["task"] / p["kind"]


def do_phase(des, phases, p, shard_dir, out, a):
    """Run one phase into out; returns its record (status done or failed, as 039's failures count) and rows."""
    if p["id"].startswith("data/"):
        return do_dataset(des["datasets"][p["dataset"]], shard_dir, out, a)
    t = next(t for t in des["tasks"] if t["id"] == p["task"])
    data = shard_dir / "data" / t["dataset"]
    task_dir = shard_dir / "tasks" / t["id"]
    cap_s = p["cap_h"] * 3600 if p["cap_h"] else None
    gate = Path(a.gate).expanduser().resolve()
    if p["kind"] in ("P", "fit"):
        init = start_matrix(t, shard_dir) if p["kind"] == "P" and t["init"] else None
        seeded = t.get("gtr20_start")   # decision 046's seeded arm (sim2_gtr20.py): never start_matrix or init_check
        gtr20 = {"gtr20pi": "pi", "gtr20f": "counted"}.get(t["method"], False)
        # the smoke run prints each linked fit (-v, printing only) so that F's first start can be checked
        relaxed = t["relaxed"] and not des["smoke"]
        rec = pl.fit_into(out, data, gate, min(t["threads"], a.cores), cap_s, level2=t["level2"], relaxed=relaxed,
                          gtr20=gtr20,
                          rounds=None if relaxed else t["rounds"],
                          init_model=("seed.paml", pl.gtr20_seed(seeded["p"], common.target_text())) if seeded
                          else init,
                          extra=["-v"] if des["smoke"] and init else [])
        worse = "Tree logl gets worse" in read(out / "fit" / "stdout.txt") if (out / "fit" / "stdout.txt").exists() \
            else False
        failed = (rec["exit_code"] != 0 or rec["at_cap"] or rec.get("assertion_or_abort") or worse
                  or not rec["data_match"])
        if init and not failed:
            try:    # a check of the start, not part of the fit: its failure must not fail the fit
                fits = FIT.findall(read(out / "fit" / "stdout.txt"))
                rec.update(init_check(t, data, out, gate, init, shard_dir, float(fits[0][1]) if fits else None))
            except Exception:  # noqa: BLE001
                rec["init_check_error"] = traceback.format_exc()[-2000:]
        # also after a fit stopped at its cap, so that a rate at the bound is flagged whenever a checkpoint exists
        rows = seeded_rows(out / "fit", seeded["p"], len(json.loads(read(data / "data.json"))["genes"])) \
            if seeded and (out / "fit" / f"{PREFIX}.ckp.gz").exists() else []
        return {"status": "failed" if failed else "done", "record": {**rec, "logl_worse": worse}}, rows
    tree_text = read(data / pl.TREE_RATE1).strip()
    level2 = pl.level2_state(task_dir / "P" / "fit", pl.leaf_order(tree_text)) if t["level2"] else None
    if p["kind"] == "B":
        rec, rows = stage_b(data, task_dir / "P" / "fit", out, min(t["workers"], a.cores), cap_s, tree_text, level2)
        return {"status": "done", "record": rec}, rows
    if p["kind"] == "cert":
        return {"status": "done", "record": certify(data, task_dir / "P" / "fit", task_dir / "B",
                                                    min(p["cores"], a.cores), cap_s, tree_text, level2)}, []
    return {"status": "done", "record": metrics(t, data, shard_dir / "data" / t["heldout"], task_dir, out, gate)}, []


SEED_READ = "Reading model parameters from file"            # ModelMarkov::readParameters, at every verbosity
BOUNDARY = "Estimated model parameters are at boundary"     # IQTree::optimizeModelParameters' warning


def seeded_rows(fitdir, p, genes):
    """Decision 046's rows for a seeded GTR20 fit (021): from the checkpoint, the fixed Tyr-Val's 10-digit text
    against p's and no free rate above 99 (asserted), the free rates at the floor (reported); from stdout,
    reported only (printed output), the seed's reads against the gene count and IQ-TREE's boundary warning. A
    failed row, or a checkpoint that cannot be read, fails in the shard's report; the phase stays done."""
    try:
        b = pl.rate_bounds(pl.gtr20_state(read_checkpoint(Path(fitdir) / f"{PREFIX}.ckp.gz")))
    except Exception:  # noqa: BLE001
        return [check("seeded GTR20: the checkpoint's rates can be read (046)", False,
                      error=traceback.format_exc()[-1000:])]
    text = read(Path(fitdir) / "stdout.txt") if (Path(fitdir) / "stdout.txt").exists() else ""
    reads = sum(1 for line in text.splitlines() if line.startswith(SEED_READ))
    return [check("seeded GTR20: the fixed Tyr-Val's checkpoint text equals p's (046, 021)",
                  b["fixed_text"] == f"{p:.10g}", fixed_text=b["fixed_text"], p=p),
            check("seeded GTR20: no free rate above 99, IQ-TREE's 0.99 x MAX_RATE (046, 021)", b["above_upper"] == 0,
                  above_upper=b["above_upper"], largest_free=b["largest_free"], headroom=b["headroom"]),
            check("seeded GTR20: free rates at the floor, <= 1.0001e-4", None, at_floor=b["at_floor"]),
            check("seeded GTR20: printed seed reads against the genes, and boundary warnings (printed output)", None,
                  seed_reads=reads, genes=genes, boundary_warnings=text.count(BOUNDARY))]


def do_dataset(ds, shard_dir, out, a):
    if ds["kind"] == "simulate":
        tree, loci = pl.read_scaffold(SCAFFOLD)
        manifest, rows = pl.simulate_into(out, Path(a.frozen).expanduser().resolve(), tree, loci, ds["truth"],
                                          ds["genes"], BASE, ds["replicate"])
        failed = [r["check"] for r in rows if r["passed"] is False]
        return {"status": "failed" if failed else "done", "record": {"failed": failed,
                                                                     "sha256": [g["sha256"] for g in manifest["genes"]]}}, rows
    src = shard_dir / "data" / ds["source"]
    if ds["kind"] == "mask":
        curve = Path(a.curve).expanduser().resolve()
        man, same = pl.mask_dataset(src, out, ds["s"], *pl.read_curve(curve), curve_sha256=sha256(curve))
        return {"status": "done", "record": man["masking"]}, []
    keep = pl.pruned_taxa(pl.read_scaffold(SCAFFOLD)[0], 26, TAXA26_SEED)
    man, same = pl.prune_dataset(src, out, keep)
    return {"status": "done" if same else "failed", "record": {"charsets_kept": same, "taxa": man["taxa"]}}, \
        [check("pruning keeps the charsets", same)]


def start_matrix(t, shard_dir):
    """Cell F's --init-model file: the GTR20+F{pi*} fit's checkpoint carried as PAML (decision 025), or NQ.bird's
    flux transposed, a full -m FILE matrix."""
    if t["init"] == "gtr20":
        ck = read_checkpoint(shard_dir / "tasks" / t["init_from"] / "fit" / "fit" / f"{PREFIX}.ckp.gz")
        return "start.paml", pl.gtr20_paml(pl.gtr20_state(ck), common.target_text())
    pi, truths = simulate.truths()
    path = shard_dir / "tasks" / t["id"] / "start.qmat"
    path.parent.mkdir(parents=True, exist_ok=True)
    write_matrix(path, simulate.transposed(truths["NQ.bird"], pi), pi)
    return "start.qmat", read(path)


def start_q(t, shard_dir, pi):
    """The matrix cell F's start file stands for, as NQC builds it at pi: reversibleSeed of the carried
    exchangeabilities, or the transposed flux's generator (already at pi)."""
    if t["init"] == "gtr20":
        ck = read_checkpoint(shard_dir / "tasks" / t["init_from"] / "fit" / "fit" / f"{PREFIX}.ckp.gz")
        return chart.reversible_seed(pl.gtr20_state(ck)["R"], pi)
    pi_s, truths = simulate.truths()
    return simulate.transposed(truths["NQ.bird"], pi_s)


def init_check(t, data, out, gate, init, shard_dir, start_score, timeout=3600):
    """Whether --init-model reached the linked -q fit (CODE_PLAN 3.8, F): IQ-TREE's --show-lh total of the start
    matrix on the training data (a PAML file with +F{pi*}, a full matrix as itself; 3 decimals) against the start
    lnL of NQC's first linked fit, printed only at -v (the smoke run), decided only when no start coordinate lies
    outside the default domain (else the box clamps the start, which is reported)."""
    man = json.loads(read(data / "data.json"))
    pi = np.array([float(v) for v in man["target"].split(",")])
    z0 = chart.inverse(start_q(t, shard_dir, pi), nqc.lg_seed(pi)[1])
    outside = int(np.sum((z0 < nqc.LOWER) | (z0 > nqc.UPPER)))
    model = f"{init[0]}+F{{{man['target']}}}" if init[0].endswith(".paml") else init[0]
    files = {"parts.nex": pl.q_file(read(data / "sets.nex"), man["genes"], man["rates"], model=model), init[0]: init[1]}
    run = pl.Runner(gate, out, timeout).run("init-lnl", ["-s", data / "alignment.phy", "-q", "parts.nex", "-te",
                                                         data / "tree.nwk", "--show-lh", "-keep-ident"], files)
    printed = pl.TOTAL.findall(run.stdout())
    lnl = float(printed[-1]) if run.code == 0 and printed else None
    reached = None
    if outside == 0 and lnl is not None and start_score is not None:
        reached = abs(lnl - start_score) <= max(TOL_ACCEPTED * abs(lnl), 1e-3)   # the total prints 3 decimals
    return {"init_lnl": lnl, "start_outside_domain": outside, "init_reached": reached}


def metrics(t, data, held, task_dir, out, gate):
    """The estimands at P and B (NQC) or of the GTR20 fit as fitted and with pi* in place of its frequencies."""
    pi = common.target_values()
    if t["method"] == "nqc":
        qp = task_dir / "P" / "fit" / f"{PREFIX}.NQC.qmat"
        rec = {"P": pl.matrix_metrics(readers.read_qmat(read(qp))[0], data, held, gate, rundir=out, name="P",
                                      qmat_text=read(qp))}
        brec = task_dir / "B" / "record.json"
        b = json.loads(read(brec)) if brec.exists() else {}
        if b.get("status") != "done":   # P's estimands survive a failed B
            rec["B"] = None
            return rec
        b = b["record"]
        rec["B"] = pl.matrix_metrics(readers.read_qmat(read(task_dir / "B" / "B.qmat"))[0], data, held, gate,
                                     rundir=out, name="B", qmat_text=read(task_dir / "B" / "B.qmat"))
        if b["truth_lnl"] is not None:
            rec["P"]["training_excess"] = b["lnl_p_oracle"] - b["truth_lnl"]
            rec["B"]["training_excess"] = b["lnl_b"] - b["truth_lnl"]
        return rec
    st = pl.gtr20_state(read_checkpoint(task_dir / "fit" / "fit" / f"{PREFIX}.ckp.gz"))
    pooled = all(np.allclose(f, st["freq"], rtol=1e-9, atol=1e-12) for f in st["part_freqs"].values())
    return {"fitted": pl.matrix_metrics(chart.reversible_seed(st["R"], st["freq"]), data, held, gate, rundir=out,
                                        name="fitted", pi_fit=st["freq"]),
            "replaced": pl.matrix_metrics(chart.reversible_seed(st["R"], pi), data, held, gate, rundir=out,
                                          name="replaced"),
            "exchangeabilities": st["R"][np.triu_indices(20, 1)].tolist(), "frequencies": st["freq"].tolist(),
            "pooled": pooled, "partitions": len(st["part_freqs"])}


def cmd_phase(a, design_fn=design):
    """One phase in its own process: works in a .partial directory, writes record.json, renames it. design_fn is
    the study's design (sim2_gtr20.py passes its own)."""
    des = design_fn(a.smoke)
    shard_dir = Path(a.out).expanduser().resolve() / f"shard-{a.shard}"
    phases = shard_phases(des, a.shard)
    p = phases[a.id]
    final = phase_dir(shard_dir, p)
    partial = final.with_name(final.name + ".partial")
    partial.mkdir(parents=True)
    t0 = time.monotonic()
    try:
        rec, rows = do_phase(des, phases, p, shard_dir, partial, a)
    except Exception:  # noqa: BLE001  a phase that raises is a counted failure (039)
        rec, rows = {"status": "failed", "error": traceback.format_exc()[-3000:]}, []
    rec.update({"phase": p["id"], "seconds": round(time.monotonic() - t0, 1), "cores": p["cores"],
                "est_h": p["est_h"], "rows": rows})
    write(partial / "record.json", json.dumps(rec, indent=1, default=common._json))
    partial.rename(final)


# ---------------------------------------------------------------- checks

def checks(out, frozen, gate, git, allow_dirty, write_fixture=False, prov=None):
    """Rows of decision 021's job-start checks; with write_fixture, the cross-machine genes are written to FIXTURE
    instead of compared with it, with prov's commit, cleanliness, date and run directory."""
    tree, loci = pl.read_scaffold(SCAFFOLD)
    rows = []
    cmd = [sys.executable, str(HERE / "alisim_check.py"), "run", "--frozen", str(frozen), "--gate", str(gate), "--out",
           str(out / "alisim"), "--git", git] + (["--allow-dirty"] if allow_dirty else [])
    code = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT).returncode
    res = json.loads(read(out / "alisim" / "results.json")) if (out / "alisim" / "results.json").exists() else {}
    failed = [r["check"] for r in res.get("results", []) if r.get("passed") is False]
    rows.append(check("AliSim exact (021's three AliSim rows)", code == 0 and bool(res) and not failed, exit_code=code,
                      failed=failed))
    for name in TRUTHS:
        d = out / f"pipeline-{name}"
        d.mkdir()
        for r in pl.checks_rows(d, frozen, gate, tree, loci, TRUTHS[name], 20, pl.SEED, 3600):
            rows.append({**r, "check": f"{name}: {r['check']}"})
    genes = {}
    for name, truth in (("R0", TRUTHS["R0"]), ("profiles", PROFILES)):
        d = out / f"crossmachine-{name}"
        d.mkdir()
        man, found = pl.simulate_into(d, frozen, tree, loci, truth, 30, BASE, 1)
        genes[name] = [g["sha256"] for g in man["genes"]]
        rows += found
    if write_fixture:
        FIXTURE.parent.mkdir(parents=True, exist_ok=True)
        made = {k: (prov or {}).get(k) for k in ("git_commit", "git_clean", "started_utc", "out")}
        write(FIXTURE, json.dumps({"provenance": {"script": SCRIPT, **made, "base": BASE, "genes": 30, "replicate": 1,
                                                  "frozen_sha256": sha256(frozen)}, "genes": genes}, indent=1) + "\n")
        rows.append(check("cross-machine fixture written", None, path=str(FIXTURE)))
    else:
        want = json.loads(read(FIXTURE))["genes"] if FIXTURE.exists() else None
        rows.append(check("cross-machine: genes byte-identical to the committed fixture (021)", want == genes,
                          fixture=str(FIXTURE), present=want is not None))
    return rows


def binaries(a):
    frozen, gate = (Path(p).expanduser().resolve() for p in (a.frozen, a.gate))
    return frozen, gate, {"frozen": {"path": str(frozen), "sha256": sha256(frozen)},
                          "gate": {"path": str(gate), "sha256": sha256(gate)}}


def cmd_checks(a):
    prov = common.provenance(a, SCRIPT, "checks", "confirmatory", write_fixture=a.write_fixture)
    out = common.new_out(a.out, prov, runs=RUNS)
    frozen, gate, prov["binaries"] = binaries(a)
    common.finish(out, prov, "SIM-2 checks", checks(out, frozen, gate, a.git, a.allow_dirty, a.write_fixture, prov))


# ---------------------------------------------------------------- design, run, collect

def design_body(des):
    return json.dumps(des, indent=1, sort_keys=True, default=common._json)


def cmd_design(a):
    prov = common.provenance(a, SCRIPT, "design-smoke" if a.smoke else "design", "design", smoke=a.smoke)
    out = common.new_out(a.out, prov, runs=RUNS)
    des = design(a.smoke)
    write(out / "design.json", design_body(des))
    proj = projection(des, cores=a.cores, deadline_h=a.deadline_hours)
    counts = {}
    for t in des["tasks"]:
        key = f"{t['cell']} {t['method']}"
        counts[key] = counts.get(key, 0) + 1
    rows = [check("tasks by cell and method", None, **counts),
            check("datasets", None, n=len(des["datasets"])),
            check("projection: charged core-hours by shard at 1 and 1.5 times the estimates", None, **proj),
            check("projection within the 40,000 charged core-hour ceiling (PLAN.md risk 26)",
                  all(v["charged_core_h"] + v["round2_core_h"] <= 40000 for v in proj.values()),
                  **{k: v["charged_core_h"] + v["round2_core_h"] for k, v in proj.items()})]
    common.finish(out, prov, "SIM-2 design", rows, {"projection": proj})


def measured(rec, p):
    """A finished phase's measured-to-estimated duration, for the estimates' correction: phases that ran to
    the end, including fits stopped at their cap; None for a phase that died or raised."""
    if p["est_h"] <= 0 or "seconds" not in rec or "error" in rec or "exit_code" in rec and "record" not in rec:
        return None
    if rec["status"] == "failed" and not rec.get("record", {}).get("at_cap"):
        return None
    return rec["seconds"] / 3600 / p["est_h"]


def scan(shard_dir, phases):
    """Finished phases' states from their records, with their measured ratios by phase kind and size (so a later
    round starts from them); a leftover .partial renamed .interrupted-N."""
    state, interrupted, ratios = {}, [], {}
    for p in phases.values():
        final = phase_dir(shard_dir, p)
        partial = final.with_name(final.name + ".partial")
        if partial.exists():
            n = 1
            while final.with_name(f"{final.name}.interrupted-{n}").exists():
                n += 1
            partial.rename(final.with_name(f"{final.name}.interrupted-{n}"))
            interrupted.append(p["id"])
        if (final / "record.json").exists():
            rec = json.loads(read(final / "record.json"))
            state[p["id"]] = rec["status"]
            m = measured(rec, p)
            if m is not None:
                ratios.setdefault((p["kind"], p["size"]), []).append(m)
    return state, interrupted, ratios


def cmd_run(a, design_fn=design, script=Path(__file__).resolve(), title="SIM-2 run"):
    """One shard of the study whose design is design_fn and whose driver is script (each phase runs as
    `script phase`); sim2_gtr20.py passes its own."""
    run_dir = Path(a.out).expanduser().resolve()
    if run_dir.is_relative_to(REPO):
        sys.exit(f"{run_dir}: the run directory must lie outside the repository")
    prov = common.provenance(a, Path(script).relative_to(REPO).as_posix(), f"run-shard{a.shard}",
                             "confirmatory" if not a.smoke else "development",
                             shard=a.shard, cores=a.cores, deadline_hours=a.deadline_hours, smoke=a.smoke,
                             skip_checks=a.skip_checks, curve=str(a.curve))
    frozen, gate, prov["binaries"] = binaries(a)
    shard_dir = run_dir / f"shard-{a.shard}"
    shard_dir.mkdir(parents=True, exist_ok=True)
    t_start = time.monotonic()
    des = design_fn(a.smoke)
    if des["shards"] != a.shards:
        sys.exit(f"the design has {des['shards']} shards, not {a.shards}")
    body = design_body(des)
    if (run_dir / "design.json").exists():
        if read(run_dir / "design.json") != body:
            sys.exit(f"{run_dir / 'design.json'} differs from this commit's design")
    else:   # five jobs on different nodes may start together: a temporary name per host and process
        tmp = run_dir / f"design.json.{socket.gethostname()}.{os.getpid()}"
        write(tmp, body)
        tmp.replace(run_dir / "design.json")
    n = 1
    while (shard_dir / f"run-{n}").exists():
        n += 1
    log_dir = shard_dir / f"run-{n}"
    log_dir.mkdir()
    rows = []
    phases = shard_phases(des, a.shard)
    if any(p["kind"] == "mask" for p in phases.values()) and not Path(a.curve).expanduser().exists():
        sys.exit(f"shard {a.shard} holds cell D, and CLOAK's curve {a.curve} is missing (CLOAK contingency)")
    if not a.skip_checks:
        (log_dir / "checks").mkdir()
        rows = checks(log_dir / "checks", frozen, gate, a.git, a.allow_dirty)
        if any(r["passed"] is False for r in rows):
            common.finish(log_dir, prov, f"{title}, shard {a.shard}: checks failed", rows)
    state, interrupted, ratios = scan(shard_dir, phases)
    deadline_h = a.deadline_hours - (time.monotonic() - t_start) / 3600
    t0 = time.monotonic()
    running = {}

    def ratio(p):
        r = ratios.get((p["kind"], p["size"]))
        return float(np.median(r)) if r else 1.0

    base = [sys.executable, str(script), "phase", "--out", str(run_dir), "--shard", str(a.shard),
            "--frozen", str(frozen), "--gate", str(gate), "--curve", str(a.curve), "--cores", str(a.cores)] \
        + (["--smoke"] if a.smoke else [])
    print(f"shard {a.shard}: {len(phases)} phases, {len(state)} already finished, {len(interrupted)} interrupted",
          flush=True)
    while True:
        for pid, (proc, started, fh) in list(running.items()):
            if proc.poll() is None:
                continue
            fh.close()
            del running[pid]
            p = phases[pid]
            final = phase_dir(shard_dir, p)
            if (final / "record.json").exists():
                state[pid] = json.loads(read(final / "record.json"))["status"]
            else:   # the phase process died before its record: counted as failed, kept for inspection
                partial = final.with_name(final.name + ".partial")
                if partial.exists():
                    partial.rename(final)
                final.mkdir(parents=True, exist_ok=True)
                write(final / "record.json", json.dumps({"status": "failed", "exit_code": proc.returncode,
                                                         "phase": pid}))
                state[pid] = "failed"
            hours = (time.monotonic() - started) / 3600
            m = measured(json.loads(read(final / "record.json")), p)
            if m is not None:
                ratios.setdefault((p["kind"], p["size"]), []).append(m)
            print(f"{datetime.datetime.now():%H:%M:%S} {state[pid]} {pid} ({hours:.2f} h)", flush=True)
        free = a.cores - sum(min(phases[i]["cores"], a.cores) for i in running)
        now_h = (time.monotonic() - t0) / 3600
        ends = {i: (s - t0) / 3600 + phases[i]["est_h"] * ratio(phases[i]) for i, (_, s, _) in running.items()}
        for p in startable(phases, state, ends, free, now_h, deadline_h, a.cores, ratio):
            fh = open(log_dir / (p["id"].replace("/", "__") + ".log"), "w", encoding="utf-8")
            proc = subprocess.Popen(base + ["--id", p["id"]], stdout=fh, stderr=subprocess.STDOUT)
            running[p["id"]] = (proc, time.monotonic(), fh)
            print(f"{datetime.datetime.now():%H:%M:%S} start {p['id']} ({min(p['cores'], a.cores)} cores, "
                  f"est {p['est_h'] * ratio(p):.2f} h)", flush=True)
        if not running:
            break
        time.sleep(2 if a.smoke else 20)
    counts = {s: sum(1 for v in state.values() if v == s) for s in ("done", "failed", "blocked")}
    counts["not_started"] = len(phases) - len(state)
    rows.append(check(f"shard {a.shard}: phases", None, interrupted=interrupted, **counts))
    for p in phases.values():
        final = phase_dir(shard_dir, p)
        if (final / "record.json").exists():
            for r in json.loads(read(final / "record.json")).get("rows", []):
                rows.append({**r, "check": f"{p['id']}: {r['check']}"})
    common.finish(log_dir, prov, f"{title}, shard {a.shard}", rows, {"state": state})


COLLECT = ("record.json", "results.json", "report.md", "data.json", "design.json", "B.qmat", f"{PREFIX}.NQC.info",
           f"{PREFIX}.NQC.qmat", "start.paml", "start.qmat", "seed.paml")


def cmd_collect(a, prefix="collect"):
    """One tarball of RUNDIR's records, named prefix-<stamp>.tar.gz: no alignments, checkpoints or logs, but a
    failed phase's stdout."""
    run_dir = Path(a.out).expanduser().resolve()
    stamp = datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    target = run_dir / f"{prefix}-{stamp}.tar.gz"
    with tarfile.open(target, "w:gz") as tar:
        for path in sorted(run_dir.rglob("*")):
            if not path.is_file() or path.name.startswith("collect-"):
                continue
            keep = path.name in COLLECT or (path.suffix == ".log" and path.parent.name.startswith("run-"))
            if path.name == "stdout.txt":
                rec = next((q / "record.json" for q in path.parents if (q / "record.json").exists()), None)
                keep = rec is not None and json.loads(read(rec)).get("status") == "failed"
            if keep:
                tar.add(path, arcname=str(path.relative_to(run_dir)))
    print(f"wrote {target} ({target.stat().st_size / 1e6:.1f} MB)")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="step", required=True)
    d = sub.add_parser("design")
    d.add_argument("--smoke", action="store_true")
    d.add_argument("--cores", type=int, default=94)
    d.add_argument("--deadline-hours", type=float, default=46.0)
    c = sub.add_parser("checks")
    c.add_argument("--write-fixture", action="store_true")
    r = sub.add_parser("run")
    r.add_argument("--shard", type=int, required=True)
    r.add_argument("--shards", type=int, default=SHARDS)
    r.add_argument("--cores", type=int, required=True)
    r.add_argument("--deadline-hours", type=float, default=46.0)
    r.add_argument("--smoke", action="store_true")
    r.add_argument("--skip-checks", action="store_true", help="development only (the smoke run)")
    k = sub.add_parser("collect")
    ph = sub.add_parser("phase", help="internal: one phase, started by run")
    ph.add_argument("--shard", type=int, required=True)
    ph.add_argument("--id", required=True)
    ph.add_argument("--smoke", action="store_true")
    ph.add_argument("--cores", type=int, required=True, help="the shard's cores, a cap on threads and workers")
    for p in (c, r, ph):
        p.add_argument("--frozen", required=True, help="the frozen binary (AliSim)")
        p.add_argument("--gate", required=True, help="the gate binary (IQ-TREE fits)")
    for p in (r, ph):
        p.add_argument("--curve", default=str(CURVE), help="CLOAK's binned removal curve")
    for p in (d, c, r):
        p.add_argument("--git", default="git", help="git executable (git.exe under WSL)")
        p.add_argument("--allow-dirty", action="store_true")
    for p in (d, c):
        p.add_argument("--out", help="new directory outside the repository (default under ~/iqtree3-runs/sim2/)")
    for p in (r, k, ph):
        p.add_argument("--out", required=True, help="the run directory")
    a = ap.parse_args()
    {"design": cmd_design, "checks": cmd_checks, "run": cmd_run, "collect": cmd_collect, "phase": cmd_phase}[a.step](a)


if __name__ == "__main__":
    main()
