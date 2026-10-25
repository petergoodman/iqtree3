"""Exploratory: S3c's first part, the stopping study and the converged evidence for PLAN.md risks 22
and 24 (CODE_PLAN.md section 4, S3c). Checks S-1 to S-3 and each case's oracle match are asserted;
everything else is reported.

  python test_scripts/fixedeq/explore/convergence.py run --binary BIN --out DIR --git git.exe \
      [--cases A,B,C] [--rounds 300] [--max-iter 2000] [--jobs 7] [--workers 11] [--allow-dirty]
  python test_scripts/fixedeq/explore/convergence.py reread --gate DIR --probe DIR --out DIR --git git.exe

run: three Level 1 cases with every nuisance fixed: A aa_example (+G4{0.5}, baseline run02's rooted
tree); B turtle -S (+G4{0.5}, run09's trees); C turtle -q with the tree, shapes and partition rates
of an NQC -p Level 2 +G4 fit in the default box (decision 024). In each, NONREV, NQC in the default
box and NQC at 1e-8,1e7 run at IQ-TREE's default settings and with the round rules relaxed (-me
1e-6, --loop-model and -nparam for ROUNDS rounds), and NQC relaxed at 1e-5,1e5 and 1e-6,1e6. The
oracle's reference climber (oracle/optimize.py lbfgsb_fit) then starts at NQC's relaxed box
endpoint and climbs to the peak of each domain in turn, P0, P1, P1b (added after case A's pilot),
P2 and P3, each from the last, and climbs P2's domain again from NQC's relaxed endpoint at
1e-8,1e7 (P2w), which may lie by another peak. The six built-in NQ matrices are placed against
the domains. IQ-TREE runs at -T 1, several at once, so its times are
indicative only; the oracle spreads each gradient's probes over worker processes, which changes no
value. Writes results.json and results.md in the new directory DIR.

reread: recomputes from S3b's gate outputs (GATE is the gate's run directory) and S0 probe c2 (PROBE
is the probe record directory) the numbers that correct the S3b record: the turtle -q
continuations' per-round gains, evaluations per round and dfpmin's test at their first gradient;
c2's per-round gains; final coordinates outside the default box, by row, in the gate's NQC runs.
Writes reread.json and reread.md in the new directory DIR.
"""
import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse  # noqa: E402
import datetime  # noqa: E402
import json  # noqa: E402
import math  # noqa: E402
import multiprocessing  # noqa: E402
import platform  # noqa: E402
import re  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
import traceback  # noqa: E402
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
import scipy  # noqa: E402

HERE = Path(__file__).resolve().parent
for _p in ("", "nqc", "probes", "regression", "differential"):
    sys.path.insert(0, str(HERE.parent / _p))
from oracle import builtin, chart, nqc, optimize, readers, target  # noqa: E402
from probes import AA, NEX, TURTLE, Run, Runner, baseline_runs, tree_text, turtle_target  # noqa: E402
from regress import REPO, capture, read, sha256, write  # noqa: E402
from nqcchecks import (ALPHA, BUG_ABORT, FALL, FIT, INITIAL, SWITCH, TOL_LNL, WIDE, check, final_struct,  # noqa: E402
                       first_gradient, floats, joint_args, part_alignment, q_charpartition, q_route, turtle_files)

SCRIPT = Path(__file__).resolve().relative_to(REPO).as_posix()
ME_DEFAULT, LOOP_DEFAULT, NPARAM_DEFAULT = 0.01, 10, 100   # utils/tools.cpp's defaults
ME_RELAXED = 1e-6
PGTOL = 1e-3        # the reference climber's line: largest projected-gradient component, -lnL per unit
MATERIAL = 0.1      # Peter, 2026-10-24: a converged box cost above this is material
STALL = 0.01        # S-3
FLAT = 1e-3         # P3: below this the coordinates leaving P2's domain add nothing
CONVERGED = -1e-6   # a later peak further below an earlier one than this was not reached
COLLAPSE = 3        # coordinates of one row at the upper bound that mark its reference as collapsed
RATE_FLOOR = 1e-3   # rates compared between peaks
GTOL = 1e-4         # dfpmin's gtol, max(gradient_epsilon, TOL_RATE)
# (peak, log lower, log upper, --nqc-domain value). P1b was added after case A's pilot (Peter,
# 2026-10-24): the smallest symmetric box holding every within-row rate ratio NONREV's bounds
# [1e-4, 100] allow.
DOMAINS = (("P0", nqc.LOWER, nqc.UPPER, None), ("P1", math.log(1e-5), math.log(1e5), "1e-5,1e5"),
           ("P1b", math.log(1e-6), math.log(1e6), "1e-6,1e6"), ("P2", math.log(1e-8), math.log(1e7), WIDE),
           ("P3", -50.0, 50.0, None))
DOMAIN = {d[0]: d for d in DOMAINS}
CASES = {"A": ("aa_example Level 1", "single", 0.01), "B": ("turtle Level 1 -S", "-S", 0.1),
         "C": ("turtle Level 1 -q", "-q", 0.01)}   # label, path, the default round threshold
MODELS = (("nonrev", False, None), ("nqc-box", True, None), ("nqc-wide", True, WIDE))
NQ = ("NQ.PFAM", "NQ.INSECT", "NQ.YEAST", "NQ.PLANT", "NQ.BIRD", "NQ.MAMMAL")
WIDTH = 18          # coordinates per row (n - 2)

CURRENT = re.compile(r"^\d+\. Current log-likelihood: (\S+)", re.M)
TOOK = re.compile(r"^Parameters optimization took (\d+) rounds", re.M)
ROUND_S = re.compile(r"^\d+\. Log-likelihood: (\S+)$", re.M)
INITIAL_Q = re.compile(r"^Initial log-likelihood: (\S+)$", re.M)
STEP_Q = re.compile(r"^Current log-likelihood at step \d+: (\S+)$", re.M)
BEST = re.compile(r"^BEST SCORE FOUND : (\S+)", re.M)


class Context:
    def __init__(self):
        self.target_text = turtle_target()
        self.pi = target.parse_target(self.target_text)[0]
        self.ref = nqc.lg_seed(self.pi)[1]
        self.t02 = tree_text(baseline_runs(), "run02")
        self.model = f"NQC+F{{{self.target_text}}}"


# ---------------------------------------------------------------- the oracle's side

def objective(spec, pi, ref):
    """The oracle's -lnL for a case: under -S each partition on its own tree, which holds only its
    taxa with data; under -q one tree, missing sequences as unknowns, partition rates normalized as
    IQ-TREE normalizes them."""
    if spec["mode"] == "single":
        return nqc.Objective(readers.read_phylip(AA), readers.parse_newick(spec["tree"]), pi, ref, shape=ALPHA)
    aln, charsets = readers.read_fasta(TURTLE), readers.read_charsets(NEX)
    trees = [t for t in spec["trees"].splitlines() if t.strip()] if spec["mode"] == "-S" else None
    rates = (nqc.normalize_part_rates(spec["rates"], [len(s) for s in charsets.values()])
             if spec["mode"] == "-q" else None)
    parts = []
    for k, sites in enumerate(charsets.values()):
        tree = readers.parse_newick(trees[k] if trees else spec["tree"])
        leaves = {leaf.name for leaf in tree.leaves()}
        sub = {n: s for n, s in part_alignment(aln, sites).items() if n in leaves}
        shape = ALPHA if spec["mode"] == "-S" else spec["shapes"][k]
        parts.append(nqc.Objective(sub, tree, pi, ref, shape=shape, part_rate=rates[k] if rates else 1.0))
    return nqc.LinkedObjective(parts)


_OBJECTIVE = None


def _init_worker(spec, pi, ref):
    global _OBJECTIVE
    _OBJECTIVE = objective(spec, np.array(pi), np.array(ref))


def _evaluate(x):
    return _OBJECTIVE(np.asarray(x))


def climb_chain(job):
    """The reference climber through P0 to P3, each climb starting at the last one's peak; then P2w,
    P2's domain climbed from IQ-TREE's relaxed endpoint at 1e-8,1e7, in case that domain holds a
    higher peak than the one the chain reaches. Each gradient's probes are spread over worker
    processes."""
    key, spec, pi, ref, z0, z_wide, max_iter, workers = job
    f = objective(spec, np.array(pi), np.array(ref))
    climbs = [(name, lo, hi, None) for name, lo, hi, _ in DOMAINS]
    if z_wide is not None:
        climbs.append(("P2w", DOMAIN["P2"][1], DOMAIN["P2"][2], z_wide))
    z, peaks = np.array(z0, dtype=float), []
    with ProcessPoolExecutor(max_workers=workers, mp_context=multiprocessing.get_context("spawn"),
                             initializer=_init_worker, initargs=(spec, pi, ref)) as pool:
        def evaluate(pts):
            return list(pool.map(_evaluate, pts, chunksize=max(1, len(pts) // (4 * workers))))

        for name, lo, hi, start in climbs:
            t0 = time.perf_counter()
            x0 = z if start is None else np.array(start, dtype=float)
            r = optimize.lbfgsb_fit(f, x0, np.full(len(x0), lo), np.full(len(x0), hi), PGTOL, max_iter,
                                    evaluate=evaluate)
            peaks.append({"peak": name, "lower": lo, "upper": hi, "lnl": -r.f, "x": r.x.tolist(),
                          "evaluations": r.evaluations, "phases": r.phases, "pg_max": r.pg_max,
                          "gradient": r.gradient.tolist(), "at_lower": r.at_lower.tolist(),
                          "at_upper": r.at_upper.tolist(), "wall_s": round(time.perf_counter() - t0, 1)})
            print(f"  oracle {key} {name}: lnL {-r.f:.6f}, pg {r.pg_max:.2e}, {r.phases}, "
                  f"{r.evaluations} evaluations, {peaks[-1]['wall_s']} s", flush=True)
            if start is None:
                z = r.x
    return key, peaks


# ---------------------------------------------------------------- reading IQ-TREE's runs

def by_row(indices):
    """Coordinate indices -> {amino acid of the row: count}."""
    out = {}
    for i in indices:
        row = builtin.STATES[int(i) // WIDTH]
        out[row] = out.get(row, 0) + 1
    return out


def outside(z, lo, hi):
    z = np.asarray(z)
    return {"below": by_row(np.flatnonzero(z < lo)), "above": by_row(np.flatnonzero(z > hi)),
            "min": float(z.min()), "max": float(z.max())}


def trajectory(run, mode):
    """The log-likelihood after each round, the initial one first, and the round count. The single
    path prints a round only when it gains more than -me; its last round's gain reaches the model
    but not the printed score."""
    text = run.stdout()
    if mode == "single":
        m, took = INITIAL.search(text), TOOK.search(text)
        vals = ([float(m.group(1))] if m else []) + [float(v) for v in CURRENT.findall(text)]
        return vals, int(took.group(1)) if took else None
    if mode == "-S":
        pairs = [float(v) for v in ROUND_S.findall(text)]
        return pairs[:1] + pairs[1::2], len(pairs) // 2
    m = INITIAL_Q.search(text)
    vals = ([float(m.group(1))] if m else []) + [float(v) for v in STEP_Q.findall(text)]
    return vals, len(vals) - 1


def ended_by(mode, rounds, gains, me, loop, nparam):
    """Which rule ended IQ-TREE's round loop (model/modelfactory.cpp, partitionmodel.cpp,
    partitionmodelplen.cpp)."""
    if rounds is None:
        return None
    if mode == "single":
        cap = rounds >= nparam - 1
    elif mode == "-S":
        cap = rounds >= loop and bool(gains) and gains[-1] >= 10 * me
    else:
        cap = rounds >= nparam - 1 and bool(gains) and abs(gains[-1]) >= me
    return "round cap" if cap else "round gain below the threshold"


def run_state(run, mode, nqc_model, setting, rounds_relaxed):
    """One IQ-TREE run's ending, final log-likelihood and trajectory. NQC's final is the class's
    accepted score (17 digits); NONREV's is BEST SCORE FOUND, which on the single path omits the last
    round's gain (below -me)."""
    if run.code != 0:
        return {"exit_code": run.code}
    me, loop, nparam = ((ME_DEFAULT, LOOP_DEFAULT, NPARAM_DEFAULT) if setting == "default"
                        else (ME_RELAXED, rounds_relaxed, rounds_relaxed + 2))
    vals, rounds = trajectory(run, mode)
    text = run.stdout()
    ck = final_struct(run) if nqc_model else None
    best = BEST.search(text)
    final = float(ck["accepted_score"]) if ck else (float(best.group(1)) if best else None)
    if final is not None and vals and abs(final - vals[-1]) > 1e-9:   # printed at 10 decimals
        vals = vals + [final]
    gains = list(np.diff(vals)) if len(vals) > 1 else []
    out = {"exit_code": 0, "final_lnl": final, "rounds": rounds,
           "ended_by": ended_by(mode, rounds, gains, me, loop, nparam), "wall_s": run.wall,
           "initial_lnl": vals[0] if vals else None, "gains": gains,
           "switches": text.count(SWITCH), "falls": len(FALL.findall(text)),
           "assertion_or_abort": "Assertion" in text or BUG_ABORT in text}
    if ck:
        z = floats(ck["coordinates"])
        out.update({"evaluations": int(ck["evaluations"]), "fits": int(ck["fits"]), "bounds": ck["bounds"],
                    "outside_default_box": outside(z, nqc.LOWER, nqc.UPPER), "coordinates": z.tolist(),
                    "fit_gains": [float(f[3]) - float(f[1]) for f in FIT.findall(text)]})
    return out


# ---------------------------------------------------------------- the analysis

def peak_report(peaks, pi, ref):
    """Per peak: the gain over the last, bound sets by row and side, the bound coordinates whose
    central gradient pushes outward (binding), collapsed-reference rows, and the largest change in
    log rate from P0 over rates above RATE_FLOOR."""
    Q0 = chart.build(np.array(peaks[0]["x"]), pi, ref)[0]
    off = ~np.eye(len(pi), dtype=bool)
    named = {p["peak"]: p for p in peaks}
    rows, prev = [], None
    for p in peaks:
        if p["peak"] == "P2w":
            rows.append(check("peak P2w: P2's domain climbed from IQ-TREE's relaxed endpoint there (report)", None,
                              lnl=p["lnl"], above_P2=p["lnl"] - named["P2"]["lnl"], pg_max=p["pg_max"],
                              phases=p["phases"], evaluations=p["evaluations"], wall_s=p["wall_s"],
                              outside_default_box=outside(p["x"], nqc.LOWER, nqc.UPPER)))
            continue
        x, g = np.array(p["x"]), np.array(p["gradient"])
        lower, upper = np.array(p["at_lower"], dtype=int), np.array(p["at_upper"], dtype=int)
        bind_lo, bind_hi = lower[g[lower] > 0], upper[g[upper] < 0]
        at_upper_rows = by_row(upper)
        Q = chart.build(x, pi, ref)[0]
        keep = off & ((Q0 > RATE_FLOOR) | (Q > RATE_FLOOR))
        dlog = np.zeros_like(Q)
        dlog[keep] = np.abs(np.log(Q[keep]) - np.log(Q0[keep]))
        row_change = {builtin.STATES[i]: float(dlog[i].max()) for i in np.argsort(-dlog.max(axis=1))[:5]}
        gain = None if prev is None else p["lnl"] - prev["lnl"]
        values = {"lnl": p["lnl"], "gain_over_previous": gain, "pg_max": p["pg_max"], "phases": p["phases"],
                  "evaluations": p["evaluations"], "wall_s": p["wall_s"],
                  "at_lower": by_row(lower), "at_upper": at_upper_rows,
                  "binding_lower": by_row(bind_lo), "binding_upper": by_row(bind_hi),
                  "largest_push": float(max(np.abs(g[np.concatenate([bind_lo, bind_hi])]), default=0.0)),
                  "collapsed_reference_rows": sorted(r for r, k in at_upper_rows.items() if k >= COLLAPSE),
                  "largest_log_rate_change_from_P0_by_row": row_change}
        if prev is not None and p["peak"] == "P3":
            values["past_P2"] = outside(x, prev["lower"], prev["upper"])
            values["flat"] = gain < FLAT
        rows.append(check(f"peak {p['peak']} [{p['lower']:.4f}, {p['upper']:.4f}] (report)", None, **values))
        if gain is not None:
            rows.append(check(f"peak {p['peak']}: no lower than the last peak by more than {-CONVERGED:g}",
                              gain >= CONVERGED, gain=gain))
        prev = p
    return rows


def case_checks(key, states, peaks, oracle_rows):
    label, mode, threshold = CASES[key]
    rows = list(oracle_rows)
    for name, s in states.items():
        values = {k: v for k, v in s.items() if k != "coordinates"}
        rows.append(check(f"{name} (report)", None, **values))
    box_d, box_r = states.get(f"{key}-nqc-box-default", {}), states.get(f"{key}-nqc-box-relaxed", {})
    if peaks:
        p0 = peaks[0]["lnl"]
        if box_d.get("final_lnl") is not None:
            rows.append(check(f"S-2: the box peak P0 above NQC's default box endpoint by more than {threshold:g}",
                              p0 - box_d["final_lnl"] > threshold, p0=p0, default=box_d["final_lnl"],
                              difference=p0 - box_d["final_lnl"]))
        if box_r.get("final_lnl") is not None:
            stalled = box_r["ended_by"] == "round gain below the threshold" and p0 - box_r["final_lnl"] > STALL
            rows.append(check(f"S-3: NQC's relaxed box run does not stall more than {STALL:g} below P0", not stalled,
                              p0=p0, relaxed=box_r["final_lnl"], ended_by=box_r["ended_by"],
                              difference=p0 - box_r["final_lnl"]))
        named = {q["peak"]: q for q in peaks}
        for peak, run in (("P0", f"{key}-nqc-box-relaxed"), ("P1", f"{key}-nqc-p1-relaxed"),
                          ("P1b", f"{key}-nqc-p1b-relaxed"),
                          ("P2w" if "P2w" in named else "P2", f"{key}-nqc-wide-relaxed")):
            if peak not in named:
                continue
            final = states.get(run, {}).get("final_lnl")
            rows.append(check(f"{run}: below the oracle's {peak} by (report)", None,
                              shortfall=None if final is None else named[peak]["lnl"] - final))
        p2 = max(named[k]["lnl"] for k in ("P2", "P2w") if k in named)
        rows.append(check(f"box cost, the higher of P2 and P2w minus P0, against the material line {MATERIAL:g} "
                          "(report; Peter, 2026-10-24)", None, p0=p0, p2=named["P2"]["lnl"],
                          p2w=named.get("P2w", {}).get("lnl"), cost=p2 - p0, material=p2 - p0 > MATERIAL))
    return rows


def nq_placement(c):
    """The built-in NQ matrices' jump chains in the chart (diagonals rebuilt as minus the row sums),
    against P0 to P2, with references from LG at each matrix's own pi and at the turtle target."""
    rows = []
    for name in NQ:
        Q0 = builtin.full_matrix(name)[0].copy()
        np.fill_diagonal(Q0, 0.0)
        np.fill_diagonal(Q0, -Q0.sum(axis=1))
        own = chart.solve_stationary(Q0)
        for label, ref in (("references from LG at its own pi", nqc.lg_seed(own)[1]),
                           ("references from LG at the turtle target", c.ref)):
            try:
                z = chart.inverse(Q0, ref)
            except chart.BuildError as e:
                rows.append(check(f"{name}, {label} (report)", None, error=str(e)))
                continue
            per = {d: {"below": int(np.sum(z < lo)), "above": int(np.sum(z > hi))}
                   for d, lo, hi, _ in DOMAINS if d != "P3"}
            rows.append(check(f"{name}, {label}: coordinates outside each domain (report)", None, domains=per,
                              rows_above_P0=by_row(np.flatnonzero(z > nqc.UPPER)), min=float(z.min()),
                              max=float(z.max())))
    return rows


# ---------------------------------------------------------------- run

def case_setup(key, c, r):
    """Files and argument builder for a case; C first fits NQC under -p in the default box for its
    nuisances."""
    if key == "A":
        return ({"input.tree": c.t02},
                lambda m: ["-s", AA, "-m", m + f"+G4{{{ALPHA}}}", "-te", "input.tree", "-blfix", "-v"],
                {"mode": "single", "tree": c.t02})
    if key == "B":
        files = turtle_files("-S")
        return (files, lambda m: joint_args(TURTLE, "-S", NEX, m, f"+G4{{{ALPHA}}}", "partition.trees", "-blfix", "-v"),
                {"mode": "-S", "trees": files["partition.trees"]})
    source = r.run("C-source-p-level2-g4", joint_args(TURTLE, "-p", NEX, c.model, "+G4", "input.tree", "-v"),
                   turtle_files("-p"))
    names = list(readers.read_charsets(NEX))
    route = q_route(source, names)
    if route is None:
        raise RuntimeError("the -p fit left out its tree, shapes or partition rates")
    return ({"parts.nex": q_charpartition(read(NEX), names, route), "input.tree": route["tree"]},
            lambda m: ["-s", TURTLE, "-q", "parts.nex", "--model-joint", m, "-te", "input.tree", "-blfix", "-v"],
            {"mode": "-q", "tree": route["tree"], "shapes": [float(a) for a in route["shapes"]],
             "rates": [float(x) for x in route["rates"]]})


def oracle_match(key, spec, c, run):
    """Asserted: the oracle's objective at an NQC endpoint equals IQ-TREE's accepted score (decision 021)."""
    ck = final_struct(run) if run.code == 0 else None
    name = f"{run.name}: oracle -lnL at the endpoint against the accepted score (decision 021)"
    if not ck:
        return check(name, False, reason="no NQC endpoint", exit_code=run.code)
    oracle = -objective(spec, c.pi, c.ref)(floats(ck["coordinates"]))
    accepted = float(ck["accepted_score"])
    return check(name, abs(accepted - oracle) <= TOL_LNL * abs(oracle), accepted=accepted, oracle=oracle)


def new_out(path):
    out = Path(path).expanduser().resolve()
    if out.exists() or out.is_relative_to(REPO):
        sys.exit(f"{out}: needs a new directory outside the repository")
    return out


def provenance(a):
    return {"script": SCRIPT, "command": sys.argv,
            "started_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
            "git_commit": capture(a.git, "rev-parse", "HEAD", cwd=REPO),
            "git_clean": capture(a.git, "status", "--porcelain", cwd=REPO) == "",
            "python": sys.version, "numpy": np.__version__, "scipy": scipy.__version__,
            "platform": platform.platform(), "label": "exploratory"}


def markdown(title, prov, results):
    lines = [f"# {title}", "", f"Exploratory. Generated by `{SCRIPT}` on {prov['started_utc']} at commit "
             f"`{prov['git_commit']}` (clean: {prov['git_clean']}).", ""]
    for case, rows in results.items():
        lines += [f"## {case}", "", "| Check | Result | Values |", "|---|---|---|"]
        for row in rows:
            verdict = {True: "pass", False: "FAIL", None: "report"}[row["passed"]]
            values = {k: v for k, v in row.items() if k not in ("check", "passed")}
            lines.append(f"| {row['check']} | {verdict} | `{json.dumps(values, default=str)[:600]}` |")
        lines.append("")
    return "\n".join(lines) + "\n"


def cmd_run(a):
    out = new_out(a.out)
    binary = Path(a.binary).expanduser().resolve()
    prov = provenance(a)
    if not prov["git_clean"] and not a.allow_dirty:
        sys.exit(f"the repository is not clean ({prov['git_commit']}); commit first or pass --allow-dirty")
    keys = a.cases.split(",")
    if set(keys) - set(CASES):
        sys.exit(f"unknown cases: {sorted(set(keys) - set(CASES))}")
    relaxed = ["-me", f"{ME_RELAXED:g}", "--loop-model", str(a.rounds), "-nparam", str(a.rounds + 2)]
    prov.update({"binary": {"path": str(binary), "sha256": sha256(binary)}, "cases": keys, "seed": "IQ-TREE -seed 1",
                 "settings": {"relaxed": relaxed, "rounds": a.rounds, "max_iter": a.max_iter, "pgtol": PGTOL,
                              "material": MATERIAL, "stall": STALL, "flat": FLAT, "converged": CONVERGED,
                              "collapse": COLLAPSE, "rate_floor": RATE_FLOOR, "domains": DOMAINS,
                              "jobs": a.jobs, "workers": a.workers, "threads": "-T 1"}})
    out.mkdir(parents=True)
    r, c = Runner(binary, out, a.timeout), Context()
    setups, jobs = {}, []
    for key in keys:
        files, args, spec = case_setup(key, c, r)
        setups[key] = spec
        for config, is_nqc, domain in MODELS + (("nqc-p1", True, DOMAIN["P1"][3]), ("nqc-p1b", True, DOMAIN["P1b"][3])):
            model = c.model if is_nqc else "NONREV"
            extra = ["--nqc-domain", domain] if domain else []
            for setting in (("relaxed",) if config in ("nqc-p1", "nqc-p1b") else ("default", "relaxed")):
                jobs.append((f"{key}-{config}-{setting}", args(model) + extra + (relaxed if setting == "relaxed" else []),
                             files, key, is_nqc, setting))
    print(f"{len(jobs)} IQ-TREE runs, {a.jobs} at a time", flush=True)
    with ThreadPoolExecutor(max_workers=a.jobs) as pool:
        runs = dict(zip([j[0] for j in jobs], pool.map(lambda j: r.run(j[0], j[1], j[2]), jobs)))
    states = {key: {} for key in keys}
    for name, _, _, key, is_nqc, setting in jobs:
        states[key][name] = run_state(runs[name], CASES[key][1], is_nqc, setting, a.rounds)
    results, chains, oracle_rows = {}, [], {key: [] for key in keys}
    for key in keys:
        for setting in ("default", "relaxed"):
            oracle_rows[key].append(oracle_match(key, setups[key], c, runs[f"{key}-nqc-box-{setting}"]))
        start = states[key][f"{key}-nqc-box-relaxed"].get("coordinates")
        wide = states[key][f"{key}-nqc-wide-relaxed"].get("coordinates")
        if start is not None:
            chains.append([key, setups[key], c.pi.tolist(), [int(v) for v in c.ref], start, wide, a.max_iter])
    for job in chains:
        job.append(max(1, a.workers // len(chains)))
    print(f"{len(chains)} oracle chains, {a.workers // max(1, len(chains))} worker processes each", flush=True)
    with ThreadPoolExecutor(max_workers=max(1, len(chains))) as pool:
        peaks = dict(pool.map(climb_chain, chains))
    for key in keys:
        try:
            rows = case_checks(key, states[key], peaks.get(key), oracle_rows[key])
            if key in peaks:
                rows += peak_report(peaks[key], c.pi, c.ref)
        except Exception:  # one case's crash is recorded and the others still report
            rows = [check("case analysed", False, exception=traceback.format_exc())]
        results[f"{key}: {CASES[key][0]}"] = rows
    s1 = {key: (states[key].get(f"{key}-nonrev-relaxed", {}).get("final_lnl"),
                states[key].get(f"{key}-nonrev-default", {}).get("final_lnl")) for key in keys}
    gains = {key: (v[0] - v[1] if None not in v else None) for key, v in s1.items()}
    results["S-1: NONREV stops short at default settings"] = [check(
        "S-1: NONREV's relaxed run ends above its default run by more than the case's round threshold in "
        "at least one case", any(g is not None and g > CASES[k][2] for k, g in gains.items()),
        relaxed_minus_default=gains, thresholds={k: CASES[k][2] for k in keys},
        note="single path: BEST SCORE FOUND omits the last round's gain, below -me")]
    results["NQ matrices placed against the domains"] = nq_placement(c)
    prov["finished_utc"] = datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")
    write(out / "results.json", json.dumps({"provenance": prov, "runs": r.runs, "states": states, "peaks": peaks,
                                            "results": results}, indent=1, default=str))
    write(out / "results.md", markdown("Stopping study and converged risk 22 evidence", prov, results))
    failed = [row["check"] for rows in results.values() for row in rows if row["passed"] is False]
    for f in failed:
        print(f"  FAIL {f}")
    print(f"{len(failed)} asserted checks failed; wrote {out / 'results.md'}")
    sys.exit(1 if failed else 0)


# ---------------------------------------------------------------- reread

def dfpmin_test(text):
    """dfpmin's gtol test (utils/optimization.cpp:849-855) at the first printed gradient."""
    got = first_gradient(text)
    if not got:
        return None
    f, x, _, g = got
    return float(np.max(np.abs(g) * np.maximum(np.abs(x), 1.0)) / max(abs(f), 1.0))


def cmd_reread(a):
    gate, probe, out = Path(a.gate).expanduser().resolve(), Path(a.probe).expanduser().resolve(), new_out(a.out)
    prov = provenance(a)
    if not prov["git_clean"] and not a.allow_dirty:
        sys.exit(f"the repository is not clean ({prov['git_commit']}); commit first or pass --allow-dirty")
    prov.update({"gate": str(gate), "probe": str(probe)})
    used, results = {}, {}
    rows, finals = [], {}
    for side in ("box", "wide"):
        d = gate / "nqc" / f"warm-turtle-level-2--p-g4-(-q-route)-{side}"
        run = Run(d, d.name, [], 0, 0.0)
        used[d.name] = sha256(d / "stdout.txt")
        vals, rounds = trajectory(run, "-q")
        ck = final_struct(run)
        finals[side] = float(ck["accepted_score"])
        rows.append(check(f"turtle -q continuation, {side} (report)", None, rounds=rounds,
                          gains=list(np.diff(vals)), total_gain=vals[-1] - vals[0], final=finals[side],
                          evaluations=int(ck["evaluations"]), fits=int(ck["fits"]),
                          evaluations_per_fit=int(ck["evaluations"]) / int(ck["fits"]),
                          dfpmin_test_at_first_gradient=dfpmin_test(run.stdout()), gtol=GTOL))
    rows.append(check("turtle -q continuations: wide final minus box final (report)", None,
                      difference=finals["wide"] - finals["box"]))
    results["turtle -q continuations (S3b gate)"] = rows
    run = Run(probe / "c2", "c2", [], 0, 0.0)
    used["c2"] = sha256(probe / "c2" / "stdout.txt")
    vals, rounds = trajectory(run, "-S")
    results["S0 probe c2: NONREV, turtle -S Level 1"] = [check("c2 per-round gains (report)", None, rounds=rounds,
                                                               gains=list(np.diff(vals)), final=vals[-1])]
    rows = []
    for d in sorted(p for p in (gate / "nqc").iterdir() if p.is_dir() and (p / "probe.ckp.gz").exists()):
        ck = final_struct(Run(d, d.name, [], 0, 0.0))
        if ck and "coordinates" in ck:
            dom = floats(ck["domain"])
            rows.append(check(f"{d.name} (report)", None, domain=dom.tolist(),
                              **outside(floats(ck["coordinates"]), nqc.LOWER, nqc.UPPER)))
    results["final coordinates outside the default box, by row (S3b gate)"] = rows
    prov["inputs_sha256"] = used
    prov["finished_utc"] = datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")
    out.mkdir(parents=True)
    write(out / "reread.json", json.dumps({"provenance": prov, "results": results}, indent=1, default=str))
    write(out / "reread.md", markdown("S3b's risk 22 numbers, re-read", prov, results))
    print(f"wrote {out / 'reread.md'}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="command", required=True)
    p = sub.add_parser("run")
    p.add_argument("--binary", required=True)
    p.add_argument("--out", required=True, help="new directory outside the repository")
    p.add_argument("--cases", default="A,B,C")
    p.add_argument("--rounds", type=int, default=300, help="the relaxed runs' round cap")
    p.add_argument("--max-iter", type=int, default=2000, help="L-BFGS-B iterations per pass")
    p.add_argument("--jobs", type=int, default=7, help="IQ-TREE runs at a time")
    p.add_argument("--workers", type=int, default=11, help="oracle worker processes, shared among the chains")
    p.add_argument("--timeout", type=float, default=7200, help="seconds per IQ-TREE run")
    q = sub.add_parser("reread")
    q.add_argument("--gate", required=True, help="the S3b gate's run directory")
    q.add_argument("--probe", required=True, help="the S0 probe record directory holding c2")
    q.add_argument("--out", required=True, help="new directory outside the repository")
    for x in (p, q):
        x.add_argument("--git", default="git", help="git executable (git.exe under WSL)")
        x.add_argument("--allow-dirty", action="store_true")
    a = ap.parse_args()
    cmd_run(a) if a.command == "run" else cmd_reread(a)


if __name__ == "__main__":
    main()
