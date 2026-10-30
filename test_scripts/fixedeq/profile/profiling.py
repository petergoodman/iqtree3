"""S3c's profiling: NQC against --model-joint NONREV (CODE_PLAN.md section 4, S3c; PLAN.md risks 5, 18
and 24). The thresholds were approved on 2026-10-22 before any measurement; the scaled design, the
measurement and the simulation model by Peter on 2026-10-25.

  python test_scripts/fixedeq/profile/profiling.py simulate --alisim FROZEN --out DATA --git git.exe
  python test_scripts/fixedeq/profile/profiling.py windows --binary BIN --data DATA --out DIR --git git.exe
  python test_scripts/fixedeq/profile/profiling.py fits --binary BIN --data DATA --windows DIR --out DIR \
      --git git.exe [--only KEYS]

simulate: two collections of 100 alignments, 20 taxa and 500 sites each, simulated by the frozen
binary's AliSim from NQ.pfam's jump chain at decision 023's target with +G4{0.5}. Under -S each
alignment has its own seeded rooted tree; under -p all share one tree T, each with its branch
lengths times a partition rate r_k (log-uniform on [0.5, 2], mean 1 over the 100). The 10-, 30- and
100-sets are the first alignments of each collection, written concatenated with their charsets,
their trees, and the -q partition file holding the true shape and rates. Since IQ-TREE rescales
given rates to mean 1 without rescaling the tree (tree/phylosupertreeplen.cpp, normalizePartRate),
a -p set of n alignments gets the tree T * m_n and the rates r_k / m_n, m_n the mean of its r_k.

windows: time per linked evaluation. gdb stops at PartitionModel::targetFunk, called once per linked
evaluation for both models, reads /proc/uptime at evaluation 1 and 1001 and ends the run. Level 1
+G4 under -S and -q on turtle and the three sets at -T 1, and on the 10- and 100-sets at -T 2 and -T
4. First the counter's checks: NONREV's per-partition calls (a second breakpoint, on
ModelMarkov::targetFunk) against its linked ones, and the cost of a gdb stop from plain and counted
NQC runs.

fits: full fits at default settings, each run under gdb counting its evaluations: turtle Level 1 +G4
under -S and -q and Level 2 under -S and -p with +G4, +I+G4 and free +R4; the 10-set Level 1 +G4
under -S and -q and Level 2 +G4 under -S and -p. Per-fit ratios, split into the per-evaluation ratio
(from the windows) and the evaluation-count ratio, against the thresholds; the labelled
extrapolation to real training.

Every IQ-TREE run uses -seed 1, -v and NQC at its default domain, and runs alone, NONREV and NQC in
turn; turtle's -q inputs come from NQC's -p Level 2 +G4 fit (decision 024), run first. Writes
results.json and results.md (data.json for simulate) in the new directory, outside the repository.
"""
import argparse
import datetime
import json
import platform
import re
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
for _p in ("", "nqc", "probes", "regression", "differential", "explore"):
    sys.path.insert(0, str(HERE.parent / _p))
from oracle import readers  # noqa: E402
from probes import NEX, PREFIX, TURTLE, Run, Runner, turtle_target, with_charpartition  # noqa: E402
from regress import REPO, SCALARS, capture, read, sha256, write  # noqa: E402
from nqcchecks import (ALPHA, BUG_ABORT, CPU, FALL, SWITCH, WALL, check, final_struct, floats,  # noqa: E402
                       joint_args, one_round, q_charpartition, q_route, turtle_files)
from convergence import LOOP_DEFAULT, ME_DEFAULT, NPARAM_DEFAULT, by_row, ended_by, new_out, trajectory  # noqa: E402
import toy  # noqa: E402
import twopart  # noqa: E402

SCRIPT = Path(__file__).resolve().relative_to(REPO).as_posix()
PER_EVAL_LINE, PER_FIT_LINE, BLOCK_LINE = 1.25, 2.0, 3.0   # Peter, 2026-10-22
WINDOW, PAIRS = 1000, 3   # evaluations timed per window; plain and counted overhead pairs
IGNORE_ALL = 2_000_000_000
TAXA, SITES, COLLECTION, SETS = 20, 500, 100, (10, 30, 100)
RATE_RANGE = (0.5, 2.0)
SEED_S_TREE, SEED_S_ALN = 2026102600, 2026102700   # plus the alignment's number, 1 to 100
SEED_P_TREE, SEED_P_ALN, SEED_RATES = 20261028, 2026102900, 20261030
REPEATS = {"turtle": 3, "10": 3, "30": 1, "100": 1}
THREADED = ((10, 2), (10, 4), (100, 2), (100, 4))
MIN_RATE, MAX_RATE = 1e-4, 100.0   # model/modelmarkov.h, NONREV's bounds
# Normal runs keep 0.92 to 0.93 of wall time as CPU under gdb; a development run paused for 55 s kept 0.44
CPU_FLOOR, CLOCK_TICK = 0.75, 100   # CLOCK_TICK: getconf CLK_TCK on this Debian
STATES = "ARNDCQEGHILKMFPSTWYV"
COUNTED = "PartitionModel::targetFunk(double*)"
PER_PARTITION = "ModelMarkov::targetFunk(double*)"
FIT_START = "PartitionModel::optimizeLinkedModel(bool, double)"   # once per linked fit, so once a round
GDB_OPTIONS = ["-nx", "-batch", "-iex", "set auto-load python-scripts off", "-iex", "set debuginfod enabled off"]
UPTIME = "uptime.txt"
RAM = re.compile(r"^NOTE: (\d+) MB RAM", re.M)
KERNEL = re.compile(r"^Kernel: .* - (\d+) threads", re.M)


# ---------------------------------------------------------------- gdb

def gdb_script(mode, extra=(), window=WINDOW):
    """gdb's own command language (decision 013 keeps project Python out of gdb's embedded system
    Python). The shell start stays on: with it off, gdb keeps the backslashes it puts before braces.
    In count mode each linked fit's start prints the evaluations so far, for evaluations per round."""
    head = ["set pagination off", "set confirm off", "set breakpoint pending off",
            "handle SIGABRT SIGSEGV SIGFPE SIGBUS SIGILL nostop noprint pass"]
    if mode == "window":
        # each reading: seconds of uptime, then IQ-TREE's user and system CPU ticks (runs go one at a time)
        reading = (f'  shell echo $(cut -d" " -f1 /proc/uptime) '
                   f'$(cut -d" " -f14,15 /proc/$(pgrep -n -x iqtree3)/stat) >> {UPTIME}')
        body = ["set $profile_n = 0", f"break '{COUNTED}'", "commands 1", "  silent",
                "  set $profile_n = $profile_n + 1", reading,
                "  if $profile_n < 2", f"    ignore 1 {window - 1}", "    continue", "  end", "end", "run"]
    else:
        body = []
        for i, name in enumerate((COUNTED, *extra), start=1):
            body += [f"break '{name}'", f"ignore {i} {IGNORE_ALL}"]
        body += [f"break '{FIT_START}'", f"commands {len(extra) + 2}", "  silent", "  echo PROFILE-FIT\\n",
                 "  info breakpoints 1", "  continue", "end", "run"]
    tail = ["echo PROFILE-EXITCODE\\n", "print $_exitcode", "echo PROFILE-EXITSIGNAL\\n", "print $_exitsignal",
            "echo PROFILE-BREAKPOINTS\\n", "info breakpoints"] + (["kill"] if mode == "window" else [])
    return "\n".join(head + body + tail) + "\n"


HIT = re.compile(r"^\s+breakpoint already hit (\d+) times?")


def parse_gdb(text):
    """The exit code and signal gdb printed (None when void), each breakpoint's hit count, and the
    evaluations counted at the start of each linked fit."""
    def value(marker):
        m = re.search(re.escape(marker) + r"\n\$\d+ = (\S+)", text)
        return int(m.group(1)) if m and m.group(1).lstrip("-").isdigit() else None
    lines = text.splitlines()
    starts = []
    for i, line in enumerate(lines):
        if line.endswith("PROFILE-FIT"):
            m = next((HIT.match(x) for x in lines[i + 1:i + 5] if HIT.match(x)), None)
            starts.append(int(m.group(1)) if m else 0)
    hits, current = {}, None
    section = text.split("PROFILE-BREAKPOINTS\n", 1)[1] if "PROFILE-BREAKPOINTS\n" in text else ""
    for line in section.splitlines():
        m = re.match(r"^\d+\s+breakpoint\s+keep\s+\S+\s+\S+\s+<(.+?)(\+\d+)?>", line)
        if m:
            current = m.group(1)
            hits[current] = 0
            continue
        m = HIT.match(line)
        if m and current:
            hits[current] = int(m.group(1))
    return {"exit_code": value("PROFILE-EXITCODE"), "exit_signal": value("PROFILE-EXITSIGNAL"), "hits": hits,
            "fit_starts": starts}


def evaluations_per_fit(starts, total):
    """Each linked fit's evaluations, from the counts at the fits' starts and the run's total."""
    return [b - a for a, b in zip(starts, starts[1:] + [total])] if starts and total is not None else []


def window_seconds(uptime_text):
    """Seconds between the two /proc/uptime readings (10 ms resolution), or None without both."""
    vals = [float(line.split()[0]) for line in uptime_text.splitlines() if line.strip()]
    return vals[1] - vals[0] if len(vals) == 2 else None


def window_cpu(uptime_text, tick=CLOCK_TICK):
    """IQ-TREE's CPU seconds between the two readings, from its /proc stat ticks, or None."""
    rows = [line.split() for line in uptime_text.splitlines() if line.strip()]
    if len(rows) != 2 or any(len(r) != 3 for r in rows):
        return None
    return (int(rows[1][1]) + int(rows[1][2]) - int(rows[0][1]) - int(rows[0][2])) / tick


def cpu_row(name, cpu, wall, threads):
    """Asserted at -T 1, set before any gate run: a timed run whose CPU time is below CPU_FLOOR of its
    wall time was paused or shared the CPU, so its time is not a measurement; it is rerun, and the
    rerun recorded. Reported above -T 1."""
    frac = cpu / wall if cpu is not None and wall else None
    return check(f"{name}: CPU time at least {CPU_FLOOR} of wall time (not paused or sharing the CPU)"
                 + ("" if threads == 1 else " (report above -T 1)"),
                 bool(frac is not None and frac >= CPU_FLOOR) if threads == 1 else None, cpu_s=cpu, wall_s=wall,
                 fraction=frac)


def run_code(mode, g, window, timed_out):
    """IQ-TREE's outcome as gdb saw it: 0 for a complete window or a normal exit, else the reason."""
    if timed_out:
        return "timeout"
    if mode == "window":
        return 0 if g.get("window_s") is not None and g["hits"].get(COUNTED) == window + 1 else "window incomplete"
    if g["exit_signal"] is not None:
        return f"signal {g['exit_signal']}"
    return g["exit_code"] if g["exit_code"] is not None else "no exit code"


class GdbRun(Run):
    """A run under gdb: IQ-TREE's own output is its log, since gdb writes into stdout too."""

    def stdout(self):
        log = self.file(".log")
        return read(log) if log.exists() else read(self.dir / "stdout.txt")

    def raw(self):
        return read(self.dir / "stdout.txt")


class GdbRunner(Runner):
    """probes.Runner under gdb. gdb gets its own session so that a timeout, an interrupt or an error
    kills it and IQ-TREE together, and the run's code is IQ-TREE's exit code as gdb reports it."""

    def __init__(self, binary, out, timeout, gdb="gdb", window=WINDOW):
        super().__init__(binary, out, timeout)
        self.gdb, self.window = gdb, window

    def plain(self, name, args, files=None, threads=1):
        return super().run(name, args, files, threads)

    def run(self, name, args, files=None, threads=1, mode="count", extra=()):
        rundir = self.out / name
        rundir.mkdir()
        for fname, content in (files or {}).items():
            write(rundir / fname, content)
        write(rundir / "profile.gdb", gdb_script(mode, extra, self.window))
        argv = [self.gdb, *GDB_OPTIONS, "-x", "profile.gdb", "--args", str(self.binary), *[str(a) for a in args],
                "-seed", "1", "-T", str(threads), "--prefix", PREFIX]
        t0, timed_out = time.perf_counter(), False
        with open(rundir / "stdout.txt", "w", encoding="utf-8") as fh:
            proc = subprocess.Popen(argv, cwd=rundir, stdout=fh, stderr=subprocess.STDOUT, start_new_session=True)
            try:
                proc.wait(timeout=self.timeout)
            except subprocess.TimeoutExpired:
                timed_out = True
            finally:
                if proc.poll() is None:
                    subprocess.run(["pkill", "-KILL", "-s", str(proc.pid)])
                    proc.wait()
        wall = round(time.perf_counter() - t0, 1)
        g = parse_gdb(read(rundir / "stdout.txt"))
        g["mode"] = mode
        if mode == "window":
            up = rundir / UPTIME
            g["window_s"] = window_seconds(read(up)) if up.exists() else None
            g["window_cpu_s"] = window_cpu(read(up)) if up.exists() else None
        code = run_code(mode, g, self.window, timed_out)
        r = GdbRun(rundir, name, [str(a) for a in args], code, wall, threads)
        r.gdb = g
        self.runs.append({**r.record(), "gdb": g})
        print(f"  {name}: {code}, {wall} s, hits {g['hits'].get(COUNTED)}", flush=True)
        return r


# ---------------------------------------------------------------- simulated data

def scale_tree(text, factor):
    """Every branch length times factor, at 17 significant digits."""
    return re.sub(r":([0-9.eE+-]+)", lambda m: f":{float(m.group(1)) * factor:.17g}", text)


def partition_rates():
    """The -p collection's rates: seeded log-uniform on RATE_RANGE, normalized to mean 1."""
    rng = np.random.default_rng(SEED_RATES)
    r = np.exp(rng.uniform(np.log(RATE_RANGE[0]), np.log(RATE_RANGE[1]), COLLECTION))
    return r / r.mean()


def set_rates(rates, n):
    """A -p set's tree factor m_n and its -q rates r_k / m_n, which IQ-TREE's renormalization leaves
    unchanged and whose products with T * m_n are the simulated branch lengths."""
    m = float(np.mean(rates[:n]))
    return m, [float(x) / m for x in rates[:n]]


def nexus_sets(n, sites=SITES):
    lines = [f"charset p{k:03d} = {(k - 1) * sites + 1}-{k * sites};" for k in range(1, n + 1)]
    return "#nexus\nbegin sets;\n" + "\n".join(lines) + "\nend;\n"


def concatenate(alignments):
    """Concatenated PHYLIP text of alignments sharing the taxa T1..Tn."""
    taxa = sorted(alignments[0], key=lambda name: int(name[1:]))
    length = sum(len(a[taxa[0]]) for a in alignments)
    return f"{len(taxa)} {length}\n" + "".join(f"{t} {''.join(a[t] for a in alignments)}\n" for t in taxa)


def cmd_simulate(a):
    out = new_out(a.out)
    prov = provenance(a)
    if not prov["git_clean"] and not a.allow_dirty:
        sys.exit(f"the repository is not clean ({prov['git_commit']}); commit first or pass --allow-dirty")
    alisim = Path(a.alisim).expanduser().resolve()
    prov["alisim"] = {"path": str(alisim), "sha256": sha256(alisim)}
    out.mkdir(parents=True)
    pi = np.array([float(v) for v in turtle_target().split(",")])
    Q, rate = toy.truth(pi), f"+G4{{{ALPHA}}}"
    rates, shared = partition_rates(), toy.tree_text(SEED_P_TREE, TAXA)
    taxa = {f"T{i + 1}" for i in range(TAXA)}
    made, found, alns = {"S": [], "P": []}, [], {"S": [], "P": []}
    for coll in ("S", "P"):
        d = out / coll
        d.mkdir()
        for k in range(1, COLLECTION + 1):
            name = f"{coll.lower()}{k:03d}"
            tree = (toy.tree_text(SEED_S_TREE + k, TAXA) if coll == "S" else scale_tree(shared, rates[k - 1]))
            seed = (SEED_S_ALN if coll == "S" else SEED_P_ALN) + k
            made[coll].append({"name": name, "seed": seed, **toy.simulate(d, alisim, name, tree, Q, pi, SITES, seed, rate)})
            aln = readers.read_phylip(d / f"{name}.phy")
            if set(aln) != taxa:
                found.append(f"{name}: taxa {sorted(taxa ^ set(aln))} missing or extra")
            found += [f"{name}: {p}" for p in twopart.problems(aln)]
            alns[coll].append(aln)
            print(f"  {name}", end="\r", flush=True)
    if found:
        sys.exit("simulated data unusable: " + "; ".join(found))
    write(out / "P" / "tree.nwk", shared)
    files, factors = {}, {}
    for n in SETS:
        for coll in ("S", "P"):
            stem = f"{coll}{n}"
            write(out / f"{stem}.phy", concatenate(alns[coll][:n]))
            write(out / f"{stem}.nex", nexus_sets(n))
            files[stem] = [f"{stem}.phy", f"{stem}.nex"]
        trees = [read(out / "S" / f"s{k:03d}.nwk").strip() for k in range(1, n + 1)]
        write(out / f"S{n}.trees", "\n".join(trees) + "\n")
        m, q = set_rates(rates, n)
        factors[n] = m
        write(out / f"P{n}.tree", scale_tree(shared, m))
        write(out / f"P{n}-q.nex", with_charpartition(nexus_sets(n), [(f"LG+G4{{{ALPHA}}}", f"p{k:03d}", x)
                                                                      for k, x in enumerate(q, start=1)]))
        files[f"S{n}"].append(f"S{n}.trees")
        files[f"P{n}"] += [f"P{n}.tree", f"P{n}-q.nex"]
    prov["finished_utc"] = now()
    manifest = {"provenance": prov, "settings": {"taxa": TAXA, "sites": SITES, "collection": COLLECTION, "sets": SETS,
                                                 "rate": rate, "rate_range": RATE_RANGE,
                                                 "seeds": {"S_tree": SEED_S_TREE, "S_alignment": SEED_S_ALN,
                                                           "P_tree": SEED_P_TREE, "P_alignment": SEED_P_ALN,
                                                           "rates": SEED_RATES}},
                "rates": rates.tolist(), "tree_factors": factors, "alignments": made,
                "files": {stem: {f: sha256(out / f) for f in fs} for stem, fs in files.items()}}
    write(out / "data.json", json.dumps(manifest, indent=1, default=str))
    print(f"\nwrote {out / 'data.json'}")


# ---------------------------------------------------------------- configurations

class Data:
    """Inputs of the turtle and simulated configurations; the simulated files checked against
    simulate's SHA-256 record."""

    def __init__(self, data_dir):
        self.dir = Path(data_dir).expanduser().resolve()
        self.manifest = json.loads(read(self.dir / "data.json"))
        self.target_text = turtle_target()
        self.model = f"NQC+F{{{self.target_text}}}"
        self.names = list(readers.read_charsets(NEX))
        self.route = None

    def verify(self):
        bad = [f for files in self.manifest["files"].values() for f, h in files.items() if sha256(self.dir / f) != h]
        return check("simulated data match simulate's SHA-256 record", not bad, mismatched=bad,
                     data_json_sha256=sha256(self.dir / "data.json"))

    def level1(self, dataset, arr, model):
        """(args, files) of a Level 1 +G4 run, -S or -q."""
        if dataset == "turtle":
            if arr == "-S":
                return (joint_args(TURTLE, "-S", NEX, model, f"+G4{{{ALPHA}}}", "partition.trees", "-blfix", "-v"),
                        turtle_files("-S"))
            return (["-s", TURTLE, "-q", "parts.nex", "--model-joint", model, "-te", "input.tree", "-blfix", "-v"],
                    {"parts.nex": q_charpartition(read(NEX), self.names, self.route), "input.tree": self.route["tree"]})
        d, n = self.dir, int(dataset)
        if arr == "-S":
            return joint_args(d / f"S{n}.phy", "-S", d / f"S{n}.nex", model, f"+G4{{{ALPHA}}}", d / f"S{n}.trees",
                              "-blfix", "-v"), {}
        return ["-s", d / f"P{n}.phy", "-q", d / f"P{n}-q.nex", "--model-joint", model, "-te", d / f"P{n}.tree",
                "-blfix", "-v"], {}

    def level2(self, dataset, arr, rate, model):
        """(args, files) of a Level 2 run from the given trees, -S or -p."""
        if dataset == "turtle":
            return joint_args(TURTLE, arr, NEX, model, rate, "partition.trees" if arr == "-S" else "input.tree",
                              "-v"), turtle_files(arr)
        d, n = self.dir, int(dataset)
        stem = f"S{n}" if arr == "-S" else f"P{n}"
        trees = d / (f"S{n}.trees" if arr == "-S" else f"P{n}.tree")
        return joint_args(d / f"{stem}.phy", arr, d / f"{stem}.nex", model, rate, trees, "-v"), {}


def source_fit(r, data):
    """NQC's turtle -p Level 2 +G4 fit, counted; its tree, shapes and partition rates are turtle's -q
    Level 1 inputs (decision 024)."""
    args, files = data.level2("turtle", "-p", "+G4", data.model)
    run = r.run("source-turtle-L2-p-G4-nqc", args, files)
    data.route = q_route(run, data.names)
    if data.route is None:
        sys.exit("the turtle -p Level 2 fit left out its tree, shapes or partition rates")
    return run


def models(data):
    return (("nonrev", "NONREV"), ("nqc", data.model))


# ---------------------------------------------------------------- reading a counted fit

def wall_cpu(run):
    rep = run.report()
    w, c = WALL.search(rep), CPU.search(rep)
    return (float(w.group(1)) if w else None), (float(c.group(1)) if c else None)


NONREV_RATES = re.compile(r"^PartitionModel(Plen)?!.+!ModelProtein!rates$")   # -p and -q name the struct ...Plen


def nonrev_bounds(ck):
    """NONREV's linked rates at its bounds, by row (rates are row-major off-diagonals,
    ModelMarkov::decomposeRateMatrixNonrev), from a flat checkpoint's 10-digit values."""
    vals = next((v for k, v in ck.items() if NONREV_RATES.match(k)), None)
    if vals is None:
        return None
    r = np.array([float(v) for v in vals.split(",")])
    low, high = np.flatnonzero(r <= MIN_RATE * (1 + 1e-9)), np.flatnonzero(r >= MAX_RATE * (1 - 1e-9))
    rows = lambda idx: {STATES[int(i) // 19]: int(np.sum(idx // 19 == i // 19)) for i in idx}  # noqa: E731
    return {"at_lower": rows(low), "at_upper": rows(high), "lower": len(low), "upper": len(high)}


def fit_state(run, arr, is_nqc):
    """One counted full fit: time, rounds and ending rule, evaluations in total and per linked fit,
    scores, parameters, memory, bounds and the monitors."""
    text = run.stdout()
    wall, cpu = wall_cpu(run)
    tmode = "-S" if arr == "-S" else "-q"
    vals, rounds = trajectory(run, tmode)
    gains = list(np.diff(vals)) if len(vals) > 1 else []
    params = re.search(SCALARS["free_parameters"], run.report(), re.M)
    ram = RAM.search(text)
    hits = run.gdb["hits"].get(COUNTED)
    out = {"exit_code": run.code, "complete": run.code == 0 and hits is not None, "wall_s": wall, "cpu_s": cpu,
           "driver_wall_s": run.wall, "hits": hits,
           "per_fit_evaluations": evaluations_per_fit(run.gdb["fit_starts"], hits),
           "rounds": rounds, "ended_by": ended_by(tmode, rounds, gains, ME_DEFAULT, LOOP_DEFAULT, NPARAM_DEFAULT),
           "final_lnl": vals[-1] if vals else None, "initial_lnl": vals[0] if vals else None,
           "per_round_s": wall / rounds if wall and rounds else None,
           "per_hit_s": wall / hits if wall and hits else None,
           "free_parameters": int(params.group(1)) if params else None, "ram_mb": int(ram.group(1)) if ram else None,
           "switches": text.count(SWITCH), "falls": len(FALL.findall(text)),
           "assertion_or_abort": any(s in t for s in ("Assertion", BUG_ABORT) for t in (text, run.raw())),
           "gains": gains}
    if is_nqc:
        ck = final_struct(run) or {}
        out["evaluations_key"] = int(ck["evaluations"]) if "evaluations" in ck else None
        out["accepted_score"] = ck.get("accepted_score")
        if "coordinates" in ck:
            z, (lo, hi) = floats(ck["coordinates"]), floats(ck["domain"])
            out["bounds"] = {"at_lower": by_row(np.flatnonzero(z <= lo)), "at_upper": by_row(np.flatnonzero(z >= hi)),
                             "key": ck.get("bounds")}
    else:
        out["bounds"] = nonrev_bounds(run.checkpoint())
    return out


# ---------------------------------------------------------------- ratios and verdicts

def median(xs):
    xs = [x for x in xs if x is not None]
    return float(np.median(xs)) if xs else None


def ratio(a, b):
    return a / b if a is not None and b else None


def without_stop(x, stop):
    """A time with the measured cost of its gdb stops removed (reported beside the raw figure)."""
    return x - stop if x is not None and stop is not None else None


def encode_key(key):
    return "|".join(map(str, key))


def decode_key(text):
    """(dataset, arrangement, threads, model): the dataset stays a string, threads an int."""
    dataset, arr, threads, model = text.split("|")
    return dataset, arr, int(threads), model


def per_eval_rows(per_eval, stop=None):
    """NQC to NONREV per evaluation for each set, arrangement and thread count, raw and with the cost
    of a gdb stop (measured at -T 1) removed from both; verdicts at -T 1 only, on the worse of the two."""
    rows = []
    for (dataset, arr, threads) in sorted({k[:3] for k in per_eval}, key=lambda k: (k[2], str(k[0]), k[1])):
        nq, nr = per_eval.get((dataset, arr, threads, "nqc"), []), per_eval.get((dataset, arr, threads, "nonrev"), [])
        x = ratio(median(nq), median(nr))
        y = ratio(without_stop(median(nq), stop), without_stop(median(nr), stop)) if threads == 1 else None
        worst = max(v for v in (x, y) if v is not None) if x is not None else None
        flag = (None if threads != 1 else "missing" if worst is None else
                "investigate (above 1.25)" if worst > PER_EVAL_LINE else "within")
        rows.append(check(f"{dataset} Level 1 {arr} -T {threads}: NQC to NONREV per evaluation (report; "
                          f"{'line 1.25' if threads == 1 else 'no verdict above -T 1'})", None,
                          ratio=x, ratio_without_stop_cost=y, nqc_s=nq, nonrev_s=nr, verdict=flag))
    return rows


def split(nqc, nonrev, per_eval_nqc, per_eval_nonrev, stop=None):
    """A per-fit ratio as the per-evaluation ratio times the evaluation-count ratio; the per-fit ratio
    with each run's gdb stops removed; each model's time outside its counted evaluations (nuisance
    updates, setup, report). Nothing is computed from an incomplete run."""
    if not (nqc.get("complete") and nonrev.get("complete")):
        return {"per_fit": None, "incomplete": {"nqc": nqc.get("exit_code"), "nonrev": nonrev.get("exit_code")}}

    def rest(state, e):
        return state["wall_s"] - state["hits"] * e if state["wall_s"] and state["hits"] and e else None

    def clean(state):
        return state["wall_s"] - state["hits"] * stop if state["wall_s"] and state["hits"] and stop else None
    return {"per_fit": ratio(nqc["wall_s"], nonrev["wall_s"]),
            "per_fit_without_stop_cost": ratio(clean(nqc), clean(nonrev)),
            "per_evaluation": ratio(per_eval_nqc, per_eval_nonrev),
            "evaluation_count": ratio(nqc["hits"], nonrev["hits"]),
            "outside_evaluations_s": {"nqc": rest(nqc, per_eval_nqc), "nonrev": rest(nonrev, per_eval_nonrev)}}


def per_fit_row(key, parts):
    """Asserted: a per-fit ratio of at most 3, raw and with the gdb stops removed (above, the gate is
    blocked until fixed or accepted); above 2 it is flagged to investigate. A missing ratio fails."""
    vals = [v for v in (parts["per_fit"], parts.get("per_fit_without_stop_cost")) if v is not None]
    x = max(vals) if vals else None
    flag = ("missing" if x is None else "blocks the gate (above 3)" if x > BLOCK_LINE else
            "investigate (above 2)" if x > PER_FIT_LINE else "within")
    return check(f"{key}: NQC to NONREV per fit at most 3 (line 2 to investigate)", bool(x is not None and x <= BLOCK_LINE),
                 verdict=flag, **parts)


def linear_fit(ns, ys):
    """Least squares y = a + b n, with the largest relative residual."""
    b, a = np.polyfit(np.asarray(ns, float), np.asarray(ys, float), 1)
    resid = max(abs(a + b * n - y) / y for n, y in zip(ns, ys))
    return float(a), float(b), float(resid)


def at_cap(per_fit, cap):
    """Evaluations of a fit run to its round cap: its first linked fit, then cap - 1 rounds at the mean
    of its later ones; None with no later round to take the mean from."""
    return per_fit[0] + (cap - 1) * float(np.mean(per_fit[1:])) if len(per_fit) > 1 else None


def extrapolate(per_eval, fits, ram):
    """Labelled extrapolation, evaluation time only (nuisance updates excluded): 2000 alignments' worth
    of 500 columns, about a million columns, at 20 taxa (scale by the real taxon count); evaluations
    per fit at each round cap from the 10-set's fits; 32 cores at the 100-set's -T 4 efficiency."""
    out = {"taxa": TAXA, "columns": 2000 * SITES, "note": "evaluation time only; nuisance updates excluded"}
    # window arrangement -> (round cap, the 10-set's Level 1 and Level 2 fit keys)
    caps = {"-S": (LOOP_DEFAULT, "10-L1-S-G4", "10-L2-S-G4"), "-q": (NPARAM_DEFAULT - 1, "10-L1-q-G4", "10-L2-p-G4")}
    for model in ("nonrev", "nqc"):
        for arr, (cap, *keys) in caps.items():
            ys = [median(per_eval.get((str(n), arr, 1, model), [])) for n in SETS]
            if None in ys:
                continue
            a, b, resid = linear_fit(SETS, ys)
            e = a + 2000 * b
            t1, t4 = ys[-1], median(per_eval.get(("100", arr, 4, model), []))
            eff = t1 / t4 / 4 if t4 else None
            row = {"per_evaluation_fit": {"a_s": a, "b_s_per_alignment": b, "largest_relative_residual": resid},
                   "per_evaluation_s_million_columns": e, "efficiency_T4": eff}
            for level, key in zip(("L1", "L2"), keys):
                f = fits.get(key, {}).get(model)
                n = at_cap(f["per_fit_evaluations"], cap) if f and f.get("complete") else None
                if n is not None:
                    row[f"{level}_at_cap"] = {"evaluations": n, "hours_T1": n * e / 3600,
                                              "hours_32_cores": n * e / 3600 / (32 * eff) if eff else None}
            out[f"{model} {arr}"] = row
    rams = [ram.get(n) for n in SETS]
    if None not in rams:
        a, b, _ = linear_fit(SETS, rams)
        out["memory_mb_million_columns"] = a + 2000 * b
    return out


# ---------------------------------------------------------------- windows

def cmd_windows(a):
    out, binary, prov = start(a)
    data = Data(a.data)
    r = GdbRunner(binary, out, a.timeout, a.gdb, a.window)
    results = {"data": [data.verify()]}
    only = set(a.only.split(","))
    if only & {"checks", "turtle"}:
        src = source_fit(r, data)
        state = fit_state(src, "-p", True)
        results["source"] = [key_row(src.name, state), cpu_row(src.name, state["cpu_s"], state["wall_s"], 1),
                             check("turtle -p Level 2 +G4 source fit (report)", None,
                                   **{k: v for k, v in state.items() if k != "gains"})]
    checks = []
    if "checks" in only:
        checks += counter_checks(r, data, a.pairs)
    results["counter"] = checks
    per_eval, ram, rows = {}, {}, []
    plan = []
    for dataset in ("turtle", "10", "30", "100"):
        if dataset in only:
            plan += [(dataset, arr, 1, rep) for rep in range(REPEATS[dataset]) for arr in ("-S", "-q")]
    if "threads" in only:
        plan += [(str(n), arr, t, 0) for n, t in THREADED for arr in ("-S", "-q")]
    for i, (dataset, arr, threads, rep) in enumerate(plan):
        order = models(data) if i % 2 == 0 else tuple(reversed(models(data)))
        for label, model in order:
            args, files = data.level1(dataset, arr, model)
            run = r.run(f"w-{dataset}-{arr[1:]}-T{threads}-{label}-r{rep}", args, files, threads, mode="window")
            secs = run.gdb.get("window_s")
            per_eval.setdefault((dataset, arr, threads, label), []).append(secs / r.window if secs else None)
            text = run.stdout()
            m, k = RAM.search(text), KERNEL.search(text)
            if m and dataset != "turtle" and threads == 1:
                ram[int(dataset)] = max(ram.get(int(dataset), 0), int(m.group(1)))
            rows.append(check(f"{run.name}: window of {r.window} evaluations complete", run.code == 0,
                              per_evaluation_s=secs / r.window if secs else None, kernel_threads=k and int(k.group(1)),
                              ram_mb=m and int(m.group(1)), exit=run.code))
            rows.append(cpu_row(run.name, run.gdb.get("window_cpu_s"), secs, threads))
    results["windows"] = rows
    stop = next((row.get("mean_s") for row in checks if "mean_s" in row), None)
    results["per-evaluation ratios"] = per_eval_rows(per_eval, stop)
    spread = {f"{d} {arr}": ratio(median(per_eval.get((d, arr, 1, "nqc"), [])),
                                  median(per_eval.get((d, arr, 1, "nonrev"), [])))
              for d in ("turtle", *map(str, SETS)) for arr in ("-S", "-q")}
    results["expected"] = [check("per-evaluation ratio about equal at 10, 30 and 100 alignments, largest on turtle "
                                 "(report; stated before the runs)", None, ratios=spread)]
    finish(out, prov, r, results, {"per_eval": {encode_key(k): v for k, v in per_eval.items()}, "ram_mb": ram,
                                   "stop_cost_s": stop},
           "S3c profiling, Stage A: counter checks and per-evaluation windows")


def counter_checks(r, data, pairs):
    """NONREV's per-partition calls against its linked ones (3 turtle partitions), one round each; the
    cost of a gdb stop from plain and counted pairs of NQC's turtle -S Level 1 round."""
    rows = []
    parts = len(data.names)
    for arr in ("-S", "-q"):
        args, files = data.level1("turtle", arr, "NONREV")
        run = r.run(f"check-nonrev-turtle-{arr[1:]}", args + one_round(arr), files, extra=(PER_PARTITION,))
        h = run.gdb["hits"]
        rows.append(check(f"NONREV turtle Level 1 {arr}, one round: {PER_PARTITION} hits equal {parts} x "
                          f"{COUNTED} hits", bool(run.code == 0 and bool(h.get(COUNTED))
                          and h.get(PER_PARTITION) == parts * h[COUNTED]), hits=h, exit=run.code))
    args, files = data.level1("turtle", "-S", data.model)
    args = args + one_round("-S")
    costs, base = [], None
    for i in range(pairs):
        pair = {}
        for kind in (("plain", "counted") if i % 2 == 0 else ("counted", "plain")):
            name = f"overhead-{i}-{kind}"
            pair[kind] = r.plain(name, args, files) if kind == "plain" else r.run(name, args, files)
        p, c = pair["plain"], pair["counted"]
        pk, ck = final_struct(p) or {}, final_struct(c) or {}
        (pw, pc), (cw, cc) = wall_cpu(p), wall_cpu(c)
        rows += [cpu_row(p.name, pc, pw, 1), cpu_row(c.name, cc, cw, 1)]
        hits = c.gdb["hits"].get(COUNTED)
        # 17 significant digits identify a double, so equal strings are equal bits
        same = (p.code == 0 and c.code == 0 and pk.get("evaluations") == ck.get("evaluations")
                and pk.get("accepted_score") is not None and pk["accepted_score"] == ck["accepted_score"])
        rows.append(check(f"overhead pair {i}: plain and counted runs make the same evaluations and the same "
                          "accepted score, bitwise", same, plain=pk.get("accepted_score"), counted=ck.get("accepted_score"),
                          evaluations=(pk.get("evaluations"), ck.get("evaluations")), hits=hits))
        rows.append(check(f"counted NQC run {i}: hits equal the evaluations key", bool(hits is not None
                          and str(hits) == ck.get("evaluations")), hits=hits, key=ck.get("evaluations")))
        if same and pw and cw and hits:
            costs.append((cw - pw) / hits)
            base = pw / hits
    rows.append(check("cost of a gdb stop, from IQ-TREE's wall-clock times (report; ratios are also shown with it "
                      "removed, and verdicts take the worse)",
                      None, per_stop_s=costs, mean_s=float(np.mean(costs)) if costs else None,
                      plain_per_evaluation_s=base,
                      fraction=float(np.mean(costs)) / base if costs and base else None,
                      note="at -T 2 and -T 4 gdb stops every thread; thread results include a cost measured at -T 1"))
    return rows


# ---------------------------------------------------------------- fits

def fit_plan(only):
    """(key, dataset, level, arrangement, rate) of Stage B, filtered by key prefixes."""
    plan = [("turtle", "L1", "-S", "+G4"), ("turtle", "L1", "-q", "+G4")]
    plan += [("turtle", "L2", arr, rate) for arr in ("-S", "-p") for rate in ("+G4", "+I+G4", "+R4")]
    plan += [("10", "L1", "-S", "+G4"), ("10", "L1", "-q", "+G4"), ("10", "L2", "-S", "+G4"), ("10", "L2", "-p", "+G4")]
    keyed = [(f"{d}-{lv}-{arr[1:]}-{rate.replace('+', '')}", d, lv, arr, rate) for d, lv, arr, rate in plan]
    return [k for k in keyed if any(k[0].startswith(p) for p in only)]


def windows_problems(windows, binary_sha, data_sha):
    """Why a windows run cannot supply this fits run's per-evaluation times: another binary, other
    data, or a failed asserted row."""
    found = []
    if windows["provenance"]["binary"]["sha256"] != binary_sha:
        found.append("another binary")
    if windows["results"]["data"][0].get("data_json_sha256") != data_sha:
        found.append("other data")
    failed = [row["check"] for rows in windows["results"].values() for row in rows if row["passed"] is False]
    if failed:
        found.append(f"failed rows: {failed}")
    return found


def key_row(name, state):
    """Asserted for every counted NQC run: gdb's hits equal NQC's evaluations key."""
    return check(f"{name}: hits equal the evaluations key", bool(state["hits"] is not None
                 and state["hits"] == state["evaluations_key"]), hits=state["hits"], key=state["evaluations_key"])


def cmd_fits(a):
    out, binary, prov = start(a)
    data = Data(a.data)
    wdir = Path(a.windows).expanduser().resolve()
    windows = json.loads(read(wdir / "results.json"))
    problems = windows_problems(windows, prov["binary"]["sha256"], sha256(data.dir / "data.json"))
    if problems:
        sys.exit(f"{wdir}: unusable windows run ({'; '.join(problems)})")
    per_eval = {decode_key(k): v for k, v in windows["extra"]["per_eval"].items()}
    ram = {int(k): v for k, v in windows["extra"]["ram_mb"].items()}
    stop = windows["extra"].get("stop_cost_s")
    prov["windows"] = {"path": str(wdir), "results_sha256": sha256(wdir / "results.json")}
    r = GdbRunner(binary, out, a.timeout, a.gdb, a.window)
    results = {"data": [data.verify()]}
    plan = fit_plan(a.only.split(","))
    # turtle's source fit runs first and is also the NQC row of turtle-L2-p-G4
    src = source_fit(r, data) if any(k[1] == "turtle" for k in plan) else None
    fits, rows = {}, []
    for i, (key, dataset, level, arr, rate) in enumerate(plan):
        order = models(data) if i % 2 == 0 else tuple(reversed(models(data)))
        for label, model in order:
            if key == "turtle-L2-p-G4" and label == "nqc":
                run = src
            else:
                args, files = (data.level1(dataset, arr, model) if level == "L1" else
                               data.level2(dataset, arr, rate, model))
                run = r.run(f"f-{key}-{label}", args, files)
            state = fit_state(run, arr, label == "nqc")
            fits.setdefault(key, {})[label] = state
            rows.append(check(f"{run.name} completes, with its evaluations counted", state["complete"], exit=run.code,
                              hits=state["hits"]))
            if label == "nqc":
                rows.append(key_row(run.name, state))
            rows.append(cpu_row(run.name, state["cpu_s"], state["wall_s"], 1))
            rows.append(check(f"{run.name} (report)", None, **{k: v for k, v in state.items() if k != "gains"}))
    results["fits"] = rows
    verdicts = []
    for key, dataset, level, arr, _ in plan:
        f = fits.get(key, {})
        if "nqc" not in f or "nonrev" not in f:
            continue
        w = "-S" if arr == "-S" else "-q"
        pe = {m: median(per_eval.get((dataset, w, 1, m), [])) for m in ("nqc", "nonrev")}
        parts = split(f["nqc"], f["nonrev"], pe["nqc"], pe["nonrev"], stop)
        diff = (f["nonrev"]["free_parameters"] - f["nqc"]["free_parameters"]
                if f["nqc"]["free_parameters"] and f["nonrev"]["free_parameters"] else None)
        parts.update({"rounds": (f["nqc"]["rounds"], f["nonrev"]["rounds"]),
                      "ended_by": (f["nqc"]["ended_by"], f["nonrev"]["ended_by"]),
                      "final_lnl": (f["nqc"]["final_lnl"], f["nonrev"]["final_lnl"]),
                      "per_fit_evaluations": (f["nqc"]["per_fit_evaluations"], f["nonrev"]["per_fit_evaluations"]),
                      "parameters_nonrev_minus_nqc": diff, "parameters_as_expected_19": diff == 19})
        if level == "L1":
            parts["window_against_fit_per_evaluation"] = {m: (pe[m], f[m]["per_hit_s"]) for m in ("nqc", "nonrev")}
        verdicts.append(per_fit_row(key, parts))
    results["per-fit ratios"] = verdicts
    results["extrapolation"] = [check("training cost extrapolated (report; labelled extrapolation: about a million "
                                      "columns at 20 taxa, evaluation time only, round caps, -T 4 efficiency)", None,
                                      **extrapolate(per_eval, fits, ram))]
    finish(out, prov, r, results, {"fits": fits}, "S3c profiling, Stage B: full fits")


# ---------------------------------------------------------------- common

def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")


def provenance(a):
    return {"script": SCRIPT, "command": sys.argv, "started_utc": now(),
            "git_commit": capture(a.git, "rev-parse", "HEAD", cwd=REPO),
            "git_clean": capture(a.git, "status", "--porcelain", cwd=REPO) == "",
            "python": sys.version, "numpy": np.__version__, "platform": platform.platform(),
            "label": "gate measurement; thresholds set 2026-10-22 before any measurement"}


def start(a):
    out = new_out(a.out)
    binary = Path(a.binary).expanduser().resolve()
    prov = provenance(a)
    if not prov["git_clean"] and not a.allow_dirty:
        sys.exit(f"the repository is not clean ({prov['git_commit']}); commit first or pass --allow-dirty")
    beside = binary.parent / (binary.name + ".sha256")
    digest = sha256(binary)
    if beside.exists() and read(beside).split()[0] != digest:
        sys.exit(f"SHA-256 mismatch for {binary}: {digest}")
    prov.update({"binary": {"path": str(binary), "sha256": digest}, "gdb": capture(a.gdb, "--version").splitlines()[0],
                 "data": str(Path(a.data).expanduser().resolve()), "settings": {
                     "window": a.window, "pairs": getattr(a, "pairs", None), "repeats": REPEATS,
                     "threaded": THREADED, "lines": {
                         "per_evaluation": PER_EVAL_LINE, "per_fit": PER_FIT_LINE, "block": BLOCK_LINE},
                     "iqtree": "-seed 1 -v, NQC at its default domain, one run at a time"}})
    out.mkdir(parents=True)
    return out, binary, prov


def markdown(title, prov, results):
    lines = [f"# {title}", "", f"Generated by `{SCRIPT}` on {prov['started_utc']} at commit `{prov['git_commit']}` "
             f"(clean: {prov['git_clean']}); binary SHA-256 `{prov['binary']['sha256']}`.", ""]
    for section, rows in results.items():
        lines += [f"## {section}", "", "| Check | Result | Values |", "|---|---|---|"]
        for row in rows:
            verdict = {True: "pass", False: "FAIL", None: "report"}[row["passed"]]
            values = {k: v for k, v in row.items() if k not in ("check", "passed")}
            lines.append(f"| {row['check']} | {verdict} | `{json.dumps(values, default=str)[:700]}` |")
        lines.append("")
    return "\n".join(lines) + "\n"


def finish(out, prov, r, results, extra, title):
    prov["finished_utc"] = now()
    write(out / "results.json", json.dumps({"provenance": prov, "runs": r.runs, "results": results, "extra": extra},
                                           indent=1, default=str))
    write(out / "results.md", markdown(title, prov, results))
    failed = [row["check"] for rows in results.values() for row in rows if row["passed"] is False]
    for f in failed:
        print(f"  FAIL {f}")
    print(f"{len(failed)} asserted checks failed; wrote {out / 'results.md'}")
    sys.exit(1 if failed else 0)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="command", required=True)
    s = sub.add_parser("simulate")
    s.add_argument("--alisim", required=True, help="the frozen binary, whose AliSim simulates")
    w = sub.add_parser("windows")
    w.add_argument("--only", default="checks,turtle,10,30,100,threads")
    w.add_argument("--timeout", type=float, default=3600, help="seconds per IQ-TREE run")
    f = sub.add_parser("fits")
    f.add_argument("--windows", required=True, help="the windows run's directory")
    f.add_argument("--only", default="turtle,10", help="configuration key prefixes, such as turtle,10-L1-S,10-L2-S")
    f.add_argument("--timeout", type=float, default=14400, help="seconds per IQ-TREE run")
    w.add_argument("--pairs", type=int, default=PAIRS, help="plain and counted overhead pairs (fewer for development)")
    for x in (w, f):
        x.add_argument("--binary", required=True)
        x.add_argument("--data", required=True, help="simulate's directory")
        x.add_argument("--gdb", default="gdb")
        x.add_argument("--window", type=int, default=WINDOW, help="evaluations per window (fewer for development)")
    for x in (s, w, f):
        x.add_argument("--out", required=True, help="new directory outside the repository")
        x.add_argument("--git", default="git", help="git executable (git.exe under WSL)")
        x.add_argument("--allow-dirty", action="store_true")
    a = ap.parse_args()
    {"simulate": cmd_simulate, "windows": cmd_windows, "fits": cmd_fits}[a.command](a)


if __name__ == "__main__":
    main()
