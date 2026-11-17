"""SIM-1's exact fits (CODE_PLAN.md 3.8, "Exact fits"): recover, confirmatory (decision 021's SIM-1 recovery
rows, with the convergence, certificate and residual rows of decisions 039 and 042); screen, exploratory.

  python test_scripts/fixedeq/sim/exact.py recover --git git.exe [--problems R1,R2,R3,R4,R6] [--workers 10] \
      [--out DIR] [--allow-dirty]
  python test_scripts/fixedeq/sim/exact.py screen --git git.exe [--scenarios S1a,S1b,S2,S3,S4a,S4b,S5] \
      [--taxa 3|4] [--loci sim/scaffold/compositions.tsv] [--workers 10] [--out DIR] [--allow-dirty]

screen: NQC at pi* from LG at pi* by the flux reference, on exact weights of misspecified data, each scenario
on both base truths (Q.bird at pi*, C = 0; NQ.bird): S1a and S1b a global shift toward LG's composition (s =
0.25, 1), S2 per-locus compositions, S3 C20's site profiles, S4a residue-dependent and S4b column-context
masking, S5 gene-tree discordance; truths moved by T3 transport (decision 044). Reported per fit: status
("not evaluable" unless converged with a passing certificate), C, dC against the base truth and against the
generating reference, the per-site log ratio from data to fit, the largest net-flux errors and bound activity;
a fit with |dC against the generating reference| >= 0.1 x 0.0577 is flagged and refitted from its reference.
Interpretations (Peter to confirm): S5's alternative topologies keep each taxon's pendant length and the
internal branch; under masking the log ratio is not a KL (gaps are scored as unknown); the 4-taxon repeats are
the two flagged reversible-truth scenarios with the largest |dC| (repeat_selection), each on both truths, and
S5 has no 4-taxon form.

Each problem fits NQC at pi* (sim/target.txt) by the flux reference (decision 042) in the matched default
box with references from LG at pi* (decision 011), to exact pattern weights over every column (N = 1e6):
R1 NQ.bird, 3 taxa; R2 Q.bird at pi*; R3 LG at pi*, started at NQ.bird; R4 NQ.bird with +G4{0.5} in data
and fit; R5 NQ.bird, 4 taxa; R6 NQ.bird with every residue masked with probability 0.133, gaps scored as
unknown. Every other start is LG at pi*. Recorded before each fit: the truth inside the box (asserted);
the Hessian H of -lnL in z at the truth by central differences of the analytic gradient; each rate's
predicted error p_r = 1e-6 sum_k |(J H^-1)_rk|, J = d log q / dz. Asserted: converged; flux.certificate;
residuals; the lnL deficit against the truth; every rate within max(1e-3, 10 p_r) of the truth.
The default --problems leaves out R5 (4 taxa, hours), which runs in its own invocation. Each problem
writes NAME-before-fit.json before its fit, and results.json is rewritten after each problem.
"""
import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse  # noqa: E402
import json  # noqa: E402
import multiprocessing  # noqa: E402
import resource  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
import traceback  # noqa: E402
from dataclasses import dataclass  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
for _p in ("", "regression", "reference"):
    sys.path.insert(0, str(HERE.parent / _p))
sys.path.insert(0, str(HERE))
from oracle import builtin, chart, flux, likelihood as lk, nqc, readers, simulate, t3  # noqa: E402
from regress import REPO, read, sha256, write  # noqa: E402
from reference import Pool, fit_record  # noqa: E402
from nqcchecks import check  # noqa: E402
import common  # noqa: E402

SCRIPT = Path(__file__).resolve().relative_to(REPO).as_posix()
SITES = 1e6
THREE = "((a:0.2,b:0.3):0.1,c:0.4);"
FOUR = "((a:0.2,b:0.3):0.1,(c:0.25,d:0.35):0.15);"
ALPHA = 0.5
MASK = 0.133
CAP_3, CAP_4 = 4 * 3600.0, 8 * 3600.0
MAX_STEPS = 50_000
HESSIAN_STEP = 1e-5
STOP_GRADIENT = flux.PG_STOP      # 1e-6: the gradient size the declared stop leaves, behind p_r
RATE_FLOOR, RATE_FACTOR = 1e-3, 10.0
DEFICIT = 1e-3
RESIDUAL = 1e-10
CERTIFICATE = 1e-9                # times max(1, |lnL|)
MEMORY_SHARE = 0.8                # of MemAvailable given to the pool


@dataclass
class Problem:
    name: str
    truth: np.ndarray
    tree: str
    start: np.ndarray
    shape: float = None
    mask: float = None
    cap: float = CAP_3


def problems():
    """R1 to R6 at pi*, from the oracle's truths."""
    t = simulate.truths()[1]
    nq, lg = t["NQ.bird"], t["LG"]
    return {p.name: p for p in (Problem("R1", nq, THREE, lg), Problem("R2", t["Q.bird"], THREE, lg),
                                 Problem("R3", lg, THREE, nq), Problem("R4", nq, THREE, lg, shape=ALPHA),
                                 Problem("R5", nq, FOUR, lg, cap=CAP_4), Problem("R6", nq, THREE, lg, mask=MASK))}


def objective(p, pi, ref):
    """The NQC objective at pi on p's exact weights: every column of 20 residues (21 letters, the last
    scored as unknown, when masked), with p's rate categories in data and fit."""
    tree = readers.parse_newick(p.tree)
    names = tuple(sorted(leaf.name for leaf in tree.leaves()))
    rates, weights = lk.rate_categories(shape=p.shape)
    w = simulate.exact_weights([simulate.Component(1.0, tree, p.truth, pi, rates, weights)], names, SITES)
    if p.mask:
        w = simulate.mask(w, np.full((len(w), len(names)), p.mask))
    return nqc.Objective.from_patterns(tree, simulate.columns(names, gaps=bool(p.mask)), w, pi, ref, shape=p.shape)


# ---------------------------------------------------------------- identifiability

def central_hessian(gradients, z, step=HESSIAN_STEP):
    """The Hessian of -lnL in z, symmetrized, from central differences of the analytic gradient;
    gradients(list of z) returns [(value, gradient)] in order. Returns H and the gradients taken."""
    z = np.asarray(z, dtype=float)
    probes = []
    for k in range(len(z)):
        for s in (step, -step):
            x = z.copy()
            x[k] += s
            probes.append(x)
    g = [np.asarray(gk, dtype=float) for _, gk in gradients(probes)]
    H = np.empty((len(z), len(z)))
    for k in range(len(z)):
        H[:, k] = (g[2 * k] - g[2 * k + 1]) / (probes[2 * k][k] - probes[2 * k + 1][k])
    return (H + H.T) / 2, len(probes)


def log_rate_jacobian(z, pi, ref):
    """J[r, k] = d log q_r / d z_k, r over the off-diagonal rates in row-major order (chart.derivative)."""
    Q = chart.build(z, pi, ref)[0]
    off = ~np.eye(len(pi), dtype=bool)
    return np.array([chart.derivative(z, pi, ref, e)[off] / Q[off] for e in np.eye(len(z))]).T


def predicted_errors(H, J, gradient=STOP_GRADIENT):
    """p_r = gradient sum_k |(J H^-1)_rk|: each log rate's error when every z-gradient component is at
    most gradient (021's recovery row)."""
    return gradient * np.abs(np.linalg.solve(H, J.T).T).sum(axis=1)


# ---------------------------------------------------------------- workers

def gradient_peak(obj, z):
    """Resident memory a forked process adds to take one gradient (bytes): its peak less its own
    resident set at the start, shared with the parent; and the gradient's seconds. A child that dies
    without reporting raises RuntimeError."""
    ctx = multiprocessing.get_context("fork")
    parent, child = ctx.Pipe()

    def work():
        base = _rss()
        t0 = time.perf_counter()
        obj.gradient(z)
        child.send((base, resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024, time.perf_counter() - t0))

    proc = ctx.Process(target=work)
    proc.start()
    child.close()
    try:
        base, peak, seconds = parent.recv()
    except EOFError:
        proc.join()
        raise RuntimeError(f"the memory probe exited with {proc.exitcode} before reporting") from None
    proc.join()
    return {"child_start_rss": base, "child_peak": peak, "added": max(peak - base, 1), "gradient_seconds": seconds}


def _rss():
    with open("/proc/self/statm") as fh:
        return int(fh.read().split()[1]) * os.sysconf("SC_PAGE_SIZE")


def mem_available():
    with open("/proc/meminfo") as fh:
        return next(int(line.split()[1]) * 1024 for line in fh if line.startswith("MemAvailable:"))


def workers_for(requested, added, available):
    return max(1, min(requested, int(MEMORY_SHARE * available // added)))


# ---------------------------------------------------------------- one recovery

def _log(text):
    print(text, flush=True)


def fit_and_check(name, obj, truth, start, pi, ref, cap, workers, log=_log, before=None):
    """Rows and a record for one recovery: the identifiability measures at the truth, the flux
    reference from start, and 021's rows against the truth. before, if given, is called with the
    record made before the fit (the truth's z, H's eigenvalues, every rate's p_r and line)."""
    pi = np.asarray(pi, dtype=float)
    dom = flux.matched(ref)
    zt, ft = chart.inverse(truth, ref), flux.flux_of(truth, pi)
    rows = [check(f"{name}: the truth inside the box", dom.contains(ft), smallest_slack=float(np.min(dom.slack(ft))),
                  z_min=float(zt.min()), z_max=float(zt.max()))]
    pool = Pool(obj, workers)
    try:
        t0 = time.monotonic()
        H, taken = central_hessian(pool.gradients, zt)
        hessian_seconds = time.monotonic() - t0
        log(f"  {name}: Hessian at the truth, {taken} gradients, {hessian_seconds:.0f} s")
        J = log_rate_jacobian(zt, pi, ref)
        p_r = predicted_errors(H, J)
        rate_line = np.maximum(RATE_FLOOR, RATE_FACTOR * p_r)
        if before is not None:
            before({"z_truth": zt, "hessian_eigenvalues": np.linalg.eigvalsh(H), "predicted_errors": p_r,
                    "rate_lines": rate_line, "hessian_seconds": hessian_seconds})
        f0, moved = flux.inward(start, pi, dom)
        last = [time.monotonic()]

        def progress(state):
            if time.monotonic() - last[0] >= 300:
                last[0] = time.monotonic()
                log(f"  {name}: step {state['step']}, mu {state['mu']:.1e}, lnL {state['value']:.6f}")

        r = flux.fit(flux.LikelihoodFlux(obj, pi), f0, dom, objective_map=pool.fluxes, max_iter=MAX_STEPS,
                     max_seconds=cap, progress=progress)
    finally:
        pool.close()
        if pool.pool is not None:
            pool.pool.join()
    record ={"fit": fit_record(r), "start_moved": moved, "hessian_seconds": round(hessian_seconds, 1),
              "hessian_gradients": taken, "hessian_smallest_eigenvalue": float(np.linalg.eigvalsh(H).min()),
              "predicted_error_max": float(p_r.max()), "predicted_error_median": float(np.median(p_r)),
              "workers": workers}
    rows.append(check(f"{name}: converged under the declared stop", r.status == "converged", status=r.status,
                      message=r.message, projected_z_gradient=r.pg, mu=r.mu, seconds=round(r.seconds, 1)))
    if r.status == "failed":
        return rows, record
    value, gf = flux.LikelihoodFlux(obj, pi)(r.f)
    cert = flux.certificate(r.f, gf, dom)
    line = CERTIFICATE * max(1.0, abs(value))
    res = flux.residuals(r.f, pi)
    Q = flux.generator(r.f, pi)
    lnl_t, lnl = obj.log_likelihood(truth), obj.log_likelihood(Q)
    off = ~np.eye(len(pi), dtype=bool)
    err = np.abs(Q[off] - truth[off]) / truth[off]
    net = 2 * (simulate.flux_parts(Q, pi)[1] - simulate.flux_parts(truth, pi)[1])
    rows += [check(f"{name}: certificate <= 1e-9 max(1, |lnL|)", cert["residual"] <= line, line=line, **cert),
             check(f"{name}: residuals <= 1e-10", max(res.values()) <= RESIDUAL, **res),
             check(f"{name}: lnL deficit against the truth <= 1e-3", lnl_t - lnl <= DEFICIT,
                   deficit=lnl_t - lnl, truth_lnl=lnl_t, lnl=lnl),
             check(f"{name}: every rate within max(1e-3, 10 p_r) of the truth", bool(np.all(err <= rate_line)),
                   outside=int(np.sum(err > rate_line)), worst_ratio=float(np.max(err / rate_line)),
                   largest_relative_error=float(err.max()), relative_errors=err, lines=rate_line),
             check(f"{name}: circulation and net flux against the truth", None,
                   C_truth=simulate.circulation(truth, pi), C_fit=simulate.circulation(Q, pi),
                   C_error=abs(simulate.circulation(Q, pi) - simulate.circulation(truth, pi)),
                   net_flux_error_max=float(np.max(np.abs(net))), net_flux_error_rms=float(np.sqrt(np.mean(net ** 2))))]
    return rows, record


def cmd_recover(a):
    names = a.problems.split(",")
    every = problems()
    if not names or any(n not in every for n in names):
        sys.exit(f"--problems takes a comma-separated subset of {', '.join(every)}")
    prov = common.provenance(a, SCRIPT, "recover-" + "-".join(names), "confirmatory", problems=names,
                             sites=SITES, max_steps=MAX_STEPS, caps_seconds={n: every[n].cap for n in names},
                             hessian_step=HESSIAN_STEP)
    out = common.new_out(a.out, prov)
    pi = common.target_values()
    prov["target_sha256"] = sha256(common.TARGET)
    ref = nqc.lg_seed(pi)[1]
    rows = [check("pi* file equals the oracle's pi*", np.array_equal(pi, simulate.truths()[0]))]
    records, title = {}, f"SIM-1 exact recovery ({', '.join(names)})"
    for name in names:
        p = every[name]
        t0 = time.monotonic()
        try:
            obj = objective(p, pi, ref)
            peak = gradient_peak(obj, chart.inverse(p.truth, ref))
            available = mem_available()
            workers = workers_for(a.workers, peak["added"], available)
            print(f"{name}: {len(obj.counts)} columns, workers {workers} (gradient adds "
                  f"{peak['added'] / 2**20:.0f} MiB, {peak['gradient_seconds']:.2f} s; available "
                  f"{available / 2**30:.1f} GiB)", flush=True)
            records[name] = {"columns": len(obj.counts), "memory": peak, "mem_available": available}
            write(out / f"{name}-before-fit.json", json.dumps({"provenance": prov, name: records[name]}, indent=1,
                                                               default=common._json))

            def before(pre, name=name):
                records[name]["before_fit"] = pre
                write(out / f"{name}-before-fit.json", json.dumps({"provenance": prov, name: records[name]},
                                                                   indent=1, default=common._json))

            found, record = fit_and_check(name, obj, p.truth, p.start, pi, ref, p.cap, workers, before=before)
            records[name].update(record)
        except Exception:  # noqa: BLE001  a problem that raises is a failed row; the others still run
            found = [check(f"{name}: completes", False, error=traceback.format_exc()[-2000:])]
        records.setdefault(name, {})["seconds"] = round(time.monotonic() - t0, 1)
        rows += found
        failed = [row["check"] for row in found if row["passed"] is False]
        print(f"{name}: {records[name].get('fit', {}).get('status')}, {records[name]['seconds']:.0f} s, "
              f"failed: {failed or 'none'}", flush=True)
        common.save(out, prov, title + ", in progress", rows, {"records": records})
    common.finish(out, prov, title, rows, {"records": records})


# ---------------------------------------------------------------- screen (exploratory)

SCENARIOS = ("S1a", "S1b", "S2", "S3", "S4a", "S4b", "S5")
SHIFTS = {"S1a": 0.25, "S1b": 1.0}          # global shift toward LG's composition
MASK_RATES = {**{c: 0.2 for c in "ARGQSPEK"}, **{c: 0.05 for c in "WCFIYVLN"}}   # S4a; the rest 0.1
MASK_REST = 0.1
MASK_TARGET = 0.133                          # S4b's overall masked fraction, S4a's measured one
LOCUS_FLOOR = 1e-4                           # S2: every recentred locus entry at least this
DISCORDANCE = (0.69, 0.155, 0.155)           # S5: decision 043's 31% discordant, split evenly
TOPOLOGIES = (THREE, "((a:0.2,c:0.4):0.1,b:0.3);", "((b:0.3,c:0.4):0.1,a:0.2);")   # each taxon keeps its length
FLAG = 0.1 * 0.0577                          # |dC against the generating reference| (Peter)
BASES = ("Q.bird", "NQ.bird")
NEAR_BOUND = 1e-6


def residue_rates():
    return np.array([MASK_RATES.get(c, MASK_REST) for c in lk.AMINO])


def read_loci(path):
    """S2's per-locus table (sim/birddata.py inventory, sim/scaffold/compositions.tsv): '#' provenance lines,
    a header, then one line per locus: name, length and the 20 residue counts in IQ-TREE's order, tab-separated."""
    lines = [line for line in read(path).splitlines() if line.strip() and not line.startswith("#")]
    if not lines or lines[0].split("\t")[2:] != list(lk.AMINO):
        raise ValueError(f"{path}: the header's residues are not in IQ-TREE's order {lk.AMINO}")
    rows = [line.split("\t") for line in lines[1:]]
    if not rows or any(len(r) != 22 for r in rows):
        raise ValueError(f"{path}: expected lines of a name, a length and 20 counts")
    return np.array([float(r[1]) for r in rows]), np.array([[float(v) for v in r[2:]] for r in rows])


def composition(weights, k):
    """The residue composition of masked weights over 21^k columns, gaps left out."""
    digits = simulate.column_digits(k, simulate.STATES + 1)
    counts = np.zeros(simulate.STATES + 1)
    for leaf in range(k):
        counts += np.bincount(digits[:, leaf], weights=weights, minlength=simulate.STATES + 1)
    return counts[:simulate.STATES] / counts[:simulate.STATES].sum()


def screen_data(scenario, base, pi, tree_text, loci=None):
    """Exact weights for one scenario on one base truth at pi*, over the columns of the fitted tree's leaves;
    returns the weights, whether they hold gaps, the generating reference (Q, its pi) and a record."""
    tree = readers.parse_newick(tree_text)
    names = tuple(sorted(leaf.name for leaf in tree.leaves()))
    k = len(names)
    if scenario in SHIFTS:
        pd = simulate.shift(pi, builtin.exchangeabilities("LG")[1], SHIFTS[scenario])
        Qd = simulate.t3_transport(base, pi, pd)
        w = simulate.exact_weights([simulate.Component(1.0, tree, Qd, pd)], names, SITES)
        return w, False, (Qd, pd), {"s": SHIFTS[scenario], "total_variation": simulate.total_variation(pi, pd)}
    if scenario in ("S2", "S3"):
        if scenario == "S3":
            profiles, weights = builtin.frequency_profiles("C20")
            lam = simulate.recentre_limit(pi, profiles, weights) / 2
        else:
            if loci is None:
                raise ValueError("S2 needs the per-locus table (--loci)")
            weights, counts = loci
            profiles = (counts + 1) / (counts + 1).sum(axis=1, keepdims=True)
            lam = min(1.0, simulate.recentre_limit(pi, profiles, weights, LOCUS_FLOOR))
        moved = simulate.recentre(pi, profiles, weights, lam)
        comps = [simulate.Component(wk, tree, simulate.t3_transport(base, pi, pk), pk) for wk, pk in zip(weights, moved)]
        tv = np.array([simulate.total_variation(pk, pi) for pk in moved])
        w = simulate.exact_weights(comps, names, SITES)
        return w, False, (flux.generator(simulate.mean_flux(comps), pi), pi), {
            "lambda": lam, "components": len(comps), "smallest_entry": float(moved.min()),
            "weighted_mean_tv": float(np.asarray(weights) @ tv / np.sum(weights))}
    if scenario in ("S4a", "S4b"):
        w0 = simulate.exact_weights([simulate.Component(1.0, tree, base, pi)], names, SITES)
        if scenario == "S4a":
            prob = simulate.residue_mask(k, residue_rates())
            p = None
        else:
            p = MASK_TARGET / simulate.masked_fraction(w0, simulate.context_mask(k, 1.0))
            if p > 1:
                raise ValueError(f"S4b: masking {MASK_TARGET} needs probability {p:.3g} > 1")
            prob = simulate.context_mask(k, p)
        w = simulate.mask(w0, prob)
        return w, True, (base, pi), {"masked_fraction": simulate.masked_fraction(w0, prob), "probability": p,
                                     "composition_shift_tv": simulate.total_variation(composition(w, k), pi)}
    if scenario == "S5":
        if k != 3:
            raise ValueError("S5's three rooted topologies are defined at three taxa only")
        comps = [simulate.Component(wt, readers.parse_newick(t), base, pi) for wt, t in zip(DISCORDANCE, TOPOLOGIES)]
        return simulate.exact_weights(comps, names, SITES), False, (base, pi), {"weights": DISCORDANCE}
    raise ValueError(f"unknown scenario {scenario}")


def screen_fit(label, obj, start_flux, pi, ref, base, reference, workers, cap, log=_log):
    """One screen fit: NQC at pi from start_flux by the flux reference, and the screen's report (CODE_PLAN 3.8)."""
    dom = flux.matched(ref)
    f0, moved = flux.inward_flux(start_flux, dom, pi)
    pool = Pool(obj, workers)
    try:
        r = flux.fit(flux.LikelihoodFlux(obj, pi), f0, dom, objective_map=pool.fluxes, max_iter=MAX_STEPS,
                     max_seconds=cap)
    finally:
        pool.close()
        if pool.pool is not None:
            pool.pool.join()
    out = {"fit": fit_record(r), "start_moved": moved}
    if r.status == "failed":
        out.update({"evaluable": False, "flag": False})
        return out
    value, gf = flux.LikelihoodFlux(obj, pi)(r.f)
    cert = flux.certificate(r.f, gf, dom)
    evaluable = r.status == "converged" and cert["residual"] <= CERTIFICATE * max(1.0, abs(value))
    Q = flux.generator(r.f, pi)
    Qr, pr = reference
    c_hat, c_ref, c_base = (simulate.circulation(Q, pi), simulate.circulation(Qr, pr), simulate.circulation(base, pi))
    lnl = obj.log_likelihood(Q)
    w = obj.counts[obj.counts > 0]
    net = 2 * (simulate.flux_parts(Q, pi)[1] - simulate.flux_parts(Qr, pr)[1])
    i, j, _ = t3._edges(len(pi))
    upper = np.flatnonzero(i < j)       # net flux is antisymmetric: each pair once
    top = upper[np.argsort(-np.abs(net[upper]))[:3]]
    z = chart.inverse(Q, ref)
    out.update({"evaluable": bool(evaluable), "certificate": cert["residual"], "lnl": lnl, "C_fit": c_hat,
                "C_reference": c_ref, "C_base": c_base, "dC_reference": c_hat - c_ref, "dC_base": c_hat - c_base,
                "log_ratio_per_site": float((w @ np.log(w / SITES) - lnl) / SITES),
                "largest_net_flux_errors": [[lk.AMINO[i[e]] + lk.AMINO[j[e]], float(net[e])] for e in top],
                "at_lower": int(np.sum(z <= dom.lower + NEAR_BOUND)), "at_upper": int(np.sum(z >= dom.upper - NEAR_BOUND)),
                "m_active": r.m_active, "flag": bool(evaluable and abs(c_hat - c_ref) >= FLAG)})
    log(f"  {label}: {r.status}, evaluable {out['evaluable']}, C {c_hat:.5f} (reference {c_ref:.5f}, "
        f"dC {c_hat - c_ref:+.5f}), flag {out['flag']}, {r.seconds:.0f} s")
    return out


def best_endpoint(fit):
    """A flagged fit's endpoint with the higher lnL, its own or its restart's from the generating reference
    (an evaluable restart only): the flag and dC the 4-taxon selection reads."""
    again = fit.get("from_reference") or {}
    if again.get("evaluable") and again.get("lnl") is not None and again["lnl"] > fit["lnl"]:
        return again
    return fit


def repeat_selection(fits):
    """CODE_PLAN 3.8's 4-taxon repeats: among the reversible truth's evaluable fits, the two flagged with the
    largest |dC against the generating reference|, or the single largest if none flags; each fit read at
    its higher-lnL endpoint (best_endpoint)."""
    rev = [best_endpoint(f) | {"scenario": f["scenario"]} for f in fits if f["base"] == "Q.bird" and f.get("evaluable")]
    flagged = sorted((f for f in rev if f["flag"]), key=lambda f: -abs(f["dC_reference"]))
    chosen = flagged[:2] if flagged else sorted(rev, key=lambda f: -abs(f["dC_reference"]))[:1]
    return [f["scenario"] for f in chosen]


def cmd_screen(a):
    scenarios = a.scenarios.split(",")
    if any(s not in SCENARIOS for s in scenarios):
        sys.exit(f"--scenarios takes a comma-separated subset of {', '.join(SCENARIOS)}")
    tree_text, cap = (THREE, CAP_3) if a.taxa == 3 else (FOUR, CAP_4)
    prov = common.provenance(a, SCRIPT, f"screen-{a.taxa}taxa-" + "-".join(scenarios), "exploratory",
                             scenarios=scenarios, taxa=a.taxa, tree=tree_text, sites=SITES, flag=FLAG,
                             cap_seconds=cap,
                             loci={"path": a.loci, "sha256": sha256(a.loci)} if "S2" in scenarios else None)
    out = common.new_out(a.out, prov)
    pi = common.target_values()
    prov["target_sha256"] = sha256(common.TARGET)
    _, truths = simulate.truths()
    ref = nqc.lg_seed(pi)[1]
    start = flux.flux_of(truths["LG"], pi)
    loci = read_loci(a.loci) if "S2" in scenarios else None
    tree = readers.parse_newick(tree_text)
    names = tuple(sorted(leaf.name for leaf in tree.leaves()))
    rows, fits, title = [], [], f"SIM-1 exact screen, {a.taxa} taxa ({', '.join(scenarios)})"
    for scenario in scenarios:
        if scenario == "S5" and a.taxa != 3:
            rows.append(check("S5: skipped (its three rooted topologies are defined at three taxa)", None))
            continue
        for base in BASES:
            label = f"{scenario} {base}"
            t0 = time.monotonic()
            try:
                w, gaps, reference, info = screen_data(scenario, truths[base], pi, tree_text, loci)
                obj = nqc.Objective.from_patterns(tree, simulate.columns(names, gaps=gaps), w, pi, ref)
                peak = gradient_peak(obj, chart.inverse(truths[base], ref))
                workers = workers_for(a.workers, peak["added"], mem_available())
                fit = screen_fit(label, obj, start, pi, ref, truths[base], reference, workers, cap)
            except Exception:  # noqa: BLE001  logged as not evaluable; the other fits still run
                fits.append({"scenario": scenario, "base": base, "evaluable": False, "flag": False,
                             "seconds": round(time.monotonic() - t0, 1)})
                rows.append(check(f"{label}: raised", None, error=traceback.format_exc()[-2000:]))
                common.save(out, prov, title + ", in progress", rows, {"fits": fits})
                continue
            fit.update({"scenario": scenario, "base": base, "data": info, "columns": len(w), "workers": workers,
                        "seconds": round(time.monotonic() - t0, 1)})
            fits.append(fit)
            common.save(out, prov, title + ", in progress", rows, {"fits": fits})
            if fit["flag"]:     # the restart from the generating reference, its failure kept apart from the fit
                try:
                    again = screen_fit(label + " from its reference", obj, flux.flux_of(*reference), pi, ref,
                                       truths[base], reference, workers, cap)
                    fit["from_reference"] = {k: again.get(k) for k in ("evaluable", "lnl", "C_fit", "dC_reference",
                                                                       "flag", "certificate")}
                    fit["from_reference"]["status"] = again["fit"]["status"]
                except Exception:  # noqa: BLE001
                    fit["from_reference"] = {"error": traceback.format_exc()[-2000:], "evaluable": False}
                fit["higher_endpoint"] = "restart" if best_endpoint(fit) is not fit else "first"
            rows.append(check(f"{label}", None, **{k: v for k, v in fit.items() if k != "fit"},
                              status=fit["fit"]["status"]))
            common.save(out, prov, title + ", in progress", rows, {"fits": fits})
    counts = {"fits": len(fits), "not_evaluable": sum(not f.get("evaluable") for f in fits),
              "flagged": [f"{f['scenario']} {f['base']}" for f in fits if f.get("flag")]}
    rows.append(check("screen summary", None, **counts, repeat_selection=repeat_selection(fits)))
    common.finish(out, prov, title, rows, {"fits": fits})


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="step", required=True)
    r = sub.add_parser("recover")
    r.add_argument("--problems", default="R1,R2,R3,R4,R6")
    r.add_argument("--workers", type=int, default=10)
    r.add_argument("--out", help="new directory outside the repository (default under ~/iqtree3-runs/sim1/)")
    r.add_argument("--git", default="git", help="git executable (git.exe under WSL)")
    r.add_argument("--allow-dirty", action="store_true")
    s = sub.add_parser("screen")
    s.add_argument("--scenarios", default=",".join(SCENARIOS))
    s.add_argument("--taxa", type=int, choices=(3, 4), default=3)
    s.add_argument("--loci", default=str(HERE / "scaffold" / "compositions.tsv"), help="S2's per-locus table")
    s.add_argument("--workers", type=int, default=10)
    s.add_argument("--out", help="new directory outside the repository (default under ~/iqtree3-runs/sim1/)")
    s.add_argument("--git", default="git", help="git executable (git.exe under WSL)")
    s.add_argument("--allow-dirty", action="store_true")
    a = ap.parse_args()
    {"recover": cmd_recover, "screen": cmd_screen}[a.step](a)


if __name__ == "__main__":
    main()
