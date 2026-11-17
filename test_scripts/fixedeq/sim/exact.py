"""Confirmatory: SIM-1's exact recovery (CODE_PLAN.md 3.8, "Exact fits"; decision 021's SIM-1 recovery rows,
with the convergence, certificate and residual rows of decisions 039 and 042). The screen step comes with
the bird data (SIM-1 commit 5).

  python test_scripts/fixedeq/sim/exact.py recover --git git.exe [--problems R1,R2,R3,R4,R6] [--workers 10] \
      [--out DIR] [--allow-dirty]

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
from oracle import chart, flux, likelihood as lk, nqc, readers, simulate  # noqa: E402
from regress import REPO, sha256, write  # noqa: E402
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


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="step", required=True)
    r = sub.add_parser("recover")
    r.add_argument("--problems", default="R1,R2,R3,R4,R6")
    r.add_argument("--workers", type=int, default=10)
    r.add_argument("--out", help="new directory outside the repository (default under ~/iqtree3-runs/sim1/)")
    r.add_argument("--git", default="git", help="git executable (git.exe under WSL)")
    r.add_argument("--allow-dirty", action="store_true")
    a = ap.parse_args()
    cmd_recover(a)


if __name__ == "__main__":
    main()
