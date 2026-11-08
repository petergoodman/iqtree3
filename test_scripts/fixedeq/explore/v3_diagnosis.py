"""Exploratory: why S5's V3 ends at a projected z-gradient of 2.2e-5 against its declared 1e-6
(CHANGELOG.md, 2026-11-08). Not a test; nothing here is a pass line, and no remedy is chosen here.

  python test_scripts/fixedeq/explore/v3_diagnosis.py \
      --recorded ~/iqtree3-runs/s5-gate/20261108T062232Z-775dea60-validate/results.json \
      --out DIR --git git.exe [--workers 10] [--allow-dirty]

[0] V3 rebuilt as reference.py's validate builds it, its random draws replayed in order; the
    truth's lnL against the recorded run's.
[1] Floor: at the exact truth, whose exact gradient is zero (exact pattern frequencies), the
    analytic projected z-gradient by the flux route and by the z route; V2's truth likewise.
[2] V3's fit rerun with every step logged, its endpoint against the recorded one, and the 360
    coordinates there: truth on a wall or not, held or not, gradient, distance to the truth.
[3] Repairs from that endpoint (skipped when [1] finds a floor at or above the stop): active-set
    Newton steps in the reference's flux basis with the finish's step (flux._run.finish, copied so
    flux.py stays unchanged), faces held {the endpoint's, the truth's} by Hessian {forward, central}
    at step 1e-6; and L-BFGS-B's analytic pass in z at pgtol 1e-6.
"""
import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse  # noqa: E402
import copy  # noqa: E402
import json  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
from scipy.linalg import LinAlgError, cho_factor, cho_solve, lstsq, null_space  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "reference"))
import reference as R  # noqa: E402
from oracle import cases, chart, flux, nqc, optimize, target  # noqa: E402
from probes import turtle_target  # noqa: E402
from regress import REPO, sha256, write  # noqa: E402

SCRIPT = Path(__file__).resolve().relative_to(REPO).as_posix()
LO, HI = R.LO, R.HI
N_LO, N_HI = 160, 5          # reference.extreme_truth's defaults
PG_STOP = flux.PG_STOP
STEP = flux.HESSIAN_STEP
REPAIR_STEPS = 20
PATIENCE = 5                 # a repair ends after this many steps without a new smallest gradient
LBFGS_SECONDS = 600


def say(text):
    print(text, flush=True)


# ---------------------------------------------------------------- [0] V3 as validate builds it

def rebuild():
    """cmd_validate's draws in order: V1's lognormal, kkt_truth, V2's cycle flux, extreme_truth. The
    generator is copied before extreme_truth so its coordinates and picks can be redrawn."""
    pi = target.parse_target(turtle_target())[0]
    Qlg, ref = nqc.lg_seed(pi)
    rng = np.random.default_rng(R.SEED)
    n = len(pi)
    rng.lognormal(0, 1, n * (n - 1))
    R.kkt_truth(pi, ref, rng)
    Q2 = cases.generator_from_flux(cases.cycle_flux(n, rng), pi)
    probe = copy.deepcopy(rng)
    Q3 = R.extreme_truth(pi, Qlg, ref, rng)
    z3 = chart.inverse(Qlg, ref)
    pick = probe.choice(len(z3), N_LO + N_HI, replace=False)
    z3[pick[:N_LO]], z3[pick[N_LO:]] = LO, HI
    same = bool(np.array_equal(chart.build(z3, pi, ref)[0], Q3))
    return pi, Qlg, ref, Q2, Q3, z3, np.sort(pick[:N_LO]), np.sort(pick[N_LO:]), same


def recorded(path):
    d = json.loads(Path(path).expanduser().read_text())
    rows = {r["check"]: r for r in d["results"]}
    deficit = rows["V3 lnL deficit against the truth <= 1e-3"]
    rec = d["records"]["V3"]
    return {"commit": d["provenance"]["git_commit"], "truth_lnl": deficit["truth_lnl"], "lnl": deficit["lnl"],
            "flux_value": rec["lnl"], "projected_z_gradient": rec["projected_z_gradient"],
            "m_active": rec["m_active"], "iterations": rec["iterations"]}


# ---------------------------------------------------------------- [1] floor at the truth

def floor(obj, pi, ref, dom, Qt, zt, wall=None):
    lf = flux.LikelihoodFlux(obj, pi)
    ft = flux.flux_of(Qt, pi)
    lnl, gf = lf(ft)
    zg = flux.z_gradient(ft, gf, ref)
    fz, gz = obj.gradient(zt)
    out = {"lnl_flux_route": float(lnl), "lnl_z_route": float(-fz),
           "projected_flux_route": flux.projected_z_gradient(ft, gf, dom),
           "projected_z_route": float(np.max(np.abs(optimize.projected_gradient(zt, gz, LO, HI)))),
           "max_abs_gradient_flux_route": float(np.max(np.abs(zg))),
           "max_abs_gradient_z_route": float(np.max(np.abs(gz))),
           "smallest_flux": float(np.min(ft)), "smallest_slack": float(np.min(dom.slack(ft)))}
    if wall is not None:
        inside = np.setdiff1d(np.arange(len(zt)), wall)
        out["max_abs_gradient_z_route_wall"] = float(np.max(np.abs(gz[wall])))
        out["max_abs_gradient_z_route_interior"] = float(np.max(np.abs(gz[inside])))
    return out


# ---------------------------------------------------------------- [2] the endpoint, coordinate by coordinate

def coordinates(obj, pi, ref, dom, f, Qt, zt, lo_k, hi_k):
    lnl, gf = flux.LikelihoodFlux(obj, pi)(f)
    z = flux.z_of_flux(f, ref)
    g = -flux.z_gradient(f, gf, ref)               # of -lnL, as the declared stop's
    pg = optimize.projected_gradient(z, g, LO, HI)
    m = len(z)
    slack = dom.slack(f)
    held = np.where(slack[:m] <= flux.ACTIVE, "lower", np.where(slack[m:] <= flux.ACTIVE, "upper", ""))
    truth = np.full(m, "interior", dtype=object)
    truth[lo_k], truth[hi_k] = "lower", "upper"
    pairs = [(i, j) for i, free in enumerate(chart.free_destinations(ref)) for j in free]
    Qe = flux.generator(f, pi)
    rows = [{"k": k, "i": i, "j": j, "truth": truth[k], "held": str(held[k]), "z_truth": zt[k], "z_end": z[k],
             "dz": z[k] - zt[k], "to_lower": z[k] - LO, "to_upper": HI - z[k], "slack_lower": slack[k],
             "slack_upper": slack[m + k], "gradient": g[k], "projected": pg[k], "q_truth": Qt[i, j],
             "q_end": Qe[i, j]} for k, (i, j) in enumerate(pairs)]
    classes = {}
    for r in rows:
        key = f"truth {r['truth']}, " + (f"held at {r['held']}" if r["held"] else "not held")
        c = classes.setdefault(key, {"count": 0, "max_abs_projected": 0.0, "max_abs_gradient": 0.0, "max_abs_dz": 0.0,
                                     "largest_q_truth": 0.0})
        c["count"] += 1
        c["max_abs_projected"] = max(c["max_abs_projected"], abs(float(r["projected"])))
        c["max_abs_gradient"] = max(c["max_abs_gradient"], abs(float(r["gradient"])))
        c["max_abs_dz"] = max(c["max_abs_dz"], abs(float(r["dz"])))
        c["largest_q_truth"] = max(c["largest_q_truth"], float(r["q_truth"]))
    top = sorted(rows, key=lambda r: -abs(r["projected"]))[:10]
    return rows, classes, top, float(np.max(np.abs(pg)))


# ---------------------------------------------------------------- [3] repairs from the endpoint

def hessian(pool, basis, u, g0, central):
    """Of -lnL in u, by forward or central differences of the analytic gradient at STEP, each probe's
    step shrunk by 10 until its fluxes stay positive (as flux._run.hessian), symmetrized."""
    n = len(u)
    probes, steps = [], []
    for k in range(n):
        h = STEP
        while True:
            e = np.zeros(n)
            e[k] = h
            plus, minus = basis.flux(u + e), basis.flux(u - e)
            if np.all(plus > 0) and (not central or np.all(minus > 0)):
                break
            h *= 0.1
            if h < 1e-15:
                raise RuntimeError("no Hessian probe keeps the fluxes positive")
        probes += [plus, minus] if central else [plus]
        steps.append(h)
    grads = [-basis.chain(np.asarray(gf, dtype=float)) for _, gf in pool.fluxes(probes)]
    H = np.empty((n, n))
    for k in range(n):
        H[:, k] = (grads[2 * k] - grads[2 * k + 1]) / (2 * steps[k]) if central else (grads[k] - g0) / steps[k]
    return (H + H.T) / 2


def newton_step(H, g):
    tau, scale = 0.0, max(1.0, float(np.max(np.abs(np.diag(H))))) if len(H) else 1.0
    while True:
        try:
            c = cho_factor(H + tau * np.eye(len(H)))
            break
        except LinAlgError:
            tau = 1e-12 * scale if tau == 0.0 else 10 * tau
    return -cho_solve(c, g), tau


def repair(pool, lf, dom, f_end, faces, central, truth_value):
    """The finish's step (working rows as equalities, their slacks closed in the same step; a blocking
    row joins), taken whole up to the first blocking row, with no value test and no release, so that
    only the faces and the Hessian differ between variants."""
    basis = flux.Basis(f_end, dom)
    u = np.zeros(basis.N.shape[1])
    W = faces.copy()
    trace, best, since, ending = [], np.inf, 0, "step cap"
    t0 = time.monotonic()
    for step in range(REPAIR_STEPS + 1):
        f = basis.flux(u)
        value, gf = lf(f)
        g = -basis.chain(gf)
        pg = flux.projected_z_gradient(f, gf, dom)
        s = basis.slack(u)
        entry = {"step": step, "deficit": truth_value - float(value), "pg": pg, "faces": int(W.sum()),
                 "max_abs_working_slack": float(np.max(np.abs(s[W]))) if W.any() else 0.0,
                 "min_other_slack": float(np.min(s[~W]))}
        trace.append(entry)
        say(f"    step {step}: pg {pg:.3g}, deficit {entry['deficit']:.3g}, faces {entry['faces']}")
        if pg <= PG_STOP:
            ending = "reached the stop"
            break
        best, since = (pg, 0) if pg < best else (best, since + 1)
        if since >= PATIENCE:
            ending = f"no smaller gradient in {PATIENCE} steps"
            break
        if step == REPAIR_STEPS:
            break
        H = hessian(pool, basis, u, g, central)
        GW = basis.G[W]
        dp = lstsq(GW, -s[W])[0] if W.any() else np.zeros(len(u))
        P = null_space(GW) if W.any() else np.eye(len(u))
        dn, tau = newton_step(P.T @ H @ P, P.T @ (g + H @ dp))
        d = dp + P @ dn
        Gd = basis.G @ d
        block = ~W & (Gd < 0)
        a_block = float(np.min(-s[block] / Gd[block])) if block.any() else np.inf
        alpha = min(1.0, a_block)
        if block.any() and alpha >= a_block:
            W[np.flatnonzero(block)[int(np.argmin(-s[block] / Gd[block]))]] = True
        while not np.all(basis.flux(u + alpha * d) > 0):
            alpha *= 0.5
        entry.update({"alpha": alpha, "damping": tau, "predicted_gain": float(-(g @ d) / 2)})
        u = u + alpha * d
    return {"ending": ending, "steps": len(trace) - 1, "first_pg": trace[0]["pg"],
            "smallest_pg": min(e["pg"] for e in trace), "final_pg": trace[-1]["pg"],
            "final_deficit": trace[-1]["deficit"], "final_faces": trace[-1]["faces"],
            "seconds": round(time.monotonic() - t0, 1), "trace": trace}


def climb_z(obj, pool, z0, truth_value, wall):
    t0 = time.monotonic()
    start_pg = float(np.max(np.abs(optimize.projected_gradient(z0, obj.gradient(z0)[1], LO, HI))))
    c = optimize.lbfgsb_fit(obj, z0, LO, HI, pgtol=PG_STOP, max_iter=20000, central_step=R.CENTRAL_STEP,
                            evaluate=pool.values, gradient=obj.gradient, maxcor=20, max_seconds=LBFGS_SECONDS)
    at_bound = np.union1d(c.at_lower, c.at_upper)
    return {"start_pg": start_pg, "analytic_pg": c.analytic_pg, "central_pg": c.pg_max,
            "deficit": truth_value - (-c.f), "phases": c.phases, "gradients": c.gradients,
            "at_bound": int(len(at_bound)), "truth_wall_at_bound": int(len(np.intersect1d(at_bound, wall))),
            "seconds": round(time.monotonic() - t0, 1)}


# ---------------------------------------------------------------- report

def report(out, prov, res):
    s0, s1, s2, s3 = res["step0"], res["step1"], res["step2"], res.get("step3")
    g = lambda v: f"{v:.3g}" if isinstance(v, float) else str(v)  # noqa: E731
    L = ["# V3 diagnosis", "",
         f"Exploratory. Generated by `{SCRIPT}` (sha256 {prov['script_sha256'][:12]}) on {prov['started_utc']} at "
         f"commit `{prov['git_commit']}` (clean: {prov['git_clean']}). Not a test.", "",
         "## [0] V3 rebuilt", "",
         f"Truth lnL {s0['truth_lnl']!r} against the recorded {s0['recorded']['truth_lnl']!r}: "
         f"{'equal' if s0['truth_equal'] else 'DIFFERENT'}; redrawn coordinates rebuild the truth bitwise: "
         f"{s0['coordinates_rebuild_truth']}.", "",
         "## [1] Gradient at the exact truth (exactly zero in exact arithmetic)", "",
         "| Truth | projected, flux route | projected, z route | max \\|g\\|, z route | wall | interior | smallest flux |",
         "|---|---|---|---|---|---|---|"]
    for name in ("V2", "V3"):
        r = s1[name]
        L.append(f"| {name} | {g(r['projected_flux_route'])} | {g(r['projected_z_route'])} | "
                 f"{g(r['max_abs_gradient_z_route'])} | {g(r.get('max_abs_gradient_z_route_wall', '-'))} | "
                 f"{g(r.get('max_abs_gradient_z_route_interior', '-'))} | {g(r['smallest_flux'])} |")
    L += ["", f"Floor at or above the stop {PG_STOP:g}: {s1['floor_at_stop']}.", "",
          "## [2] V3's fit rerun", "",
          f"{s2['status']}: {s2['message']}. lnL {s2['flux_value']!r} (recorded {s2['recorded_flux_value']!r}, "
          f"equal: {s2['value_equal']}); projected gradient {s2['pg']!r} (recorded {s2['recorded_pg']!r}, "
          f"equal: {s2['pg_equal']}); {s2['seconds']} s.", "",
          "Last steps logged:", "", "| step | mu | decrement | alpha | pure | faces | pg |", "|---|---|---|---|---|---|---|"]
    for e in s2["log"][-12:]:
        L.append(f"| {e['step']} | {g(e['mu'])} | {g(e['decrement'])} | {g(e['alpha'])} | {e['pure']} | "
                 f"{e.get('faces', '')} | {g(e['pg'])} |")
    L += ["", "By class at the endpoint:", "",
          "| class | count | max \\|projected\\| | max \\|g\\| | max \\|dz\\| | largest true rate |", "|---|---|---|---|---|---|"]
    for k, c in sorted(s2["classes"].items()):
        L.append(f"| {k} | {c['count']} | {g(c['max_abs_projected'])} | {g(c['max_abs_gradient'])} | "
                 f"{g(c['max_abs_dz'])} | {g(c['largest_q_truth'])} |")
    L += ["", "Ten largest projected components:", "",
          "| k | (i, j) | truth | held | projected | gradient | dz | to lower | true rate |", "|---|---|---|---|---|---|---|---|---|"]
    for r in s2["top"]:
        L.append(f"| {r['k']} | ({r['i']}, {r['j']}) | {r['truth']} | {r['held'] or '-'} | {g(float(r['projected']))} | "
                 f"{g(float(r['gradient']))} | {g(float(r['dz']))} | {g(float(r['to_lower']))} | {g(float(r['q_truth']))} |")
    L += ["", "## [3] Repairs from the endpoint", ""]
    if s3 is None:
        L.append("Skipped: the floor in [1] is at or above the stop.")
    else:
        L += ["| faces held | Hessian | ending | steps | first pg | smallest pg | final pg | final deficit | final faces |",
              "|---|---|---|---|---|---|---|---|---|"]
        for key, r in s3["repairs"].items():
            faces, kind = key.split(" / ")
            L.append(f"| {faces} | {kind} | {r['ending']} | {r['steps']} | {g(r['first_pg'])} | {g(r['smallest_pg'])} | "
                     f"{g(r['final_pg'])} | {g(r['final_deficit'])} | {r['final_faces']} |")
        z = s3["lbfgsb_z"]
        L += ["", f"L-BFGS-B in z (analytic, pgtol {PG_STOP:g}, memory 20, {LBFGS_SECONDS} s cap): start pg "
              f"{g(z['start_pg'])}, end analytic pg {g(z['analytic_pg'])}, central pg {g(z['central_pg'])}, deficit "
              f"{g(z['deficit'])}, {z['at_bound']} coordinates at a bound ({z['truth_wall_at_bound']} of the truth's "
              f"{N_LO + N_HI}), {z['gradients']} gradients, {z['seconds']} s; {z['phases']}."]
    write(out / "report.md", "\n".join(L) + "\n")


# ---------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--recorded", required=True, help="the failed validate run's results.json")
    ap.add_argument("--out", required=True)
    ap.add_argument("--git", default="git")
    ap.add_argument("--workers", type=int, default=10)
    ap.add_argument("--allow-dirty", action="store_true")
    a = ap.parse_args()
    prov = R.provenance(a, "v3-diagnosis")
    prov.update({"script": SCRIPT, "label": "exploratory", "script_sha256": sha256(Path(__file__).resolve())})
    out = R.new_out(a.out)
    res = {}

    pi, Qlg, ref, Q2, Q3, z3, lo_k, hi_k, same = rebuild()
    rec = recorded(a.recorded)
    obj3 = R.exact_objective(pi, ref, Q3)
    truth_lnl = obj3.log_likelihood(Q3)
    res["step0"] = {"recorded": rec, "truth_lnl": truth_lnl, "truth_equal": truth_lnl == rec["truth_lnl"],
                    "coordinates_rebuild_truth": same, "lower": lo_k, "upper": hi_k}
    say(f"[0] truth lnL {truth_lnl!r}, recorded {rec['truth_lnl']!r}; coordinates rebuild the truth: {same}")
    if not (res["step0"]["truth_equal"] and same):
        write(out / "results.json", json.dumps({"provenance": prov, **res}, indent=1, default=R._json))
        sys.exit("[0] V3 is not rebuilt exactly; stopping")

    dom = flux.matched(ref)
    wall = np.union1d(lo_k, hi_k)
    obj2 = R.exact_objective(pi, ref, Q2)
    s1 = {"V2": floor(obj2, pi, ref, dom, Q2, chart.inverse(Q2, ref)),
          "V3": floor(obj3, pi, ref, dom, Q3, z3, wall)}
    s1["floor_at_stop"] = bool(max(s1["V3"]["projected_flux_route"], s1["V3"]["projected_z_route"]) >= PG_STOP)
    res["step1"] = s1
    say(f"[1] V3 truth: projected {s1['V3']['projected_flux_route']:.3g} (flux), {s1['V3']['projected_z_route']:.3g} (z); "
        f"V2 truth: {s1['V2']['projected_flux_route']:.3g}, {s1['V2']['projected_z_route']:.3g}")

    pool = R.Pool(obj3, a.workers)
    lf = flux.LikelihoodFlux(obj3, pi)
    truth_value = float(lf(flux.flux_of(Q3, pi))[0])
    f0, moved = flux.inward(Qlg, pi, dom)
    log = []
    r = flux.fit(lf, f0, dom, objective_map=pool.fluxes, progress=log.append)
    rows, classes, top, pg_end = coordinates(obj3, pi, ref, dom, r.f, Q3, z3, lo_k, hi_k)
    res["step2"] = {"status": r.status, "message": r.message, "flux_value": r.value, "pg": r.pg,
                    "recorded_flux_value": rec["flux_value"], "recorded_pg": rec["projected_z_gradient"],
                    "value_equal": r.value == rec["flux_value"], "pg_equal": r.pg == rec["projected_z_gradient"],
                    "deficit": truth_lnl - obj3.log_likelihood(flux.generator(r.f, pi)), "m_active": r.m_active,
                    "start_moved": moved, "seconds": round(r.seconds, 1), "totals": r.totals(), "log": log,
                    "classes": classes, "top": top, "pg_from_table": pg_end, "coordinates": rows,
                    "endpoint_flux": r.f}
    say(f"[2] {r.status}: {r.message}; pg {r.pg!r} (recorded {rec['projected_z_gradient']!r}); {r.seconds:.0f} s")

    if not s1["floor_at_stop"]:
        m = len(z3)
        endpoint_faces = dom.slack(r.f) <= flux.ACTIVE
        truth_faces = np.zeros(2 * m, dtype=bool)
        truth_faces[lo_k], truth_faces[m + hi_k] = True, True
        repairs = {}
        for fname, faces in (("endpoint's", endpoint_faces), ("truth's", truth_faces)):
            for kind in ("forward", "central"):
                key = f"{fname} {int(faces.sum())} / {kind}"
                repairs[key] = repair(pool, lf, dom, r.f, faces, kind == "central", truth_value)
                say(f"[3] {key}: {repairs[key]['ending']}, pg {repairs[key]['first_pg']:.3g} -> "
                    f"{repairs[key]['final_pg']:.3g} in {repairs[key]['steps']} steps")
        z0 = np.clip(flux.z_of_flux(r.f, ref), LO, HI)
        res["step3"] = {"repairs": repairs, "lbfgsb_z": climb_z(obj3, pool, z0, truth_value, wall)}
        say(f"[3] L-BFGS-B in z: analytic pg {res['step3']['lbfgsb_z']['analytic_pg']:.3g}, "
            f"deficit {res['step3']['lbfgsb_z']['deficit']:.3g}")
    pool.close()

    prov["finished_utc"] = time.strftime("%Y-%m-%dT%H:%M:%S+00:00", time.gmtime())
    write(out / "results.json", json.dumps({"provenance": prov, **res}, indent=1, default=R._json))
    report(out, prov, res)
    say(f"wrote {out / 'report.md'}")


if __name__ == "__main__":
    main()
