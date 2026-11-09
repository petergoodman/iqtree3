"""Exploratory: why S5's z climbs Z_e and Z_e5 end above the declared analytic projected z-gradient
1e-6 (A5 i) on toy (CHANGELOG.md, 2026-11-09). Not a test; nothing here is a pass line, and no remedy
is chosen here.

  python test_scripts/fixedeq/explore/z_climb_diagnosis.py \
      --calibrated ~/iqtree3-runs/s5-gate/20261109T011054Z-feec3cd5-toy-calibrate \
      --out DIR --git git.exe [--workers 10] [--allow-dirty]

[0] The problem rebuilt from the calibration's data (SHA-256 against its record); f and the analytic
    projected gradient at the recorded Z_e and Z_e5 against their records.
[1] At each endpoint: the analytic projected gradient by coordinate class (at a bound, within 1e-3
    of one, interior) and its ten largest components; f's rounding, as the spread of f over 20
    points moved 1e-10 at random (the slope moves f by about 1e-14 there).
[2] Repairs from each endpoint: (a) L-BFGS-B's analytic pass restarted with its own memory, 120 s
    cap; (b) projected Newton in the box, the Hessian by central differences of the analytic
    gradient (step 1e-6 max(1, |z|)) on the free coordinates, those not at a bound the gradient
    pushes against, recomputed every step; a free coordinate at a bound that the step would push out
    is held for that step; the step is cut at the first bound reached; acceptance as flux.py's
    finish: Armijo above 1e-12 |f| (flux.ROUNDING), below it no rise above that and a decrement at
    least halving each step; 20 steps. Step 0 records the Newton decrement and the free Hessian's
    extreme eigenvalues.
"""
import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse  # noqa: E402
import json  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
from scipy.linalg import LinAlgError, cho_factor, cho_solve  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "reference"))
import reference as R  # noqa: E402
from oracle import chart, flux, optimize  # noqa: E402
from regress import REPO, read, sha256, write  # noqa: E402

AA = R.AA

SCRIPT = Path(__file__).resolve().relative_to(REPO).as_posix()
LO, HI = R.LO, R.HI
PG_STOP = R.PG_ANALYTIC
STEP = flux.HESSIAN_STEP
ROUNDING, ARMIJO, ALPHA_MIN = flux.ROUNDING, flux.ARMIJO, flux.ALPHA_MIN
NEAR = 1e-3
NEWTON_STEPS = 20
RESTART_SECONDS = 120
NOISE_POINTS, NOISE_MOVE = 20, 1e-10


def say(text):
    print(text, flush=True)


def _grad(z):
    return R._OBJ.gradient(np.asarray(z))


def gradients(pool, zs):
    if pool.pool is None:
        return [_grad(z) for z in zs]
    return pool.pool.map(_grad, zs, chunksize=max(1, len(zs) // (4 * pool.workers)))


def pgrad(z, g):
    return optimize.projected_gradient(z, g, LO, HI)


# ---------------------------------------------------------------- [1] the endpoint

def endpoint(obj, z, ref, rng):
    f, g = obj.gradient(z)
    pg = pgrad(z, g)
    dist = np.minimum(z - LO, HI - z)
    cls = np.where(dist == 0, "at a bound", np.where(dist <= NEAR, f"within {NEAR:g}", "interior"))
    classes = {c: {"count": int(np.sum(cls == c)), "max_abs_projected": float(np.max(np.abs(pg[cls == c]))),
                   "max_abs_gradient": float(np.max(np.abs(g[cls == c])))} for c in np.unique(cls)}
    pairs = [(i, int(j)) for i, free in enumerate(chart.free_destinations(ref)) for j in free]
    top = [{"k": int(k), "i": pairs[k][0], "j": pairs[k][1], "z": float(z[k]), "to_bound": float(dist[k]),
            "gradient": float(g[k]), "projected": float(pg[k]), "class": str(cls[k])}
           for k in np.argsort(-np.abs(pg))[:10]]
    values = [obj.gradient(z + NOISE_MOVE * rng.choice([-1.0, 1.0], len(z)))[0] for _ in range(NOISE_POINTS)]
    return {"f": f, "pg": float(np.max(np.abs(pg))), "classes": classes, "top": top,
            "noise_spread": float(np.max(values) - np.min(values)), "rounding_line": ROUNDING * abs(f),
            "spacing": float(np.spacing(abs(f)))}


# ---------------------------------------------------------------- [2] repairs

def restart(obj, pool, z0, memory):
    t0 = time.monotonic()
    c = optimize.lbfgsb_fit(obj, z0, LO, HI, pgtol=PG_STOP, max_iter=50_000, central_step=R.CENTRAL_STEP,
                            evaluate=pool.values, gradient=obj.gradient, maxcor=memory, max_seconds=RESTART_SECONDS)
    return {"memory": memory, "analytic_pg": c.analytic_pg, "central_pg": c.pg_max, "f": c.f, "phases": c.phases,
            "gradients": c.gradients, "at_lower": len(c.at_lower), "at_upper": len(c.at_upper), "z": c.x,
            "seconds": round(time.monotonic() - t0, 1)}


def hessian(pool, z, free):
    h = STEP * np.maximum(1.0, np.abs(z[free]))
    probes = []
    for k, step in zip(free, h):
        for s in (step, -step):
            v = z.copy()
            v[k] += s
            probes.append(v)
    g = np.array([gk for _, gk in gradients(pool, probes)])
    H = ((g[0::2] - g[1::2]) / (2 * h)[:, None])[:, free].T
    return (H + H.T) / 2, len(probes)


def newton_step(H, g):
    tau, scale = 0.0, max(1.0, float(np.max(np.abs(np.diag(H))))) if len(H) else 1.0
    while True:
        try:
            c = cho_factor(H + tau * np.eye(len(H)))
            break
        except LinAlgError:
            tau = 1e-12 * scale if tau == 0.0 else 10 * tau
    return -cho_solve(c, g), tau


def projected_newton(obj, pool, z0):
    t0 = time.monotonic()
    z = np.array(z0, dtype=float)
    f, g = obj.gradient(z)
    f0, trace, last_pure, ending, hgrads = f, [], None, "step cap", 0
    for step in range(NEWTON_STEPS + 1):
        pg = float(np.max(np.abs(pgrad(z, g))))
        held = ((z <= LO) & (g > 0)) | ((z >= HI) & (g < 0))
        entry = {"step": step, "f": f, "gain": f0 - f, "pg": pg, "held": int(held.sum())}
        trace.append(entry)
        say(f"    step {step}: pg {pg:.3g}, gain {f0 - f:.3g}, held {int(held.sum())}")
        if pg <= PG_STOP:
            ending = "reached the stop"
            break
        if step == NEWTON_STEPS:
            break
        free = np.flatnonzero(~held)
        H, count = hessian(pool, z, free)
        hgrads += count
        if step == 0:
            ev = np.linalg.eigvalsh(H)
            entry.update({"eig_min": float(ev[0]), "eig_max": float(ev[-1])})
        keep = np.ones(len(free), dtype=bool)
        while True:   # hold a free coordinate at a bound that the step would push out
            sub = free[keep]
            d_sub, tau = newton_step(H[np.ix_(keep, keep)], g[sub])
            out = ((z[sub] <= LO) & (d_sub < 0)) | ((z[sub] >= HI) & (d_sub > 0))
            if not out.any():
                break
            keep[np.flatnonzero(keep)[out]] = False
        d = np.zeros_like(z)
        d[sub] = d_sub
        ratio = np.full(len(z), np.inf)
        ratio[d < 0] = (LO - z[d < 0]) / d[d < 0]
        ratio[d > 0] = (HI - z[d > 0]) / d[d > 0]
        a_block = float(np.min(ratio))
        slope = float(g @ d)
        decrement = -slope / 2
        rounding = ROUNDING * max(1.0, abs(f))
        pure = decrement <= rounding
        entry.update({"decrement": decrement, "pure": pure, "damping": tau, "a_block": a_block,
                      "held_for_step": int(len(free) - keep.sum())})
        if pure and last_pure is not None and decrement > 0.5 * last_pure:
            ending = "stalled below rounding"
            break
        alpha = min(1.0, a_block)
        while alpha >= ALPHA_MIN:
            zn = np.clip(z + alpha * d, LO, HI)
            fn, gn = obj.gradient(zn)
            if (fn <= f + rounding) if pure else (fn <= f + ARMIJO * alpha * slope):
                break
            alpha *= 0.5
        else:
            ending = "no acceptable step"
            break
        entry["alpha"] = alpha
        last_pure = decrement if pure and alpha == 1.0 else None
        z, f, g = zn, fn, gn
    return {"ending": ending, "steps": len(trace) - 1, "first_pg": trace[0]["pg"], "final_pg": trace[-1]["pg"],
            "gain": f0 - f, "hessian_gradients": hgrads, "seconds": round(time.monotonic() - t0, 1),
            "trace": trace, "z": z}


# ---------------------------------------------------------------- report

def report(out, prov, res):
    g = lambda v: f"{v:.3g}" if isinstance(v, float) else str(v)  # noqa: E731
    L = ["# Z climb diagnosis", "",
         f"Exploratory. Generated by `{SCRIPT}` (sha256 {prov['script_sha256'][:12]}) on {prov['started_utc']} at "
         f"commit `{prov['git_commit']}` (clean: {prov['git_clean']}). Not a test.", "",
         f"Calibration: `{prov['calibrated']['path']}` (commit {prov['calibrated']['commit'][:8]}).", "",
         "## [0] Rebuilt", "", f"Data identical: {res['data_identical']}.", "",
         "| Climb | recorded lnL | lnL here | equal | recorded analytic pg | here |", "|---|---|---|---|---|---|"]
    for name, r in res["climbs"].items():
        L.append(f"| {name} | {r['recorded_lnl']!r} | {r['lnl']!r} | {r['lnl_equal']} | {g(r['recorded_pg'])} | "
                 f"{g(r['pg'])} |")
    L += ["", f"Z_e against Z_e5: max |dz| {g(res['ze_vs_ze5']['max_abs_dz'])}, lnL difference "
          f"{g(res['ze_vs_ze5']['lnl_difference'])}.", ""]
    for name, r in res["climbs"].items():
        e = r["endpoint"]
        L += [f"## [1] {name}'s endpoint", "",
              f"f {e['f']!r}; rounding: spread of f over {NOISE_POINTS} points moved {NOISE_MOVE:g} is "
              f"{g(e['noise_spread'])}; spacing of f {g(e['spacing'])}; flux.ROUNDING line {g(e['rounding_line'])}.", "",
              "| class | count | max \\|projected\\| | max \\|g\\| |", "|---|---|---|---|"]
        for c, v in sorted(e["classes"].items()):
            L.append(f"| {c} | {v['count']} | {g(v['max_abs_projected'])} | {g(v['max_abs_gradient'])} |")
        L += ["", "| k | (i, j) | z | to bound | gradient | projected | class |", "|---|---|---|---|---|---|---|"]
        for t in e["top"]:
            L.append(f"| {t['k']} | ({t['i']}, {t['j']}) | {g(t['z'])} | {g(t['to_bound'])} | {g(t['gradient'])} | "
                     f"{g(t['projected'])} | {t['class']} |")
        rs, pn = r["restart"], r["newton"]
        s0 = pn["trace"][0]
        L += ["", f"## [2] Repairs from {name}", "",
              f"(a) L-BFGS-B restarted (memory {rs['memory']}, {RESTART_SECONDS} s cap): analytic pg {g(rs['analytic_pg'])}, "
              f"central pg {g(rs['central_pg'])}, gain {g(r['endpoint']['f'] - rs['f'])}, {rs['gradients']} gradients, "
              f"{rs['seconds']} s; {rs['phases']}.", "",
              f"(b) projected Newton: {pn['ending']} after {pn['steps']} steps, pg {g(pn['first_pg'])} -> "
              f"{g(pn['final_pg'])}, gain {g(pn['gain'])}, {pn['hessian_gradients']} Hessian gradients, {pn['seconds']} s. "
              f"Step 0: decrement {g(s0.get('decrement'))} (pure: {s0.get('pure')}), free Hessian eigenvalues "
              f"{g(s0.get('eig_min'))} to {g(s0.get('eig_max'))}, damping {g(s0.get('damping'))}.", "",
              "| step | pg | gain | held | decrement | pure | alpha |", "|---|---|---|---|---|---|---|"]
        for t in pn["trace"]:
            L.append(f"| {t['step']} | {g(t['pg'])} | {g(t['gain'])} | {t['held']} | {g(t.get('decrement', '-'))} | "
                     f"{t.get('pure', '-')} | {g(t.get('alpha', '-'))} |")
        L.append("")
    if "newton_endpoints" in res:
        L += [f"Projected Newton's endpoints from Z_e and Z_e5: max |dz| {g(res['newton_endpoints']['max_abs_dz'])}, "
              f"lnL difference {g(res['newton_endpoints']['lnl_difference'])}."]
    write(out / "report.md", "\n".join(L) + "\n")


# ---------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--calibrated", required=True, help="the failed calibrate run's directory")
    ap.add_argument("--out", required=True)
    ap.add_argument("--git", default="git")
    ap.add_argument("--workers", type=int, default=10)
    ap.add_argument("--allow-dirty", action="store_true")
    a = ap.parse_args()
    prov = R.provenance(a, "z-climb-diagnosis")
    cal_dir = Path(a.calibrated).expanduser().resolve()
    cal = json.loads(read(cal_dir / "results.json"))
    prov.update({"script": SCRIPT, "label": "exploratory", "script_sha256": sha256(Path(__file__).resolve()),
                 "calibrated": {"path": str(cal_dir), "sha256": sha256(cal_dir / "results.json"),
                                "commit": cal["provenance"]["git_commit"]}})
    out = R.new_out(a.out)
    p = R.Problem(cal["case"], cal_dir / "data")
    made = cal["made"]["sha256"]
    res = {"data_identical": {k: sha256(AA if k == "aa_example.phy" else p.data / k) for k in made} == made,
           "climbs": {}}
    say(f"[0] data identical: {res['data_identical']}")
    if not res["data_identical"]:
        sys.exit("[0] the calibration's data changed; stopping")
    obj = p.objective()
    pool = R.Pool(obj, a.workers)
    rng = np.random.default_rng(R.SEED)
    for name, memory in (("Z_e", 20), ("Z_e5", 5)):
        c = cal["climbs"][name]
        z = np.array(c["z"])
        f, gz = obj.gradient(z)
        r = {"recorded_lnl": c["lnl"], "lnl": -f, "lnl_equal": -f == c["lnl"], "recorded_pg": c["analytic_pg"],
             "pg": float(np.max(np.abs(pgrad(z, gz))))}
        say(f"[0] {name}: lnL {-f!r} (recorded {c['lnl']!r}), pg {r['pg']:.3g}")
        r["endpoint"] = endpoint(obj, z, p.ref, rng)
        e = r["endpoint"]
        say(f"[1] {name}: f noise spread {e['noise_spread']:.3g}, rounding line {e['rounding_line']:.3g}; "
            + ", ".join(f"{k} {v['count']} (pg {v['max_abs_projected']:.3g})" for k, v in e["classes"].items()))
        r["restart"] = restart(obj, pool, z, memory)
        say(f"[2a] {name} restarted: pg {r['restart']['analytic_pg']:.3g}, {r['restart']['phases']}")
        say(f"[2b] {name} projected Newton:")
        r["newton"] = projected_newton(obj, pool, z)
        say(f"[2b] {name}: {r['newton']['ending']}, pg {r['newton']['final_pg']:.3g}, gain {r['newton']['gain']:.3g}")
        res["climbs"][name] = r
    za, zb = (np.array(cal["climbs"][k]["z"]) for k in ("Z_e", "Z_e5"))
    res["ze_vs_ze5"] = {"max_abs_dz": float(np.max(np.abs(za - zb))),
                        "lnl_difference": cal["climbs"]["Z_e"]["lnl"] - cal["climbs"]["Z_e5"]["lnl"]}
    na, nb = (res["climbs"][k]["newton"] for k in ("Z_e", "Z_e5"))
    res["newton_endpoints"] = {"max_abs_dz": float(np.max(np.abs(na["z"] - nb["z"]))),
                               "lnl_difference": obj.gradient(nb["z"])[0] - obj.gradient(na["z"])[0]}
    pool.close()
    prov["finished_utc"] = time.strftime("%Y-%m-%dT%H:%M:%S+00:00", time.gmtime())
    write(out / "results.json", json.dumps({"provenance": prov, **res}, indent=1, default=R._json))
    report(out, prov, res)
    say(f"wrote {out / 'report.md'}")


if __name__ == "__main__":
    main()
