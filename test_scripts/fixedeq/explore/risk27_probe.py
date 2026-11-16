"""Exploratory (risk 27): why aa_example's z climb Z_e met its first-order test 0.299 below the box's
peak, and whether the z climb can certify a peak in SIM-1-type problems. Not a test; no pass lines.

  python risk27_probe.py --calibrated AA_CAL_DIR --validate VALIDATE_DIR --out FILE.json [--parts ...]

Parts (in order): ze (flux-space KKT at Z_e, N_r, N_d), exact (z climbs from LG at pi* on S5's V2 and
V3 exact problems), seed (z climb from LG at pi* on aa), tight (z climb from Z_e at pgtol 1e-9, capped),
peak (flux fit from Z_e, then Z_e against the peak).
"""
import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse  # noqa: E402
import datetime  # noqa: E402
import hashlib  # noqa: E402
import json  # noqa: E402
import subprocess  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

REPO = Path("/mnt/c/Users/peter/BioInformatics/iqtree3")
sys.path.insert(0, str(REPO / "test_scripts/fixedeq/reference"))
import reference as R  # noqa: E402
from oracle import cases, chart, flux, nqc, optimize, target  # noqa: E402

STATES = "ARNDCQEGHILKMFPSTWYV"
OUT = {}


def log(*a):
    print(*a, flush=True)


def pairs(ref):
    return [(i, int(j)) for i, free in enumerate(chart.free_destinations(ref)) for j in free]


THRESHOLDS = (0.0, 1e-12, 1e-9, 1e-6)


def fixed_kkt(f, gf, dom):
    """KKT residual with the active set fixed by normalized slack <= t (flux.kkt minimizes over t,
    "all rows" included, which drops complementarity), with nu . slack reported."""
    out = {}
    for t in THRESHOLDS:
        act, eta, mult, fitted = flux._multipliers(f, gf, dom, t)
        out[f"{t:g}"] = {"residual": float(np.max(np.abs(gf - fitted))), "active": int(len(act)),
                         "complementarity": float(mult @ (dom.A[act] @ f)) if len(act) else 0.0}
    return out


def flux_state(obj, pi, ref, z):
    """lnL, flux, d lnL/df and the flux-space KKT record (fixed active thresholds) at the chart point z."""
    Q, nu = chart.build(np.asarray(z), pi, ref)
    f = flux.flux_of(Q, pi)
    lnl, gf = flux.LikelihoodFlux(obj, pi)(f)
    dom = flux.matched(ref)
    k = fixed_kkt(f, gf, dom)
    k["gradient_norm"] = float(np.max(np.abs(gf)))
    k["residual"] = k["1e-12"]["residual"]
    k["threshold"] = 1e-12
    return lnl, f, gf, k, nu, Q


def kkt_detail(f, gf, ref, pi, k, top=8):
    dom = flux.matched(ref)
    act, eta, mult, fitted = flux._multipliers(f, gf, dom, k["threshold"])
    r = gf - fitted
    i, j, _ = flux.t3._edges(len(ref))
    order = np.argsort(-np.abs(r))[:top]
    rows = [{"edge": f"{STATES[i[e]]}->{STATES[j[e]]}", "flux": float(f[e]), "residual": float(r[e]),
             "flux_times_residual": float(f[e] * abs(r[e]))} for e in order]
    return rows


def row_detail(name, z, f, ref, rows=(4, 17)):
    """For rows C and W: the reference destination's flux, then each destination's z and flux."""
    idx = flux.edge_index(len(ref))
    prs = pairs(ref)
    out = {}
    for i in rows:
        r = int(ref[i])
        items = [(STATES[j], float(z[k]), float(f[idx[i, j]])) for k, (ii, j) in enumerate(prs) if ii == i]
        out[STATES[i]] = {"reference": STATES[r], "reference_flux": float(f[idx[i, r]]), "destinations": items}
        log(f"    {name} row {STATES[i]}: reference {STATES[r]} flux {f[idx[i, r]]:.3g}; "
            + ", ".join(f"{d} z {zz:.2f} f {ff:.2g}" for d, zz, ff in items))
    return out


def climb(obj, z0, pgtol, workers, max_seconds=R.CAP_SECONDS):
    pool = R.Pool(obj, workers)
    t0 = time.monotonic()
    c = optimize.lbfgsb_fit(obj, np.asarray(z0), R.LO, R.HI, pgtol=pgtol, max_iter=50_000,
                            central_step=R.CENTRAL_STEP, evaluate=pool.values, gradient=obj.gradient, maxcor=20,
                            max_seconds=max_seconds, gradient_map=pool.gradients)
    pool.close()
    return c, time.monotonic() - t0


def climb_record(c, secs):
    return {"lnl": -c.f, "analytic_pg": c.analytic_pg, "phases": c.phases, "gradients": c.gradients,
            "at_lower": int(len(c.at_lower)), "at_upper": int(len(c.at_upper)), "seconds": round(secs, 1)}


def part_ze(a, cal, p, obj):
    log("== ze: flux-space first-order check at Z_e, N_r and N_d (aa_example)")
    pts = {"Z_e": cal["climbs"]["Z_e"]["z"], "Z_e5": cal["climbs"]["Z_e5"]["z"],
           "N_r": cal["endpoints"]["relaxed"]["z"], "N_d": cal["endpoints"]["default"]["z"]}
    res = {}
    for name, z in pts.items():
        z = np.asarray(z)
        fz, gz = obj.gradient(z)
        pgz = float(np.max(np.abs(optimize.projected_gradient(z, gz, R.LO, R.HI))))
        lnl, f, gf, k, nu, Q = flux_state(obj, p.pi, p.ref, z)
        out = {"lnl": lnl, "z_pg": pgz, "kkt": k, "relative_kkt": k["residual"] / k["gradient_norm"],
               "fluxes_below_1e-10": int(np.sum(f < 1e-10)), "smallest_flux": float(f.min()),
               "at_lower": int(np.sum(z <= R.LO)), "at_upper": int(np.sum(z >= R.HI)),
               "top_residuals": kkt_detail(f, gf, p.ref, p.pi, k)}
        F = np.zeros((20, 20))
        i, j, _ = flux.t3._edges(20)
        F[i, j] = f
        out["state_outflow"] = {STATES[s]: float(F[s].sum()) for s in range(20)}
        out["jump_chain_nu"] = {STATES[s]: float(nu[s]) for s in range(20)}
        res[name] = out
        log(f"{name}: lnL {lnl:.10f}, z pg {pgz:.3g}, gradient norm {k['gradient_norm']:.4g}, fluxes < 1e-10: "
            f"{out['fluxes_below_1e-10']}, smallest {f.min():.3g}, at lower {out['at_lower']}, upper {out['at_upper']}")
        log("    flux KKT by active threshold: " + "; ".join(
            f"{t}: residual {k[t]['residual']:.4g}, active {k[t]['active']}, compl {k[t]['complementarity']:.3g}"
            for t in (f"{x:g}" for x in THRESHOLDS)))
        for row in out["top_residuals"]:
            log(f"    {row['edge']}: flux {row['flux']:.3g}, residual {row['residual']:.4g}, "
                f"flux x |residual| {row['flux_times_residual']:.3g}")
        if name in ("Z_e", "N_r"):
            out["rows"] = row_detail(name, z, f, p.ref)
        small = sorted(out["state_outflow"].items(), key=lambda kv: kv[1])[:4]
        log("    smallest state outflows: " + ", ".join(f"{s} {v:.3g}" for s, v in small))
    OUT["ze"] = res


def exact_truths(pi, ref, Qlg):
    """V2 and V3's truths, by cmd_validate's draws from SEED in order."""
    rng = np.random.default_rng(R.SEED)
    n = len(pi)
    rng.lognormal(0, 1, n * (n - 1))
    R.kkt_truth(pi, ref, rng)
    return {"V2": cases.generator_from_flux(cases.cycle_flux(n, rng), pi), "V3": R.extreme_truth(pi, Qlg, ref, rng)}


def part_exact(a):
    log("== exact: z climbs from LG at pi* on S5's exact problems V2 (interior truth) and V3 (extreme truth)")
    pi = target.parse_target(R.turtle_target())[0]
    Qlg, ref = nqc.lg_seed(pi)
    val = json.loads(Path(a.validate, "results.json").read_text())
    recorded = {row["check"].split()[0]: row for row in val["results"] if "deficit against the truth" in row["check"]}
    res = {}
    for name, Qt in exact_truths(pi, ref, Qlg).items():
        obj = R.exact_objective(pi, ref, Qt)
        truth_lnl = obj.log_likelihood(Qt)
        zt = chart.inverse(Qt, ref)
        log(f"{name}: truth lnL {truth_lnl:.10f} (validate recorded {recorded[name]['truth_lnl']:.10f}), truth at lower "
            f"{int(np.sum(zt <= R.LO + 1e-12))}, upper {int(np.sum(zt >= R.HI - 1e-12))}; flux fit's deficit then "
            f"{recorded[name]['deficit']:.3g}")
        z0 = np.clip(chart.inverse(Qlg, ref), R.LO + 1e-9, R.HI - 1e-9)
        c, secs = climb(obj, z0, R.PG_ANALYTIC, a.workers)
        lnl, f, gf, k, _, _ = flux_state(obj, pi, ref, c.x)
        rec = climb_record(c, secs)
        rec.update({"truth_lnl": truth_lnl, "deficit": truth_lnl - (-c.f), "kkt": k,
                    "relative_kkt": k["residual"] / k["gradient_norm"], "smallest_flux": float(f.min()),
                    "max_abs_dz_to_truth": float(np.max(np.abs(c.x - np.clip(zt, R.LO, R.HI)))),
                    "top_residuals": kkt_detail(f, gf, ref, pi, k, top=5)})
        lt, ft, gft, kt, _, _ = flux_state(obj, pi, ref, np.clip(zt, R.LO, R.HI))
        rec["truth_kkt"] = kt
        rec["truth_relative_kkt"] = kt["residual"] / kt["gradient_norm"]
        res[name] = rec
        log(f"{name} z climb: lnL {-c.f:.10f}, deficit {rec['deficit']:.4g}, pg {c.analytic_pg:.3g}, phases "
            f"{c.phases}, at lower {rec['at_lower']}, upper {rec['at_upper']}, {secs:.0f} s; flux KKT residual "
            f"{k['residual']:.4g} (relative {rec['relative_kkt']:.3g}); at the truth {kt['residual']:.4g} "
            f"(relative {rec['truth_relative_kkt']:.3g}); max |dz| to truth {rec['max_abs_dz_to_truth']:.3g}")
        OUT.setdefault("exact_z", {})[name] = c.x.tolist()
    OUT["exact"] = res


def part_seed(a, p, obj, peak):
    log("== seed: z climb from LG at pi* on aa_example")
    z0 = np.clip(chart.inverse(p.seed, p.ref), R.LO + 1e-9, R.HI - 1e-9)
    c, secs = climb(obj, z0, R.PG_ANALYTIC, a.workers)
    lnl, f, gf, k, _, _ = flux_state(obj, p.pi, p.ref, c.x)
    rec = climb_record(c, secs)
    rec.update({"below_peak": peak - (-c.f), "kkt": k, "relative_kkt": k["residual"] / k["gradient_norm"]})
    log(f"seed z climb: lnL {-c.f:.10f}, below the peak {rec['below_peak']:.4g}, pg {c.analytic_pg:.3g}, phases "
        f"{c.phases}, at lower {rec['at_lower']}, upper {rec['at_upper']}, {secs:.0f} s; flux KKT residual "
        f"{k['residual']:.4g} (relative {rec['relative_kkt']:.3g})")
    OUT["seed"] = rec
    OUT["seed_z"] = c.x.tolist()


def part_tight(a, cal, p, obj, peak):
    log("== tight: z climb from Z_e at pgtol 1e-9, capped at 1800 s")
    c, secs = climb(obj, np.asarray(cal["climbs"]["Z_e"]["z"]), 1e-9, a.workers, max_seconds=1800.0)
    rec = climb_record(c, secs)
    rec["below_peak"] = peak - (-c.f)
    log(f"tight z climb: lnL {-c.f:.10f}, below the peak {rec['below_peak']:.4g}, pg {c.analytic_pg:.3g}, phases "
        f"{c.phases}, at lower {rec['at_lower']}, upper {rec['at_upper']}, {secs:.0f} s")
    OUT["tight"] = rec


def part_peak(a, cal, p, obj):
    log("== peak: flux fit from Z_e (S5's cross-start), then Z_e against the peak")
    pool = R.Pool(obj, a.workers)
    dom = flux.matched(p.ref)
    lf = flux.LikelihoodFlux(obj, p.pi)
    ze = np.asarray(cal["climbs"]["Z_e"]["z"])
    Qz = chart.build(ze, p.pi, p.ref)[0]
    fz = flux.flux_of(Qz, p.pi)
    r = flux.fit(lf, flux.inward(Qz, p.pi, dom)[0], dom, objective_map=pool.fluxes)
    pool.close()
    lnl = obj.log_likelihood(flux.generator(r.f, p.pi))
    zp = flux.z_of_flux(r.f, p.ref)
    log(f"flux fit from Z_e: {r.status}, lnL {lnl:.10f}, {r.seconds:.0f} s, KKT {r.kkt['residual']:.3g}")
    i, j, _ = flux.t3._edges(20)
    ratio = r.f / fz
    order = np.argsort(-np.abs(np.log(ratio)))[:12]
    moves = [{"edge": f"{STATES[i[e]]}->{STATES[j[e]]}", "flux_Z_e": float(fz[e]), "flux_peak": float(r.f[e])}
             for e in order]
    log("largest relative flux changes Z_e -> peak:")
    for m in moves:
        log(f"    {m['edge']}: {m['flux_Z_e']:.3g} -> {m['flux_peak']:.3g}")
    prs = pairs(p.ref)
    dz = zp - ze
    big = np.argsort(-np.abs(dz))[:12]
    log("largest z changes Z_e -> peak: " + ", ".join(f"{STATES[prs[k][0]]}->{STATES[prs[k][1]]} {ze[k]:.2f}->{zp[k]:.2f}"
                                                        for k in big))
    row_detail("peak", zp, r.f, p.ref)
    lnl_p, f_p, gf_p, k_p, _, _ = flux_state(obj, p.pi, p.ref, np.clip(zp, R.LO, R.HI))
    log(f"flux KKT at the peak (fixed thresholds): " + "; ".join(
        f"{t}: residual {k_p[t]['residual']:.3g}, active {k_p[t]['active']}" for t in (f"{x:g}" for x in THRESHOLDS)))
    zprof = [-obj(np.clip((1 - t) * ze + t * zp, R.LO, R.HI)) for t in np.linspace(0, 1, 21)]
    log("lnL along the straight z segment Z_e -> peak: " + ", ".join(f"{v:.4f}" for v in zprof))
    OUT["peak"] = {"status": r.status, "lnl": lnl, "seconds": round(r.seconds, 1), "kkt": r.kkt, "flux_moves": moves,
                   "z_profile": zprof, "z_peak": zp.tolist()}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--calibrated", required=True)
    ap.add_argument("--validate", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--parts", default="ze,exact,seed,tight,peak")
    ap.add_argument("--workers", type=int, default=10)
    a = ap.parse_args()
    me = Path(__file__).resolve()
    commit = subprocess.run(["git.exe", "rev-parse", "HEAD"], cwd=REPO, capture_output=True, text=True).stdout.strip()
    OUT["provenance"] = {"script": str(me), "sha256": hashlib.sha256(me.read_bytes()).hexdigest(), "commit": commit,
                         "started_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
                         "argv": sys.argv, "cpus": os.cpu_count(), "label": "exploratory", "seed": R.SEED}
    log(json.dumps(OUT["provenance"]))
    cal_dir = Path(a.calibrated).expanduser().resolve()
    cal = json.loads((cal_dir / "results.json").read_text())
    p = R.Problem(cal["case"], cal_dir / "data")
    obj = p.objective()
    peak = -6958.732388341286   # F_s = F_e at S5's gate (aa compare results.json)
    parts = a.parts.split(",")
    for part in parts:
        t0 = time.monotonic()
        if part == "ze":
            part_ze(a, cal, p, obj)
        elif part == "exact":
            part_exact(a)
        elif part == "seed":
            part_seed(a, p, obj, peak)
        elif part == "tight":
            part_tight(a, cal, p, obj, peak)
        elif part == "peak":
            part_peak(a, cal, p, obj)
        log(f"-- {part}: {time.monotonic() - t0:.0f} s")
        OUT["provenance"]["finished_utc"] = datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")
        Path(a.out).expanduser().write_text(json.dumps(OUT, indent=1, default=R._json))


if __name__ == "__main__":
    main()
