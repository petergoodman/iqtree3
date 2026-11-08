"""Confirmatory: S5's independent reference (decision 042; CODE_PLAN.md 3.1 and 3.3; lines in decision
021's S5 rows). One step per invocation, each writing results.json and report.md with provenance into
a new directory outside the repository, and exiting 1 on any asserted failure.

  python test_scripts/fixedeq/reference/reference.py validate --out DIR --git git.exe [--workers 10]
  python test_scripts/fixedeq/reference/reference.py calibrate --case CASE --binary BIN --alisim ALISIM \
      [--calibration TOYCAL/results.json] --out DIR --git git.exe [--workers 10]
  python test_scripts/fixedeq/reference/reference.py compare --calibrated DIR --out DIR --git git.exe [--workers 10]

validate: known answers at n = 20, never compared with NQC: V1 (sum N_e log f_e, open, against the
technical plan's dual Newton), V1b (the same with ratio constraints binding at a known optimum),
V2 and V3 (exact pattern frequencies, N = 1e6, on a rooted three-taxon tree, from an interior
cycle-flux truth and from LG at pi* with 160 coordinates at the lower bound and 5 at the upper),
and the analytic gradients (AG) at V2's start.

calibrate: one problem (toy, toy300, aa, twopart): NQC Level 1 at IQ-TREE's defaults (N_d, -vv) and
relaxed (N_r: -me 1e-6, 300-round caps), both with --nqc-domain 1e-5,1e5; the same-objective checks
(A4); the z climbs Z_e and Z_e5 from N_r's endpoint (analytic L-BFGS-B, memory 20 and 5, A5); and
the same-peak line L = max(1e-3, 10 |lnL(Z_e) - lnL(Z_e5)|), all written before any flux fit.

compare: from a calibrate directory: the reference from LG at pi* (F_s) and from N_r's endpoint
(F_e); reported only, two seeded random starts on toy and toy300 and the open domain from F_e's
endpoint; checks A3, A5 to A8, the outcome table with cross-starts on a disagreement, and the
reported differences.
"""
import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse  # noqa: E402
import datetime  # noqa: E402
import itertools  # noqa: E402
import json  # noqa: E402
import math  # noqa: E402
import multiprocessing  # noqa: E402
import platform  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
import traceback  # noqa: E402
from concurrent.futures import ThreadPoolExecutor  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
import scipy  # noqa: E402

HERE = Path(__file__).resolve().parent
for _p in ("", "nqc", "probes", "regression", "differential", "explore"):
    sys.path.insert(0, str(HERE.parent / _p))
from oracle import cases, chart, flux, likelihood as lk, nqc, optimize, readers, target  # noqa: E402
from probes import AA, Runner, baseline_runs, tree_text, turtle_target  # noqa: E402
from regress import REPO, capture, read, sha256, write  # noqa: E402
from nqcchecks import ALPHA, TOL_LNL, check, final_struct, first_gradient, floats, joint_args  # noqa: E402
from convergence import ME_RELAXED, run_state  # noqa: E402
import toy  # noqa: E402
import twopart  # noqa: E402

SCRIPT = Path(__file__).resolve().relative_to(REPO).as_posix()
DOMAIN = "1e-5,1e5"
LO, HI = nqc.LOWER, nqc.UPPER
ROUNDS = 300
RESIDUAL = 1e-10        # 021: residuals at every accepted reference point
PG_ANALYTIC = 1e-6      # A5 (i)
PG_CENTRAL = 1e-3       # A5 (ii), central differences at step 1e-4 max(1, |z|)
CENTRAL_STEP = 1e-4
M_MU = 1e-6             # A5 (iii)
L_FLOOR = 1e-3          # 021's L = max(1e-3, 10 s)
DEFICIT = 1e-3          # V2, V3
GAP = 1e-9              # V1, V1b: duality gap over |F|
RATE_FLOOR = 1e-3       # relative Q differences are reported above this rate
SEED = 20261107
SITES_EXACT = 1e6
THREE = "((a:0.2,b:0.3):0.1,c:0.4);"
CASES = ("toy", "toy300", "aa", "twopart")
CAP_SECONDS = 4 * 3600.0


# ---------------------------------------------------------------- worker pool

_OBJ = None


def _init(obj):
    global _OBJ
    _OBJ = obj


def _flux_work(f):
    return flux.LikelihoodFlux(_OBJ, _OBJ.pi)(np.asarray(f))


def _z_value(z):
    return _OBJ(np.asarray(z))


class Pool:
    """Worker processes holding one objective; spreads the Newton Hessian's probes and the central
    differences, which changes no value."""

    def __init__(self, obj, workers):
        self.workers = workers
        self.pool = multiprocessing.get_context("fork").Pool(workers, _init, (obj,)) if workers > 1 else None
        if self.pool is None:
            _init(obj)

    def _map(self, fn, xs):
        if self.pool is None:
            return [fn(x) for x in xs]
        return self.pool.map(fn, xs, chunksize=max(1, len(xs) // (4 * self.workers)))

    def fluxes(self, fs):
        return self._map(_flux_work, fs)

    def values(self, zs):
        return self._map(_z_value, zs)

    def close(self):
        if self.pool is not None:
            self.pool.close()


# ---------------------------------------------------------------- shared helpers

def new_out(path):
    out = Path(path).expanduser().resolve()
    if out.exists() or out.is_relative_to(REPO):
        sys.exit(f"{out}: needs a new directory outside the repository")
    out.mkdir(parents=True)
    return out


def provenance(a, step):
    prov = {"script": SCRIPT, "step": step, "command": sys.argv, "label": "confirmatory",
            "started_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
            "git_commit": capture(a.git, "rev-parse", "HEAD", cwd=REPO),
            "git_clean": capture(a.git, "status", "--porcelain", cwd=REPO) == "",
            "python": sys.version, "numpy": np.__version__, "scipy": scipy.__version__,
            "platform": platform.platform(), "workers": a.workers, "seed": SEED}
    if not prov["git_clean"] and not a.allow_dirty:
        sys.exit(f"the repository is not clean ({prov['git_commit']}); commit first or pass --allow-dirty")
    return prov


def finish(out, prov, title, results, extra=None):
    prov["finished_utc"] = datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")
    write(out / "results.json", json.dumps({"provenance": prov, "results": results, **(extra or {})}, indent=1,
                                           default=_json))
    lines = [f"# {title}", "", f"Confirmatory. Generated by `{SCRIPT}` ({prov['step']}) on {prov['started_utc']} at "
             f"commit `{prov['git_commit']}` (clean: {prov['git_clean']}).", "",
             "| Check | Result | Values |", "|---|---|---|"]
    for row in results:
        verdict = {True: "pass", False: "FAIL", None: "report"}[row["passed"]]
        values = {k: v for k, v in row.items() if k not in ("check", "passed")}
        lines.append(f"| {row['check']} | {verdict} | `{json.dumps(values, default=_json)[:700]}` |")
    write(out / "report.md", "\n".join(lines) + "\n")
    failed = [row["check"] for row in results if row["passed"] is False]
    for f in failed:
        print(f"  FAIL {f}")
    print(f"{len(failed)} asserted checks failed; wrote {out / 'report.md'}")
    sys.exit(1 if failed else 0)


def _json(v):
    if isinstance(v, np.ndarray):
        return v.tolist()
    if isinstance(v, (np.floating, np.integer, np.bool_)):
        return v.item()
    return str(v)


def fit_record(r):
    """A flux fit's record: status, value, the declared stop's quantities and the counts."""
    return {"status": r.status, "message": r.message, "lnl": r.value, "start_lnl": r.start_value,
            "projected_z_gradient": r.pg, "mu": r.mu, "m_active": r.m_active, "m_mu": r.m_active * r.mu,
            "kkt": r.kkt, "runs": len(r.runs), "seconds": round(r.seconds, 1), **r.totals(),
            "run_endings": [run["ending"] for run in r.runs], "run_gains": [run["gain"] for run in r.runs],
            "finish_steps": [run.get("finish_steps") for run in r.runs],
            "finish_released": [run.get("released") for run in r.runs],
            "smallest_flux": float(np.min(r.f))}


def same_peak_line(lnl_e, lnl_e5):
    """021's L = max(1e-3, 10 s), s the spread of the two z climbs."""
    return max(L_FLOOR, 10 * abs(lnl_e - lnl_e5))


def label_pair(lnl_h, lnl_w, from_w, from_h, profile_min, L):
    """The outcome table's diagnosis of a failed same-peak row (CODE_PLAN.md 3.3), H the higher and W
    the lower endpoint. from_w: (lnL, converged) of H's fitter started at W's endpoint; from_h: W's
    fitter started at H's endpoint; profile_min: the lowest lnL along the feasible segment between
    them. The failed row stays failed whatever the label."""
    (hw, hw_ok), (wh, wh_ok) = from_w, from_h
    if not (hw_ok and wh_ok) or hw < lnl_w - L or wh < lnl_h - L:
        return "cross-start unconverged"
    if abs(hw - lnl_h) <= L and abs(wh - lnl_h) <= L:
        return "W stopped short"
    if wh > lnl_h + L:
        return "H not at a peak"
    if abs(hw - lnl_w) <= L and abs(wh - lnl_h) <= L:
        return "distinct stationary points, one ridge" if profile_min >= lnl_w - L else \
            "distinct stationary points, distinct peaks"
    return "undetermined"


def flux_lnl(obj, p, Q):
    """The value the flux fit climbs: log_likelihood_gradient's lnL at Q's flux (A4)."""
    return float(flux.LikelihoodFlux(obj, p.pi)(flux.flux_of(Q, p.pi))[0])


def segment_profile(obj, pi, fa, fb, points=21):
    """lnL along the feasible segment between two fluxes (P_box is convex)."""
    return [obj.log_likelihood(flux.generator((1 - t) * fa + t * fb, pi)) for t in np.linspace(0, 1, points)]


def z_central(pool, z):
    """Central-difference gradient of -lnL in z at step CENTRAL_STEP max(1, |z|)."""
    return optimize.central_gradient(pool.values, np.asarray(z), CENTRAL_STEP)


def a5_rows(name, z, pg_analytic, pool, m_mu=None):
    g = z_central(pool, z)
    pg_c = float(np.max(np.abs(optimize.projected_gradient(z, g, LO, HI))))
    rows = [check(f"A5 {name}: analytic projected z-gradient <= {PG_ANALYTIC:g}", pg_analytic <= PG_ANALYTIC,
                  value=pg_analytic),
            check(f"A5 {name}: central-difference projected z-gradient <= {PG_CENTRAL:g}", pg_c <= PG_CENTRAL,
                  value=pg_c)]
    if m_mu is not None:
        rows.append(check(f"A5 {name}: m mu <= {M_MU:g}", m_mu <= M_MU, value=m_mu))
    return rows, all(r["passed"] for r in rows)


# ---------------------------------------------------------------- validate

def exact_objective(pi, ref, Qt):
    tree = readers.parse_newick(THREE)
    E = np.eye(len(pi))
    cols = [{"a": E[x], "b": E[y], "c": E[w]} for x, y, w in itertools.product(range(len(pi)), repeat=3)]
    weights = SITES_EXACT * np.exp(lk.site_log_likelihoods(tree, cols, Qt, pi))
    return nqc.Objective.from_patterns(tree, cols, weights, pi, ref)


def kkt_truth(pi, ref, rng, n_lo=40, n_hi=10):
    """V1b's known optimum: f* with n_lo coordinates at the lower bound and n_hi at the upper, and
    N = f* (C^T eta - A^T nu), nu > 0 on exactly those rows."""
    dom = flux.matched(ref)
    z = flux.z_of_flux(flux.flux_of(cases.generator_from_flux(cases.cycle_flux(len(pi), rng), pi), pi), ref)
    m = len(z)
    pick = rng.choice(m, n_lo + n_hi, replace=False)
    z[pick[:n_lo]], z[pick[n_lo:]] = LO, HI
    fstar = flux.flux_of(chart.build(z, pi, ref)[0], pi)
    eta = np.r_[rng.normal(0, 50, len(pi) - 1), 1000.0]
    nu = np.zeros(2 * m)
    nu[pick[:n_lo]] = rng.uniform(100, 300, n_lo)
    nu[m + pick[n_lo:]] = rng.uniform(1e-3, 3e-3, n_hi)
    w = flux.geometry(len(pi))[2].T @ eta - dom.A.T @ nu
    if not np.all(w > 0):
        sys.exit("V1b construction: a non-positive dual denominator")
    return fstar, fstar * w


def extreme_truth(pi, Qlg, ref, rng, n_lo=160, n_hi=5):
    z = chart.inverse(Qlg, ref)
    pick = rng.choice(len(z), n_lo + n_hi, replace=False)
    z[pick[:n_lo]], z[pick[n_lo:]] = LO, HI
    return chart.build(z, pi, ref)[0]


def cmd_validate(a):
    prov = provenance(a, "validate")
    out = new_out(a.out)
    pi = target.parse_target(turtle_target())[0]
    Qlg, ref = nqc.lg_seed(pi)
    rows, records = [], {}
    rng = np.random.default_rng(SEED)
    n = 20

    # V1
    N = rng.lognormal(0, 1, n * (n - 1))
    dom = flux.open_domain(ref)
    r = flux.fit(flux.Concave(N), flux.inward(Qlg, pi, dom)[0], dom)
    gap = flux.duality_gap(r.f, N, dom, r.mu)
    fd = flux.dual_newton(N, n)
    records["V1"] = fit_record(r)
    rows += [check("V1 converged", r.status == "converged", **records["V1"]),
             check("V1 duality gap <= 1e-9 |F|", gap <= GAP * abs(r.value), gap=gap, line=GAP * abs(r.value)),
             check("V1 against the dual Newton optimum (report)", None, flux_max_abs=float(np.max(np.abs(r.f - fd))),
                   objective_difference=float(N @ np.log(fd)) - r.value),
             check("V1 residuals <= 1e-10", max(flux.residuals(r.f, pi).values()) <= RESIDUAL,
                   **flux.residuals(r.f, pi))]
    print(f"V1: {r.status}, gap {gap:.3g}, {r.seconds:.0f} s", flush=True)

    # V1b
    fstar, Nb = kkt_truth(pi, ref, rng)
    dom = flux.matched(ref)
    f0, moved = flux.inward(Qlg, pi, dom)
    r = flux.fit(flux.Concave(Nb), f0, dom)
    gap = flux.duality_gap(r.f, Nb, dom, r.mu)
    records["V1b"] = fit_record(r)
    rows += [check("V1b converged", r.status == "converged", **records["V1b"]),
             check("V1b duality gap <= 1e-9 |F|", gap <= GAP * abs(r.value), gap=gap, line=GAP * abs(r.value)),
             check("V1b against the known optimum (report)", None, flux_max_abs=float(np.max(np.abs(r.f - fstar))),
                   objective_difference=float(Nb @ np.log(fstar)) - r.value, binding_rows_found=r.m_active,
                   start_moved=moved),
             check("V1b residuals <= 1e-10", max(flux.residuals(r.f, pi).values()) <= RESIDUAL,
                   **flux.residuals(r.f, pi))]
    print(f"V1b: {r.status}, gap {gap:.3g}, {r.seconds:.0f} s", flush=True)

    # V2, V3 and AG at V2's start
    truths = {"V2": cases.generator_from_flux(cases.cycle_flux(n, rng), pi),
              "V3": extreme_truth(pi, Qlg, ref, rng)}
    for name, Qt in truths.items():
        obj = exact_objective(pi, ref, Qt)
        dom = flux.matched(ref)
        if name == "V2":
            rows += ag_rows(obj, pi, ref, Qlg, rng)
        pool = Pool(obj, a.workers)
        f0, moved = flux.inward(Qlg, pi, dom)
        r = flux.fit(flux.LikelihoodFlux(obj, pi), f0, dom, objective_map=pool.fluxes)
        pool.close()
        truth_lnl = obj.log_likelihood(Qt)
        lnl = obj.log_likelihood(flux.generator(r.f, pi))
        records[name] = fit_record(r)
        rows += [check(f"{name} truth inside the box", dom.contains(flux.flux_of(Qt, pi)) if name == "V2" else None,
                       smallest_slack=float(np.min(dom.slack(flux.flux_of(Qt, pi))))),
                 check(f"{name} converged", r.status == "converged", **records[name]),
                 check(f"{name} lnL deficit against the truth <= 1e-3", truth_lnl - lnl <= DEFICIT,
                       deficit=truth_lnl - lnl, truth_lnl=truth_lnl, lnl=lnl, start_moved=moved),
                 check(f"{name} residuals <= 1e-10", max(flux.residuals(r.f, pi).values()) <= RESIDUAL,
                       **flux.residuals(r.f, pi))]
        print(f"{name}: {r.status}, deficit {truth_lnl - lnl:.3g}, {r.seconds:.0f} s", flush=True)
    finish(out, prov, "S5 validation", rows, {"records": records})


def ag_rows(obj, pi, ref, Q, rng):
    """AG: analytic gradients in u (flux) and z against central differences at h and 2h, h = 1e-4,
    on five coordinates and five random directions, at Q."""
    z = chart.inverse(Q, ref)
    f = flux.flux_of(Q, pi)
    basis = flux.Basis(f, flux.matched(ref))
    lf = flux.LikelihoodFlux(obj, pi)
    dirs = list(np.eye(len(z))[rng.choice(len(z), 5, replace=False)])
    dirs += [d / np.linalg.norm(d) for d in rng.normal(size=(5, len(z)))]
    u = flux.gradient_check(lambda v: -lf(basis.flux(v))[0], np.zeros(len(z)), -basis.chain(lf(f)[1]), dirs)
    zz = flux.gradient_check(obj, z, obj.gradient(z)[1], dirs)
    return [check(f"AG in {name}: |g.d - D(h)| within 021's line on 10 directions", all(e <= l for e, l in res),
                  worst_ratio=max(e / l for e, l in res)) for name, res in (("u", u), ("z", zz))]


# ---------------------------------------------------------------- problems

class Problem:
    """One comparison problem: its data (made once by calibrate, read again by compare), the NQC
    command and the oracle's objective."""

    def __init__(self, case, data):
        self.case, self.data = case, Path(data)
        self.target_text = turtle_target()
        self.pi = target.parse_target(self.target_text)[0]
        self.seed, self.ref = nqc.lg_seed(self.pi)
        self.model = f"NQC+F{{{self.target_text}}}"
        self.mode = "-S" if case == "twopart" else "single"

    def make(self, alisim):
        """Write the data into the new directory self.data; returns what to record."""
        self.data.mkdir()
        if self.case in ("toy", "toy300"):
            made = toy.make(self.data, alisim, self.pi)
            if self.case == "toy300":
                aln = readers.read_phylip(self.data / "toy.phy")
                write(self.data / "toy300.phy", f"{len(aln)} 300\n" + "".join(f"{k} {s[:300]}\n" for k, s in aln.items()))
            return made
        if self.case == "aa":
            write(self.data / "input.tree", tree_text(baseline_runs(), "run02"))
            return {"sha256": {"input.tree": sha256(self.data / "input.tree"), "aa_example.phy": sha256(AA)}}
        made = twopart.make(self.data, alisim, self.pi)
        write(self.data / "partition.trees", read(self.data / "twopart.nwk") * 2)
        return made

    def command(self, relaxed):
        extra = ["--nqc-domain", DOMAIN] + (["-me", f"{ME_RELAXED:g}", "--loop-model", str(ROUNDS), "-nparam",
                                             str(ROUNDS + 2)] if relaxed else [])
        d = self.data
        if self.case in ("toy", "toy300"):
            return ["-s", d / f"{self.case}.phy", "-m", self.model, "-te", d / "toy.nwk", "-blfix", "-vv", *extra]
        if self.case == "aa":
            return ["-s", AA, "-m", self.model + f"+G4{{{ALPHA}}}", "-te", d / "input.tree", "-blfix", "-vv", *extra]
        return joint_args(d / "twopart.phy", "-S", d / "twopart.nex", self.model, f"+G4{{{ALPHA}}}",
                          d / "partition.trees", "-blfix", "-vv", *extra)

    def objective(self):
        d = self.data
        if self.case in ("toy", "toy300"):
            return nqc.Objective(readers.read_phylip(d / f"{self.case}.phy"), readers.read_newick(d / "toy.nwk"),
                                 self.pi, self.ref)
        if self.case == "aa":
            return nqc.Objective(readers.read_phylip(AA), readers.read_newick(d / "input.tree"), self.pi, self.ref,
                                 shape=ALPHA)
        aln, tree = readers.read_phylip(d / "twopart.phy"), readers.read_newick(d / "twopart.nwk")
        return nqc.LinkedObjective([nqc.Objective({k: "".join(s[i] for i in sites) for k, s in aln.items()}, tree,
                                                  self.pi, self.ref, shape=ALPHA)
                                    for sites in readers.read_charsets(d / "twopart.nex").values()])


# ---------------------------------------------------------------- calibrate

def cmd_calibrate(a):
    prov = provenance(a, "calibrate")
    if a.case in ("toy", "toy300") and not a.calibration:
        sys.exit("toy and toy300 need --calibration (explore/toy_calibration.py's results.json)")
    out = new_out(a.out)
    binary = Path(a.binary).expanduser().resolve()
    alisim = Path(a.alisim).expanduser().resolve()
    prov.update({"case": a.case, "binary": {"path": str(binary), "sha256": sha256(binary)},
                 "alisim": {"path": str(alisim), "sha256": sha256(alisim)}, "domain": DOMAIN, "rounds": ROUNDS})
    p = Problem(a.case, out / "data")
    made = p.make(alisim)
    rows = []
    if a.case in ("toy", "toy300"):
        cal = json.loads(read(a.calibration))
        rows.append(check("toy data identical to the toy calibration's", made["sha256"] == cal["provenance"]["toy"]["sha256"]))
    obj = p.objective()
    runner = Runner(binary, out, a.timeout)
    with ThreadPoolExecutor(max_workers=2) as ex:
        runs = dict(zip(("default", "relaxed"), ex.map(lambda s: runner.run(f"nqc-{s}", p.command(s == "relaxed")),
                                                         ("default", "relaxed"))))
    states, ends = {}, {}
    for setting, run in runs.items():
        states[setting] = run_state(run, p.mode, True, setting, ROUNDS)
        ck = final_struct(run) if run.code == 0 else None
        rows.append(check(f"N_{setting[0]} exits 0 with an NQC checkpoint", ck is not None, exit_code=run.code))
        if ck is None:
            continue
        z = floats(ck["coordinates"])
        ends[setting] = {"z": z.tolist(), "accepted_score": float(ck["accepted_score"]),
                         "ended_by": states[setting].get("ended_by"), "rounds": states[setting].get("rounds")}
        lnl_flux = flux_lnl(obj, p, chart.build(z, p.pi, p.ref)[0])
        rows.append(check(f"A4: the reference's lnL (flux map) at N_{setting[0]}'s endpoint against accepted_score",
                          abs(lnl_flux - float(ck["accepted_score"])) <= TOL_LNL * abs(float(ck["accepted_score"])),
                          reference=lnl_flux, chart_objective=-obj(z), iqtree=float(ck["accepted_score"])))
        if ck.get("references") is not None:
            rows.append(check(f"N_{setting[0]} references equal the oracle's", list(map(int, ck["references"].split(",")))
                              == [int(v) for v in p.ref]))
        if ck.get("target") is not None:
            rows.append(check(f"N_{setting[0]} target bitwise equal to the parsed one",
                              np.array_equal(floats(ck["target"]), p.pi)))
    got = first_gradient(runs["default"].stdout()) if runs["default"].code == 0 else None
    start = -flux_lnl(obj, p, p.seed)   # -lnL, as the -vv head prints it
    rows.append(check("A4: the reference's -lnL (flux map) at LG at pi* against N_d's first -lnL (-vv)",
                      got is not None and abs(got[0] - start) <= TOL_LNL * abs(start),
                      reference=start, chart_objective=obj(chart.inverse(p.seed, p.ref)), iqtree=got[0] if got else None))
    climbs = {}
    if "relaxed" in ends:
        pool = Pool(obj, a.workers)
        for name, memory in (("Z_e", 20), ("Z_e5", 5)):
            t0 = time.monotonic()
            c = optimize.lbfgsb_fit(obj, np.array(ends["relaxed"]["z"]), LO, HI, pgtol=PG_ANALYTIC, max_iter=50_000,
                                    central_step=CENTRAL_STEP, evaluate=pool.values, gradient=obj.gradient,
                                    maxcor=memory, max_seconds=CAP_SECONDS)
            climbs[name] = {"lnl": -c.f, "z": c.x.tolist(), "analytic_pg": c.analytic_pg, "central_pg": c.pg_max,
                            "phases": c.phases, "gradients": c.gradients, "evaluations": c.evaluations,
                            "at_lower": len(c.at_lower), "at_upper": len(c.at_upper),
                            "seconds": round(time.monotonic() - t0, 1)}
            rows += [check(f"A5 {name}: analytic projected z-gradient <= {PG_ANALYTIC:g}", c.analytic_pg <= PG_ANALYTIC,
                           value=c.analytic_pg, phases=c.phases),
                     check(f"A5 {name}: central-difference projected z-gradient <= {PG_CENTRAL:g}",
                           c.pg_max <= PG_CENTRAL, value=c.pg_max)]
            print(f"{name}: lnL {-c.f:.10f}, analytic pg {c.analytic_pg:.3g}, {climbs[name]['seconds']} s", flush=True)
        pool.close()
    L = same_peak_line(climbs["Z_e"]["lnl"], climbs["Z_e5"]["lnl"]) if len(climbs) == 2 else None
    rows.append(check("L recorded before any flux fit (report)", None, L=L,
                      spread=abs(climbs["Z_e"]["lnl"] - climbs["Z_e5"]["lnl"]) if L else None))
    finish(out, prov, f"S5 calibration: {a.case}", rows,
           {"case": a.case, "made": made, "states": states, "endpoints": ends, "climbs": climbs, "L": L,
            "runs": runner.runs})


# ---------------------------------------------------------------- compare

def cmd_compare(a):
    prov = provenance(a, "compare")
    cal_dir = Path(a.calibrated).expanduser().resolve()
    cal = json.loads(read(cal_dir / "results.json"))
    prov["calibrated"] = {"path": str(cal_dir), "sha256": sha256(cal_dir / "results.json"),
                          "commit": cal["provenance"]["git_commit"]}
    out = new_out(a.out)
    case, L = cal["case"], cal["L"]
    if L is None or "relaxed" not in cal["endpoints"] or "default" not in cal["endpoints"]:
        sys.exit("the calibration has no L or no NQC endpoints; nothing to compare")
    p = Problem(case, cal_dir / "data")
    rows, fits = [], {}
    made = cal["made"]["sha256"]
    found = {k: sha256(AA if k == "aa_example.phy" else p.data / k) for k in made}
    rows += [check("data identical to the calibration's", found == made, found=found),
             check("calibrated at this commit", cal["provenance"]["git_commit"] == prov["git_commit"],
                   calibrated=cal["provenance"]["git_commit"])]
    obj = p.objective()
    pi, ref = p.pi, p.ref
    pool = Pool(obj, a.workers)
    dom = flux.matched(ref)
    lf = flux.LikelihoodFlux(obj, pi)

    def save_partial():
        write(out / "partial.json", json.dumps({"rows": rows, "fits": {k: {x: y for x, y in v.items() if x != "f"}
                                                                       for k, v in fits.items()}},
                                               indent=1, default=_json))

    def run_fit(name, Q, domain=dom):
        f0, moved = flux.inward(Q, pi, domain)
        shift = obj.log_likelihood(flux.generator(f0, pi)) - obj.log_likelihood(Q)
        r = flux.fit(lf, f0, domain, objective_map=pool.fluxes)
        rec = fit_record(r)
        rec.update({"start_moved": moved, "start_shift": shift,
                    "lnl_documented": obj.log_likelihood(flux.generator(r.f, pi)), "f": r.f})
        res = flux.residuals(r.f, pi)
        rows.append(check(f"A3 {name}: residuals <= 1e-10", max(res.values()) <= RESIDUAL, **res))
        print(f"{name}: {r.status}, lnL {rec['lnl_documented']:.10f}, {r.seconds:.0f} s", flush=True)
        fits[name] = rec
        save_partial()
        return r, rec

    def converged_row(name, r, rec):
        z = flux.z_of_flux(r.f, ref)
        a5, ok = a5_rows(name, z, r.pg, pool, rec["m_mu"])
        rows.extend(a5)
        rows.append(check(f"{name} ended by the declared stop", r.status == "converged", message=r.message))
        return ok and r.status == "converged"

    q_of = lambda z: chart.build(np.asarray(z), pi, ref)[0]  # noqa: E731
    z_r = np.array(cal["endpoints"]["relaxed"]["z"])
    r_s, F_s = run_fit("F_s", p.seed)
    ok_s = converged_row("F_s", r_s, F_s)
    r_e, F_e = run_fit("F_e", q_of(z_r))
    ok_e = converged_row("F_e", r_e, F_e)
    ze, ze5 = cal["climbs"]["Z_e"], cal["climbs"]["Z_e5"]
    ok_z = all(c["analytic_pg"] <= PG_ANALYTIC and c["central_pg"] <= PG_CENTRAL for c in (ze, ze5))
    nd, nr = cal["endpoints"]["default"]["accepted_score"], cal["endpoints"]["relaxed"]["accepted_score"]
    lnl_s, lnl_e = F_s["lnl_documented"], F_e["lnl_documented"]
    evaluable = ok_s and ok_e and ok_z
    if not evaluable:
        rows += [check(f"{name} (not evaluable: a compared fit is unconverged)", False, not_evaluable=True)
                 for name in ("A6", "A7", "A8")]
    else:
        a6 = abs(lnl_e - ze["lnl"]) <= L
        a7 = abs(lnl_s - lnl_e) <= L
        a8 = nd <= lnl_s + L and nr <= lnl_s + L
        rows += [check("A6: |lnL(F_e) - lnL(Z_e)| <= L", a6, difference=lnl_e - ze["lnl"], L=L),
                 check("A7: |lnL(F_s) - lnL(F_e)| <= L", a7, difference=lnl_s - lnl_e, L=L),
                 check("A8: NQC's default and relaxed endpoints not above lnL(F_s) + L", a8, n_d=nd, n_r=nr,
                       f_s=lnl_s, L=L)]
        if not a6:
            rows.extend(cross_starts("A6", (r_e.f, lnl_e, "flux"), (flux.flux_of(q_of(ze["z"]), pi), ze["lnl"], "z"),
                                     obj, pi, ref, dom, lf, pool, L))
            save_partial()
        if not a7:
            rows.extend(cross_starts("A7", (r_s.f, lnl_s, "flux"), (r_e.f, lnl_e, "flux"), obj, pi, ref, dom, lf,
                                     pool, L))
            save_partial()
        if not a8:
            rows.append(check("A8 failed: where the climbs part (report)", None, f_e=lnl_e, f_s=lnl_s, n_d=nd, n_r=nr))
    rows.append(check("NQC's shortfall below the reference's peak, with its ending rule (report)", None,
                      default=lnl_s - nd, relaxed=lnl_s - nr, ended_by={k: v["ended_by"] for k, v in
                                                                        cal["endpoints"].items()}))
    rows.append(differences("F_s against F_e", r_s.f, r_e.f, pi, ref))
    rows.append(differences("F_e against Z_e", r_e.f, flux.flux_of(q_of(ze["z"]), pi), pi, ref))
    if case in ("toy", "toy300"):
        rng = np.random.default_rng(SEED + CASES.index(case))
        for name in ("F_m1", "F_m2"):
            z0 = np.clip(rng.normal(0.0, 2.0, len(z_r)), LO + 1e-9, HI - 1e-9)
            r, rec = run_fit(name, q_of(z0))
            rows.append(check(f"{name}: seeded random start (report)", None, status=r.status,
                              lnl=rec["lnl_documented"], below_f_s=lnl_s - rec["lnl_documented"]))
    open_dom = flux.open_domain(ref)
    r_o, F_o = run_fit("F_o", flux.generator(r_e.f, pi), open_dom)
    rows.append(check("F_o: the open domain from F_e's endpoint, the box's cost (report)", None, status=r_o.status,
                      gain=F_o["lnl_documented"] - lnl_e, smallest_flux=F_o["smallest_flux"]))
    for name, z in (("N_d", cal["endpoints"]["default"]["z"]), ("N_r", z_r), ("Z_e", ze["z"]),
                    ("F_s", flux.z_of_flux(r_s.f, ref)), ("F_e", flux.z_of_flux(r_e.f, ref))):
        rows.append(export_label(name, obj, np.asarray(z), pool))
    pool.close()
    for rec in fits.values():
        rec.pop("f")
    finish(out, prov, f"S5 comparison: {case}", rows, {"case": case, "L": L, "fits": fits})


def cross_starts(row, H, W, obj, pi, ref, dom, lf, pool, L):
    """A failed same-peak row: each fitter from the other's endpoint, moved 1e-9 inside. H and W are
    (flux vector, lnL, fitter) with fitter "flux" (the reference) or "z" (the analytic z climb).
    Returns each cross-start's A5 rows and the diagnosis."""
    if H[1] < W[1]:
        H, W = W, H
    rows = []

    def refit(label, fitter, start_flux):
        Q = flux.generator(start_flux, pi)
        if fitter == "flux":
            r = flux.fit(lf, flux.inward(Q, pi, dom)[0], dom, objective_map=pool.fluxes)
            a5, ok = a5_rows(label, flux.z_of_flux(r.f, ref), r.pg, pool, r.m_active * r.mu)
            rows.extend(a5)
            return obj.log_likelihood(flux.generator(r.f, pi)), ok and r.status == "converged"
        z0 = np.clip(chart.inverse(Q, ref), LO + 1e-9, HI - 1e-9)
        c = optimize.lbfgsb_fit(obj, z0, LO, HI, pgtol=PG_ANALYTIC, max_iter=50_000, central_step=CENTRAL_STEP,
                                evaluate=pool.values, gradient=obj.gradient, maxcor=20, max_seconds=CAP_SECONDS)
        a5 = [check(f"A5 {label}: analytic projected z-gradient <= {PG_ANALYTIC:g}", c.analytic_pg <= PG_ANALYTIC,
                    value=c.analytic_pg),
              check(f"A5 {label}: central-difference projected z-gradient <= {PG_CENTRAL:g}", c.pg_max <= PG_CENTRAL,
                    value=c.pg_max)]
        rows.extend(a5)
        return -c.f, all(x["passed"] for x in a5)

    hw, hw_ok = refit(f"{row} cross-start, H's fitter from W", H[2], W[0])
    wh, wh_ok = refit(f"{row} cross-start, W's fitter from H", W[2], H[0])
    profile = segment_profile(obj, pi, H[0], W[0])
    label = label_pair(H[1], W[1], (hw, hw_ok), (wh, wh_ok), min(profile), L)
    return rows + [check(f"{row} failed: cross-start diagnosis (report)", None, label=label, h=H[1], w=W[1],
                         h_fitter_from_w=hw, w_fitter_from_h=wh, converged=(hw_ok, wh_ok), profile=profile)]


def differences(name, fa, fb, pi, ref):
    za, zb = flux.z_of_flux(fa, ref), flux.z_of_flux(fb, ref)
    Qa, Qb = flux.generator(fa, pi), flux.generator(fb, pi)
    off = ~np.eye(len(pi), dtype=bool) & (np.minimum(Qa, Qb) > RATE_FLOOR)
    return check(f"{name}: z and Q differences (report)", None, max_abs_dz=float(np.max(np.abs(za - zb))),
                 max_relative_q_above_floor=float(np.max(np.abs(Qa - Qb)[off] / np.abs(Qb)[off])) if off.any() else None,
                 at_lower=[int(np.sum(za <= LO + 1e-6)), int(np.sum(zb <= LO + 1e-6))],
                 at_upper=[int(np.sum(za >= HI - 1e-6)), int(np.sum(zb >= HI - 1e-6))])


def export_label(name, obj, z, pool):
    """Decision 009's label at an endpoint, from the export's central differences. A reference endpoint
    lies strictly inside the box, so a coordinate at a binding face counts its outward gradient and the
    label can read "failed" there; it is reported, not asserted."""
    x = np.clip(z, LO, HI)
    h = chart.FIRST_ORDER_ETA * np.maximum(1.0, np.abs(x))
    up, down = x + h, x - h
    values = pool.values(list(optimize.probes(x, h)) + list(optimize.probes(x, -h)))
    n = len(x)
    g = (np.array(values[:n]) - np.array(values[n:])) / ((up - x) + (x - down))
    res = chart.first_order_test(x, g, LO, HI, float(obj(x)))
    return check(f"{name}: decision 009's export label (report)", None, label="passed" if res["passed"] else "failed",
                 statistic=res["statistic"], at_lower=res["at_lower"], at_upper=res["at_upper"],
                 within_1e6_of_a_bound=int(np.sum((x - LO < 1e-6) | (HI - x < 1e-6))))


# ---------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="step", required=True)
    for name in ("validate", "calibrate", "compare"):
        s = sub.add_parser(name)
        s.add_argument("--out", required=True)
        s.add_argument("--git", default="git")
        s.add_argument("--workers", type=int, default=10)
        s.add_argument("--allow-dirty", action="store_true")
        if name == "calibrate":
            s.add_argument("--case", required=True, choices=CASES)
            s.add_argument("--binary", required=True)
            s.add_argument("--alisim", required=True)
            s.add_argument("--calibration")
            s.add_argument("--timeout", type=int, default=6 * 3600)
        if name == "compare":
            s.add_argument("--calibrated", required=True)
    a = ap.parse_args()
    try:
        {"validate": cmd_validate, "calibrate": cmd_calibrate, "compare": cmd_compare}[a.step](a)
    except SystemExit:
        raise
    except Exception:
        traceback.print_exc()
        sys.exit(2)


if __name__ == "__main__":
    main()
