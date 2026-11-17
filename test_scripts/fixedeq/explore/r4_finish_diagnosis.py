"""Exploratory: why the flux reference stopped short on SIM-1's R4 (NQ.bird with +G4{0.5}, three taxa,
exact weights N = 1e6; CHANGELOG 2026-11-16). Not a rerun of the gate: R4's confirmatory result stays
failed and counted.

  python test_scripts/fixedeq/explore/r4_finish_diagnosis.py --confirmatory RESULTS_JSON --git git.exe \
      [--workers 10] [--steps 8] [--out DIR] [--allow-dirty]

1. Refit R4 exactly as sim/exact.py recover did (same objective, start, domain, caps, pool), with every
   Newton step logged, and compare the endpoint with the confirmatory record (deterministic: same lnL
   and projected z-gradient expected).
2. At the endpoint, the active-set finish's own quantities with no face held: the central-difference
   Hessian in the re-centred basis, the Newton step, its decrement against the finish's rounding
   (1e-12 max(1, |lnL|)) and its stationarity threshold (1e-6 times that rounding); the gradient's
   share along the Hessian's eigenvectors.
3. Plain Newton steps from the endpoint (re-centred each step, a step cut only to stay inside the box),
   recording after each the lnL, the projected z-gradient against the declared 1e-6 and flux.certificate
   against 1e-9 max(1, |lnL|).

Run once, 2026-11-17 UTC, on the solver before decision 042's finish revision (outputs in
~/iqtree3-runs/sim1/explore-r4-finish-20261117T025318Z/): the refit repeated the confirmatory endpoint
bitwise; there the finish's stationary exit fired (decrement 2e-12 against 5.5e-12) with the projected
z-gradient at 2.4e-6, 97% of the gradient along the stiffer half of the Hessian; one Newton step gave
4e-11 and a certificate of 5e-9 (line 5.5e-3). Rerun on the revised solver, step 1's refit no longer
stops short, so the script documents the old behaviour.
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

HERE = Path(__file__).resolve().parent
for _p in ("", "regression", "reference", "sim"):
    sys.path.insert(0, str(HERE.parent / _p))
from oracle import flux, nqc  # noqa: E402
from regress import REPO, read, write  # noqa: E402
from reference import Pool, fit_record  # noqa: E402
import common  # noqa: E402
import exact as ex  # noqa: E402

SCRIPT = Path(__file__).resolve().relative_to(REPO).as_posix()


def state(obj, pi, dom, f):
    value, gf = flux.LikelihoodFlux(obj, pi)(f)
    cert = flux.certificate(f, gf, dom)
    return {"lnl": float(value), "projected_z_gradient": flux.projected_z_gradient(f, gf, dom),
            "certificate": cert["residual"], "certificate_line": ex.CERTIFICATE * max(1.0, abs(value)),
            "active": cert["active"], "smallest_flux": float(f.min()),
            "smallest_slack": float(np.min(dom.slack(f)))}, gf


def newton(pool, obj, pi, dom, f):
    """The finish's step with no face held: Hessian by central differences in the basis centred at f."""
    basis = flux.Basis(f, dom)
    u0 = np.zeros(basis.N.shape[1])
    value, gf = flux.LikelihoodFlux(obj, pi)(f)
    g = -basis.chain(gf)
    H, count = flux.difference_hessian(pool.fluxes, basis, u0, g, central=True)
    w, V = np.linalg.eigh(H)
    d = -V @ ((V.T @ g) / w) if w.min() > 0 else None
    rounding = flux.ROUNDING * max(1.0, abs(value))
    out = {"hessian_gradients": count, "eig_min": float(w.min()), "eig_max": float(w.max()),
           "gradient_norm_u": float(np.max(np.abs(g))), "rounding": rounding,
           "stationary_threshold": 1e-6 * rounding}
    if d is not None:
        dec = float(-(g @ d) / 2)
        share = (V.T @ g) ** 2 / w / 2       # each eigendirection's part of the decrement
        stiff = w > np.median(w)
        out.update({"decrement": dec, "pure": dec <= rounding, "stationary_branch": dec <= 1e-6 * rounding,
                    "decrement_share_stiff_half": float(share[stiff].sum() / max(share.sum(), 1e-300)),
                    "gradient_share_stiff_half": float(np.sum((V.T @ g)[stiff] ** 2) / max(np.sum(g ** 2), 1e-300))})
    return basis, d, out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--confirmatory", required=True, help="results.json of the recover run that holds R4")
    ap.add_argument("--workers", type=int, default=10)
    ap.add_argument("--steps", type=int, default=8)
    ap.add_argument("--out")
    ap.add_argument("--git", default="git")
    ap.add_argument("--allow-dirty", action="store_true")
    a = ap.parse_args()
    prov = common.provenance(a, SCRIPT, "r4-finish-diagnosis", "exploratory", workers=a.workers)
    out = common.new_out(a.out, prov)
    confirmed = json.loads(read(Path(a.confirmatory).expanduser()))["records"]["R4"]["fit"]
    pi = common.target_values()
    ref = nqc.lg_seed(pi)[1]
    dom = flux.matched(ref)
    p = ex.problems()["R4"]
    obj = ex.objective(p, pi, ref)
    pool = Pool(obj, a.workers)
    steps = []
    t0 = time.monotonic()
    r = flux.fit(flux.LikelihoodFlux(obj, pi), flux.inward(p.start, pi, dom)[0], dom, objective_map=pool.fluxes,
                 max_iter=ex.MAX_STEPS, max_seconds=p.cap, progress=steps.append)
    refit = fit_record(r)
    result = {"refit": refit, "refit_seconds": round(time.monotonic() - t0, 1), "steps": steps,
              "same_as_confirmatory": {"lnl": refit["lnl"] == confirmed["lnl"],
                                       "projected_z_gradient": refit["projected_z_gradient"] ==
                                       confirmed["projected_z_gradient"], "message": refit["message"] ==
                                       confirmed["message"]}}
    print(f"refit: {r.status}, lnL {r.value:.9f}, pg {r.pg:.3g}; same as confirmatory "
          f"{result['same_as_confirmatory']}", flush=True)
    f = r.f
    newton_steps = []
    for k in range(a.steps + 1):
        st, _ = state(obj, pi, dom, f)
        basis, d, quant = newton(pool, obj, pi, dom, f)
        st.update(quant)
        newton_steps.append(st)
        print(f"step {k}: lnL {st['lnl']:.9f} pg {st['projected_z_gradient']:.3g} cert {st['certificate']:.3g} "
              f"(line {st['certificate_line']:.3g}) decrement {st.get('decrement', float('nan')):.3g} "
              f"stationary branch {st.get('stationary_branch')}", flush=True)
        if d is None or k == a.steps:
            break
        s, Gd = basis.slack(np.zeros(len(d))), basis.G @ d
        shrink = Gd < 0
        alpha = min(1.0, 0.995 * float(np.min(-s[shrink] / Gd[shrink]))) if shrink.any() else 1.0
        f = basis.flux(alpha * d)
        newton_steps[-1]["alpha"] = alpha
    pool.close()
    result["newton_from_endpoint"] = newton_steps
    write(out / "results.json", json.dumps({"provenance": prov, **result}, indent=1, default=common._json))
    print(f"wrote {out / 'results.json'}")


if __name__ == "__main__":
    main()
