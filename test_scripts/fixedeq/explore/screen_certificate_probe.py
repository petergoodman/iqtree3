"""Exploratory: why most of SIM-1's 3-taxon screen fits were "not evaluable" (converged under the declared
stop, flux.certificate failed). Not a rerun of the screen's record: those fits stay counted as recorded.

  python test_scripts/fixedeq/explore/screen_certificate_probe.py --screen RESULTS_JSON --git git.exe \
      [--workers 6] [--out DIR] [--allow-dirty]

For each screen fit recorded as not evaluable: the same objective, start, domain and caps refitted (flux.fit
is deterministic; the endpoint's lnL is compared with the record's), then at the endpoint the rows'
normalized slacks counted at 1e-12, 1e-9, 1e-8 and 1e-6, and flux.certificate's residual at each of those
activity thresholds against decision 021's line 1e-9 max(1, |lnL|), with the analytic projected z-gradient.
The question: are these endpoints first-order points with rows on the box's bounds that the barrier leaves
at slacks above the certificate's 1e-12 (an endpoint the active-set finish never reached), or short of a
peak? The endpoint fluxes are written out for reuse.
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
from oracle import flux, nqc, readers, simulate  # noqa: E402
from regress import REPO, read, write  # noqa: E402
from reference import Pool, fit_record  # noqa: E402
import common  # noqa: E402
import exact as ex  # noqa: E402

SCRIPT = Path(__file__).resolve().relative_to(REPO).as_posix()
THRESHOLDS = (1e-12, 1e-9, 1e-8, 1e-6)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--screen", required=True, help="the screen run's results.json")
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--out")
    ap.add_argument("--git", default="git")
    ap.add_argument("--allow-dirty", action="store_true")
    a = ap.parse_args()
    screen = json.loads(read(Path(a.screen).expanduser()))
    prov = common.provenance(a, SCRIPT, "screen-certificate-probe", "exploratory", screen=a.screen,
                             screen_commit=screen["provenance"]["git_commit"], thresholds=THRESHOLDS)
    out = common.new_out(a.out, prov)
    pi = common.target_values()
    _, truths = simulate.truths()
    ref = nqc.lg_seed(pi)[1]
    dom = flux.matched(ref)
    start = flux.flux_of(truths["LG"], pi)
    tree = readers.parse_newick(ex.THREE)
    names = tuple(sorted(leaf.name for leaf in tree.leaves()))
    loci = ex.read_loci(screen["provenance"]["loci"]["path"]) if screen["provenance"].get("loci") else None
    found = []
    for rec in screen["fits"]:
        if rec.get("evaluable") or "fit" not in rec:
            continue
        label = f"{rec['scenario']} {rec['base']}"
        t0 = time.monotonic()
        w, gaps, reference, _ = ex.screen_data(rec["scenario"], truths[rec["base"]], pi, ex.THREE, loci)
        obj = nqc.Objective.from_patterns(tree, simulate.columns(names, gaps=gaps), w, pi, ref)
        pool = Pool(obj, a.workers)
        try:
            r = flux.fit(flux.LikelihoodFlux(obj, pi), flux.inward_flux(start, dom, pi)[0], dom,
                         objective_map=pool.fluxes, max_iter=ex.MAX_STEPS, max_seconds=ex.CAP_3)
        finally:
            pool.close()
            if pool.pool is not None:
                pool.pool.join()
        value, gf = flux.LikelihoodFlux(obj, pi)(r.f)
        slack = dom.slack(r.f)
        line = ex.CERTIFICATE * max(1.0, abs(value))
        row = {"label": label, "record_lnl": rec.get("lnl"), "refit_lnl": float(value),
               "same_endpoint": rec.get("lnl") == float(value), "fit": fit_record(r),
               "projected_z_gradient": r.pg, "line": line, "C_fit": simulate.circulation(flux.generator(r.f, pi), pi),
               "rows_at_slack": {f"{t:g}": int(np.sum(slack <= t)) for t in THRESHOLDS},
               "smallest_slacks": np.sort(slack)[:6].tolist(),
               "certificate": {f"{t:g}": flux.certificate(r.f, gf, dom, threshold=t)["residual"] for t in THRESHOLDS},
               "seconds": round(time.monotonic() - t0, 1)}
        row["passes_at"] = [t for t, v in row["certificate"].items() if v <= line]
        found.append(row)
        np.save(out / f"endpoint-{rec['scenario']}-{rec['base']}.npy", r.f)
        write(out / "results.json", json.dumps({"provenance": prov, "fits": found}, indent=1, default=common._json))
        print(f"{label}: same endpoint {row['same_endpoint']}, rows at slack {row['rows_at_slack']}, "
              f"certificate {row['certificate']} (line {line:.3g}), passes at {row['passes_at']}", flush=True)
    write(out / "results.json", json.dumps({"provenance": prov, "fits": found}, indent=1, default=common._json))
    print(f"wrote {out / 'results.json'}")


if __name__ == "__main__":
    main()
