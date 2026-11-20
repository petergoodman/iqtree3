"""Engineering measurement (CODE_PLAN.md 3.8, "Laptop work before submission", step 1): SIM-2's stage B timed
on the SIM-1 pilot's capped NQ.bird endpoints, which sets B's caps and the scheduler's estimates, and a flux.fit
probe that decides the certified subset (Peter's rule, 2026-11-19: each of its 6 fits within 8 h on Puma and all
6 within 500 charged core-hours).

  python test_scripts/fixedeq/explore/stage_b_calibration.py --pilot DIR --git git.exe [--allow-dirty] [--out DIR]
      [--steps speedup,climb30,probe,climb100]

speedup: three gradients at W = 1, 2, 4 and 8 gene-parallel workers on the 100-gene data, each W's values checked
equal to W = 1's. climb30: sim2.stage_b at W = 1, capped at 1.54 h (2 h on Puma at 1.3 times the laptop). probe:
flux.fit from climb30's endpoint on 8 workers, capped at 1 h. climb100: stage_b at W = 4, capped at 4 h; a climb
not stalled by then is a stop for Peter. Results are saved after each step.
"""
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "sim"))
sys.path.insert(0, str(HERE.parent / "reference"))
import sim2  # noqa: E402  (BLAS threads at 1; puts the drivers on sys.path)
import common  # noqa: E402
import pipeline as pl  # noqa: E402
from oracle import chart, flux, nqc  # noqa: E402
from nqcchecks import check  # noqa: E402
from reference import Pool, fit_record  # noqa: E402
from regress import REPO, read, sha256  # noqa: E402

SCRIPT = Path(__file__).resolve().relative_to(REPO).as_posix()
STEPS = ("speedup", "climb30", "probe", "climb100")
CELL30, CELL100 = "main-NQ.bird-r1", "genes100-NQ.bird"
PUMA = 1.3                                  # Puma's assumed slowdown per core against the laptop (CODE_PLAN 3.8)
CAP30, CAP100, CAP_PROBE = 2 * 3600 / PUMA, 4 * 3600.0, 3600.0
PROBE_WORKERS = 8
SUBSET_FITS, SUBSET_HOURS, SUBSET_CORE_HOURS = 6, 8.0, 500.0


def target(data):
    return np.array([float(v) for v in json.loads(read(Path(data) / "data.json"))["target"].split(",")])


def speedup(pilot, tree):
    data, fit = pilot / CELL100 / "data", pilot / CELL100 / "fit"
    pi, z = target(data), sim2.nqc_state(fit)["z"]
    obj = sim2.objective_for(data, tree, pi)[0]
    out, first = {}, None
    for w in (1, 2, 4, 8):
        times = []
        with nqc.GeneParallel(obj.parts, w) as par:
            for _ in range(3):
                t0 = time.monotonic()
                f, g = par.gradient(z)
                times.append(round(time.monotonic() - t0, 2))
        first = first or (f, g)
        out[f"W{w}"] = {"seconds": times, "equal_to_w1": bool(f == first[0] and np.array_equal(g, first[1]))}
        print(f"speedup W={w}: {times}", flush=True)
    return out


def probe(pilot, tree, z_b):
    data, fit = pilot / CELL30 / "data", pilot / CELL30 / "fit"
    pi, state = target(data), sim2.nqc_state(fit)
    obj = sim2.objective_for(data, tree, pi)[0]
    dom = flux.matched(obj.ref, *state["domain"])
    f0, moved = flux.inward(chart.build(np.asarray(z_b), pi, obj.ref)[0], pi, dom)
    pool = Pool(obj, PROBE_WORKERS)
    try:
        r = flux.fit(flux.LikelihoodFlux(obj, pi), f0, dom, objective_map=pool.fluxes, max_seconds=CAP_PROBE)
    finally:
        pool.close()
        if pool.pool is not None:
            pool.pool.join()
    rec = {"fit": fit_record(r), "start_moved": moved, "workers": PROBE_WORKERS}
    if r.status == "converged":
        hours = r.seconds * PUMA / 3600
        charged = PROBE_WORKERS * hours
        rec.update({"puma_hours_per_fit": hours, "charged_core_hours_6": SUBSET_FITS * charged,
                    "include_subset": hours <= SUBSET_HOURS and SUBSET_FITS * charged <= SUBSET_CORE_HOURS})
    else:
        rec.update({"include_subset": False, "reason": f"the probe ended {r.status} within its {CAP_PROBE:.0f} s cap"})
    return rec


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pilot", required=True, help="the SIM-1 pilot's run directory")
    ap.add_argument("--scaffold", default=str(HERE.parent / "sim" / "scaffold"))
    ap.add_argument("--steps", default=",".join(STEPS))
    ap.add_argument("--out", help="new directory outside the repository (default under ~/iqtree3-runs/sim2/)")
    ap.add_argument("--git", default="git")
    ap.add_argument("--allow-dirty", action="store_true")
    a = ap.parse_args()
    steps = a.steps.split(",")
    if any(s not in STEPS for s in steps):
        sys.exit(f"--steps takes a comma-separated subset of {', '.join(STEPS)}")
    pilot = Path(a.pilot).expanduser().resolve()
    prov = common.provenance(a, SCRIPT, "calibration", "engineering measurement", pilot=str(pilot), steps=steps,
                             caps={"climb30": CAP30, "climb100": CAP100, "probe": CAP_PROBE},
                             scaffold=pl.scaffold_record(a.scaffold))
    out = common.new_out(a.out, prov, runs="~/iqtree3-runs/sim2")
    tree = pl.read_scaffold(a.scaffold)[0]
    pilot_scaffold = json.loads(read(pilot / "results.json"))["provenance"]["scaffold"]["sha256"]
    rows = [check("the scaffold is the pilot's", pilot_scaffold == prov["scaffold"]["sha256"],
                  pilot=pilot_scaffold, here=prov["scaffold"]["sha256"])]
    records, title = {}, "SIM-2 stage B calibration"

    def save(name, rec, new_rows=()):
        records[name] = rec
        rows.extend(new_rows)
        rows.append(check(name, None, **{k: v for k, v in rec.items() if k != "trace"}))
        common.save(out, prov, title + ", in progress", rows, {"records": records})

    if "speedup" in steps:
        save("speedup", speedup(pilot, tree))
    if "climb30" in steps or "probe" in steps:
        (out / "climb30").mkdir()
        rec, found = sim2.stage_b(pilot / CELL30 / "data", pilot / CELL30 / "fit", out / "climb30", 1, CAP30, tree)
        save("climb30", rec, found)
        print(f"climb30: {rec['stop']}, {rec['iterations']} iterations, {rec['seconds']} s, gain {rec['gain']:.3f}",
              flush=True)
        if "probe" in steps:
            save("probe", probe(pilot, tree, rec["z_b"]))
    if "climb100" in steps:
        (out / "climb100").mkdir()
        rec, found = sim2.stage_b(pilot / CELL100 / "data", pilot / CELL100 / "fit", out / "climb100", 4, CAP100,
                                  tree)
        found.append(check("the 100-gene climb stalled within 4 h (else a stop for Peter)", rec["stop"] == "stalled",
                           stop=rec["stop"], seconds=rec["seconds"]))
        save("climb100", rec, found)
    common.finish(out, prov, title, rows, {"records": records,
                                           "inputs": {c: sha256(pilot / c / "fit" / "probe.ckp.gz")
                                                      for c in (CELL30, CELL100)}})


if __name__ == "__main__":
    main()
