"""Exploratory: where the SIM-1 pilot's per-rate log errors sit. At 100 genes the RMS log error rose (1.11
against 0.72 to 0.97 at 30 genes) while fewer rates were off twofold. Per fit: the log error e = log(q_fit /
q_true) over the 380 off-diagonal rates against each rate's true flux pi*_i q_ij, as the RMS within flux
quintiles, a flux-weighted RMS, counts off twofold and tenfold, the largest errors with their flux rank, and
the coordinates within 1 of the domain's bounds (the chart's z, from the export record).

  python test_scripts/fixedeq/explore/pilot_rate_errors.py --git git.exe [--out DIR] PILOT_DIR
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "sim"))
import pipeline as pl  # noqa: E402  (also puts the drivers on sys.path)
import common  # noqa: E402
from oracle import builtin, readers  # noqa: E402
from nqcchecks import check  # noqa: E402
from regress import REPO, read  # noqa: E402

SCRIPT = Path(__file__).resolve().relative_to(REPO).as_posix()


def errors(fitdir, datadir):
    pi = np.array([float(v) for v in json.loads(read(datadir / "data.json"))["target"].split(",")])
    Qt = readers.read_qmat(read(datadir / "truth.qmat"))[0]
    Qh = readers.read_qmat(read(fitdir / f"{pl.PREFIX}.NQC.qmat"))[0]
    info = readers.read_info(read(fitdir / f"{pl.PREFIX}.NQC.info"))
    off = ~np.eye(20, dtype=bool)
    e, flux = np.log(Qh[off] / Qt[off]), (pi[:, None] * Qt)[off]
    rank = np.argsort(np.argsort(-flux))                 # 0 is the largest true flux
    quint = np.array_split(np.argsort(-flux), 5)
    pairs = [f"{builtin.STATES[i]}{builtin.STATES[j]}" for i in range(20) for j in range(20) if i != j]
    z = np.array([float(v) for v in info["coordinates"].split(",")])
    lo, hi = (float(v) for v in info["domain"].split(","))
    top = np.argsort(-np.abs(e))[:8]
    return {"rms": float(np.sqrt(np.mean(e ** 2))),
            "rms_by_flux_quintile": [round(float(np.sqrt(np.mean(e[q] ** 2))), 3) for q in quint],
            "flux_weighted_rms": float(np.sqrt(np.sum(flux * e ** 2) / np.sum(flux))),
            "twofold": int(np.sum(np.abs(e) > np.log(2))), "tenfold": int(np.sum(np.abs(e) > np.log(10))),
            "largest": [{"pair": pairs[k], "e": round(float(e[k]), 3), "flux_rank": int(rank[k]),
                         "q_true": float(Qt[off][k]), "q_fit": float(Qh[off][k])} for k in top],
            "z_within_1_of_lower": int(np.sum(z < lo + 1)), "z_within_1_of_upper": int(np.sum(z > hi - 1)),
            "z_min": float(z.min())}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("pilot")
    ap.add_argument("--out")
    ap.add_argument("--git", default="git")
    ap.add_argument("--allow-dirty", action="store_true")
    a = ap.parse_args()
    prov = common.provenance(a, SCRIPT, "pilot-rate-errors", "exploratory", pilot=a.pilot)
    out = common.new_out(a.out, prov)
    d = Path(a.pilot).expanduser().resolve()
    rows = []
    for name, rec in json.loads(read(d / "results.json"))["records"].items():
        fit = d / name / "fit"
        if not (fit / f"{pl.PREFIX}.NQC.qmat").exists():
            continue
        r = errors(fit, d / name / "data")
        rows.append(check(name, None, columns=rec.get("columns"), taxa=rec.get("taxa"), **r))
        print(f"{name}: rms {r['rms']:.3f}, by flux quintile {r['rms_by_flux_quintile']}, flux-weighted "
              f"{r['flux_weighted_rms']:.3f}, twofold {r['twofold']}, tenfold {r['tenfold']}, z near lower "
              f"{r['z_within_1_of_lower']} (min {r['z_min']:.2f}); largest {[(x['pair'], x['e'], x['flux_rank']) for x in r['largest'][:4]]}",
              flush=True)
    common.finish(out, prov, "SIM-1 pilot: per-rate log errors against true flux", rows)


if __name__ == "__main__":
    main()
