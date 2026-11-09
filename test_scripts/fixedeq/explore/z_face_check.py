"""Exploratory: is a z climb's endpoint stuck on the lower face of NQC's box, where moving one
coordinate inward would still raise lnL? Asked after S5's gate labelled aa_example's A6 failure "W
stopped short" (CHANGELOG.md, S5 close). Not a test; nothing here is a pass line.

  python test_scripts/fixedeq/explore/z_face_check.py --calibrated CALIBRATION_DIR [--climb Z_e] [--workers 10]

Reads a reference.py calibration; at the climb's endpoint, prints the coordinates at each bound, the
largest projected-gradient component and how many lower-bound coordinates have a gradient pulling
inward; then, for every lower-bound coordinate alone, lnL on a grid from the lower bound to the upper
(step 0.5), and the largest gain over the endpoint.
"""
import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse  # noqa: E402
import json  # noqa: E402
import sys  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "reference"))
import reference as R  # noqa: E402
from oracle import chart, optimize  # noqa: E402

STEP = 0.5


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--calibrated", required=True)
    ap.add_argument("--climb", default="Z_e", choices=("Z_e", "Z_e5"))
    ap.add_argument("--workers", type=int, default=10)
    a = ap.parse_args()
    cal = Path(a.calibrated).expanduser().resolve()
    res = json.loads((cal / "results.json").read_text())
    p = R.Problem(res["case"], cal / "data")
    obj = p.objective()
    z = np.array(res["climbs"][a.climb]["z"])
    f0, g = obj.gradient(z)
    print(f"{res['case']} {a.climb}: lnL {-f0!r} (recorded {res['climbs'][a.climb]['lnl']!r})")
    pg = optimize.projected_gradient(z, g, R.LO, R.HI)
    lo, hi = np.flatnonzero(z <= R.LO), np.flatnonzero(z >= R.HI)
    inward = g[lo] < 0
    print(f"at lower {len(lo)}, at upper {len(hi)}, max |projected| {np.max(np.abs(pg)):.3g}; at lower pulling "
          f"inward: {int(inward.sum())}, largest |g| there {np.max(np.abs(g[lo][inward]), initial=0.0):.3g}")
    grid = np.arange(R.LO + STEP, R.HI + 1e-9, STEP)
    pts = []
    for k in lo:
        for t in grid:
            v = z.copy()
            v[k] = t
            pts.append(v)
    pool = R.Pool(obj, a.workers)
    gains = f0 - np.array(pool.values(pts)).reshape(len(lo), len(grid))   # lnL gain over the endpoint
    pool.close()
    best = gains.max(axis=1)
    order = np.argsort(-best)
    pairs = [(i, int(j)) for i, free in enumerate(chart.free_destinations(p.ref)) for j in free]
    print(f"single-coordinate moves from the lower face: {int(np.sum(best > 1e-3))} gain more than 1e-3, "
          f"{int(np.sum(best > 0))} gain anything; largest {best[order[0]]:.4g}")
    print(f"k (i, j): best gain, at z, gradient at the endpoint, gain at lower + {STEP:g}, + 2, + 5")
    for o in order[:12]:
        k = lo[o]
        print(f"  {k} {pairs[k]}: {best[o]:.4g} at z {grid[gains[o].argmax()]:.1f}, g {g[k]:.3g}, "
              f"{gains[o, 0]:.3g}, {gains[o, 3]:.3g}, {gains[o, 9]:.3g}")


if __name__ == "__main__":
    main()
