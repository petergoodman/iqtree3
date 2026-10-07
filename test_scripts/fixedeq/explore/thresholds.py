"""Exploratory, read-only: the measurements behind decisions 026 and 027 (the compiled pass lines,
the coordinate box and the start policy). Not a test; nothing here is a pass line.

  python test_scripts/fixedeq/explore/thresholds.py --runs ~/iqtree3-runs --git git.exe > OUT.txt

Reads the probe (g) record and baseline runs 2 and 8 under --runs, the G0 manifest's target, and
the oracle. Writes to standard output only, with a provenance header.

[A] where starts and fitted matrices sit in the default box, in jump-chain log-ratio coordinates
    against the references of LG rebuilt at the target (decision 011), after checking that the
    probe g0 checkpoint parses to the exchangeabilities its report prints;
[B] how far two correct floating-point routes to the same Q differ, by case and conditioning;
[C] the sum of the target when written at fewer decimals (decision 006's 1e-6 line).
"""
import argparse
import datetime
import json
import math
import subprocess
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / "test_scripts/fixedeq"))
sys.path.insert(0, str(REPO / "test_scripts/fixedeq/probes"))
from oracle import builtin, cases, chart  # noqa: E402
from probes import ckp_values, read_checkpoint  # noqa: E402

SEED = 20261005
LO, HI = math.log(1e-5), math.log(10.0)
N = 20


def git(exe, *args):
    try:
        return subprocess.run([exe, *args], cwd=REPO, capture_output=True, text=True, check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError) as e:
        return f"unavailable ({e})"


def full_from_offdiag(rates):
    """IQ-TREE's non-reversible rates[]: row-major, diagonal skipped (model/modelmarkov.cpp:122-131)."""
    Q = np.zeros((N, N))
    Q[~np.eye(N, dtype=bool)] = rates
    np.fill_diagonal(Q, -Q.sum(axis=1))
    return Q


def sym_from_upper(rates):
    """IQ-TREE's reversible rates[]: upper triangle, row-major (model/modelmarkov.cpp:114-118)."""
    R = np.zeros((N, N))
    R[np.triu_indices(N, 1)] = rates
    return R + R.T


def paml_lower(text):
    """The 19 lower-triangle rows of the PAML block in an .iqtree report."""
    lines = text[text.index("PAML format"):].splitlines()[2:N + 1]
    R = np.zeros((N, N))
    for i, line in enumerate(lines, start=1):
        vals = [float(x) for x in line.split()]
        assert len(vals) == i, (i, len(vals))
        R[i, :i] = vals
    return R + R.T


def placement(label, z, note=""):
    z = np.asarray(z)
    lo_d, hi_d = z.min() - LO, HI - z.max()
    near = [int(np.sum((z - LO < m) | (HI - z < m))) for m in (0.5, 1.0, 2.0)]
    outside = int(np.sum((z < LO) | (z > HI)))
    print(f"  {label}: z in [{z.min():.3f}, {z.max():.3f}]; distance to lower {lo_d:.3f}, to upper {hi_d:.3f}; "
          f"outside {outside}/360; within 0.5/1/2 of an edge {near[0]}/{near[1]}/{near[2]} {note}")


def build_alt(z, pi, ref):
    """The same mathematics by a different floating-point path: softmax summed with fsum in reverse
    order, nu by LU on the square system with one balance equation replaced by sum(nu) = 1, and the
    diagonal by fsum. Returns Q, nu and the condition number of the augmented nu system."""
    n = len(ref)
    z = np.asarray(z, dtype=float).reshape(n, n - 2)
    K = np.zeros((n, n))
    for i, free in enumerate(chart.free_destinations(ref)):
        logit = np.full(n, -np.inf)
        logit[free] = z[i]
        logit[ref[i]] = 0.0
        e = np.exp(logit - logit.max())
        K[i] = e / math.fsum(e[::-1])
    A = (K - np.eye(n)).T.copy()
    A[-1, :] = 1.0
    b = np.zeros(n)
    b[-1] = 1.0
    nu = np.linalg.solve(A, b)
    Q = (nu / pi)[:, None] * K
    np.fill_diagonal(Q, 0.0)
    for i in range(n):
        Q[i, i] = -math.fsum(Q[i])
    return Q, nu, np.linalg.cond(np.vstack([np.ones(n), (K - np.eye(n)).T]))


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--runs", required=True, help="the run directory, ~/iqtree3-runs")
    ap.add_argument("--git", default="git", help="git executable (git.exe under WSL)")
    a = ap.parse_args()
    runs = Path(a.runs).expanduser()

    print(f"script test_scripts/fixedeq/explore/thresholds.py; commit {git(a.git, 'rev-parse', 'HEAD')}; "
          f"clean {git(a.git, 'status', '--porcelain') == ''}; "
          f"{datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds')}; "
          f"seed {SEED}; numpy {np.__version__}; box [{LO:.6f}, {HI:.6f}]")

    man = json.loads((REPO / "test_scripts/fixedeq/manifest/g0_manifest.json").read_text())["manifest"]
    pi = np.array([float(v) for v in man["target"]["values_17_digits"]])
    lg = chart.reversible_seed(builtin.exchangeabilities("LG")[0], pi)
    ref = chart.references(lg)

    print("\n[A] Placement in the default box (references: LG at the decision 023 target)")
    g0 = runs / "probes/record-20261005T174330Z/g0"
    g_rates = np.array(ckp_values(read_checkpoint(g0 / "probe.ckp.gz")["ModelProtein!rates"]))
    R_g = sym_from_upper(g_rates)
    print(f"  probe g0 parse check: max |checkpoint - report| over 190 exchangeabilities "
          f"{np.max(np.abs(R_g - paml_lower((g0 / 'probe.iqtree').read_text()))):.2e} (report: 6 decimals)")
    print(f"  probe g0 GTR20 exchangeabilities at 1e-4: {int(np.sum(np.abs(g_rates - 1e-4) < 1e-9))}/190, "
          f"at 100: {int(np.sum(np.abs(g_rates - 100) < 1e-7))}/190")
    Q_g = chart.reversible_seed(R_g, pi)
    nq = builtin.full_matrix("NQ.PFAM")[0]
    placement("LG at pi* (default start)", chart.inverse(lg, ref))
    placement("NQ.pfam", chart.inverse(nq, ref), "[jump chain only]")
    placement("GTR20+F{pi*} fit on aa_example (probe g0; decision 025 incumbent)", chart.inverse(Q_g, ref))
    for label, path, key in [
            ("NONREV fit, aa_example (baseline run 2)", "baseline/rep1/run02/run02.ckp.gz", "ModelProtein!rates"),
            ("joint NONREV fit, turtle -p (baseline run 8)", "baseline/rep1/run08/run08.ckp.gz",
             "PartitionModelPlen!NONREV+FO!ModelProtein!rates")]:
        r = np.array(ckp_values(read_checkpoint(runs / path)[key]))
        placement(label, chart.inverse(full_from_offdiag(r), ref),
                  "[own pi, jump chain only; checkpoint from before the final model optimization]")
        print(f"    raw rates at IQ-TREE's MIN_RATE 1e-4: {int(np.sum(np.abs(r - 1e-4) < 1e-9))}/380, "
              f"at MAX_RATE 100: {int(np.sum(np.abs(r - 100) < 1e-7))}/380")

    print("\n[B] Two correct routes to the same Q from the same z")
    rng = np.random.default_rng(SEED)
    skewed = np.full(N, 1e-4)
    skewed[0] = 1 - skewed[1:].sum()
    blocks_ref = np.array([1] + [0] * 9 + [11] + [10] * 9)
    blocks = np.array([[0.0 if (j < 10) == (i < 10) else -11.5 for j in free]
                       for i, free in enumerate(chart.free_destinations(blocks_ref))])
    trials = []
    for _ in range(20):
        p = cases.random_target(N, rng)
        Qf = cases.generator_from_flux(cases.cycle_flux(N, rng), p)
        r_ = chart.references(Qf)
        trials.append(("ordinary: random balanced flux, random target", p, chart.inverse(Qf, r_), r_))
    trials.append(("LG at the turtle pi*", pi, chart.inverse(lg, ref), ref))
    trials.append(("GTR20 incumbent at the turtle pi*", pi, chart.inverse(Q_g, ref), ref))
    for _ in range(5):
        trials.append(("hard: skewed target, 19 entries at 1e-4", skewed, rng.normal(0, 1, N * (N - 2)), ref))
        trials.append(("hard: wide coordinates over the whole box", cases.random_target(N, rng),
                       rng.uniform(LO, HI, N * (N - 2)), ref))
    trials.append(("hard: near-reducible, two blocks joined at the lower edge", cases.random_target(N, rng),
                   blocks.ravel(), blocks_ref))
    summary = {}
    for label, p, z, r_ in trials:
        Q1, nu1 = chart.build(z, p, r_)
        Q2, nu2, cond = build_alt(z, p, r_)
        d1, d2 = chart.diagnostics(Q1, p, nu1), chart.diagnostics(Q2, p, nu2)
        s = summary.setdefault(label, {"n": 0, "rel": 0.0, "res": 0.0, "cond": 0.0, "minq": 1.0})
        s["n"] += 1
        s["rel"] = max(s["rel"], float(np.max(np.abs(Q1 - Q2) / np.abs(Q1))))
        s["res"] = max(s["res"], *(d[k] for d in (d1, d2) for k in ("stationarity", "mean_rate_error", "row_sums")))
        s["cond"] = max(s["cond"], cond)
        s["minq"] = min(s["minq"], d1["min_rate"])
    for label, s in summary.items():
        print(f"  {label} (x{s['n']}): Q relative difference {s['rel']:.1e}; worst residual {s['res']:.1e}; "
              f"condition number up to {s['cond']:.1e}; smallest rate {s['minq']:.1e}")

    print("\n[C] Sum of the target when written at fewer decimals")
    for d in (8, 6, 5, 4):
        print(f"  {d} decimals: sum - 1 = {math.fsum(float(f'{x:.{d}f}') for x in pi) - 1:+.2e}")
    print(f"  17 significant digits: sum - 1 = {math.fsum(float(s) for s in man['target']['values_17_digits']) - 1:+.2e}")


if __name__ == "__main__":
    main()
