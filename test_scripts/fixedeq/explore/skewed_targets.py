"""Exploratory, read-only: whether LG rebuilt at a skewed target lies in NQC's default box (PLAN.md
risk 22; decisions 023 and 025). Not a test; nothing here is a pass line.

  python test_scripts/fixedeq/explore/skewed_targets.py --git git.exe > OUT.txt

For each target: LG rebuilt at it (D05), its row-maximum references (decision 011) and its jump-chain
log-ratio coordinates (D01), with the number outside the box [log 1e-5, log 10] (decision 007).
Targets: decision 023's turtle target raised to powers and renormalized, LG's own frequencies, the
oracle's extreme fixture (19 entries at 1e-4), and seeded Dirichlet draws. Writes to standard output
only, with a provenance header.
"""
import argparse
import datetime
import math
import subprocess
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / "test_scripts/fixedeq"))
from oracle import builtin, cases, chart, nqc, readers, target  # noqa: E402

SEED = 20261022
DRAWS = 1000
MIN_TARGET = 1e-4  # NQC refuses smaller entries (D07, decision 006)


def git(exe, *args):
    try:
        return subprocess.run([exe, *args], cwd=REPO, capture_output=True, text=True, check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError) as e:
        return f"unavailable ({e})"


def outside(pi):
    Q, ref = nqc.lg_seed(pi)
    z = chart.inverse(Q, ref)
    return int(np.sum(z < nqc.LOWER)) + int(np.sum(z > nqc.UPPER)), float(z.min()), float(z.max())


def row(label, pi):
    n, zmin, zmax = outside(pi)
    note = "" if pi.min() >= MIN_TARGET else "  [refused by NQC: an entry below 1e-4]"
    print(f"  {label:<26} min {pi.min():.2e}  max {pi.max():.3f}  max/min {pi.max() / pi.min():9.1f}  "
          f"z [{zmin:7.2f}, {zmax:5.2f}]  outside {n:3d}  smallest ratio {math.exp(zmin):.1e}{note}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--git", default="git", help="git executable (git.exe under WSL)")
    a = ap.parse_args()
    print(f"script test_scripts/fixedeq/explore/skewed_targets.py; commit {git(a.git, 'rev-parse', 'HEAD')}; "
          f"clean {git(a.git, 'status', '--porcelain') == ''}; "
          f"{datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds')}; "
          f"seed {SEED}; draws {DRAWS}; numpy {np.__version__}; box [{nqc.LOWER:.4f}, {nqc.UPPER:.4f}]")

    data = REPO / "test_scripts/test_data"
    turtle = np.asarray(target.pooled_frequencies(readers.read_fasta(data / "turtle_aa.fasta"),
                                                  readers.read_charsets(data / "turtle_aa.nex"), "p"))
    print("\nTurtle target (decision 023) raised to a power and renormalized")
    for alpha in (1, 1.5, 2, 2.5, 3, 4, 5):
        p = turtle ** alpha
        row(f"power {alpha}", p / p.sum())

    print("\nLG's own frequencies, and the oracle's extreme fixture")
    row("LG frequencies", builtin.exchangeabilities("LG")[1])
    extreme = np.full(20, MIN_TARGET)
    extreme[0] = 1 - extreme[1:].sum()
    row("19 entries at 1e-4", extreme)

    print(f"\nDirichlet targets, {DRAWS} draws per concentration (draws NQC would refuse are skipped)")
    rng = np.random.default_rng(SEED)
    for conc in (1.0, 2.0, 5.0, 20.0):
        counts, ratios, skipped = [], [], 0
        for _ in range(DRAWS):
            pi = cases.random_target(20, rng, conc)
            if pi.min() < MIN_TARGET:
                skipped += 1
                continue
            counts.append(outside(pi)[0])
            ratios.append(pi.max() / pi.min())
        c = np.array(counts)
        print(f"  concentration {conc:>4}: kept {len(c):4d}, skipped {skipped:3d}; median max/min "
              f"{np.median(ratios):7.1f}; with any outside {np.mean(c > 0):6.1%}; outside median "
              f"{np.median(c):.0f}, max {c.max()}")


if __name__ == "__main__":
    main()
