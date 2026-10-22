import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "nqc"))
import fitrecords  # noqa: E402


def block(n, start, rescored, upper, fall=False):
    lines = [f"NQC fit {n}: start lnL {start!r}, optimizer lnL {rescored!r}, re-scored lnL {rescored!r}"]
    if fall:
        lines.append(f"NQC fit {n}: the re-scored lnL is below the start by {start - rescored!r} (reported, not restored)")
    lines += [f"NQC coordinates at the lower bound 0, at the upper bound {upper}, outside the domain 0",
              "NQC residuals: |pi*Q| 1e-17, row sums 2e-16, |mean rate - 1| 0, nu 3e-17, minimum rate 1e-05"]
    return "\n".join(lines) + "\n"


def test_a_record_printed_twice_counts_once_and_a_fall_is_kept():
    a, b = block(1, -10.0, -9.0, 2), block(2, -9.0, -9.5, 3, fall=True)
    lnl, bounds, _ = fitrecords.records("Gamma shape\n" + a + "BEST\n" + a + b + "x\n" + b)
    assert bounds == [(0, 2, 0), (0, 3, 0)]
    assert [r < s for s, _, r in lnl] == [False, True]


def test_distinct_fits_sharing_a_number_are_both_kept():
    # the +I+G restore puts the class's fit counter back (model/modelnonrevfixedeq.cpp:422)
    lnl, _, _ = fitrecords.records(block(3, -10.0, -9.0, 1) + block(3, -9.0, -8.5, 1))
    assert len(lnl) == 2


def test_bound_counts_use_exact_equality_as_the_class_does():
    z = np.array([-2.0, 1.0, np.nextafter(1.0, 0.0), np.nextafter(-2.0, -3.0), 0.0])
    assert fitrecords.bound_counts(z, -2.0, 1.0) == (1, 1, 1)
