"""Numbers the design documents report, recomputed by the oracle (decisions 019, 020, 021)."""
import sys
from pathlib import Path

import numpy as np
from scipy.linalg import expm

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import design_rerun  # noqa: E402
from oracle import chart, t3  # noqa: E402


def test_dimensions_360_189_171():
    n = 20
    pairs = [(i, j) for i in range(n) for j in range(n) if i != j]
    B = np.zeros((n, len(pairs)))
    for e, (i, j) in enumerate(pairs):
        B[j, e], B[i, e] = 1.0, -1.0
    C = np.vstack([B[:-1], np.ones(len(pairs))])
    assert np.linalg.matrix_rank(C) == 20
    assert len(pairs) - 20 == 360 == n * (n - 2)
    assert len(t3.strength_pairs(n)) == 189 and len(t3.tilt_pairs(n)) == 171


def test_triangle_bound_20_72():
    """P-log section 4.5 and D09 (1): forward rates 100, reverse 1e-4, all else 1, uniform target."""
    n = 20
    q = np.ones((n, n))
    q[0, 1] = q[1, 2] = q[2, 0] = 100.0
    q[1, 0] = q[2, 1] = q[0, 2] = 1e-4
    np.fill_diagonal(q, 0.0)
    np.fill_diagonal(q, -q.sum(1))
    pi = np.full(n, 1.0 / n)
    assert chart.diagnostics(q / (pi @ -np.diag(q)), pi)["stationarity"] <= 1e-10
    a, h = t3.inverse(q, pi)
    cycle = h[0, 1] + h[1, 2] + h[2, 0]
    star = h[1, 2] + h[0, 1] - h[0, 2]
    for value in (cycle, star):
        assert design_rerun.reproduced("digits", "20.72", value)


def test_non_concavity_curvature_0_015378():
    """P-positive Theorem 8 and check 5: a rooted two-tip tree under a circulant at uniform pi."""
    def circulant(x):
        return np.array([[-1, (1 + x) / 2, (1 - x) / 2],
                         [(1 - x) / 2, -1, (1 + x) / 2],
                         [(1 + x) / 2, (1 - x) / 2, -1.0]])
    t1, t2 = 0.1, 5.1
    k = np.sqrt(3) * (t2 - t1) / 2
    amp = 2 * np.exp(-1.5 * (t1 + t2))
    x0, step = np.pi / k, 1e-3

    def log_l00(x):
        return np.log(expm(t1 * circulant(x))[:, 0] @ expm(t2 * circulant(x))[:, 0] / 3)
    numeric = (log_l00(x0 + step) - 2 * log_l00(x0) + log_l00(x0 - step)) / step ** 2
    for value in (numeric, amp * k * k / (1 - amp)):
        assert design_rerun.reproduced("digits", "0.015378", value)
