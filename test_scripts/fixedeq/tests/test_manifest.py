import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "manifest"))
import manifest  # noqa: E402
from oracle import builtin, chart  # noqa: E402


def test_state_order_read_from_source_is_the_oracle_order():
    assert manifest.state_order() == builtin.STATES


def test_lg_seed_is_stationary_at_the_target():
    pi = manifest.test_target()
    d = chart.diagnostics(manifest.lg_seed(pi), pi)
    assert d["stationarity"] <= 1e-10 and d["mean_rate_error"] <= 1e-10   # decision 021


def test_references_are_the_row_maxima():
    Q = manifest.lg_seed(manifest.test_target())
    ref = chart.references(Q)
    off = Q.copy()
    np.fill_diagonal(off, -np.inf)
    assert all(off[i, ref[i]] == off[i].max() for i in range(20))
