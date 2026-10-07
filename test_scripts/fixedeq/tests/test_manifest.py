import sys
from pathlib import Path

import numpy as np
import pytest

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


@pytest.mark.skipif(not manifest.INCUMBENT_CKP.exists(), reason="the probe (g) record is not on this machine")
def test_incumbent_checkpoint_unpacks_to_the_reported_exchangeabilities():
    R = manifest.incumbent_exchangeabilities()
    text = (manifest.INCUMBENT_CKP.parent / "probe.iqtree").read_text()
    rows = text[text.index("PAML format"):].splitlines()[2:21]   # the 19 lower-triangle rows
    # the report rounds to 6 decimals (5e-7) and the checkpoint to 10 significant digits (5e-8 at 100)
    for i, row in enumerate(rows, start=1):
        assert np.allclose([float(x) for x in row.split()], R[i, :i], rtol=0, atol=6e-7)
