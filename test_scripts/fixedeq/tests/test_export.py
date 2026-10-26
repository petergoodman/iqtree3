"""S4's oracle side of the export: the first-order label, the export gradient and the file readers
(decision 009, decision 010, CODE_PLAN.md sections 2.1 and 2.2)."""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "differential"))
from oracle import chart, nqc, readers  # noqa: E402
import differential  # noqa: E402

LO, HI = -2.0, 3.0


def label(x, g, f=100.0):
    return chart.first_order_test(np.array(x, dtype=float), np.array(g, dtype=float), LO, HI, f)


def test_interior_coordinates_pass_below_gtol_and_fail_above_it():
    # s = |g| max(|x|, 1) / max(|f|, 1); with f = 100 and |x| <= 1, s = |g| / 100
    assert label([0.5, -1.0], [0.009, -0.0099])["passed"]
    above = label([0.5, -1.0], [0.0, 0.02])
    assert not above["passed"] and above["statistic"] == 0.02 / 100 and above["coordinate"] == 1
    # |x| > 1 scales the component up
    assert label([2.5, 0.0], [0.005, 0.0])["statistic"] == 0.005 * 2.5 / 100


def test_a_bound_counts_only_the_part_pointing_into_the_box():
    # g is the derivative of -lnL: at the lower bound g < 0 asks to move up, into the box
    lower_out = label([LO, 0.0], [5.0, 0.0])
    assert lower_out["passed"] and lower_out["at_lower"] == 1 and lower_out["statistic"] == 0.0
    lower_in = label([LO, 0.0], [-5.0, 0.0])
    assert not lower_in["passed"] and lower_in["statistic"] == 5.0 * 2.0 / 100
    upper_out = label([0.0, HI], [0.0, -5.0])
    assert upper_out["passed"] and upper_out["at_upper"] == 1
    upper_in = label([0.0, HI], [0.0, 5.0])
    assert not upper_in["passed"] and upper_in["coordinate"] == 1 and upper_in["statistic"] == 5.0 * 3.0 / 100


def test_a_coordinate_outside_the_box_fails_with_a_zero_gradient():
    r = label([HI + 1.0, 0.0], [0.0, 0.0])
    assert not r["passed"] and r["outside"] == 1 and r["statistic"] == 0.0


def test_ties_go_to_the_lower_index_and_f_below_one_is_not_scaled():
    r = chart.first_order_test(np.array([0.0, 0.0, 0.0]), np.array([1e-5, -2e-5, 2e-5]), LO, HI, 0.5)
    assert r["coordinate"] == 1 and r["statistic"] == 2e-5 and r["passed"]


@pytest.mark.parametrize("x, g, f, gtol", [
    ([], [], 1.0, 1e-4),
    ([0.0, 1.0], [0.0], 1.0, 1e-4),
    ([np.nan], [0.0], 1.0, 1e-4),
    ([0.0], [np.inf], 1.0, 1e-4),
    ([0.0], [0.0], np.nan, 1e-4),
    ([0.0], [0.0], 1.0, 0.0),
])
def test_invalid_input_is_refused(x, g, f, gtol):
    with pytest.raises(ValueError):
        chart.first_order_test(np.array(x, dtype=float), np.array(g, dtype=float), LO, HI, f, gtol)


def test_export_gradient_is_central_over_the_probes_displacements():
    c = np.array([1.0, -3.0, 0.25])

    def f(z):
        return float(c @ z ** 2 + z[0] * z[1])

    z = np.array([0.3, -2.0, 7.0])
    f0, g, span = nqc.export_gradient(f, z)
    exact = 2 * c * z + np.array([z[1], z[0], 0.0])
    for k in range(3):
        h = chart.FIRST_ORDER_ETA * max(1.0, abs(z[k]))
        assert span[k] == ((z[k] + h) - z[k]) + (z[k] - (z[k] - h))
    # a quadratic's central difference is exact up to rounding
    assert np.allclose(g, exact, rtol=1e-9, atol=1e-9)
    assert f0 == f(z)


def test_qmat_reader_returns_every_number_written(tmp_path):
    rng = np.random.default_rng(3)
    Q = rng.uniform(0.01, 2.0, (20, 20))
    np.fill_diagonal(Q, 0.0)
    np.fill_diagonal(Q, -Q.sum(axis=1))
    pi = rng.dirichlet(np.ones(20))
    path = tmp_path / "x.qmat"
    differential.write_matrix(path, Q, pi)
    Q2, pi2 = readers.read_qmat(path.read_text(encoding="utf-8"))
    assert np.array_equal(Q2, Q) and np.array_equal(pi2, pi)
    with pytest.raises(ValueError):
        readers.read_qmat("1 2 3 4 5")


def test_info_reader_keeps_values_as_text_and_refuses_repeats():
    info = readers.read_info("# NQC export\nmodel: NQC\ncommand:  -s a.phy\n\ntarget: 0.5,0.5\n")
    assert info == {"model": "NQC", "command": " -s a.phy", "target": "0.5,0.5"}
    for bad in ("model: NQC\nmodel: NQC\n", "no separator\n"):
        with pytest.raises(ValueError):
            readers.read_info(bad)
