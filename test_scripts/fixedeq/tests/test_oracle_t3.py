"""T3 balancing chart in the D04 gauge; pass lines are decision 021's."""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from oracle import cases, chart, t3  # noqa: E402

TOL = 1e-10


def rel_err(a, b):
    return float(np.max(np.abs(a - b) / np.abs(b)))


def flux_cases():
    rng = np.random.default_rng(20261005)
    for n in (3, 4, 20):
        for _ in range(4):
            pi = cases.random_target(n, rng)
            yield n, pi, cases.generator_from_flux(cases.cycle_flux(n, rng), pi)


@pytest.mark.parametrize("n, pi, Q", list(flux_cases()))
def test_round_trip_in_the_d04_gauge(n, pi, Q):
    a, h = t3.inverse(Q, pi)
    assert a[n - 2, n - 1] == 0 and np.all(h[:, n - 1] == 0)
    theta = t3.pack(a, h)
    assert theta.size == n * (n - 1) // 2 - 1 + (n - 1) * (n - 2) // 2
    Q2 = t3.build(*t3.unpack(theta, n), pi)
    d = chart.diagnostics(Q2, pi)
    assert d["stationarity"] <= TOL and d["mean_rate_error"] <= TOL and d["row_sums"] <= TOL
    assert rel_err(Q2, Q) <= TOL


@pytest.mark.parametrize("n, pi, Q", list(flux_cases()))
def test_t3_and_the_jump_chart_give_the_same_matrix(n, pi, Q):
    ref = chart.references(Q)
    Qj = chart.build(chart.inverse(Q, ref), pi, ref)[0]
    Qt = t3.build(*t3.inverse(Q, pi), pi)
    assert rel_err(Qt, Qj) <= TOL


def test_zero_tilts_give_the_reversible_matrix():
    rng = np.random.default_rng(7)
    for n in (3, 4, 20):
        pi = cases.random_target(n, rng)
        a = rng.normal(size=(n, n))
        a = (a + a.T) / 2
        Q = t3.build(a, np.zeros((n, n)), pi)
        assert rel_err(Q, chart.reversible_seed(np.exp(a) * (1 - np.eye(n)), pi)) <= TOL


def test_builder_derivative_matches_central_difference():
    rng = np.random.default_rng(8)
    for n in (3, 4, 20):
        pi = cases.random_target(n, rng)
        a, h = t3.inverse(cases.generator_from_flux(cases.cycle_flux(n, rng), pi), pi)
        da = rng.normal(size=(n, n))
        da = (da + da.T) / 2
        dh = rng.normal(size=(n, n))
        dh = (dh - dh.T) / 2
        eps = 1e-5
        dQ = t3.derivative(a, h, pi, da, dh)
        fd = (t3.build(a + eps * da, h + eps * dh, pi) - t3.build(a - eps * da, h - eps * dh, pi)) / (2 * eps)
        assert np.max(np.abs(fd - dQ)) / max(1.0, np.max(np.abs(dQ))) <= 1e-7


def test_coordinates_from_the_star_gauge_at_state_1_convert_exactly():
    """P-positive's gauge (verify_constrained_nq.py, gauge_fix): r scaled so r_{n-1,n} = 1, and
    tilts h'_ij = h_ij + h_1i - h_1j, which are zero on the star at state 1."""
    rng = np.random.default_rng(9)
    n = 20
    pi = cases.random_target(n, rng)
    Q = cases.generator_from_flux(cases.cycle_flux(n, rng), pi)
    M = Q / pi[None, :]
    np.fill_diagonal(M, 1.0)
    r, h = np.sqrt(M * M.T), 0.5 * np.log(M / M.T)
    r_star = r / r[n - 2, n - 1]
    h_star = h + h[0][:, None] - h[0][None, :]
    assert np.max(np.abs(h_star[0])) == 0
    a, h_d04 = t3.to_d04(np.log(r_star), h_star)
    assert rel_err(t3.build(a, h_d04, pi), Q) <= TOL
