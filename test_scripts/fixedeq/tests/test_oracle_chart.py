"""Log-ratio jump chart (D01); pass lines are decision 021's."""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from oracle import builtin, cases, chart  # noqa: E402

TOL = 1e-10


def rel_err(a, b):
    return float(np.max(np.abs(a - b) / np.abs(b)))


def check_family(Q, pi):
    d = chart.diagnostics(Q, pi)
    assert d["stationarity"] <= TOL and d["mean_rate_error"] <= TOL and d["row_sums"] <= TOL, d
    assert d["min_rate"] > 0


def ordinary_cases():
    rng = np.random.default_rng(20261004)
    for n in (3, 4, 20):
        for _ in range(4):
            pi = cases.random_target(n, rng)
            yield n, pi, cases.generator_from_flux(cases.cycle_flux(n, rng), pi)
            yield n, pi, chart.reversible_seed(cases.random_exchangeabilities(n, rng), pi)


@pytest.mark.parametrize("n, pi, Q", list(ordinary_cases()))
def test_round_trip_and_family_membership(n, pi, Q):
    check_family(Q, pi)
    ref = chart.references(Q)
    z = chart.inverse(Q, ref)
    assert z.shape == (n * (n - 2),)
    Q2, nu = chart.build(z, pi, ref)
    check_family(Q2, pi)
    assert rel_err(Q2, Q) <= TOL
    assert np.max(np.abs(nu - pi * -np.diag(Q))) <= TOL


def test_random_coordinates_build_family_members():
    rng = np.random.default_rng(1)
    for n in (3, 4, 20):
        pi = cases.random_target(n, rng)
        ref = rng.integers(0, n - 1, n)
        ref = np.where(ref >= np.arange(n), ref + 1, ref)
        Q, nu = chart.build(rng.normal(0, 2, n * (n - 2)), pi, ref)
        check_family(Q, pi)
        d = chart.diagnostics(Q, pi, nu)
        assert d["nu_residual"] <= TOL and d["exit_rate_mismatch"] <= TOL


def test_positive_and_log_charts_agree():
    rng = np.random.default_rng(2)
    n = 20
    pi = cases.random_target(n, rng)
    ref = chart.references(cases.generator_from_flux(cases.cycle_flux(n, rng), pi))
    z = rng.normal(0, 1.5, n * (n - 2))
    assert rel_err(chart.build_positive(chart.PIN * np.exp(z), pi, ref)[0], chart.build(z, pi, ref)[0]) <= TOL


def test_builder_derivative_matches_central_difference():
    rng = np.random.default_rng(3)
    for n in (3, 4, 20):
        pi = cases.random_target(n, rng)
        ref = chart.references(cases.generator_from_flux(cases.cycle_flux(n, rng), pi))
        z, dz, eps = rng.normal(0, 1, n * (n - 2)), rng.normal(size=n * (n - 2)), 1e-5
        dQ = chart.derivative(z, pi, ref, dz)
        fd = (chart.build(z + eps * dz, pi, ref)[0] - chart.build(z - eps * dz, pi, ref)[0]) / (2 * eps)
        assert np.max(np.abs(fd - dQ)) / max(1.0, np.max(np.abs(dQ))) <= 1e-7


def test_references_take_row_maxima_with_ties_to_the_lower_index():
    Q = np.array([[-3.0, 1.0, 1.0, 1.0],
                  [2.0, -4.0, 1.0, 1.0],
                  [1.0, 3.0, -7.0, 3.0],
                  [0.5, 0.5, 2.0, -3.0]])
    assert list(chart.references(Q)) == [1, 0, 1, 2]


def test_seeds():
    rng = np.random.default_rng(4)
    n = 20
    pi, pi0 = cases.random_target(n, rng), cases.random_target(n, rng)
    R = cases.random_exchangeabilities(n, rng)
    Qr = chart.reversible_seed(R, pi)
    check_family(Qr, pi)
    flux = pi[:, None] * Qr
    assert np.max(np.abs(flux - flux.T)) <= TOL
    off = ~np.eye(n, dtype=bool)
    assert np.max(np.abs((Qr / pi[None, :])[off] / R[off] - (Qr[0, 1] / pi[1]) / R[0, 1])) <= TOL
    Q0 = cases.generator_from_flux(cases.cycle_flux(n, rng), pi0)
    Qt = chart.transfer_seed(Q0, pi)
    check_family(Qt, pi)
    jump = lambda Q: Q[off] / np.repeat(-np.diag(Q), n - 1)
    assert rel_err(jump(Qt), jump(Q0)) <= TOL


def test_lg_rebuilt_at_a_target_is_a_family_member():
    R, pi_lg = builtin.exchangeabilities("LG")
    assert np.allclose(R, R.T) and R[0, 1] == 0.425093 and R.max() == 10.649107
    rng = np.random.default_rng(5)
    for pi in (pi_lg, np.full(20, 0.05), cases.random_target(20, rng)):
        Q = chart.reversible_seed(R, pi)
        check_family(Q, pi)
        ref = chart.references(Q)
        assert rel_err(chart.build(chart.inverse(Q, ref), pi, ref)[0], Q) <= TOL


@pytest.mark.parametrize("bad", ["underflow", "nan"])
def test_invalid_coordinates_raise_and_return_nothing(bad):
    n, pi = 4, np.full(4, 0.25)
    ref = np.array([1, 0, 0, 0])
    z = np.zeros(n * (n - 2))
    z[0] = 1e4 if bad == "underflow" else np.nan
    with pytest.raises(chart.BuildError):
        chart.build(z, pi, ref)


def test_non_positive_target_raises():
    with pytest.raises(chart.BuildError):
        chart.build(np.zeros(8), np.array([0.5, 0.5, 0.0, 0.0]), np.array([1, 0, 0, 0]))


def test_hard_cases_are_recorded_not_asserted():
    """Decision 021: residuals on hard cases are printed (run pytest with -s) and not asserted."""
    rng = np.random.default_rng(6)
    n = 20
    skewed = np.full(n, 1e-4)
    skewed[0] = 1 - skewed[1:].sum()
    ref = np.array([1] + [0] * 9 + [11] + [10] * 9)
    blocks = np.zeros((n, n - 2))
    for i, free in enumerate(chart.free_destinations(ref)):
        blocks[i] = [0.0 if (j < 10) == (i < 10) else -11.5 for j in free]
    rows = []
    for label, pi, z in [
        ("skewed target, min 1e-4", skewed, rng.normal(0, 1, n * (n - 2))),
        ("near-reducible: two blocks joined by logits of -11.5", cases.random_target(n, rng), blocks.ravel()),
        ("wide coordinate range", cases.random_target(n, rng), rng.uniform(-11.5, 2.3, n * (n - 2))),
    ]:
        try:
            Q, nu = chart.build(z, pi, ref)
            d = chart.diagnostics(Q, pi, nu)
            rows.append((label, {k: f"{v:.2e}" for k, v in d.items()}))
        except chart.BuildError as e:
            rows.append((label, f"BuildError: {e}"))
    for label, d in rows:
        print(f"hard case: {label}: {d}")
    assert len(rows) == 3
