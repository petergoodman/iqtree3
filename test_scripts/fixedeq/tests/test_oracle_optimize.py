"""The port of IQ-TREE's optimizer (decision 017), the derivative rules (decision 005), and SIM-2's
stage B climb.

Each port test checks one branch of utils/optimization.cpp at 63c330d9 against what the C++ does.
"""
import itertools
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from oracle import optimize as opt  # noqa: E402


class Recorder:
    """An objective that records every point it is evaluated at, as targetFunk sees them."""
    def __init__(self, f):
        self.f, self.points = f, []

    def __call__(self, x):
        self.points.append(np.array(x, dtype=float))
        return self.f(x)


def quadratic(A, b):
    return lambda x: 0.5 * x @ A @ x - b @ x


def test_fix_bound_clamps_each_coordinate():
    x = np.array([-1.0, 0.5, 9.0])
    opt.fix_bound(x, np.zeros(3), np.ones(3))
    assert list(x) == [0.0, 0.5, 1.0]


def test_legacy_derivative_steps_and_evaluations():
    f = Recorder(lambda x: float(np.sum(x ** 3)))
    x = np.array([0.0, 3.0, -2.0])
    fx, g = opt.legacy_derivative(f, x.copy())
    assert fx == 19.0
    assert len(f.points) == 4 and np.array_equal(f.points[0], x)
    steps = [f.points[k + 1][k] - x[k] for k in range(3)]
    assert steps[0] == 1e-4
    assert steps[1] == (3.0 + 1e-4 * 3.0) - 3.0 and steps[2] == (-2.0 + 1e-4 * 2.0) - (-2.0)
    assert np.array_equal(g, [(f.f(f.points[k + 1]) - fx) / steps[k] for k in range(3)])


def test_legacy_derivative_ignores_bounds():
    f = Recorder(lambda x: float(x[0]))
    opt.legacy_derivative(f, np.array([100.0]))
    assert f.points[1][0] > 100.0


def test_scaled_step_is_representable_and_backward_at_the_upper_bound():
    for x in (0.0, 1e-9, -0.3, 2.302585092994046, -11.5):
        h = opt.scaled_step(x, upper=np.log(10.0))
        assert (x + h) - x == h and abs(h) >= 0.5e-4 * max(1.0, abs(x))
        assert x + h <= np.log(10.0) or x > np.log(10.0)
    assert opt.scaled_step(np.log(10.0), upper=np.log(10.0)) < 0
    assert opt.scaled_step(0.0, upper=1.0) > 0


def test_scaled_derivative_uses_the_signed_step():
    f = Recorder(lambda x: float(np.sum(x ** 2)))
    x, upper = np.array([0.0, 1.0]), np.array([1.0, 1.0])
    _, g = opt.scaled_derivative(f, x.copy(), upper)
    assert f.points[2][1] < 1.0
    assert g[1] < 2.0


def test_line_search_failure_restores_x_and_keeps_the_last_trial_value():
    """lnsrch 687-690: x is reset to xold, check = 1, and f stays at the rejected trial."""
    f = Recorder(lambda x: float(x[0] ** 2))
    xold, g, p = np.array([1.0]), np.array([-1.0]), np.array([1.0])   # claims descent where f rises
    x, fret, check = opt.lnsrch(f, xold, f(xold), g, p.copy(), 100.0, np.array([-10.0]), np.array([10.0]))
    assert check == 1 and np.array_equal(x, xold)
    assert fret == f.f(f.points[-1]) and fret != 1.0


def test_line_search_accepts_the_full_newton_step():
    A, b = np.diag([2.0, 4.0]), np.array([2.0, 4.0])
    f = quadratic(A, b)
    x0 = np.zeros(2)
    g = A @ x0 - b
    x, fret, check = opt.lnsrch(f, x0, f(x0), g, -np.linalg.solve(A, g), 100.0, np.full(2, -10.0), np.full(2, 10.0))
    assert check == 0 and np.allclose(x, [1.0, 1.0], rtol=0, atol=1e-15)


def test_hessian_update_also_runs_on_negative_curvature():
    """dfpmin 875 tests fac*fac > EPS*sumdg*sumxi, so a negative fac still updates."""
    H = np.eye(2)
    H2, updated = opt.update_inverse_hessian(H.copy(), dg=np.array([-1.0, 0.0]), xi=np.array([1.0, 0.0]))
    assert updated and not np.array_equal(H2, H)
    H3, updated = opt.update_inverse_hessian(H.copy(), dg=np.array([0.0, 1e-9]), xi=np.array([1.0, 0.0]))
    assert not updated and np.array_equal(H3, H)


def test_dfpmin_stops_on_the_gradient_or_step_test_for_a_quadratic():
    A, b = np.array([[3.0, 1.0], [1.0, 2.0]]), np.array([1.0, -1.0])
    f = quadratic(A, b)
    r = opt.dfpmin(f, np.array([2.0, 2.0]), np.full(2, -10.0), np.full(2, 10.0), gtol=1e-4)
    assert r.reason in ("gtol", "TOLX") and r.iterations < opt.ITMAX
    assert f(r.x) < f(np.array([2.0, 2.0]))


def test_dfpmin_returns_silently_after_itmax():
    f = lambda x: float((1 - x[0]) ** 2 + 100 * (x[1] - x[0] ** 2) ** 2)
    r = opt.dfpmin(f, np.array([-1.2, 1.0]), np.full(2, -10.0), np.full(2, 10.0), gtol=1e-12, itmax=2)
    assert r.reason == "ITMAX" and r.iterations == 2


def test_dfpmin_after_a_failed_search_returns_the_stale_value():
    """At the corner of |x| every trial rises, so the search fails; the zero displacement then
    meets the TOLX test, and fret is the rejected trial's value, as in the C++."""
    f = Recorder(lambda x: float(abs(x[0])))
    r = opt.dfpmin(f, np.array([0.0]), np.array([-10.0]), np.array([10.0]), gtol=1e-30)
    assert r.failed_searches == 1 and r.reason == "TOLX" and r.iterations == 1
    assert np.array_equal(r.x, [0.0]) and r.f == f.f(f.points[-1]) and r.f > 0


def test_minimize_restarts_only_when_a_checked_coordinate_ends_at_a_bound():
    f = lambda x: float(np.sum((x - 5.0) ** 2))
    lower, upper = np.zeros(2), np.full(2, 3.0)
    draws = iter([0.5, 0.5, 0.25, 0.25, 0.1, 0.1, 0.9, 0.9])
    once = opt.minimize_multi_dimen(f, np.ones(2), lower, upper, np.zeros(2, bool), 1e-4)
    assert once.runs == 1
    many = opt.minimize_multi_dimen(f, np.ones(2), lower, upper, np.ones(2, bool), 1e-4, random_double=lambda: next(draws))
    assert many.runs == 1 + opt.MAX_ITER and np.allclose(many.x, upper)


def test_model_call_recomputes_the_score_at_the_returned_point():
    """ModelMarkov::optimizeParameters 1199-1235: gtol = max(epsilon, TOL_RATE), and the score is
    recomputed when the returned variables differ from the model's last evaluated state."""
    f = Recorder(lambda x: float(np.sum((x - 1.5) ** 2)))
    r = opt.optimize_model(f, np.array([3.0, 0.5]), np.full(2, 1e-4), np.full(2, 100.0), gradient_epsilon=1e-6)
    assert r.gtol == opt.TOL_RATE
    assert r.score == -f.f(r.x)


# ---------------------------------------------------------------- SIM-2's stage B

BOX = (np.full(2, -5.0), np.full(2, 5.0))
START = np.array([-1.2, 1.0])
LOOSE = dict(gain=np.inf, c_move=np.inf, low_move=np.inf)


def rosenbrock(x):
    from scipy.optimize import rosen, rosen_der
    return rosen(x), rosen_der(x)


def watch_of(fluxes, scale=1000.0):
    """C = x[0]; log rates 0 and scale x[0], with the given fluxes."""
    return lambda x: (x[0], np.array([0.0, scale * x[0]]), np.asarray(fluxes, dtype=float))


def test_low_flux_is_the_lowest_two_fifths():
    assert sorted(opt.low_flux([5.0, 1.0, 4.0, 2.0, 3.0, 9.0, 8.0, 7.0, 6.0, 10.0])) == [1, 2, 3, 4]
    assert len(opt.low_flux(np.arange(380.0))) == 152


def test_stage_b_stalls_once_a_full_window_meets_all_three_conditions():
    r = opt.stage_b_climb(rosenbrock, START, *BOX, watch_of([1.0, 2.0]), 1e9, window=5, **LOOSE)
    assert r.stop == "stalled" and r.iterations == 5 and len(r.trace) == 6
    assert r.window["span"] == 5
    assert r.window["lnl_rise"] == r.trace[0][2] - r.trace[5][2]
    assert r.window["c_move"] == abs(r.trace[5][3] - r.trace[0][3])


def test_stage_b_does_not_stall_while_any_condition_fails():
    for tight in (dict(gain=0.0), dict(c_move=0.0), dict(low_move=-1.0)):
        r = opt.stage_b_climb(rosenbrock, START, *BOX, watch_of([1.0, 2.0]), 1e9, window=5, **(LOOSE | tight))
        assert r.stop == r.message and r.iterations > 5


def test_stage_b_watches_only_the_low_flux_rates():
    tight = LOOSE | dict(low_move=0.01)
    still = opt.stage_b_climb(rosenbrock, START, *BOX, watch_of([1.0, 2.0]), 1e9, window=5, **tight)
    assert still.stop == "stalled" and still.iterations == 5 and still.window["low_flux_move"] == 0.0
    moving = opt.stage_b_climb(rosenbrock, START, *BOX, watch_of([2.0, 1.0]), 1e9, window=5, **tight)
    assert moving.iterations > 5


def test_stage_b_time_cap_is_checked_after_each_iteration():
    ticks = itertools.count()
    r = opt.stage_b_climb(rosenbrock, START, *BOX, watch_of([1.0, 2.0]), 3.0, clock=lambda: float(next(ticks)))
    assert r.stop == "time cap" and r.iterations == 3


def test_stage_b_on_a_quadratic_ends_at_pgtol_with_a_bound_held():
    """The optimum A^-1 b = (-1.2, 4.6) lies above the box's 3 in x[1]; with x[1] at 3 the best x[0] is
    -2/3, where the gradient pushes x[1] out of the box. The start (4, 4) is clipped to (3, 3)."""
    A, b = np.array([[3.0, 1.0], [1.0, 2.0]]), np.array([1.0, 8.0])
    f = quadratic(A, b)
    r = opt.stage_b_climb(lambda x: (f(x), A @ x - b), np.array([4.0, 4.0]), np.full(2, -5.0), np.full(2, 3.0),
                          lambda x: (0.0, x, np.ones(2)), 1e9)
    assert r.stop == r.message and "CONVERGENCE" in r.message
    assert r.f_start == f(np.array([3.0, 3.0])) and r.f <= r.f_start and r.gain == r.f_start - r.f
    assert r.pg_max <= 1e-6 and list(r.at_upper) == [1] and abs(r.x[0] + 2 / 3) <= 1e-6
