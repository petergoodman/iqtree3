"""The reference climber of PLAN.md risks 22 and 24 (oracle/optimize.py lbfgsb_fit)."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from oracle import cases, chart, likelihood as lk, nqc, optimize as opt, readers  # noqa: E402

TREE = "((a:0.1,b:0.2):0.05,(c:0.3,d:0.1):0.2);"


def test_projected_gradient_is_zero_only_where_a_bound_blocks_descent():
    lower, upper = np.zeros(4), np.ones(4)
    x = np.array([0.0, 1.0, 0.5, 0.0])
    g = np.array([0.3, -0.2, 0.125, -0.375])
    assert list(opt.projected_gradient(x, g, lower, upper)) == [0.0, 0.0, 0.125, -0.375]


def test_batched_gradients_equal_the_serial_rules_bitwise():
    f = lambda x: float(np.sum(np.sin(x) * np.arange(1, len(x) + 1)) + x[0] * x[-1])  # noqa: E731
    x = np.array([0.0, 1e-9, -3.0, 2.302585092994046, 0.7])
    upper = np.full(5, 2.302585092994046)
    serial = lambda pts: [f(p) for p in pts]  # noqa: E731
    assert np.array_equal(opt.forward_gradient(serial, x, upper, f(x)), opt.scaled_derivative(f, x.copy(), upper)[1])
    assert np.array_equal(opt.central_gradient(serial, x, 1e-4), opt.central_derivative(f, x.copy(), [1e-4])[0])


def test_clipped_optimum_of_a_separable_quadratic():
    a = np.array([1.0, 10.0, 100.0, 3.0, 50.0])
    m = np.array([-2.0, 0.5, 7.0, -20.0, 1.0])
    lower, upper = np.full(5, -5.0), np.full(5, 5.0)
    r = opt.lbfgsb_fit(lambda x: float(0.5 * np.sum(a * (x - m) ** 2)), np.zeros(5), lower, upper, pgtol=1e-8)
    assert np.max(np.abs(r.x - np.clip(m, lower, upper))) <= 1e-6
    assert list(r.at_lower) == [3] and list(r.at_upper) == [2]
    assert r.pg_max <= 1e-8
    assert r.phases[-1][0] == "central"   # forward differences cannot reach 1e-8 at curvature 100


def test_ends_no_lower_than_iqtrees_optimizer_on_an_nqc_objective():
    rng = np.random.default_rng(7)
    pi = cases.random_target(20, rng)
    ref = chart.references(cases.generator_from_flux(cases.cycle_flux(20, rng), pi))
    aln = {n: "".join(rng.choice(list(lk.AMINO), 25)) + "B-X" for n in "abcd"}
    f = nqc.Objective(aln, readers.parse_newick(TREE), pi, ref)
    z0 = rng.normal(-1.0, 0.5, 360)
    lower, upper = np.full(360, nqc.LOWER), np.full(360, nqc.UPPER)
    port = nqc.fit(f, z0)
    r = opt.lbfgsb_fit(f, z0, lower, upper, max_iter=20)   # 134 iterations to converge; 20 keep the test short
    assert np.all((r.x >= lower) & (r.x <= upper))
    assert -r.f >= port.score
    assert r.f == f(r.x)
