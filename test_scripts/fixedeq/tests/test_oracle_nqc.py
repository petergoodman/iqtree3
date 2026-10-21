"""The S2 objective, derivative helpers and toy problem (decision 030)."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "nqc"))
from oracle import cases, chart, likelihood as lk, nqc, optimize, readers  # noqa: E402
import toy  # noqa: E402

TREE = "((a:0.1,b:0.2):0.05,(c:0.3,d:0.1):0.2);"


def small_problem(seed=7):
    rng = np.random.default_rng(seed)
    pi = cases.random_target(20, rng)
    ref = chart.references(cases.generator_from_flux(cases.cycle_flux(20, rng), pi))
    aln = {n: "".join(rng.choice(list(lk.AMINO), 25)) + "B-X" for n in "abcd"}
    return pi, ref, aln, readers.parse_newick(TREE), rng.normal(-1.0, 0.5, 360)


def test_objective_is_the_likelihood_of_the_built_matrix():
    pi, ref, aln, tree, z = small_problem()
    Q = chart.build(z, pi, ref)[0]
    for shape, p_invar in ((None, 0.0), (0.5, 0.0), (None, 0.2), (0.5, 0.2)):
        f = nqc.Objective(aln, tree, pi, ref, shape=shape, p_invar=p_invar)
        assert f(z) == -lk.alignment_log_likelihood(aln, tree, Q, pi, shape=shape, p_invar=p_invar)
    assert f.evaluations == 1


def test_a_failed_build_scores_1e30():
    pi, ref, aln, tree, z = small_problem()
    z[5] = np.nan
    assert nqc.Objective(aln, tree, pi, ref)(z) == 1.0e30


def test_central_differences_of_a_quadratic():
    A = np.array([[3.0, 1.0, 0.0], [1.0, 2.0, 0.5], [0.0, 0.5, 4.0]])
    b = np.array([1.0, -2.0, 0.5])
    x = np.array([0.3, -2.0, 1e-9])
    g, c = nqc.central(lambda v: 0.5 * v @ A @ v - b @ v, x, 1e-3)
    assert np.allclose(g, A @ x - b, rtol=0, atol=1e-9)
    assert np.allclose(c, np.diag(A), rtol=0, atol=1e-6)


def test_legacy_step_collapses_near_zero_and_the_scaled_step_does_not():
    z = np.array([0.0, 1e-9, -3.0, nqc.UPPER])
    legacy, scaled = nqc.legacy_steps(z), nqc.scaled_steps(z)
    assert legacy[0] == optimize.ERROR_X and abs(legacy[1] - 1e-13) < 1e-20
    assert abs(scaled[1] - 1e-4) < 1e-12 and scaled[2] == (-3.0 + 1e-4 * 3.0) - (-3.0)
    assert scaled[3] < 0   # backward at the upper bound (decision 005)


def test_fit_reaches_an_interior_minimum():
    A = np.diag([2.0, 1.0, 4.0])
    m = np.array([-1.0, -5.0, 0.5])
    r = nqc.fit(lambda v: float(0.5 * (v - m) @ A @ (v - m)), np.array([-3.0, -3.0, -3.0]))
    assert np.allclose(r.x, m, atol=1e-3) and r.gtol == nqc.GRADIENT_EPSILON


def test_toy_tree_is_rooted_binary_on_eight_taxa():
    tree = readers.parse_newick(toy.tree_text())
    assert toy.tree_text() == toy.tree_text()
    assert len(tree.children) == 2
    assert sorted(n.name for n in tree.leaves()) == [f"T{i}" for i in range(1, 9)]
    nodes = [n for n in tree.preorder() if n is not tree]
    assert len(nodes) == 14 and all(n.children == [] or len(n.children) == 2 for n in nodes)
    assert all(toy.LENGTHS[0] <= round(n.length, 4) <= toy.LENGTHS[1] for n in nodes)


def test_toy_truth_is_stationary_at_the_target():
    pi = cases.random_target(20, np.random.default_rng(3))
    d = chart.diagnostics(toy.truth(pi), pi)
    assert d["stationarity"] <= 1e-10 and d["mean_rate_error"] <= 1e-10   # decision 021
