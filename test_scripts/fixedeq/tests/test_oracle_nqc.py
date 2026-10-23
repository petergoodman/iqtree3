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


def test_linked_objective_sums_partitions_in_the_order_given():
    pi, ref, aln, tree, z = small_problem()
    half = {n: s[:14] for n, s in aln.items()}
    rest = {n: s[14:] for n, s in aln.items()}
    f1, f2 = nqc.Objective(half, tree, pi, ref, shape=0.5), nqc.Objective(rest, tree, pi, ref, part_rate=2.0)
    assert nqc.LinkedObjective([f1])(z) == f1(z)
    assert nqc.LinkedObjective([f1, f2])(z) == (0.0 + f1(z)) + f2(z)
    z[3] = np.nan
    assert nqc.LinkedObjective([f1, f2])(z) == 2.0e30   # as PartitionModel::targetFunk sums 1e30 per partition


def test_a_partition_rate_scales_every_branch():
    pi, ref, aln, tree, z = small_problem()
    scaled = readers.parse_newick(TREE)
    for node in scaled.preorder():
        if node.length is not None:
            node.length *= 1.7
    folded = nqc.Objective(aln, tree, pi, ref, shape=0.5, p_invar=0.2, part_rate=1.7)(z)
    assert abs(folded - nqc.Objective(aln, scaled, pi, ref, shape=0.5, p_invar=0.2)(z)) <= 1e-12 * abs(folded)


def test_normalize_part_rates_divides_by_the_site_weighted_mean():
    assert nqc.normalize_part_rates([1.0, 3.0], [1, 1]) == [0.5, 1.5]
    r = nqc.normalize_part_rates([0.7, 1.9, 0.2], [172, 216, 386])
    assert abs(sum(a * n for a, n in zip(r, (172, 216, 386))) / 774 - 1.0) <= 1e-15
    assert abs(r[1] / r[0] - 1.9 / 0.7) <= 1e-15


def test_twopart_problems_finds_unknown_and_identical_sequences():
    import twopart
    assert twopart.problems({"a": "ARND", "b": "ARNC"}) == []
    assert twopart.problems({"a": "--X?", "b": "ARNC"}) == ["a has no known character"]
    assert twopart.problems({"a": "ARND", "b": "ARND", "c": "AAAA"}) == ["a and b are identical"]
