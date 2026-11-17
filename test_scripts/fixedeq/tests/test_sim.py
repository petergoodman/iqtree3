"""SIM-1's drivers under test_scripts/fixedeq/sim/ (CODE_PLAN.md 3.8), on constructed inputs and committed
files only, so CI runs them without a binary. Lines as set in each commit's plan."""
import sys
from collections import Counter
from pathlib import Path

import numpy as np
from scipy import stats

FIXEDEQ = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(FIXEDEQ))
sys.path.insert(0, str(FIXEDEQ / "sim"))
from oracle import cases, likelihood as lk, simulate  # noqa: E402
import alisim_check as ac  # noqa: E402
import common  # noqa: E402


def test_committed_target_is_the_oracles():
    assert np.array_equal(common.target_values(), simulate.target())


def test_exact_law_against_the_direct_product():
    """021's enumeration line: the law through site_log_likelihoods against P(t_a)^T diag(pi) P(t_b)."""
    rng = np.random.default_rng(31)
    pi = cases.random_target(20, rng)
    Q = cases.generator_from_flux(cases.cycle_flux(20, rng), pi)
    for rates, weights in ((np.ones(1), np.ones(1)), lk.rate_categories(shape=0.5)):
        law = ac.exact_law(ac.T_JOINT, Q, pi, rates, weights)
        direct = ac.direct_law(Q, pi, *ac.LENGTHS, rates, weights)
        assert np.max(np.abs(law - direct) / direct) <= 1e-12
        assert abs(law.sum() - 1) <= 1e-12


def test_pearson_pools_small_cells():
    observed = np.array([10, 20, 30, 2, 2, 1])
    expected = np.array([8.0, 16.0, 32.0, 2.0, 1.5, 0.5])
    # the three cells below 5 make one of 4.0, still below 5, merged into the smallest other (8.0)
    o, e = np.array([15.0, 20.0, 30.0]), np.array([12.0, 16.0, 32.0])
    got = ac.pearson(observed, expected)
    x2 = float(((o - e) ** 2 / e).sum())
    assert (got["cells"], got["df"], got["pooled"], got["smallest_expected"]) == (3, 2, 3, 12.0)
    assert got["x2"] == x2 and got["p"] == float(stats.chi2.sf(x2, 2))
    unpooled = ac.pearson(o, e)
    assert (unpooled["cells"], unpooled["pooled"], unpooled["x2"]) == (3, 0, x2)


def test_counts_against_counter():
    rng = np.random.default_rng(32)
    a = "".join(rng.choice(list(lk.AMINO), 3000))
    b = "".join(rng.choice(list(lk.AMINO), 3000))
    pairs = Counter(zip(a, b))
    want = np.array([[pairs[(x, y)] for y in lk.AMINO] for x in lk.AMINO])
    assert np.array_equal(ac.pair_counts(a, b), want)
    single = Counter(a)
    assert np.array_equal(ac.residue_counts(a)[0], [single[x] for x in lk.AMINO])
    assert ac.residue_counts(a)[1] == 0
    counts, other = ac.residue_counts("AR-XNb")
    assert other == 3 and counts.sum() == 3


def test_continuous_gamma_has_mean_one():
    rates, weights = ac.continuous_gamma(ac.ALPHA)
    assert abs(rates @ weights - 1) <= 1e-12
