"""SIM-1's oracle additions (oracle/simulate.py, builtin.frequency_profiles, the tip-array cache;
CODE_PLAN.md 3.1 and 3.8; decisions 041, 044), on constructed inputs and the repository's model
source. Lines as set in the SIM-1 commit-2 plan, citing decision 021's rows where they apply."""
import itertools
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from oracle import builtin, cases, chart, flux, likelihood as lk, nqc, readers, simulate as sim, t3  # noqa: E402

THREE = "((a:0.2,b:0.3):0.1,c:0.4);"
FOUR = "((a:0.1,b:0.2):0.05,(c:0.3,d:0.1):0.2);"
NAMES = ("a", "b", "c")
N = 1e6


def rel(a, b):
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    return float(np.max(np.abs(a - b) / np.abs(b)))


def residuals(Q, pi):
    d = chart.diagnostics(Q, pi)
    return max(d["stationarity"], d["row_sums"], d["mean_rate_error"])


def setting(seed, n=20):
    rng = np.random.default_rng(seed)
    pi1, pi2 = cases.random_target(n, rng), cases.random_target(n, rng)
    R = cases.random_exchangeabilities(n, rng)
    Qc = cases.generator_from_flux(cases.cycle_flux(n, rng), pi1)
    return rng, pi1, pi2, R, Qc


# ---------------------------------------------------------------- tip-array cache

def test_cached_tips_equal_the_columns_path_bitwise():
    rng, pi, _, R, Q = setting(1)
    ref = chart.references(chart.reversible_seed(R, pi))
    tree = readers.parse_newick(FOUR)
    aln = {name: "".join(rng.choice(list(lk.AMINO), 40)) + "B-XZJ" for name in "abcd"}
    _, counts, columns = lk.compress(aln, tree)
    objectives = [nqc.Objective(aln, tree, pi, ref, shape=0.5),
                  nqc.Objective(aln, tree, pi, ref, shape=0.5, p_invar=0.1),
                  nqc.Objective.from_patterns(tree, columns, counts * 1.5, pi, ref, shape=0.5)]
    for obj in objectives:
        site = lk.site_log_likelihoods(tree, obj.columns, Q, pi, obj.rates, obj.weights, obj.invariant)
        assert obj.log_likelihood(Q) == float(obj.counts @ site)
        if obj.invariant is None:
            lnl, G = obj.log_likelihood_gradient(Q)
            lnl_c, G_c = lk.log_likelihood_gradient(tree, obj.columns, Q, pi, obj.rates, obj.weights, obj.counts)
            assert lnl == lnl_c and np.array_equal(G, G_c)


# ---------------------------------------------------------------- exact weights

def test_exact_weights_sum_to_n():
    """021's enumeration row: the sum over all patterns, here N, within 1e-12 relative."""
    _, pi1, pi2, R, Qc = setting(2)
    trees = [readers.parse_newick(t) for t in
             (THREE, "((a:0.2,c:0.4):0.1,b:0.3);", "((b:0.3,c:0.4):0.1,a:0.2);")]
    rates, rw = lk.rate_categories(shape=0.5)
    Q2 = sim.t3_transport(Qc, pi1, pi2)
    cases_ = {"plain": [sim.Component(1.0, trees[0], Qc, pi1)],
              "+G4": [sim.Component(1.0, trees[0], Qc, pi1, rates, rw)],
              "composition": [sim.Component(3.0, trees[0], Qc, pi1), sim.Component(1.0, trees[0], Q2, pi2)],
              "topology": [sim.Component(w, t, Qc, pi1) for w, t in zip((0.69, 0.155, 0.155), trees)]}
    for name, comps in cases_.items():
        w = sim.exact_weights(comps, NAMES, N)
        assert w.shape == (20 ** 3,) and np.all(w > 0), name
        assert abs(w.sum() - N) <= 1e-12 * N, name


def test_exact_weights_follow_names_not_the_tree():
    _, pi, _, _, Q = setting(3)
    ab = readers.parse_newick(THREE)
    ba = readers.parse_newick("((b:0.3,a:0.2):0.1,c:0.4);")
    w = sim.exact_weights([sim.Component(1.0, ab, Q, pi)], NAMES, N)
    assert rel(sim.exact_weights([sim.Component(1.0, ba, Q, pi)], NAMES, N), w) <= 1e-12
    assert rel(N * np.exp(lk.site_log_likelihoods(ab, sim.columns(NAMES), Q, pi)), w) <= 1e-12
    swapped = sim.exact_weights([sim.Component(1.0, ab, Q, pi)], ("b", "a", "c"), N)
    assert rel(swapped.reshape(20, 20, 20), w.reshape(20, 20, 20).transpose(1, 0, 2)) <= 1e-12
    with pytest.raises(ValueError):
        sim.exact_weights([sim.Component(1.0, ab, Q, pi)], ("a", "b", "d"), N)


# ---------------------------------------------------------------- masking

def brute_mask(weights, prob, k):
    out = {}
    for x, digits in enumerate(itertools.product(range(20), repeat=k)):
        for subset in itertools.product((False, True), repeat=k):
            factor, y = 1.0, 0
            for leaf, (d, masked) in enumerate(zip(digits, subset)):
                p = prob[x, leaf]
                factor *= p if masked else 1.0 - p
                y = 21 * y + (20 if masked else d)
            out[y] = out.get(y, 0.0) + weights[x] * factor
    return np.array([out[y] for y in range(21 ** k)])


def test_mask_against_brute_force():
    rng, pi, _, _, Q = setting(4)
    tree = readers.parse_newick("(a:0.3,b:0.5);")
    w = sim.exact_weights([sim.Component(1.0, tree, Q, pi)], ("a", "b"), N)
    for prob in (sim.residue_mask(2, rng.uniform(0.01, 0.3, 20)), sim.context_mask(2, 0.25)):
        assert rel(sim.mask(w, prob), brute_mask(w, prob, 2)) <= 1e-12


def test_uniform_masking_is_a_known_answer():
    """A residue masked at random leaves the marginal likelihood of the rest: a masked column's
    weight is N P(unmasked residues) p^g (1 - p)^(k - g), with the gaps scored as unknown."""
    _, pi, _, _, Q = setting(5)
    tree, p, k = readers.parse_newick(THREE), 0.133, 3
    w = sim.exact_weights([sim.Component(1.0, tree, Q, pi)], NAMES, N)
    masked = sim.mask(w, np.full((20 ** k, k), p))
    g = (sim.column_digits(k, 21) == sim.GAP).sum(axis=1)
    known = N * np.exp(lk.site_log_likelihoods(tree, sim.columns(NAMES, gaps=True), Q, pi)) * p ** g * (1 - p) ** (k - g)
    assert rel(masked, known) <= 1e-12


def test_masked_weights_sum_to_n_and_fractions_are_linear():
    rng, pi, _, _, Q = setting(6)
    tree, k = readers.parse_newick(THREE), 3
    w = sim.exact_weights([sim.Component(1.0, tree, Q, pi)], NAMES, N)
    g = (sim.column_digits(k, 21) == sim.GAP).sum(axis=1)
    for p in (0.05, 0.2):
        for prob in (sim.context_mask(k, p), sim.residue_mask(k, rng.uniform(0, 2 * p, 20))):
            masked = sim.mask(w, prob)
            assert abs(masked.sum() - N) <= 1e-12 * N
            assert rel(masked @ g / (k * N), sim.masked_fraction(w, prob)) <= 1e-12
        assert rel(sim.masked_fraction(w, sim.context_mask(k, p)), p * sim.masked_fraction(w, sim.context_mask(k, 1.0))) \
            <= 1e-12


def test_rarity_mask_by_hand():
    """Columns (0,0,0,1): r 2/3 for the 0s, 0 for the 1; (2,3,2,3): r 1/3; (5,5,5,5): r 1, the last bin.
    With edges (0, 0.5, 1) and p (0.4, 0.1), s = 2 masks at 0.2, 0.8 and 0.2; s = 3 caps 1.2 at 1."""
    states = np.array([[0, 2, 5], [0, 3, 5], [0, 2, 5], [1, 3, 5]])
    edges, p = [0.0, 0.5, 1.0], [0.4, 0.1]
    u = np.array([[0.15, 0.79, 0.1], [0.25, 0.81, 0.3], [0.9, 0.0, 0.19], [0.7, 0.5, 0.2]])
    assert sim.rarity_mask(states, 2.0, edges, p, u).tolist() == [[True, True, True], [False, False, False],
                                                                  [False, True, True], [True, True, False]]
    assert sim.rarity_mask(states, 3.0, edges, p, np.full((4, 3), 0.999)).tolist() == \
        [[False, True, False], [False, True, False], [False, True, False], [True, True, False]]


def test_rarity_masks_are_nested_in_strength():
    rng = np.random.default_rng(12)
    states, u = rng.integers(0, 20, (52, 300)), rng.uniform(size=(52, 300))
    edges, p = np.linspace(0, 1, 11), rng.uniform(0, 0.6, 10)
    masks = [sim.rarity_mask(states, s, edges, p, u) for s in (0.5, 1.0, 2.0)]
    for weaker, stronger in zip(masks, masks[1:]):
        assert np.all(weaker <= stronger) and weaker.sum() <= stronger.sum()
    assert 0 < masks[0].sum() < masks[2].sum()


def test_rarity_mask_refuses_bad_input():
    states, u = np.zeros((4, 3), dtype=int), np.zeros((4, 3))
    for bad in (dict(states=states + 20), dict(states=states[:1], u=u[:1]), dict(u=u[:, :2]),
                dict(edges=[0.0, 0.6, 0.5, 1.0], p=[0.1, 0.1, 0.1]), dict(edges=[0.0, 1.0]), dict(p=[0.1, -0.1])):
        kw = dict(states=states, s=1.0, edges=[0.0, 0.5, 1.0], p=[0.1, 0.1], u=u) | bad
        with pytest.raises(ValueError):
            sim.rarity_mask(**kw)


# ---------------------------------------------------------------- transports

def test_t3_transport_of_a_reversible_q_is_frequency_replacement():
    _, pi1, pi2, R, _ = setting(7)
    assert rel(sim.t3_transport(chart.reversible_seed(R, pi1), pi1, pi2), chart.reversible_seed(R, pi2)) <= 1e-10


def test_t3_transport_keeps_strengths_and_affinities():
    _, pi1, pi2, _, Qc = setting(8)
    Q2 = sim.t3_transport(Qc, pi1, pi2)
    assert np.max(np.abs(t3.pack(*t3.inverse(Q2, pi2)) - t3.pack(*t3.inverse(Qc, pi1)))) <= 1e-10
    assert rel(sim.t3_transport(Q2, pi2, pi1), Qc) <= 1e-10
    assert residuals(Q2, pi2) <= 1e-10


def test_t1_transport_keeps_the_flux():
    _, pi1, pi2, _, Qc = setting(9)
    Q2 = sim.t1_transport(Qc, pi2)
    assert residuals(Q2, pi2) <= 1e-10
    assert rel(flux.flux_of(Q2, pi2), flux.flux_of(Qc, pi1)) <= 1e-10
    assert rel(sim.circulation(Q2, pi2), sim.circulation(Qc, pi1)) <= 1e-10


def test_shift_to_tv():
    pi = sim.target()
    for name in ("LG", "Q.MAMMAL", "Q.INSECT", "WAG"):
        toward = np.asarray(builtin.exchangeabilities(name)[1], dtype=float)
        toward /= toward.sum()
        s, p = sim.shift_to_tv(pi, toward, 0.06)
        assert abs(sim.total_variation(p, pi) - 0.06) <= 1e-12
        assert np.all(p > 0) and abs(p.sum() - 1) <= 1e-14 and np.array_equal(p, sim.shift(pi, toward, s))
        s, p = sim.shift_to_tv(pi, toward, sim.total_variation(toward, pi))
        assert abs(s - 1) <= 1e-9 and rel(p, toward) <= 1e-8
        assert sim.shift_to_tv(pi, toward, 0.0)[0] == 0.0
        for tv in (-0.01, 1.0):
            with pytest.raises(ValueError):
                sim.shift_to_tv(pi, toward, tv)


# ---------------------------------------------------------------- circulation and the dial

def test_circulation():
    _, pi, _, R, _ = setting(10)
    assert sim.circulation(chart.reversible_seed(R, pi), pi) <= 1e-12
    # three states, flux 0.375 around the cycle 0 -> 1 -> 2 -> 0 and 0.125 back: C = 3 (0.25) / 1.5
    F = np.full((3, 3), 0.125)
    F[0, 1] = F[1, 2] = F[2, 0] = 0.375
    np.fill_diagonal(F, 0.0)
    pi3 = np.array([0.25, 0.25, 0.5])
    assert abs(sim.circulation(cases.generator_from_flux(F, pi3), pi3) - 0.5) <= 1e-15 * 0.5


def test_dial():
    pi, truths = sim.truths()
    Q = truths["NQ.bird"]
    S, A = sim.flux_parts(Q, pi)
    i, j, _ = t3._edges(20)
    reverse = flux.edge_index(20)[j, i]
    assert np.array_equal(S, S[reverse]) and np.array_equal(A, -A[reverse])
    assert np.max(np.abs(flux.flux_of(sim.dial(Q, pi, 1.0), pi) - flux.flux_of(Q, pi))) <= 1e-15
    assert sim.circulation(sim.dial(Q, pi, 0.0), pi) <= 1e-12
    c1 = sim.circulation(Q, pi)
    for lam in np.linspace(0, 1, 5):
        Ql = sim.dial(Q, pi, lam)
        assert residuals(Ql, pi) <= 1e-10
        if lam:
            assert rel(sim.circulation(Ql, pi), lam * c1) <= 1e-12
    with pytest.raises(chart.BuildError):
        sim.dial(Q, pi, sim.dial_limit(Q, pi) * (1 + 1e-6))


def test_transposed_is_the_time_reversal():
    pi_bird, truths = sim.truths()
    _, pi_c, _, _, Qc = setting(13)
    off = ~np.eye(20, dtype=bool)
    for Q, pi in ((truths["NQ.bird"], pi_bird), (Qc, pi_c)):
        Qt = sim.transposed(Q, pi)
        assert rel(Qt[off], (pi[None, :] * Q.T / pi[:, None])[off]) <= 1e-12
        assert residuals(Qt, pi) <= 1e-10
        assert rel(sim.circulation(Qt, pi), sim.circulation(Q, pi)) <= 1e-12
        assert rel(sim.transposed(Qt, pi)[off], Q[off]) <= 1e-12


# ---------------------------------------------------------------- recentring and mixtures

def test_recentre():
    rng, pi, _, _, _ = setting(11)
    profiles, weights = rng.dirichlet(np.ones(20), 7), rng.uniform(1, 5, 7)
    for floor in (0.0, 1e-4):
        lam = sim.recentre_limit(pi, profiles, weights, floor)
        out = sim.recentre(pi, profiles, weights, lam)
        assert np.max(np.abs(weights / weights.sum() @ out - pi)) <= 1e-14
        assert abs(out.min() - floor) <= 1e-14


def test_mean_flux_of_transported_profiles_is_balanced():
    pi, truths = sim.truths()
    profiles, weights = builtin.frequency_profiles("C20")
    lam = sim.recentre_limit(pi, profiles, weights) / 2
    comps = [sim.Component(w, None, sim.t3_transport(truths["NQ.bird"], pi, p), p)
             for w, p in zip(weights, sim.recentre(pi, profiles, weights, lam))]
    assert residuals(flux.generator(sim.mean_flux(comps), pi), pi) <= 1e-10


def test_c20_reader():
    profiles, weights = builtin.frequency_profiles("C20")
    assert profiles.shape == (20, 20) and weights.shape == (20,)
    # C20pi16's N, and the FMIX weights of C20pi1 and C20pi20, at the source's 10 printed decimals
    assert rel(profiles[15, builtin.STATES.index("N")], 0.2520448083) <= 1e-9
    assert rel(weights[[0, 19]], [0.0559910600, 0.0596471901]) <= 1e-9
    assert rel(profiles.sum(axis=1), np.ones(20)) <= 1e-9 and abs(weights.sum() - 1) <= 1e-9


# ---------------------------------------------------------------- planning numbers

def test_planning_numbers():
    """The SIM-1 plan's planning measurements (made in memory, docs/agent/archive/SIM1_plan_2026-11-16.md,
    "Measured while planning"), reproduced to their printed digits."""
    pi, truths = sim.truths()
    Q = truths["NQ.bird"]
    pi_lg = builtin.exchangeabilities("LG")[1]
    assert f"{np.max(np.abs(builtin.full_matrix('NQ.BIRD')[1] - pi)):.1e}" == "6.2e-07"
    assert round(sim.circulation(Q, pi), 4) == 0.0577
    assert round(sim.dial_limit(Q, pi), 3) == 1.444
    assert round(sim.total_variation(pi, sim.shift(pi, pi_lg, 1.0)), 3) == 0.062
    assert round(sim.total_variation(pi, sim.shift(pi, pi_lg, 0.25)), 3) == 0.016
    assert round(sim.total_variation(pi, builtin.exchangeabilities("Q.MAMMAL")[1]), 3) == 0.020
    for s, c in ((0.25, 0.0565), (1.0, 0.0527)):
        p = sim.shift(pi, pi_lg, s)
        assert round(sim.circulation(sim.t3_transport(Q, pi, p), p), 4) == c
    # C20's positivity limit 0.655, and at half of it the smallest entry and the weight-averaged TV
    profiles, weights = builtin.frequency_profiles("C20")
    limit = sim.recentre_limit(pi, profiles, weights)
    assert round(limit, 3) == 0.655
    half = sim.recentre(pi, profiles, weights, limit / 2)
    assert round(half.min(), 3) == 0.010
    assert round(float(weights @ [sim.total_variation(p, pi) for p in half]), 2) == 0.16
