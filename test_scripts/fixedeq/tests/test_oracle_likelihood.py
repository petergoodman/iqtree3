"""Rooted likelihood, readers and gamma rates; pass lines are decision 021's."""
import itertools
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy.linalg import expm
from scipy.special import gammainc, gammaincinv

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from oracle import cases, chart, gamma, likelihood as lk, readers  # noqa: E402

REPO = Path(__file__).resolve().parents[3]


def rel(a, b):
    return abs(a - b) / abs(b)


def test_phylip_fasta_and_nexus_readers():
    aln = readers.read_phylip(REPO / "example" / "aa_example.phy")
    assert len(aln) == 18 and all(len(s) == 697 for s in aln.values())
    assert aln["Human"].startswith("MVEMLPTAILLVLAVSVVAK")
    turtle = readers.read_fasta(REPO / "test_scripts" / "test_data" / "turtle_aa.fasta")
    assert len(turtle) == 16 and all(len(s) == 774 for s in turtle.values())
    sets = readers.read_charsets(REPO / "test_scripts" / "test_data" / "turtle_aa.nex")
    assert [(min(v), max(v), len(v)) for v in sets.values()] == [(0, 171, 172), (172, 387, 216), (388, 773, 386)]


def test_newick_reader_keeps_the_root_and_lengths():
    tree = readers.parse_newick("((A:0.1,B:0.2):0.05,(C:0.3,D:0.4)lab:0.06);")
    assert [c.length for c in tree.children] == [0.05, 0.06] and tree.length is None
    assert sorted(leaf.name for leaf in tree.leaves()) == ["A", "B", "C", "D"]


def test_characters_map_as_iqtree_maps_them():
    v = lk.tip_vector
    assert np.array_equal(np.nonzero(v("B"))[0], [2, 3]) and np.array_equal(np.nonzero(v("Z"))[0], [5, 6])
    assert np.array_equal(np.nonzero(v("j"))[0], [9, 10]) and v("A")[0] == 1 and v("A").sum() == 1
    for c in "X?-.~!*UO":
        assert v(c).sum() == 20
    with pytest.raises(ValueError):
        v("8")


def test_file_matrix_convention_solves_pi_and_normalizes():
    rng = np.random.default_rng(10)
    pi = cases.random_target(20, rng)
    Q = 3.7 * cases.generator_from_flux(cases.cycle_flux(20, rng), pi)
    Qn, root = lk.file_model(Q)
    assert rel(float(root @ -np.diag(Qn)), 1.0) <= 1e-12
    assert np.max(np.abs(root - pi)) <= 1e-12


def brute_force(tree, columns, Q, pi):
    """Sum over every assignment of internal states; independent of the pruning code."""
    nodes = [nd for nd in tree.preorder() if nd.children]
    n = len(pi)
    total = []
    for col in columns:
        site = 0.0
        for assign in itertools.product(range(n), repeat=len(nodes)):
            state = dict(zip(map(id, nodes), assign))
            p = pi[state[id(tree)]]
            for nd in nodes:
                for ch in nd.children:
                    P = expm(ch.length * Q)[state[id(nd)]]
                    p *= P[state[id(ch)]] if ch.children else P @ col[ch.name]
            site += p
        total.append(site)
    return np.array(total)


def test_pruning_matches_brute_force_enumeration():
    rng = np.random.default_rng(11)
    n = 4
    pi = cases.random_target(n, rng)
    Q = cases.generator_from_flux(cases.cycle_flux(n, rng), pi)
    tree = readers.parse_newick("((A:0.3,B:0.7):0.2,(C:0.5,D:0.1):0.4);")
    columns = [{name: np.eye(n)[rng.integers(n)] for name in "ABCD"} for _ in range(6)]
    columns.append({"A": np.ones(n), "B": np.eye(n)[1] + np.eye(n)[2], "C": np.eye(n)[0], "D": np.ones(n)})
    site = lk.site_log_likelihoods(tree, columns, Q, pi)
    assert np.max(np.abs(np.exp(site) - brute_force(tree, columns, Q, pi)) / brute_force(tree, columns, Q, pi)) <= 1e-12


def test_probabilities_of_all_patterns_sum_to_one_with_gamma_and_invariant_sites():
    rng = np.random.default_rng(12)
    n = 3
    pi = cases.random_target(n, rng)
    Q = cases.generator_from_flux(cases.cycle_flux(n, rng), pi)
    tree = readers.parse_newick("((A:0.3,B:0.7):0.2,C:0.5);")
    rates, weights = gamma.category_rates(0.7, 4, p_invar=0.2)
    columns = [{"A": np.eye(n)[a], "B": np.eye(n)[b], "C": np.eye(n)[c]}
               for a, b, c in itertools.product(range(n), repeat=3)]
    invar = [0.2 * pi[a] if a == b == c else 0.0 for a, b, c in itertools.product(range(n), repeat=3)]
    site = lk.site_log_likelihoods(tree, columns, Q, pi, rates, weights, np.array(invar))
    assert rel(np.exp(site).sum(), 1.0) <= 1e-12


def test_reversible_likelihood_does_not_depend_on_the_root():
    rng = np.random.default_rng(13)
    pi = cases.random_target(20, rng)
    Q = chart.reversible_seed(cases.random_exchangeabilities(20, rng), pi)
    columns = [{name: np.eye(20)[rng.integers(20)] for name in "ABCD"} for _ in range(5)]
    one = readers.parse_newick("((A:0.3,B:0.7):0.2,(C:0.5,D:0.1):0.4);")
    two = readers.parse_newick("(A:0.1,(B:0.7,(C:0.5,D:0.1):0.6):0.2);")
    a, b = lk.site_log_likelihoods(one, columns, Q, pi), lk.site_log_likelihoods(two, columns, Q, pi)
    assert np.max(np.abs(a - b) / np.abs(b)) <= 1e-10


def test_an_ambiguous_site_is_the_sum_of_its_resolutions():
    rng = np.random.default_rng(14)
    pi = cases.random_target(20, rng)
    Q = cases.generator_from_flux(cases.cycle_flux(20, rng), pi)
    tree = readers.parse_newick("((A:0.3,B:0.7):0.2,C:0.5);")
    col = lambda c: {"A": lk.tip_vector("K"), "B": lk.tip_vector(c), "C": lk.tip_vector("W")}
    b, n_, d = np.exp(lk.site_log_likelihoods(tree, [col("B"), col("N"), col("D")], Q, pi))
    assert rel(n_ + d, b) <= 1e-12


def test_invariant_term_follows_the_constant_state_set():
    pi = np.arange(1, 21) / 210.0
    assert lk.invariant_term("AA-A", pi, 0.3) == 0.3 * pi[0]
    assert lk.invariant_term("NBD?", pi, 0.3) == 0.0
    assert lk.invariant_term("NB?", pi, 0.3) == 0.3 * pi[2]
    assert lk.invariant_term("BB-", pi, 0.3) == (pi[2] + pi[3]) * 0.3
    assert lk.invariant_term("--X", pi, 0.3) == 0.3
    assert lk.invariant_term("AR", pi, 0.3) == 0.0


@pytest.mark.parametrize("shape", [0.05, 0.3, 1.0, 3.0, 20.0])
def test_ported_gamma_rates_agree_with_scipy(shape):
    k = 4
    cut = gammaincinv(shape, np.arange(1, k) / k) / shape
    freq = gammainc(shape + 1, cut * shape)
    exact = k * np.diff(np.r_[0.0, freq, 1.0])
    rates, weights = gamma.category_rates(shape, k)
    assert np.max(np.abs(rates - exact) / exact) <= 1e-5
    assert np.all(weights == 0.25)
    inv_rates, inv_weights = gamma.category_rates(shape, k, p_invar=0.25)
    assert np.max(np.abs(inv_rates - exact / 0.75) / (exact / 0.75)) <= 1e-5
    assert np.all(inv_weights == 0.75 / 4)
