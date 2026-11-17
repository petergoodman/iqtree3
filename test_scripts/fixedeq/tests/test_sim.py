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
from oracle import cases, chart, flux, likelihood as lk, nqc, readers, simulate  # noqa: E402
import alisim_check as ac  # noqa: E402
import common  # noqa: E402
import exact as ex  # noqa: E402


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


# ---------------------------------------------------------------- exact.py (commit 4)

def test_log_rate_jacobian_against_central_differences():
    """021's builder line: <= 1e-7 relative to max(1, max |column|)."""
    pi = common.target_values()
    ref = nqc.lg_seed(pi)[1]
    z = chart.inverse(simulate.truths()[1]["NQ.bird"], ref)
    J = ex.log_rate_jacobian(z, pi, ref)
    off = ~np.eye(20, dtype=bool)
    for k in np.random.default_rng(41).choice(len(z), 6, replace=False):
        h = 1e-5 * max(1.0, abs(z[k]))
        up, down = z.copy(), z.copy()
        up[k] += h
        down[k] -= h
        d = (np.log(chart.build(up, pi, ref)[0][off]) - np.log(chart.build(down, pi, ref)[0][off])) / (up[k] - down[k])
        assert np.max(np.abs(J[:, k] - d)) <= 1e-7 * max(1.0, np.max(np.abs(J[:, k])))


def test_central_hessian_of_a_quadratic():
    rng = np.random.default_rng(42)
    M = rng.normal(size=(6, 6))
    A, b = M @ M.T + np.eye(6), rng.normal(size=6)
    H, taken = ex.central_hessian(lambda zs: [(0.5 * z @ A @ z + b @ z, A @ z + b) for z in zs], rng.normal(size=6))
    assert taken == 12 and np.max(np.abs(H - A)) <= 1e-8 * np.max(np.abs(A))


def test_truths_inside_the_box():
    pi = common.target_values()
    ref = nqc.lg_seed(pi)[1]
    dom = flux.matched(ref)
    for p in ex.problems().values():
        z = chart.inverse(p.truth, ref)
        assert nqc.LOWER < z.min() and z.max() < nqc.UPPER and dom.contains(flux.flux_of(p.truth, pi)), p.name


def test_r6_masked_weights_sum_to_n():
    pi = common.target_values()
    obj = ex.objective(ex.problems()["R6"], pi, nqc.lg_seed(pi)[1])
    assert len(obj.counts) == 21 ** 3 and abs(obj.counts.sum() - ex.SITES) <= 1e-12 * ex.SITES


def test_predicted_errors_by_hand():
    """H = diag(2, 4), J = [[1, 1], [0, 2]]: J H^-1 = [[0.5, 0.25], [0, 0.5]], so p = 1e-6 [0.75, 0.5]."""
    p = ex.predicted_errors(np.diag([2.0, 4.0]), np.array([[1.0, 1.0], [0.0, 2.0]]))
    assert np.max(np.abs(p - [0.75e-6, 0.5e-6]) / [0.75e-6, 0.5e-6]) <= 1e-15


def four_state_problem():
    """V2's constructed problem at n = 4: exact weights from a cycle-flux truth on three taxa."""
    rng = np.random.default_rng(204)
    pi = cases.random_target(4, rng)
    start = chart.reversible_seed(cases.random_exchangeabilities(4, rng), pi)
    ref = chart.references(start)
    truth = cases.generator_from_flux(cases.cycle_flux(4, rng), pi)
    tree = readers.parse_newick(ex.THREE)
    E = np.eye(4)
    columns = [{"a": E[x], "b": E[y], "c": E[w]} for x in range(4) for y in range(4) for w in range(4)]
    weights = ex.SITES * np.exp(lk.site_log_likelihoods(tree, columns, truth, pi))
    return nqc.Objective.from_patterns(tree, columns, weights, pi, ref), truth, start, pi, ref


def test_recovery_end_to_end_at_four_states():
    """fit_and_check at n = 4: every asserted row passes, and the before-fit record holds every rate."""
    obj, truth, start, pi, ref = four_state_problem()
    seen = []
    rows, record = ex.fit_and_check("n4", obj, truth, start, pi, ref, 600.0, 1, log=lambda s: None, before=seen.append)
    asserted = [row for row in rows if row["passed"] is not None]
    assert len(asserted) == 6 and all(row["passed"] for row in asserted), [r for r in asserted if not r["passed"]]
    assert len(seen) == 1 and len(seen[0]["predicted_errors"]) == len(seen[0]["rate_lines"]) == 12


def test_recovery_fails_against_a_moved_truth():
    """The same data checked against a truth one coordinate (0.01 in z) away: the rate row must fail."""
    obj, truth, start, pi, ref = four_state_problem()
    z = chart.inverse(truth, ref)
    z[0] += 0.01
    moved = chart.build(z, pi, ref)[0]
    rows, _ = ex.fit_and_check("n4-moved", obj, moved, start, pi, ref, 600.0, 1, log=lambda s: None)
    rate = [row for row in rows if "every rate" in row["check"]]
    assert len(rate) == 1 and rate[0]["passed"] is False


# ---------------------------------------------------------------- exact.py screen (commit 5)

def test_repeat_selection():
    def fit(s, base, dc, flag, evaluable=True):
        return {"scenario": s, "base": base, "dC_reference": dc, "flag": flag, "evaluable": evaluable}
    fits = [fit("S1a", "Q.bird", 0.001, False), fit("S3", "Q.bird", -0.02, True), fit("S4a", "Q.bird", 0.01, True),
            fit("S5", "Q.bird", 0.03, True, evaluable=False), fit("S1b", "Q.bird", 0.008, True),
            fit("S4b", "NQ.bird", 0.5, True)]
    assert ex.repeat_selection(fits) == ["S3", "S4a"]
    quiet = [fit("S1a", "Q.bird", 0.001, False), fit("S4b", "Q.bird", -0.004, False), fit("S3", "NQ.bird", 0.9, True)]
    assert ex.repeat_selection(quiet) == ["S4b"]
    assert ex.repeat_selection([]) == []


def test_s4b_masks_the_target_fraction():
    pi, truths = simulate.truths()
    _, gaps, _, info = ex.screen_data("S4b", truths["Q.bird"], pi, ex.THREE)
    assert gaps and abs(info["masked_fraction"] - ex.MASK_TARGET) <= 1e-12 * ex.MASK_TARGET


def test_s2_recentring_respects_the_floor(tmp_path):
    pi, truths = simulate.truths()
    counts = np.ones((3, 20))
    counts[0, 4] = 0.0      # an absent residue
    counts[1, 0] = 400.0    # a dominant one
    table = "locus\tlength\t" + "\t".join(lk.AMINO) + "\n" + "".join(
        f"L{k}\t{n}\t" + "\t".join(f"{v:g}" for v in row) + "\n" for k, (n, row) in enumerate(zip((300, 500, 900), counts)))
    path = tmp_path / "loci.tsv"
    path.write_text(table, encoding="utf-8")
    loci = ex.read_loci(path)
    assert np.array_equal(loci[0], [300.0, 500.0, 900.0]) and np.array_equal(loci[1], counts)
    w, gaps, (Qr, pr), info = ex.screen_data("S2", truths["NQ.bird"], pi, ex.THREE, loci)
    assert not gaps and info["lambda"] <= 1.0 and info["smallest_entry"] >= ex.LOCUS_FLOOR - 1e-15
    assert abs(w.sum() - ex.SITES) <= 1e-12 * ex.SITES
    d = chart.diagnostics(Qr, pr)
    assert max(d["stationarity"], d["row_sums"], d["mean_rate_error"]) <= 1e-10


def test_s3_reference_is_a_generator_at_pi():
    pi, truths = simulate.truths()
    _, _, (Qr, pr), info = ex.screen_data("S3", truths["Q.bird"], pi, ex.THREE)
    d = chart.diagnostics(Qr, pr)
    assert np.array_equal(pr, pi) and max(d["stationarity"], d["row_sums"], d["mean_rate_error"]) <= 1e-10
    assert info["components"] == 20 and info["smallest_entry"] > 0


def test_unmasked_composition_is_pi():
    """At a stationary root every tip's marginal is pi, so the gap-free composition of exact weights is pi."""
    pi, truths = simulate.truths()
    tree = readers.parse_newick(ex.THREE)
    w = simulate.exact_weights([simulate.Component(1.0, tree, truths["NQ.bird"], pi)], ("a", "b", "c"), ex.SITES)
    assert np.max(np.abs(ex.composition(simulate.mask(w, np.zeros((len(w), 3))), 3) - pi)) <= 1e-12


def test_s5_topologies_keep_each_taxons_length():
    lengths = {}
    for text in ex.TOPOLOGIES:
        tree = readers.parse_newick(text)
        lengths[text] = {leaf.name: leaf.length for leaf in tree.leaves()}
        assert sorted(lengths[text]) == ["a", "b", "c"]
    assert all(v == {"a": 0.2, "b": 0.3, "c": 0.4} for v in lengths.values())
    assert abs(sum(ex.DISCORDANCE) - 1) <= 1e-15


def test_repeat_selection_reads_the_higher_endpoint():
    first = {"scenario": "S3", "base": "Q.bird", "evaluable": True, "flag": True, "dC_reference": 0.02, "lnl": -10.0,
             "from_reference": {"evaluable": True, "lnl": -9.0, "flag": False, "dC_reference": 0.001}}
    other = {"scenario": "S4a", "base": "Q.bird", "evaluable": True, "flag": True, "dC_reference": 0.01, "lnl": -5.0}
    assert ex.best_endpoint(first) is first["from_reference"]
    assert ex.repeat_selection([first, other]) == ["S4a"]


def test_s2_reaches_the_floor_when_lambda_is_below_one(tmp_path):
    pi, truths = simulate.truths()
    counts = np.ones((2, 20))
    counts[0, 4] = 0.0
    counts[1, 0] = 5000.0
    path = tmp_path / "loci.tsv"
    path.write_text("# content: constructed\nlocus\tlength\t" + "\t".join(lk.AMINO) + "\n"
                    + "".join(f"L{k}\t{n}\t" + "\t".join(f"{v:g}" for v in row) + "\n"
                              for k, (n, row) in enumerate(zip((300, 900), counts))), encoding="utf-8")
    _, _, _, info = ex.screen_data("S2", truths["NQ.bird"], pi, ex.THREE, ex.read_loci(path))
    assert info["lambda"] < 1 and abs(info["smallest_entry"] - ex.LOCUS_FLOOR) <= 1e-15


def test_root_on_a_leaf():
    import birddata as bd
    rooted = bd.root_on("((a:1,b:2)0.9:0.5,(c:1,(d:1,e:1)1:0.2):0.3);", "d")
    tree = readers.parse_newick(rooted)
    assert [c.name for c in tree.children][0] == "d" and len(tree.children) == 2
    assert sorted(leaf.name for leaf in tree.leaves()) == ["a", "b", "c", "d", "e"]
    assert all(len(n.children) in (0, 2) for n in tree.preorder())
    def clades(text):
        out = set()

        def rec(nd):
            leaves = frozenset([nd.name]) if not nd.children else frozenset().union(*(rec(c) for c in nd.children))
            out.add(leaves)
            return leaves
        rec(readers.parse_newick(text))
        return out
    assert clades(rooted) == clades("(d,(e,(c,(a,b))));")
