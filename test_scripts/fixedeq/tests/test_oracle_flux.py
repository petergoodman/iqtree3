"""S5's direct balanced-flux reference (oracle/flux.py; decision 042; CODE_PLAN.md 3.1 and 3.3), with
decision 021's S5 lines: A0 to A3, AG, and V1, V1b and V2 at n = 3 and 4 (n = 20 runs in the
driver's validate step)."""
import itertools
import math
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from oracle import cases, chart, flux, likelihood as lk, nqc, optimize as opt, readers  # noqa: E402

TREE = "((a:0.1,b:0.2):0.05,(c:0.3,d:0.1):0.2);"
THREE = "((a:0.2,b:0.3):0.1,c:0.4);"


def setting(n, rng):
    """A target, a reversible start of random exchangeabilities and its row-maximum references."""
    pi = cases.random_target(n, rng)
    Qs = chart.reversible_seed(cases.random_exchangeabilities(n, rng), pi)
    return pi, Qs, chart.references(Qs)


def problem(n, seed, columns=30, shape=0.5):
    rng = np.random.default_rng(seed)
    pi, Qs, ref = setting(n, rng)
    E = np.eye(n)
    cols = [{name: E[rng.integers(n)] for name in "abcd"} for _ in range(columns)]
    cols[0]["a"] = np.ones(n)
    obj = nqc.Objective.from_patterns(readers.parse_newick(TREE), cols, rng.integers(1, 6, columns), pi, ref,
                                      shape=shape)
    return pi, Qs, ref, obj


def kkt_problem(n, rng, n_lo, n_hi):
    """V1b: a known optimum f* with n_lo coordinates at the lower bound and n_hi at the upper, and
    N = f* (C^T eta - A^T nu) for eta and nu > 0 on exactly those rows, so that f* satisfies the KKT
    conditions of the strictly concave problem and is its unique maximizer."""
    pi, Qs, ref = setting(n, rng)
    dom = flux.matched(ref)
    z = flux.z_of_flux(flux.flux_of(cases.generator_from_flux(cases.cycle_flux(n, rng), pi), pi), ref)
    m = len(z)
    pick = rng.choice(m, n_lo + n_hi, replace=False)
    z[pick[:n_lo]], z[pick[n_lo:]] = dom.lower, dom.upper
    fstar = flux.flux_of(chart.build(z, pi, ref)[0], pi)
    eta = np.r_[rng.normal(0, 50, n - 1), 1000.0]
    nu = np.zeros(2 * m)
    nu[pick[:n_lo]] = rng.uniform(100, 300, n_lo)
    nu[m + pick[n_lo:]] = rng.uniform(1e-3, 3e-3, n_hi)
    w = flux.geometry(n)[2].T @ eta - dom.A.T @ nu
    assert np.all(w > 0)
    return pi, Qs, dom, fstar, fstar * w


@pytest.mark.parametrize("n", [3, 4, 20])
def test_a0_rank_and_kernel(n):
    assert np.linalg.matrix_rank(flux.geometry(n)[2]) == n
    assert flux.canonical(n)[1].shape[1] == n * (n - 2)


@pytest.mark.parametrize("n", [3, 4, 20])
def test_a1_round_trips(n):
    rng = np.random.default_rng(10 + n)
    pi = cases.random_target(n, rng)
    ref = chart.references(cases.generator_from_flux(cases.cycle_flux(n, rng), pi))
    dom = flux.matched(ref)
    f0, Z = flux.canonical(n)
    f = flux.flux_of(cases.generator_from_flux(cases.cycle_flux(n, rng), pi), pi)
    assert np.max(np.abs(f0 + Z @ (Z.T @ (f - f0)) - f)) <= 1e-14
    z = rng.uniform(dom.lower, dom.upper, n * (n - 2))
    fz = flux.flux_of(chart.build(z, pi, ref)[0], pi)
    assert np.max(np.abs(flux.z_of_flux(fz, ref) - z)) <= 1e-10
    basis = flux.Basis(fz, dom)
    Q = flux.generator(basis.flux(0.01 * rng.normal(size=n * (n - 2))), pi)
    back = flux.generator(basis.flux(basis.coordinates(flux.flux_of(Q, pi))), pi)
    assert np.max(np.abs(back - Q) / np.abs(Q)) <= 1e-10


@pytest.mark.parametrize("n", [3, 4, 20])
def test_a2_the_box_and_the_polytope_are_one_set(n):
    rng = np.random.default_rng(20 + n)
    pi = cases.random_target(n, rng)
    ref = chart.references(cases.generator_from_flux(cases.cycle_flux(n, rng), pi))
    dom = flux.matched(ref)
    m = n * (n - 2)
    points = []
    for k in range(40):
        z = rng.choice([dom.lower, dom.upper], m) if k < 10 else rng.uniform(dom.lower, dom.upper, m)
        f = flux.flux_of(chart.build(z, pi, ref)[0], pi)
        assert np.min(dom.slack(f)) >= -1e-12
        points.append(f)
    for _ in range(20):
        z = flux.z_of_flux(rng.dirichlet(np.full(len(points), 0.2)) @ np.array(points), ref)
        assert dom.lower - 1e-10 <= z.min() and z.max() <= dom.upper + 1e-10


@pytest.mark.parametrize("n", [4, 20])
def test_z_gradient_is_the_adjoint_of_chart_derivative(n):
    rng = np.random.default_rng(n)
    pi = cases.random_target(n, rng)
    Q = cases.generator_from_flux(cases.cycle_flux(n, rng), pi)
    ref = chart.references(Q)
    z, G = chart.inverse(Q, ref), rng.normal(size=(n, n))
    e, loop = np.zeros(len(z)), np.empty(len(z))
    for k in range(len(z)):
        e[k] = 1.0
        loop[k] = np.sum(G * chart.derivative(z, pi, ref, e))
        e[k] = 0.0
    assert np.max(np.abs(nqc.z_gradient(z, pi, ref, G) - loop)) <= 1e-10 * np.max(np.abs(loop))


def test_likelihood_gradient_value_and_refusals():
    pi, Qs, ref, obj = problem(4, 7)
    lnL, _ = obj.log_likelihood_gradient(Qs)
    assert abs(lnL - obj.log_likelihood(Qs)) <= 1e-12 * abs(lnL)
    rooted = readers.parse_newick("((a:0.1,b:0.2):0.05,(c:0.3,d:0.1):0.2):0.1;")
    with pytest.raises(ValueError):
        lk.log_likelihood_gradient(rooted, obj.columns, Qs, pi)
    rng = np.random.default_rng(3)
    pi20, Q20, ref20 = setting(20, rng)
    aln = {name: "".join(rng.choice(list(lk.AMINO), 20)) for name in "abcd"}
    with pytest.raises(ValueError):
        nqc.Objective(aln, readers.parse_newick(TREE), pi20, ref20, p_invar=0.2).log_likelihood_gradient(Q20)


@pytest.mark.parametrize("n", [4, 20])
def test_ag_analytic_gradients_in_u_and_z(n):
    pi, Qs, ref, obj = problem(n, 30 + n)
    rng = np.random.default_rng(n)
    z = chart.inverse(Qs, ref) + 0.3 * rng.normal(size=n * (n - 2))
    f = flux.flux_of(chart.build(z, pi, ref)[0], pi)
    basis = flux.Basis(f, flux.matched(ref))
    lf = flux.LikelihoodFlux(obj, pi)
    dirs = list(np.eye(len(z))[rng.choice(len(z), 5, replace=False)])
    dirs += [d / np.linalg.norm(d) for d in rng.normal(size=(5, len(z)))]
    checks = flux.gradient_check(lambda u: -lf(basis.flux(u))[0], np.zeros(len(z)), -basis.chain(lf(f)[1]), dirs)
    checks += flux.gradient_check(obj, z, obj.gradient(z)[1], dirs)
    assert all(err <= line for err, line in checks), checks


def test_from_patterns_matches_the_alignment_constructor():
    rng = np.random.default_rng(5)
    pi, Q, ref = setting(20, rng)
    tree = readers.parse_newick(TREE)
    aln = {name: "".join(rng.choice(list(lk.AMINO), 25)) + "B-X" for name in "abcd"}
    _, counts, columns = lk.compress(aln, tree)
    a = nqc.Objective(aln, tree, pi, ref, shape=0.5)
    b = nqc.Objective.from_patterns(tree, columns, counts, pi, ref, shape=0.5)
    assert a.log_likelihood(Q) == b.log_likelihood(Q)
    with pytest.raises(ValueError):
        nqc.Objective.from_patterns(tree, columns, counts, pi, ref, p_invar=0.2)
    with pytest.raises(ValueError):
        nqc.Objective.from_patterns(tree, [{"a": np.ones(20)}], [1.0], pi, ref)


def test_infeasible_fluxes_are_refused_before_any_likelihood():
    class Counting:
        calls = 0

        def log_likelihood_gradient(self, Q):
            Counting.calls += 1
            return 0.0, np.zeros_like(Q)

    rng = np.random.default_rng(9)
    pi = cases.random_target(4, rng)
    good = flux.flux_of(cases.generator_from_flux(cases.cycle_flux(4, rng), pi), pi)
    lf = flux.LikelihoodFlux(Counting(), pi)
    for bad in (0.0, -1e-300, math.nan, math.inf):
        f = good.copy()
        f[3] = bad
        with pytest.raises(chart.BuildError):
            lf(f)
    assert Counting.calls == 0


def test_difference_hessian_against_the_exact_hessian():
    """Both branches against the concave objective's exact Hessian in u, J^T diag(N / f^2) J with
    J = diag(f_c) N_c: forward differences err O(h), central ones O(h^2) plus rounding, h = 1e-6."""
    rng = np.random.default_rng(5)
    pi, Qs, ref = setting(4, rng)
    dom = flux.matched(ref)
    basis = flux.Basis(flux.inward(Qs, pi, dom)[0], dom)
    N = rng.lognormal(0, 1, 12)
    objective = flux.Concave(N)
    u = rng.normal(0, 1e-2, basis.N.shape[1])
    f = basis.flux(u)
    J = basis.fc[:, None] * basis.N
    exact = J.T @ (J * (N / f ** 2)[:, None])
    g0 = -basis.chain(objective(f)[1])
    for central, line in ((False, 1e-4), (True, 1e-8)):
        H, count = flux.difference_hessian(lambda fs: [objective(p) for p in fs], basis, u, g0, central)
        assert count == (2 if central else 1) * len(u)
        assert np.max(np.abs(H - exact)) <= line * np.max(np.abs(exact))


@pytest.mark.parametrize("n", [3, 4])
def test_v1_open_concave_problem(n):
    rng = np.random.default_rng(100 + n)
    N = rng.lognormal(0, 1, n * (n - 1))
    pi, Qs, ref = setting(n, rng)
    dom = flux.open_domain(ref)
    seen = []
    objective = flux.Concave(N)
    r = flux.fit(lambda f: seen.append(np.all(f > 0)) or objective(f), flux.inward(Qs, pi, dom)[0], dom)
    assert all(seen)
    assert r.status == "converged", r.message
    assert flux.duality_gap(r.f, N, dom, r.mu) <= 1e-9 * abs(r.value)
    assert abs(float(N @ np.log(flux.dual_newton(N, n))) - r.value) <= 1e-9 * abs(r.value)
    assert max(flux.residuals(r.f, pi).values()) <= 1e-10


def test_v1b_binding_ratio_constraints_at_a_known_optimum():
    pi, Qs, dom, fstar, N = kkt_problem(4, np.random.default_rng(11), 2, 1)
    r = flux.fit(flux.Concave(N), flux.inward(Qs, pi, dom)[0], dom)
    assert r.status == "converged", r.message
    assert flux.duality_gap(r.f, N, dom, r.mu) <= 1e-9 * abs(r.value)
    assert float(N @ np.log(fstar)) - r.value <= 1e-9 * abs(r.value)
    assert r.m_active >= 3
    assert max(flux.residuals(r.f, pi).values()) <= 1e-10


def test_certificate_passes_a_known_optimum_and_fails_elsewhere():
    rng = np.random.default_rng(11)
    pi, Qs, dom, fstar, N = kkt_problem(4, rng, 2, 1)
    F = float(N @ np.log(fstar))
    assert flux.certificate(fstar, N / fstar, dom)["residual"] <= 1e-9 * abs(F)
    other = kkt_problem(4, rng, 2, 1)[4]
    assert flux.certificate(fstar, other / fstar, dom)["residual"] >= 1.0


def test_v2_with_rate_categories_meets_the_declared_stop():
    """V2's n = 3 problem with +G4{0.5} in data and fit, the setting of SIM-1's R4: the active-set
    finish once took its 'stationary' exit, on a Newton decrement below 1e-6 times its rounding scale,
    before the declared stop held (explore/r4_finish_diagnosis.py). It must end at the declared stop,
    with the certificate (decision 039) and the truth's lnL reached."""
    n = 3
    rng = np.random.default_rng(200 + n)
    pi, Qs, ref = setting(n, rng)
    Qt = cases.generator_from_flux(cases.cycle_flux(n, rng), pi)
    tree = readers.parse_newick(THREE)
    E = np.eye(n)
    cols = [{"a": E[x], "b": E[y], "c": E[w]} for x, y, w in itertools.product(range(n), repeat=3)]
    rates, weights = lk.rate_categories(shape=0.5)
    counts = 1e6 * np.exp(lk.site_log_likelihoods(tree, cols, Qt, pi, rates, weights))
    obj = nqc.Objective.from_patterns(tree, cols, counts, pi, ref, shape=0.5)
    dom = flux.matched(ref)
    r = flux.fit(flux.LikelihoodFlux(obj, pi), flux.inward(Qs, pi, dom)[0], dom)
    assert r.status == "converged", r.message
    value, gf = flux.LikelihoodFlux(obj, pi)(r.f)
    assert flux.certificate(r.f, gf, dom)["residual"] <= 1e-9 * max(1.0, abs(value))
    assert obj.log_likelihood(Qt) - obj.log_likelihood(flux.generator(r.f, pi)) <= 1e-3


@pytest.mark.parametrize("n", [3, 4])
def test_v2_exact_frequencies_recover_the_truth(n):
    rng = np.random.default_rng(200 + n)
    pi, Qs, ref = setting(n, rng)
    Qt = cases.generator_from_flux(cases.cycle_flux(n, rng), pi)
    tree = readers.parse_newick(THREE)
    E = np.eye(n)
    cols = [{"a": E[x], "b": E[y], "c": E[w]} for x, y, w in itertools.product(range(n), repeat=3)]
    weights = 1e6 * np.exp(lk.site_log_likelihoods(tree, cols, Qt, pi))
    obj = nqc.Objective.from_patterns(tree, cols, weights, pi, ref)
    dom = flux.matched(ref)
    assert dom.contains(flux.flux_of(Qt, pi))
    r = flux.fit(flux.LikelihoodFlux(obj, pi), flux.inward(Qs, pi, dom)[0], dom)
    assert r.status == "converged", r.message
    assert obj.log_likelihood(Qt) - obj.log_likelihood(flux.generator(r.f, pi)) <= 1e-3
    assert max(flux.residuals(r.f, pi).values()) <= 1e-10


def test_matched_fit_reaches_the_peak_of_lbfgsb_at_n4():
    pi, Qs, ref, obj = problem(4, 11, columns=200, shape=None)
    dom = flux.matched(ref)
    r = flux.fit(flux.LikelihoodFlux(obj, pi), flux.inward(Qs, pi, dom)[0], dom)
    climb = opt.lbfgsb_fit(obj, chart.inverse(Qs, ref), nqc.LOWER, nqc.UPPER, pgtol=1e-6, max_iter=50_000,
                           gradient=obj.gradient)
    assert r.status == "converged", r.message
    assert climb.analytic_pg <= 1e-6
    assert abs(r.value + climb.f) <= 1e-3


def test_lbfgsb_analytic_pass_on_a_clipped_quadratic():
    a = np.array([1.0, 10.0, 100.0, 3.0])
    m = np.array([0.5, -2.0, 0.25, 3.0])
    lower, upper = np.full(4, -1.0), np.full(4, 1.0)

    def f(x):
        return float(0.5 * a @ (x - m) ** 2)

    climb = opt.lbfgsb_fit(f, np.zeros(4), lower, upper, pgtol=1e-6, maxcor=5,
                           gradient=lambda x: (f(x), a * (x - m)))
    assert np.allclose(climb.x, np.clip(m, lower, upper), rtol=0, atol=1e-6 / a.min())   # pg 1e-6 at curvature a
    assert climb.analytic_pg <= 1e-6 and climb.gradients > 0 and climb.phases[0][0] == "analytic"


def test_newton_finish_on_an_ill_conditioned_box_quadratic():
    """The analytic pass's finish reaches a known box-constrained optimum: coordinates 2 and 5 on the
    lower bound, 3 on the upper, curvature from 1e-3 to 100. From the centre its steps are cut at
    bounds; from the optimum with coordinate 0 moved onto the lower bound it must leave that bound."""
    rng = np.random.default_rng(5)
    R = np.linalg.qr(rng.normal(size=(6, 6)))[0]
    A = R @ np.diag([1e-3, 1e-2, 0.1, 1.0, 10.0, 100.0]) @ R.T
    lower, upper = np.full(6, -1.0), np.full(6, 1.0)
    x_opt = np.array([0.3, -0.2, -1.0, 1.0, 0.5, -1.0])
    c = x_opt - np.linalg.solve(A, np.array([0.0, 0.0, 0.05, -0.07, 0.0, 0.02]))   # gradient there pushes out

    def fg(x):
        return float(0.5 * (x - c) @ A @ (x - c)), A @ (x - c)

    on_face = x_opt.copy()
    on_face[0] = -1.0
    for x0 in (np.zeros(6), on_face):
        f0, g0 = fg(x0)
        x, f, g, steps, ending, probes = opt._newton_finish(fg, lambda pts: [fg(p) for p in pts], x0, f0, g0,
                                                            lower, upper, 1e-9, None)
        assert ending == "reached pgtol" and 0 < steps <= opt.NEWTON_STEPS and probes > 0
        assert f <= f0
        assert np.allclose(x, x_opt, rtol=0, atol=1e-5)   # pg 1e-9 at curvature 1e-3 bounds the error by 1e-6
        assert x[2] == -1.0 and x[3] == 1.0 and x[5] == -1.0
