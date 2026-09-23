#!/usr/bin/env python3
"""
verify_constrained_nq.py

Independent numerical verification of the mathematics used in the constrained
(fixed-equilibrium) nQ estimation plan. Every check re-derives the quantity
from the definitions in this file; nothing is imported from the earlier agents'
scripts. Run: python3 verify_constrained_nq.py [--n 20] [--seed 20260922]

Outputs a JSON record (results.json next to this file) and prints a summary.
Exit code is non-zero if any hard check fails.

Provenance: written 2026-09-22 for Peter Goodman / Masel lab, session
https://claude.ai/code/session_01Pfs5eUaekCtLdKeSE6Fh7d
"""
import argparse, json, sys, time, platform
import numpy as np
from scipy.linalg import expm, expm_frechet, null_space
from scipy.optimize import minimize

ap = argparse.ArgumentParser()
ap.add_argument("--n", type=int, default=20)
ap.add_argument("--seed", type=int, default=20260922)
args = ap.parse_args()
rng = np.random.default_rng(args.seed)
N = args.n
R = {}          # results record
FAIL = []       # hard failures

def record(key, value, ok=None, tol=None):
    R[key] = value if not isinstance(value, np.generic) else float(value)
    if ok is not None and not ok:
        FAIL.append(f"{key}={value}" + (f" (tol {tol})" if tol else ""))

# ----------------------------------------------------------------------------
# Basic objects
# ----------------------------------------------------------------------------
def offdiag_index(n):
    ii, jj = np.where(~np.eye(n, dtype=bool))
    return ii, jj

def generator_from_off(off):
    Q = off.copy(); np.fill_diagonal(Q, 0.0); np.fill_diagonal(Q, -Q.sum(1)); return Q

def stationary(Q):
    n = Q.shape[0]
    A = np.vstack([np.ones(n), Q.T]); b = np.zeros(n + 1); b[0] = 1.0
    return np.linalg.lstsq(A, b, rcond=None)[0]          # same system IQ-TREE solves by QR

def mean_rate(Q, pi):
    return -float(pi @ np.diag(Q))

def normalize(Q, pi):
    return Q / mean_rate(Q, pi)

def random_pi(n):
    return rng.dirichlet(np.full(n, 2.0))

# ----------------------------------------------------------------------------
# 1. Dimension of the fiber: rank of the stationarity map on off-diagonal rates
# ----------------------------------------------------------------------------
def stationarity_matrix(pi):
    """C such that (C q)_j = (pi Q)_j for the vector q of off-diagonal rates."""
    n = len(pi); ii, jj = offdiag_index(n); m = len(ii)
    C = np.zeros((n, m))
    for e, (i, j) in enumerate(zip(ii, jj)):
        C[j, e] += pi[i]          # inflow to j
        C[i, e] -= pi[i]          # outflow from i (through the diagonal)
    return C

pi = random_pi(N)
C = stationarity_matrix(pi)
rank_C = np.linalg.matrix_rank(C)
m_row = np.zeros(C.shape[1]); ii, jj = offdiag_index(N); m_row[:] = pi[ii]   # mean rate = sum_i pi_i q_ij
rank_Cm = np.linalg.matrix_rank(np.vstack([C, m_row]))
record("1.rank_stationarity_map", int(rank_C), rank_C == N - 1)
record("1.rank_with_mean_rate_row", int(rank_Cm), rank_Cm == N)
record("1.dim_fiber", int(N * (N - 1) - rank_Cm), N * (N - 1) - rank_Cm == (N - 1) ** 2 - 1)
record("1.dim_reversible", int(N * (N - 1) // 2 - 1))
record("1.dim_surplus", int((N - 1) ** 2 - 1 - (N * (N - 1) // 2 - 1)))
# Audit's correction: the image of the drift map is 1^perp, not pi^perp.
d = np.zeros((3, 3)); d[0, 1] = 1.0; d = generator_from_off(d); p3 = np.array([.2, .3, .5])
drift = p3 @ d
record("1.drift_sum_is_zero", float(abs(drift.sum())), abs(drift.sum()) < 1e-15)
record("1.drift_dot_pi_nonzero", float(drift @ p3), abs(drift @ p3) > 1e-3)

# ----------------------------------------------------------------------------
# 2. Chart C2: jump-chain / row-shape chart.  Coordinates = row-normalized
#    off-diagonal rates (one entry per row pinned), row scales from stationarity.
# ----------------------------------------------------------------------------
def c2_build(P, pi):
    """P: zero-diagonal row-stochastic (positive off-diagonal). Returns Q in G_pi."""
    n = len(pi)
    nu = stationary(P - np.eye(n))            # stationary vector of the jump chain
    lam = nu / pi                             # exit rates -q_ii
    Q = lam[:, None] * P
    np.fill_diagonal(Q, -lam)
    return Q, nu, lam

def c2_inverse(Q):
    lam = -np.diag(Q)
    P = Q / lam[:, None]; np.fill_diagonal(P, 0.0)
    return P

def random_P(n):
    A = np.exp(rng.normal(0, 1.5, (n, n))); np.fill_diagonal(A, 0.0)
    return A / A.sum(1, keepdims=True)

P = random_P(N)
Q2, nu, lam = c2_build(P, pi)
record("2.c2_stationarity_residual", float(np.max(np.abs(pi @ Q2))), np.max(np.abs(pi @ Q2)) < 1e-13)
record("2.c2_mean_rate_minus_1", float(abs(mean_rate(Q2, pi) - 1)), abs(mean_rate(Q2, pi) - 1) < 1e-13)
record("2.c2_all_offdiag_positive", bool(np.all(Q2[~np.eye(N, dtype=bool)] > 0)), bool(np.all(Q2[~np.eye(N, dtype=bool)] > 0)))
record("2.c2_roundtrip_P", float(np.max(np.abs(c2_inverse(Q2) - P))), np.max(np.abs(c2_inverse(Q2) - P)) < 1e-13)
# nu_i = pi_i (-q_ii) for any fiber member: take an arbitrary fiber member built by C1 later; here check on Q2
record("2.c2_nu_equals_pi_times_exit", float(np.max(np.abs(nu - pi * (-np.diag(Q2))))), np.max(np.abs(nu - pi * (-np.diag(Q2)))) < 1e-14)
# Interpretation as "NONREV with per-row pinned rates": free rates r_ij with r_{i,ref}=1;
# Q_ij = lam_i * r_ij / sum_k r_ik.  Check that rescaling any row of r leaves Q unchanged (row scale is a gauge).
r = P.copy(); scale = np.exp(rng.normal(0, 2, N)); r2 = r * scale[:, None]
Q2b, _, _ = c2_build(r2 / r2.sum(1, keepdims=True), pi)
record("2.c2_row_scale_is_gauge", float(np.max(np.abs(Q2b - Q2))), np.max(np.abs(Q2b - Q2)) < 1e-13)
# Jacobian rank of the chart (ratios to reference column) -> must be (n-1)^2-1
def c2_from_ratio_coords(z, pi):
    n = len(pi); ratios = np.ones((n, n)); mask = ~np.eye(n, dtype=bool)
    # reference column for row i: last off-diagonal column (n-1, or n-2 for i=n-1)
    zz = z.reshape(n, n - 2); k = 0
    Pm = np.zeros((n, n))
    for i in range(n):
        cols = [j for j in range(n) if j != i]
        ref = cols[-1]; free = cols[:-1]
        Pm[i, ref] = 1.0
        Pm[i, free] = zz[i]
    Pm = Pm / Pm.sum(1, keepdims=True)
    return c2_build(Pm, pi)[0]

def fd_jacobian(f, z, h=1e-6):
    f0 = f(z); J = np.zeros((f0.size, z.size))
    for k in range(z.size):
        zp = z.copy(); zp[k] += h; zm = z.copy(); zm[k] -= h
        J[:, k] = (f(zp) - f(zm)).ravel() / (2 * h)
    return J

z0 = np.exp(rng.normal(0, 1, N * (N - 2)))
J2 = fd_jacobian(lambda z: c2_from_ratio_coords(z, pi), z0)
record("2.c2_jacobian_rank", int(np.linalg.matrix_rank(J2, tol=1e-7)), np.linalg.matrix_rank(J2, tol=1e-7) == (N - 1) ** 2 - 1)
# Derivative of the stationary vector: d nu = nu dP Z, Z = (I - P + 1 nu)^{-1}
dP = rng.normal(size=(N, N)); np.fill_diagonal(dP, 0.0); dP -= dP.sum(1, keepdims=True) / (N - 1); np.fill_diagonal(dP, 0.0)
Z = np.linalg.inv(np.eye(N) - P + np.outer(np.ones(N), nu))
pred = nu @ dP @ Z; eps = 1e-6
num = (stationary(P + eps * dP - np.eye(N)) - stationary(P - eps * dP - np.eye(N))) / (2 * eps)
record("2.c2_dnu_formula_error", float(np.max(np.abs(pred - num))), np.max(np.abs(pred - num)) < 1e-7)

# ----------------------------------------------------------------------------
# 3. Chart C1: balancing (T3) chart. M_ij = r_ij e^{h_ij}; N_ij = pi_i M_ij pi_j;
#    minimize Phi(g) = sum N_ij e^{g_j - g_i} on g_n = 0; q_ij = M_ij pi_j e^{g_j-g_i}; normalize.
# ----------------------------------------------------------------------------
def c1_balance(Nmat, tol=1e-14, maxit=100):
    n = Nmat.shape[0]; g = np.zeros(n); its = 0
    for its in range(1, maxit + 1):
        F = Nmat * np.exp(g[None, :] - g[:, None]); np.fill_diagonal(F, 0.0)
        grad = F.sum(0) - F.sum(1)                          # column sum - row sum
        if np.max(np.abs(grad)) / F.sum() < tol: break
        W = F + F.T; H = np.diag(W.sum(1)) - W               # Laplacian
        step = np.zeros(n); step[:-1] = np.linalg.solve(H[:-1, :-1], -grad[:-1])
        # damped step with objective decrease test (with round-off allowance)
        Phi0 = F.sum(); s = 1.0
        while True:
            g_new = g + s * step
            F_new = Nmat * np.exp(g_new[None, :] - g_new[:, None]); np.fill_diagonal(F_new, 0.0)
            if F_new.sum() <= Phi0 + 1e-12 * Phi0 or s < 1e-10: break
            s *= 0.5
        g = g_new
    return g, its

def c1_build(rmat, hmat, pi):
    Mm = rmat * np.exp(hmat); np.fill_diagonal(Mm, 0.0)
    Nmat = pi[:, None] * Mm * pi[None, :]
    g, its = c1_balance(Nmat)
    off = Mm * pi[None, :] * np.exp(g[None, :] - g[:, None])
    Q = generator_from_off(off)
    return normalize(Q, pi), g, its

def c1_inverse(Q, pi):
    Mm = Q / pi[None, :]; np.fill_diagonal(Mm, 1.0)
    rmat = np.sqrt(Mm * Mm.T); hmat = 0.5 * np.log(Mm / Mm.T); np.fill_diagonal(hmat, 0.0)
    return rmat, hmat

def gauge_fix(rmat, hmat):
    """scale gauge: r_{n-2,n-1} = 1 ; potential gauge: h_{0,j} = 0 (star at state 0)."""
    n = rmat.shape[0]
    r_g = rmat / rmat[n - 2, n - 1]
    c = hmat[0].copy()                # h'_{ij} = h_ij + c_i - c_j with c_i = h_{0i}: then h'_{0j} = h_0j + 0 - h_0j = 0
    h_g = hmat + c[:, None] - c[None, :]
    return r_g, h_g

A = rng.normal(0, 1.0, (N, N)); rmat = np.exp((A + A.T) / 2)
Hh = rng.normal(0, 0.7, (N, N)); hmat = (Hh - Hh.T) / 2
Q1, g1, its1 = c1_build(rmat, hmat, pi)
record("3.c1_stationarity_residual", float(np.max(np.abs(pi @ Q1))), np.max(np.abs(pi @ Q1)) < 1e-13)
record("3.c1_mean_rate_minus_1", float(abs(mean_rate(Q1, pi) - 1)), abs(mean_rate(Q1, pi) - 1) < 1e-13)
record("3.c1_newton_iterations_cold", int(its1))
# uniqueness of g*: perturbed start reaches same g (mod constant)
r_g, h_g = gauge_fix(*c1_inverse(Q1, pi))
Q1b, g1b, _ = c1_build(r_g, h_g, pi)
record("3.c1_gauge_fixed_roundtrip", float(np.max(np.abs(Q1b - Q1))), np.max(np.abs(Q1b - Q1)) < 1e-11)
# gauge orbit invariance: r -> kappa r, h -> h + grad c
kappa = 3.7; cvec = rng.normal(size=N)
Q1c, _, _ = c1_build(kappa * rmat, hmat + cvec[:, None] - cvec[None, :], pi)
record("3.c1_gauge_orbit_invariance", float(np.max(np.abs(Q1c - Q1))), np.max(np.abs(Q1c - Q1)) < 1e-11)
# reversible slice h = 0 equals normalized R diag(pi) with g* = 0
Q1r, g1r, _ = c1_build(rmat, np.zeros((N, N)), pi)
Qrev = normalize(generator_from_off(rmat * pi[None, :]), pi)
record("3.c1_h0_equals_R_diag_pi", float(np.max(np.abs(Q1r - Qrev))), np.max(np.abs(Q1r - Qrev)) < 1e-13)
record("3.c1_h0_gstar_is_zero", float(np.max(np.abs(g1r - g1r[-1]))), np.max(np.abs(g1r - g1r[-1])) < 1e-12)
# Jacobian rank of the gauge-fixed chart
def c1_from_coords(theta, pi):
    n = len(pi); iu = np.triu_indices(n, 1)
    nr = n * (n - 1) // 2 - 1; logr = np.zeros(n * (n - 1) // 2); logr[:nr] = theta[:nr]   # last pair pinned at log r = 0
    rm = np.zeros((n, n)); rm[iu] = np.exp(logr); rm = rm + rm.T
    hm = np.zeros((n, n)); pairs = [(i, j) for i in range(1, n) for j in range(i + 1, n)]
    for k, (i, j) in enumerate(pairs): hm[i, j] = theta[nr + k]; hm[j, i] = -theta[nr + k]
    return c1_build(rm, hm, pi)[0]
theta0 = rng.normal(0, 0.5, (N - 1) ** 2 - 1)
J1 = fd_jacobian(lambda t: c1_from_coords(t, pi), theta0, h=1e-5)
record("3.c1_jacobian_rank", int(np.linalg.matrix_rank(J1, tol=1e-6)), np.linalg.matrix_rank(J1, tol=1e-6) == (N - 1) ** 2 - 1)
# C1 and C2 agree on a common fiber point
P_from_Q1 = c2_inverse(Q1); Q1_via_c2, _, _ = c2_build(P_from_Q1, pi)
record("3.c1_c2_agree_on_fiber_point", float(np.max(np.abs(Q1_via_c2 - Q1))), np.max(np.abs(Q1_via_c2 - Q1)) < 1e-12)

# ----------------------------------------------------------------------------
# 4. Bound check (audit's counterexample): raw rates in [1e-4, 100] but the
#    gauge-fixed h' coordinate exceeds 12.  Also the sharper rate bound q_ij <= 1/(2 pi_i).
# ----------------------------------------------------------------------------
q = np.ones((N, N)); np.fill_diagonal(q, 0.0)
for (i, j) in [(0, 1), (1, 2), (2, 0)]: q[i, j] = 100.0; q[j, i] = 1e-4
Qc = generator_from_off(q); piu = np.full(N, 1.0 / N)
record("4.counterexample_uniform_stationary", float(np.max(np.abs(piu @ Qc))), np.max(np.abs(piu @ Qc)) < 1e-12)
rc, hc = c1_inverse(Qc, piu); _, hcg = gauge_fix(rc, hc)
record("4.raw_h_max", float(np.max(np.abs(hc))))
record("4.gauge_fixed_h_max", float(np.max(np.abs(hcg))), np.max(np.abs(hcg)) > 20)
record("4.gauge_fixed_h_23_formula", float(1.5 * np.log(1e6)))
# sharper bound: for any fiber member, pi_i q_ij <= 1/2 (max over many random members)
mx = 0.0
for _ in range(200):
    Pt = random_P(6); pit = random_pi(6); Qt, _, _ = c2_build(Pt, pit)
    mx = max(mx, float(np.max((pit[:, None] * Qt)[~np.eye(6, dtype=bool)])))
record("4.max_flux_entry_over_samples", mx, mx <= 0.5)

# ----------------------------------------------------------------------------
# 5. Non-concavity of the constrained likelihood (circulant counterexample).
# ----------------------------------------------------------------------------
def circ(x): return np.array([[-1, (1 + x) / 2, (1 - x) / 2], [(1 - x) / 2, -1, (1 + x) / 2], [(1 + x) / 2, (1 - x) / 2, -1.0]])
t1, t2 = 0.1, 5.1; k = np.sqrt(3) * (t2 - t1) / 2; amp = 2 * np.exp(-1.5 * (t1 + t2)); x0 = np.pi / k
def L00(x): return (expm(t1 * circ(x))[:, 0] @ expm(t2 * circ(x))[:, 0]) / 3     # rooted 2-tip, pattern (0,0), root pi uniform
formula = lambda x: (1 + amp * np.cos(k * x)) / 9
record("5.circulant_formula_error", float(abs(L00(x0) - formula(x0))), abs(L00(x0) - formula(x0)) < 1e-14)
hh = 1e-3; curv = (np.log(L00(x0 + hh)) - 2 * np.log(L00(x0)) + np.log(L00(x0 - hh))) / hh ** 2
record("5.second_derivative_numeric", float(curv), curv > 0)
record("5.second_derivative_analytic", float(amp * k * k / (1 - amp)))
record("5.circulant_stationary_uniform", float(np.max(np.abs(np.full(3, 1 / 3) @ circ(x0)))), np.max(np.abs(np.full(3, 1 / 3) @ circ(x0))) < 1e-15)

# ----------------------------------------------------------------------------
# 6. IQ-TREE's forward-difference rule on a signed coordinate near zero, versus
#    a positive coordinate.  Objective: a smooth 3-state two-tip log-likelihood
#    of size ~1e4 (pattern counts), parameterized in C1 coordinates.
# ----------------------------------------------------------------------------
pi3 = np.array([.5, .3, .2]); counts = rng.multinomial(20000, np.full(9, 1 / 9)).reshape(3, 3).astype(float)
def loglik_c1(logr, h12, ta=0.3, tb=0.9):
    rm = np.exp(np.array([[0, logr[0], logr[1]], [logr[0], 0, logr[2]], [logr[1], logr[2], 0]]))
    hm = np.zeros((3, 3)); hm[1, 2] = h12; hm[2, 1] = -h12
    Q, _, _ = c1_build(rm, hm, pi3)
    J = expm(ta * Q).T @ np.diag(pi3) @ expm(tb * Q)
    return float(np.sum(counts * np.log(J)))
logr0 = np.array([0.2, -0.4, 0.1])
def iqtree_fd(f, x):            # IQ-TREE: h = 1e-4|x|, or 1e-4 if x == 0
    h = 1e-4 * abs(x) if x != 0 else 1e-4
    return (f(x + h) - f(x)) / h
def central(f, x, h=1e-5): return (f(x + h) - f(x - h)) / (2 * h)
fd_err = {}
for hval in [1e-2, 1e-4, 1e-6, 1e-8]:
    f = lambda h12: loglik_c1(logr0, h12)
    ref = central(f, hval)
    fd_err[f"h={hval:g}"] = {"iqtree_fd_error": float(abs(iqtree_fd(f, hval) - ref)), "shifted_by_12_error": float(abs(iqtree_fd(lambda s: f(s - 12.0), hval + 12.0) - ref)), "reference_derivative": float(ref)}
record("6.fd_rule_on_signed_coordinate", fd_err)
record("6.fd_error_grows_toward_zero", fd_err["h=1e-08"]["iqtree_fd_error"] > 10 * fd_err["h=0.01"]["iqtree_fd_error"], fd_err["h=1e-08"]["iqtree_fd_error"] > 10 * fd_err["h=0.01"]["iqtree_fd_error"])
# Positive ratio coordinate (C2 style) near its lower bound behaves like NONREV's own rates: relative step 1e-4
def loglik_c2(u, ta=0.3, tb=0.9):
    Pm = np.array([[0, u, 1.0], [0.6, 0, 1.0], [0.3, 0.9, 0]]); Pm = Pm / Pm.sum(1, keepdims=True)
    Q, _, _ = c2_build(Pm, pi3); J = expm(ta * Q).T @ np.diag(pi3) @ expm(tb * Q)
    return float(np.sum(counts * np.log(J)))
fd2 = {}
for uval in [1.0, 1e-2, 1e-4]:
    ref = central(loglik_c2, uval, h=1e-6 * uval)
    fd2[f"u={uval:g}"] = {"iqtree_fd_relative_error": float(abs(iqtree_fd(loglik_c2, uval) - ref) / max(abs(ref), 1e-300))}
record("6.fd_rule_on_positive_ratio_coordinate", fd2)

# ----------------------------------------------------------------------------
# 7. Seeding: IQ-TREE's adaptStateFrequency (q_ij <- q_ij pi_new_j / pi_old_j) keeps
#    stationarity for a reversible seed and breaks it for a nonreversible seed.
# ----------------------------------------------------------------------------
pi_old = random_pi(N); pi_new = random_pi(N)
Rsym = rmat.copy(); Qrev_old = generator_from_off(Rsym * pi_old[None, :])
Qrev_new = generator_from_off(Rsym * pi_old[None, :] * (pi_new / pi_old)[None, :])
record("7.swap_reversible_seed_stationary_at_new_pi", float(np.max(np.abs(pi_new @ Qrev_new))), np.max(np.abs(pi_new @ Qrev_new)) < 1e-13)
Qnr_old, _, _ = c2_build(random_P(N), pi_old)          # nonreversible with equilibrium pi_old
Qnr_swapped = generator_from_off(Qnr_old * (pi_new / pi_old)[None, :])
record("7.swap_nonreversible_seed_residual_at_new_pi", float(np.max(np.abs(pi_new @ Qnr_swapped))), np.max(np.abs(pi_new @ Qnr_swapped)) > 1e-3)
# C2 seeding from an arbitrary Q (any equilibrium) is immediate and exact
Qseed, _, _ = c2_build(c2_inverse(Qnr_old), pi_new)
record("7.c2_seed_from_foreign_Q_stationary", float(np.max(np.abs(pi_new @ Qseed))), np.max(np.abs(pi_new @ Qseed)) < 1e-13)

# ----------------------------------------------------------------------------
# 8. Proposition: under the constrained model with root distribution pi*, every
#    node's marginal is pi* (so tip composition is pinned).
# ----------------------------------------------------------------------------
marg = pi @ expm(0.7 * Q1) @ expm(1.9 * Q1)
record("8.tip_marginal_equals_pi", float(np.max(np.abs(marg - pi))), np.max(np.abs(marg - pi)) < 1e-13)

# ----------------------------------------------------------------------------
# 9. EM M-step: normalized (n dual variables, branch lengths fixed) versus cone
#    (n-1 dual variables, then renormalize with compensating branch rescale).
# ----------------------------------------------------------------------------
n4 = 4; pi4 = random_pi(n4); i4, j4 = offdiag_index(n4); Ncount = rng.uniform(1, 5, len(i4)); Texp = rng.uniform(2, 6, n4)
B = np.zeros((n4, len(i4)));
for e, (a, b) in enumerate(zip(i4, j4)): B[b, e] += 1; B[a, e] -= 1
dvec = Texp[i4] / pi4[i4]
def solve_dual(Amat, bvec, lam0):
    lam = lam0.copy()
    val = lambda l: (bvec @ l - np.sum(Ncount * np.log(dvec + Amat.T @ l))) if np.min(dvec + Amat.T @ l) > 0 else np.inf
    for _ in range(200):
        den = dvec + Amat.T @ lam; f = Ncount / den; gr = bvec - Amat @ f
        if np.max(np.abs(gr)) < 1e-12: break
        Hm = (Amat * (Ncount / den ** 2)) @ Amat.T; p = np.linalg.solve(Hm, -gr); s = 1.0
        while val(lam + s * p) > val(lam) + 1e-4 * s * (gr @ p) + 1e-13 * max(1, abs(val(lam))): s *= 0.5
        lam = lam + s * p
    return Ncount / (dvec + Amat.T @ lam), lam
Cfull = np.vstack([B[:-1], np.ones(len(i4))]); bfull = np.r_[np.zeros(n4 - 1), 1.0]
f_norm, _ = solve_dual(Cfull, bfull, np.r_[np.zeros(n4 - 1), Ncount.sum()])
f_cone, _ = solve_dual(B[:-1], np.zeros(n4 - 1), np.zeros(n4 - 1))
obj = lambda f: float(np.sum(Ncount * np.log(f)) - dvec @ f)
record("9.mstep_normalized_residual", float(np.max(np.abs(Cfull @ f_norm - bfull))), np.max(np.abs(Cfull @ f_norm - bfull)) < 1e-10)
record("9.mstep_cone_balance_residual", float(np.max(np.abs(B @ f_cone))), np.max(np.abs(B @ f_cone)) < 1e-10)
record("9.mstep_cone_mean_rate", float(f_cone.sum()))
record("9.mstep_normalized_beats_rescaled_cone_at_fixed_T", float(obj(f_norm) - obj(f_cone / f_cone.sum())), obj(f_norm) - obj(f_cone / f_cone.sum()) > 0)
# KKT: gradient of objective at f_norm lies in the row space of Cfull
grad = Ncount / f_norm - dvec; Zn = null_space(Cfull)
record("9.mstep_kkt_projected_gradient", float(np.max(np.abs(Zn.T @ grad))), np.max(np.abs(Zn.T @ grad)) < 1e-9)
# compensating rescale leaves transition matrices unchanged
Qc4 = generator_from_off(np.zeros((n4, n4))); off = np.zeros((n4, n4)); off[i4, j4] = f_cone / pi4[i4]; Qc4 = generator_from_off(off); mth = f_cone.sum()
record("9.compensated_rescale_transition_error", float(np.max(np.abs(expm(0.4 * Qc4) - expm(0.4 * mth * (Qc4 / mth))))), np.max(np.abs(expm(0.4 * Qc4) - expm(0.4 * mth * (Qc4 / mth)))) < 1e-14)

# ----------------------------------------------------------------------------
# 10. Frechet derivative of exp(tQ) is not t E exp(tQ); score identity check.
# ----------------------------------------------------------------------------
E = np.zeros((N, N)); E[0, 1] = 1.0; E[0, 0] = -1.0; t = 0.8
fre = expm_frechet(t * Q1, t * E, compute_expm=False)
fdm = (expm(t * (Q1 + 1e-6 * E)) - expm(t * (Q1 - 1e-6 * E))) / 2e-6
record("10.frechet_vs_central_difference", float(np.max(np.abs(fre - fdm))), np.max(np.abs(fre - fdm)) < 1e-7)
record("10.naive_tE_expm_error", float(np.max(np.abs(t * E @ expm(t * Q1) - fre))), np.max(np.abs(t * E @ expm(t * Q1) - fre)) > 1e-3)

# ----------------------------------------------------------------------------
# 11. Wrong-target circulation on population data (reproduces the audit's mechanism):
#     reversible truth, 3 states, 3-tip star; fit reversible vs full at wrong pi*.
# ----------------------------------------------------------------------------
pi_true = np.array([.6, .3, .1]); pi_wrong = np.array([.5, .3, .2]); times = [.2, .7, 1.1]
def params_to_Q(x, pi_):
    rm = np.exp(np.array([[0, x[0], x[1]], [x[0], 0, -x[0] - x[1]], [x[1], -x[0] - x[1], 0]]))
    hm = np.zeros((3, 3))
    if len(x) > 2: hm[1, 2] = x[2]; hm[2, 1] = -x[2]
    return c1_build(rm, hm, pi_)[0]
def pattern(Qm, pi_):
    ps = [expm(tt * Qm) for tt in times]; return np.einsum('r,ra,rb,rc->abc', pi_, *ps).ravel()
truth = pattern(params_to_Q([.7, -.5], pi_true), pi_true)
def KL(x, pi_):
    pp = pattern(params_to_Q(x, pi_), pi_); return float(np.sum(truth * np.log(truth / pp)))
out11 = {}
for label, pi_ in [("correct", pi_true), ("wrong", pi_wrong)]:
    rev = min((minimize(lambda x: KL(x, pi_), s, method='BFGS', options={'gtol': 1e-10}) for s in [np.zeros(2), np.array([.7, -.5]), np.array([-.5, .4])]), key=lambda z: z.fun)
    full = min((minimize(lambda x: KL(x, pi_), np.r_[rev.x, s], method='BFGS', options={'gtol': 1e-10}) for s in [.05, -.05, .3]), key=lambda z: z.fun)
    out11[label] = {"reversible_KL": float(rev.fun), "full_KL": float(full.fun), "gain": float(rev.fun - full.fun), "cycle_coordinate": float(full.x[2])}
record("11.wrong_target_population_example", out11)
record("11.correct_target_recovers_reversible", out11["correct"]["gain"] < 1e-9 and abs(out11["correct"]["cycle_coordinate"]) < 1e-3, out11["correct"]["gain"] < 1e-9)
record("11.wrong_target_induces_circulation", out11["wrong"]["gain"] > 1e-7 and abs(out11["wrong"]["cycle_coordinate"]) > 1e-3, out11["wrong"]["gain"] > 1e-7)

# ----------------------------------------------------------------------------
# 12. Nesting: reversible R diag(pi*) normalized lies in the fiber (already 3.c1_h0); the fiber lies
#     in the NONREV family (trivially); and the fiber is affinely isomorphic across pi via flux.
# ----------------------------------------------------------------------------
pi_a = random_pi(N); pi_b = random_pi(N); Qa, _, _ = c2_build(P, pi_a); Fa = pi_a[:, None] * Qa; np.fill_diagonal(Fa, 0)
Qb = generator_from_off(Fa / pi_b[:, None])
record("12.flux_transport_stationary_at_pi_b", float(np.max(np.abs(pi_b @ Qb))), np.max(np.abs(pi_b @ Qb)) < 1e-13)
record("12.flux_transport_mean_rate", float(mean_rate(Qb, pi_b)), abs(mean_rate(Qb, pi_b) - 1) < 1e-13)

# ----------------------------------------------------------------------------
R["_meta"] = {"n": N, "seed": args.seed, "numpy": np.__version__, "scipy": __import__('scipy').__version__,
              "python": platform.python_version(), "date": time.strftime("%Y-%m-%d %H:%M:%S %Z"), "failures": FAIL}
with open(__file__.replace('.py', '_results.json'), 'w') as fh:
    json.dump(R, fh, indent=1, default=float)
for kk, vv in R.items():
    if kk.startswith('_'): continue
    print(f"{kk:55s} {vv}")
print("\nHARD FAILURES:", FAIL if FAIL else "none")
sys.exit(1 if FAIL else 0)
