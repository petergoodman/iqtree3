"""Jump-chain charts of the fixed-equilibrium family (D01; P-log section 3; synthesis section 3).

The coordinates are the n(n-2) log-ratios z_ij = log(K_ij / K_{i,r(i)}), for each row i and each
destination j other than i and its reference r(i), in row-major order with destinations
increasing. K is the jump matrix, nu its stationary distribution, and q_ij = nu_i K_ij / pi_i.
"""
import numpy as np
from scipy.linalg import qr, solve_triangular

PIN = 10.0  # the positive-ratio chart's pinned weight c (synthesis section 3.1)


class BuildError(ValueError):
    """A proposal that does not give a valid generator; nothing is returned (synthesis 4.3)."""


def free_destinations(ref):
    n = len(ref)
    return [[j for j in range(n) if j != i and j != ref[i]] for i in range(n)]


def solve_stationary(M):
    """x with x M = 0 and sum(x) = 1, from [1^T; M^T] x = e_1 by column-pivoted QR (CODE_PLAN 2.1)."""
    n = M.shape[0]
    A = np.vstack([np.ones(n), M.T])
    b = np.zeros(n + 1)
    b[0] = 1.0
    q, r, perm = qr(A, mode="economic", pivoting=True)
    x = np.empty(n)
    x[perm] = solve_triangular(r, q.T @ b)
    return x


def check_target(pi):
    pi = np.asarray(pi, dtype=float)
    if not np.all(np.isfinite(pi)) or np.any(pi <= 0):
        raise BuildError("the target must be finite and strictly positive")
    return pi


def jump_matrix(z, ref):
    """Stable row softmax over the off-diagonal destinations, with the reference logit at 0."""
    n = len(ref)
    z = np.asarray(z, dtype=float).reshape(n, n - 2)
    if not np.all(np.isfinite(z)):
        raise BuildError("non-finite coordinate")
    K = np.zeros((n, n))
    for i, free in enumerate(free_destinations(ref)):
        logit = np.full(n, -np.inf)
        logit[free] = z[i]
        logit[ref[i]] = 0.0
        e = np.exp(logit - logit.max())
        if np.count_nonzero(e) != n - 1:
            raise BuildError(f"softmax underflow in row {i}")
        K[i] = e / e.sum()
    return K


def generator_from_jump(K, pi):
    """Q and nu for a jump matrix K with zero diagonal and positive off-diagonal entries."""
    pi = check_target(pi)
    n = len(pi)
    nu = solve_stationary(K - np.eye(n))
    if not np.all(np.isfinite(nu)) or np.any(nu <= 0):
        raise BuildError("the jump chain's stationary distribution is not strictly positive")
    Q = (nu / pi)[:, None] * K
    np.fill_diagonal(Q, 0.0)
    np.fill_diagonal(Q, -Q.sum(axis=1))
    if not np.all(np.isfinite(Q)):
        raise BuildError("non-finite rate")
    return Q, nu


def build(z, pi, ref):
    return generator_from_jump(jump_matrix(z, ref), pi)


def jump_of(Q):
    off = ~np.eye(len(Q), dtype=bool)
    if not np.all(Q[off] > 0) or not np.all(np.isfinite(Q)):
        raise BuildError("off-diagonal rates must be finite and positive")
    K = Q / (-np.diag(Q))[:, None]
    np.fill_diagonal(K, 0.0)
    return K


def inverse(Q, ref):
    K = jump_of(np.asarray(Q, dtype=float))
    return np.concatenate([np.log(K[i, free] / K[i, ref[i]])
                           for i, free in enumerate(free_destinations(ref))])


def build_positive(w, pi, ref):
    """The positive-ratio chart: weights w with the reference weight pinned at PIN."""
    n = len(ref)
    w = np.asarray(w, dtype=float).reshape(n, n - 2)
    if not np.all(np.isfinite(w)) or np.any(w <= 0):
        raise BuildError("weights must be finite and positive")
    W = np.zeros((n, n))
    for i, free in enumerate(free_destinations(ref)):
        W[i, free] = w[i]
        W[i, ref[i]] = PIN
    return generator_from_jump(W / W.sum(axis=1, keepdims=True), pi)


def references(Q):
    """Row maxima of the off-diagonal rates; argmax returns the first maximum, so ties go to the lower index."""
    off = np.array(Q, dtype=float)
    np.fill_diagonal(off, -np.inf)
    return off.argmax(axis=1)


def normalize(Q, pi):
    return Q / float(pi @ -np.diag(Q))


def reversible_seed(R, pi):
    """q_ij = R_ij pi_j, normalized to unit mean rate (D05)."""
    pi = check_target(pi)
    Q = np.asarray(R, dtype=float) * pi[None, :]
    np.fill_diagonal(Q, 0.0)
    np.fill_diagonal(Q, -Q.sum(axis=1))
    return normalize(Q, pi)


def transfer_seed(Q0, pi):
    """Keep Q0's jump chain at the new target: the flux-preserving transfer of synthesis 6.3."""
    return generator_from_jump(jump_of(np.asarray(Q0, dtype=float)), pi)[0]


def diagnostics(Q, pi, nu=None):
    """The residuals of CODE_PLAN.md 2.1, as plain floats."""
    pi = np.asarray(pi, dtype=float)
    n = len(pi)
    off = ~np.eye(n, dtype=bool)
    exit_rates = -np.diag(Q)
    d = {
        "stationarity": float(np.max(np.abs(pi @ Q))),
        "row_sums": float(np.max(np.abs(Q.sum(axis=1)) / exit_rates)),
        "mean_rate_error": float(abs(pi @ exit_rates - 1.0)),
        "min_rate": float(Q[off].min()),
        "min_flux": float((pi[:, None] * Q)[off].min()),
    }
    if nu is not None:
        K = jump_of(Q)
        d["nu_residual"] = float(np.max(np.abs(nu @ K - nu)))
        d["exit_rate_mismatch"] = float(np.max(np.abs(exit_rates - nu / pi) / (nu / pi)))
    return d


def derivative(z, pi, ref, dz):
    """Directional derivative of Q along dz (P-log section 3.4)."""
    n = len(ref)
    K = jump_matrix(z, ref)
    Q, nu = generator_from_jump(K, pi)
    dz = np.asarray(dz, dtype=float).reshape(n, n - 2)
    dK = np.zeros((n, n))
    for i, free in enumerate(free_destinations(ref)):
        v = np.zeros(n)
        v[free] = dz[i]
        dK[i] = K[i] * (v - K[i] @ v)
    dnu = np.linalg.solve((np.eye(n) - K + np.outer(np.ones(n), nu)).T, nu @ dK)
    dQ = (dnu[:, None] * K + nu[:, None] * dK) / np.asarray(pi)[:, None]
    np.fill_diagonal(dQ, 0.0)
    np.fill_diagonal(dQ, -dQ.sum(axis=1))
    return dQ
