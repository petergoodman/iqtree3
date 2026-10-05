"""T3 balancing chart in the D04 gauge (P-log section 4 and Appendix B; synthesis section 6.1).

Symmetric log strengths a with a_{n-1,n} = 0 and antisymmetric tilts h with h_{i,n} = 0, in the
documents' one-based notation (189 + 171 coordinates at n = 20). N_ij = pi_i pi_j exp(a_ij + h_ij)
is balanced by potentials g with g_n = 0, found by damped Newton; f = x / sum(x) and
q_ij = f_ij / pi_i.
"""
import numpy as np

from .chart import BuildError, check_target


def strength_pairs(n):
    return [(i, j) for i in range(n) for j in range(i + 1, n) if (i, j) != (n - 2, n - 1)]


def tilt_pairs(n):
    return [(i, j) for i in range(n - 1) for j in range(i + 1, n - 1)]


def pack(a, h):
    n = len(a)
    return np.array([a[i, j] for i, j in strength_pairs(n)] + [h[i, j] for i, j in tilt_pairs(n)])


def unpack(theta, n):
    a, h = np.zeros((n, n)), np.zeros((n, n))
    sp = strength_pairs(n)
    for (i, j), v in zip(sp, theta[:len(sp)]):
        a[i, j] = a[j, i] = v
    for (i, j), v in zip(tilt_pairs(n), theta[len(sp):]):
        h[i, j], h[j, i] = v, -v
    return a, h


def _edges(n):
    i, j = np.where(~np.eye(n, dtype=bool))
    D = (np.eye(n)[:, j] - np.eye(n)[:, i])[:-1]
    return i, j, D


def _generator(f, pi):
    q = f / pi[:, None]
    np.fill_diagonal(q, -q.sum(axis=1))
    return q


def _balanced_flux(a, h, pi):
    n = len(pi)
    i, j, D = _edges(n)
    z = np.log(pi[i]) + np.log(pi[j]) + a[i, j] + h[i, j]
    if not np.all(np.isfinite(z)):
        raise BuildError("non-finite coordinate")
    z -= z.max()
    g = np.zeros(n - 1)
    for _ in range(100):
        x = np.exp(z + D.T @ g)
        residual = D @ x
        if np.max(np.abs(residual)) / x.sum() < 2e-14:
            break
        step = np.linalg.solve((D * x) @ D.T, -residual)
        alpha = 1.0
        while np.exp(z + D.T @ (g + alpha * step)).sum() > x.sum() + 1e-4 * alpha * (residual @ step) + 2e-15 * x.sum():
            alpha *= 0.5
            if alpha < 1e-12:
                raise BuildError("balancing stalled")
        g += alpha * step
    else:
        raise BuildError("balancing did not converge in 100 Newton steps")
    f = np.zeros((n, n))
    f[i, j] = x / x.sum()
    return f, i, j, D


def build(a, h, pi):
    pi = check_target(pi)
    return _generator(_balanced_flux(np.asarray(a, float), np.asarray(h, float), pi)[0], pi)


def to_d04(a, h):
    """Move any gauge to D04: shift the strengths so a_{n-1,n} = 0, and add the vertex gradient
    u_j - u_i with u_i = h_{i,n}, which zeroes the star at the last state."""
    n = len(a)
    a = a - a[n - 2, n - 1]
    np.fill_diagonal(a, 0.0)
    u = h[:, n - 1].copy()
    return a, h + u[None, :] - u[:, None]


def inverse(Q, pi):
    pi = check_target(pi)
    Q = np.asarray(Q, dtype=float)
    n = len(pi)
    off = ~np.eye(n, dtype=bool)
    if not np.all(Q[off] > 0):
        raise BuildError("off-diagonal rates must be positive")
    a, h = np.zeros((n, n)), np.zeros((n, n))
    a[off] = 0.5 * np.log((Q * Q.T)[off] / np.outer(pi, pi)[off])
    h[off] = 0.5 * np.log((pi[:, None] * Q)[off] / (pi[:, None] * Q).T[off])
    return to_d04(a, h)


def derivative(a, h, pi, da, dh):
    """Directional derivative of Q along (da, dh) (P-log section 4.4)."""
    pi = check_target(pi)
    f, i, j, D = _balanced_flux(np.asarray(a, float), np.asarray(h, float), pi)
    x = f[i, j]
    v = da[i, j] + dh[i, j]
    dg = np.linalg.solve((D * x) @ D.T, -D @ (x * v))
    w = v + D.T @ dg
    df = np.zeros_like(f)
    df[i, j] = x * (w - x @ w)
    return _generator(df, pi)
