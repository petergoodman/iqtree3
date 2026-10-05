"""Test generators built independently of the charts they check (P-log Appendix B)."""
import numpy as np


def random_target(n, rng, concentration=2.0):
    return rng.dirichlet(np.full(n, concentration))


def cycle_flux(n, rng):
    """A balanced positive flux summing to 1: a symmetric part plus 3n directed cycles."""
    f = np.zeros((n, n))
    for u in range(n):
        for v in range(u + 1, n):
            f[u, v] = f[v, u] = rng.lognormal(0, 1)
    for _ in range(3 * n):
        cycle = rng.choice(n, size=min(n, 3), replace=False)
        w = rng.lognormal()
        for u, v in zip(cycle, np.roll(cycle, -1)):
            f[u, v] += w
    return f / f.sum()


def generator_from_flux(f, pi):
    q = f / pi[:, None]
    np.fill_diagonal(q, 0.0)
    np.fill_diagonal(q, -q.sum(axis=1))
    return q


def random_exchangeabilities(n, rng):
    a = rng.normal(size=(n, n))
    r = np.exp((a + a.T) / 2)
    np.fill_diagonal(r, 0.0)
    return r
