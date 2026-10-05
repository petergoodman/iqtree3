"""Rooted log-likelihood under a fixed Q, with IQ-TREE's conventions (CODE_PLAN.md 3.1).

Copied from IQ-TREE at 63c330d9 on purpose: the character map, the -m FILE convention (pi solved
from Q, Q scaled to mean rate 1), the +I term, and lengths of 0 or less raised to 1e-6.
Computed independently: transition matrices by SciPy's matrix exponential, and pruning with
per-node scaling in log space.
"""
import numpy as np
from scipy.linalg import expm
from scipy.special import logsumexp

from .chart import solve_stationary
from .gamma import category_rates, invariant_rates

AMINO = "ARNDCQEGHILKMFPSTWYV"
AMBIGUOUS = {"B": (2, 3), "Z": (5, 6), "J": (9, 10)}  # model/modelprotein.cpp:1350-1366
UNKNOWN = frozenset("X?-.~!*UO")                       # alignment/alignment.cpp:1720-1730, 1812-1843
MIN_BRANCH_LENGTH = 1e-6                               # tree/phylotree.cpp:3963-3966


def states_of(c):
    c = c.upper()
    if c in AMINO:
        return (AMINO.index(c),)
    if c in AMBIGUOUS:
        return AMBIGUOUS[c]
    if c in UNKNOWN:
        return tuple(range(20))
    raise ValueError(f"{c!r} is not a protein character")


def tip_vector(c):
    v = np.zeros(20)
    v[list(states_of(c))] = 1.0
    return v


def invariant_term(column, pi, p_invar):
    """p times the frequency of the states every character allows (phylotreesse.cpp:567-602,
    alignment.cpp:1306-1338): p for an all-unknown column, 0 for a variable one."""
    common = set(range(20))
    for c in column:
        common &= set(states_of(c))
    if len(common) == 20:
        return p_invar
    total = 0.0
    for x in sorted(common):
        total += pi[x]
    return total * p_invar if common else 0.0


def file_model(Q):
    """The -m FILE convention: the diagonal rebuilt from row sums, pi solved from Q by QR, Q
    scaled to mean rate 1 at that pi, and the root distribution pi normalized to sum 1
    (model/modelmarkov.cpp:1252-1282, 2119-2131, 875-892)."""
    Q = np.array(Q, dtype=float)
    np.fill_diagonal(Q, 0.0)
    np.fill_diagonal(Q, -Q.sum(axis=1))
    pi = solve_stationary(Q)
    return Q * (1.0 / float(pi @ -np.diag(Q))), pi / pi.sum()


def _branch(node):
    if node.length is None:
        raise ValueError("every branch needs a length")
    return node.length if node.length > 0 else MIN_BRANCH_LENGTH


def site_log_likelihoods(tree, columns, Q, pi, rates=(1.0,), weights=(1.0,), invariant=None):
    """Log-likelihood of each column; a column maps each leaf name to its tip vector."""
    per_category = []
    for rate in rates:
        cache = {}

        def partial(nd):
            if not nd.children:
                return np.array([col[nd.name] for col in columns], dtype=float), np.zeros(len(columns))
            L, scale = np.ones((len(columns), len(pi))), np.zeros(len(columns))
            for ch in nd.children:
                key = (rate, _branch(ch))
                if key not in cache:
                    cache[key] = expm(Q * rate * key[1])
                Lc, sc = partial(ch)
                L *= Lc @ cache[key].T
                scale += sc
            m = L.max(axis=1)
            return L / m[:, None], scale + np.log(m)

        L, scale = partial(tree)
        root = pi if not tree.length else pi @ expm(Q * rate * tree.length)
        per_category.append(np.log(L @ root) + scale)
    site = logsumexp(np.array(per_category), axis=0, b=np.asarray(weights)[:, None])
    if invariant is not None:
        inv = np.asarray(invariant, dtype=float)
        site = np.where(inv > 0, np.logaddexp(site, np.log(np.where(inv > 0, inv, 1.0))), site)
    return site


def alignment_log_likelihood(alignment, tree, Q, pi, shape=None, ncat=4, p_invar=0.0, sites=None):
    """Total log-likelihood of an alignment (name -> sequence), summed over distinct columns."""
    names = [leaf.name for leaf in tree.leaves()]
    if set(names) != set(alignment):
        raise ValueError("the tree's leaves and the alignment's names differ")
    length = len(next(iter(alignment.values())))
    counts = {}
    for s in (range(length) if sites is None else sites):
        col = "".join(alignment[n][s] for n in names).upper()
        counts[col] = counts.get(col, 0) + 1
    patterns = list(counts)
    columns = [{n: tip_vector(c) for n, c in zip(names, col)} for col in patterns]
    if shape is not None:
        rates, weights = category_rates(shape, ncat, p_invar)
    elif p_invar:
        rates, weights = invariant_rates(p_invar)
    else:
        rates, weights = np.array([1.0]), np.array([1.0])
    invariant = [invariant_term(col, pi, p_invar) for col in patterns] if p_invar else None
    site = site_log_likelihoods(tree, columns, Q, pi, rates, weights, invariant)
    return float(np.array([counts[p] for p in patterns]) @ site)
