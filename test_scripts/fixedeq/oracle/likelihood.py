"""Rooted log-likelihood under a fixed Q, with IQ-TREE's conventions (CODE_PLAN.md 3.1).

Copied from IQ-TREE at 63c330d9 on purpose: the character map, the -m FILE convention (pi solved
from Q, Q scaled to mean rate 1), the +I term, and lengths of 0 or less raised to 1e-6.
Computed independently: transition matrices by SciPy's matrix exponential, and pruning with
per-node scaling in log space.
"""
import numpy as np
from scipy.linalg import expm, expm_frechet
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


def log_likelihood_gradient(tree, columns, Q, pi, rates=(1.0,), weights=(1.0,), counts=None):
    """lnL = counts . site log-likelihoods and G = dlnL/dQ, every entry of Q free (the diagonal
    included), for S5's reference (decision 042). One upward and one downward pass per rate
    category give, per branch c and category k, W = sum_s counts_s B(x) I(y) / L_s, with B the
    outside partial at c's parent and I c's inside partial, so that dlnL/dP_ck = W; then
    <W, L(A, E)> = <L(A^T, W), E> for the Frechet derivative L of expm at A = r_k t_c Q. Scaling as
    in site_log_likelihoods. No +I term and no root branch."""
    if tree.length:
        raise ValueError("the analytic gradient takes no root branch")
    Q, pi = np.asarray(Q, dtype=float), np.asarray(pi, dtype=float)
    nsite, n = len(columns), len(pi)
    counts = np.ones(nsite) if counts is None else np.asarray(counts, dtype=float)
    nodes = list(tree.preorder())
    zeros = np.zeros(nsite)
    passes, log_cat = [], []
    for rate in rates:
        P, inside, up = {}, {}, {}
        for nd in reversed(nodes):          # children before parents
            if nd.children:
                L, s = np.ones((nsite, n)), zeros.copy()
                for ch in nd.children:
                    a, sa = up[id(ch)]
                    L *= a
                    s += sa
                m = L.max(axis=1)
                inside[id(nd)] = (L / m[:, None], s + np.log(m))
            else:
                inside[id(nd)] = (np.array([col[nd.name] for col in columns], dtype=float), zeros)
            if nd is not tree:
                key = (rate, _branch(nd))
                if key not in P:
                    P[key] = expm(Q * rate * key[1])
                up[id(nd)] = (inside[id(nd)][0] @ P[key].T, inside[id(nd)][1])
        I, s = inside[id(tree)]
        log_cat.append(np.log(I @ pi) + s)
        passes.append((P, inside, up))
    site = logsumexp(np.array(log_cat), axis=0, b=np.asarray(weights)[:, None])
    W = {}
    for rate, w, (P, inside, up) in zip(rates, weights, passes):
        out = {id(tree): (np.tile(pi, (nsite, 1)), zeros)}
        for nd in nodes:                    # parents before children
            if not nd.children:
                continue
            O, so = out[id(nd)]
            for ch in nd.children:
                B, sb = O.copy(), so.copy()
                for d in nd.children:
                    if d is not ch:
                        B *= up[id(d)][0]
                        sb += up[id(d)][1]
                m = B.max(axis=1)
                B /= m[:, None]
                sb += np.log(m)
                key = (rate, _branch(ch))
                I, si = inside[id(ch)]
                factor = counts * np.exp(np.log(w) + sb + si - site)
                W[key] = W.get(key, 0.0) + (B * factor[:, None]).T @ I
                if ch.children:
                    o = B @ P[key]
                    m = o.max(axis=1)
                    out[id(ch)] = (o / m[:, None], sb + np.log(m))
    G = np.zeros((n, n))
    for (rate, length), Wk in W.items():
        A = Q * rate * length
        G += rate * length * expm_frechet(A.T, Wk, compute_expm=False)
    return float(counts @ site), G


def compress(alignment, tree, sites=None):
    """The alignment's distinct columns in first-seen order (leaf characters in the tree's leaf
    order), their counts, and their tip vectors (leaf name -> vector)."""
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
    return patterns, np.array([counts[p] for p in patterns], dtype=float), columns


def free_rate_categories(params):
    """Rates and weights of FreeRate categories given as weight, rate pairs, normalized as RateFree's
    constructor does (model/ratefree.cpp:31-77): rates divided by sum(w r); when the weights' sum is
    off 1 by more than 1e-5, weights multiplied by 1/sum (utils/tools.cpp:772-795) and rates divided
    by sum(w r) again. The proportions-only form is not ported."""
    w, r = [float(v) for v in params[0::2]], [float(v) for v in params[1::2]]
    if not w or len(w) != len(r):
        raise ValueError("FreeRate parameters must be weight, rate pairs")
    total, sum_w = 0.0, 0.0
    for a, b in zip(w, r):
        total += a * b
        sum_w += a
    r = [b / total for b in r]
    if abs(sum_w - 1.0) > 1e-5:
        inv = 1 / sum_w
        w = [a * inv for a in w]
        total = 0.0
        for a, b in zip(w, r):
            total += a * b
        r = [b / total for b in r]
    return np.array(r), np.array(w)


def rate_categories(shape=None, ncat=4, p_invar=0.0, free=None):
    if free is not None:
        if shape is not None or p_invar:
            raise ValueError("FreeRate categories with +G or +I are not ported")
        return free_rate_categories(free)
    if shape is not None:
        return category_rates(shape, ncat, p_invar)
    if p_invar:
        return invariant_rates(p_invar)
    return np.array([1.0]), np.array([1.0])


def alignment_log_likelihood(alignment, tree, Q, pi, shape=None, ncat=4, p_invar=0.0, sites=None, free=None):
    """Total log-likelihood of an alignment (name -> sequence), summed over distinct columns; free
    gives FreeRate weight, rate pairs in place of shape."""
    patterns, counts, columns = compress(alignment, tree, sites)
    rates, weights = rate_categories(shape, ncat, p_invar, free)
    invariant = [invariant_term(col, pi, p_invar) for col in patterns] if p_invar else None
    site = site_log_likelihoods(tree, columns, Q, pi, rates, weights, invariant)
    return float(counts @ site)
