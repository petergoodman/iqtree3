"""SIM-1's simulation inputs (CODE_PLAN.md 3.8; decisions 040, 041, 044): the composition pi*, the
truths at it, the circulation statistic C, the dial from a truth's reversible part to its full
circulation, transports to other compositions, and exact pattern weights over every column, for a
mixture of components (one truth; a composition or profile mixture; a topology mixture), with
masking.

Fluxes f_ij = pi_i q_ij are vectors in t3._edges order (flux.flux_of). Columns run in
itertools.product order, the first leaf name most significant (reference.exact_objective's order).
A masked residue is a 21st letter scored as unknown (all ones), as IQ-TREE scores a gap.
"""
import itertools
import math
from collections import namedtuple

import numpy as np

from . import builtin, chart, flux, likelihood, t3

STATES = 20
GAP = STATES  # the masked letter in a 21-letter column

# A mixture component: weight (normalized over the mixture), rooted tree, Q, its root distribution
# pi, and rate categories as likelihood.rate_categories gives them.
Component = namedtuple("Component", "weight tree Q pi rates rate_weights", defaults=((1.0,), (1.0,)))


def truths():
    """pi* and the truths at it, each at unit mean rate: NQ.bird under the -m FILE convention (its
    stored matrix with the diagonal rebuilt from row sums, pi* solved from it by QR), and Q.bird and
    LG from their exchangeabilities at pi* (chart.reversible_seed)."""
    Q, pi = likelihood.file_model(builtin.full_matrix("NQ.BIRD")[0])
    return pi, {"NQ.bird": Q,
                "Q.bird": chart.reversible_seed(builtin.exchangeabilities("Q.BIRD")[0], pi),
                "LG": chart.reversible_seed(builtin.exchangeabilities("LG")[0], pi)}


def target():
    """pi*, NQ.bird's stationary distribution (decision 040)."""
    return truths()[0]


def flux_parts(Q, pi):
    """S and A, the symmetric and antisymmetric parts of Q's flux at pi, in t3._edges order; exactly
    symmetric and antisymmetric, since each pair's sum and difference is formed once."""
    f = flux.flux_of(Q, pi)
    n = len(pi)
    i, j, _ = t3._edges(n)
    r = f[flux.edge_index(n)[j, i]]
    return (f + r) / 2, (f - r) / 2


def circulation(Q, pi):
    """Decision 041's C = sum_{i<j} |f_ij - f_ji| / sum_{i!=j} f_ij."""
    S, A = flux_parts(Q, pi)
    return float(np.abs(A).sum() / S.sum())


def dial(Q, pi, lam):
    """The generator at pi of the flux S + lam A: Q at lam = 1, its reversible part at 0, and
    C = lam C(Q). Past dial_limit a flux is not positive, and flux.generator refuses it."""
    S, A = flux_parts(Q, pi)
    return flux.generator(S + lam * A, pi)


def transposed(Q, pi):
    """The generator at pi of Q's transposed flux f_ji, dial(Q, pi, -1) since S - A is the transpose: the
    time reversal pi_j q_ji / pi_i, with C unchanged (SIM-2's cell F start)."""
    return dial(Q, pi, -1.0)


def dial_limit(Q, pi):
    """The largest lam at which every flux of the dial stays non-negative, min S / |A|."""
    S, A = flux_parts(Q, pi)
    moving = A != 0
    return float(np.min(S[moving] / np.abs(A[moving]))) if moving.any() else math.inf


def t3_transport(Q, pi_from, pi_to):
    """Decision 044: Q's T3 strengths and cycle affinities kept, its potentials re-balanced at
    pi_to, at unit mean rate; for a reversible Q, frequency replacement."""
    return t3.build(*t3.inverse(Q, pi_from), pi_to)


def t1_transport(Q, pi_to):
    """Decision 044's alternative, not chosen: the flux of a unit-mean-rate Q kept (the jump-chain
    transfer, chart.transfer_seed), and so C."""
    return chart.transfer_seed(Q, pi_to)


def shift(pi, toward, s):
    """pi (toward / pi)^s, normalized: pi at s = 0, toward at s = 1 (the screen's global shift)."""
    pi = np.asarray(pi, dtype=float)
    v = pi * (np.asarray(toward, dtype=float) / pi) ** s
    return v / v.sum()


def total_variation(p, q):
    return float(np.abs(np.asarray(p, dtype=float) - np.asarray(q, dtype=float)).sum() / 2)


def shift_to_tv(pi, toward, tv):
    """s and shift(pi, toward, s) at total variation tv from pi (SIM-2's shifted targets): the root in
    the first of the brackets [0, 1], [1, 2], [2, 4], ... whose upper end reaches tv, since TV is not
    proven monotone in s. Raises ValueError for a tv below 0 or out of the direction's reach."""
    from scipy.optimize import brentq

    if tv < 0:
        raise ValueError("tv must be non-negative")
    if tv == 0:
        return 0.0, shift(pi, toward, 0.0)

    def gap(s):
        return total_variation(shift(pi, toward, s), pi) - tv

    lo, hi = 0.0, 1.0
    with np.errstate(over="ignore", invalid="ignore"):
        for _ in range(64):
            g = gap(hi)
            if not np.isfinite(g):
                break
            if g >= 0:
                s = brentq(gap, lo, hi, xtol=1e-15)
                return s, shift(pi, toward, s)
            lo, hi = hi, 2 * hi
    raise ValueError(f"total variation {tv} is out of reach toward this composition")


def _deviations(profiles, weights):
    P = np.asarray(profiles, dtype=float)
    w = np.asarray(weights, dtype=float)
    return P - (w / w.sum()) @ P


def recentre(pi, profiles, weights, lam):
    """pi + lam (pi_k - pibar) for each profile pi_k, pibar their weighted mean, so that the
    weighted mean of the results is pi (the screen's per-gene compositions and site profiles)."""
    return np.asarray(pi, dtype=float) + lam * _deviations(profiles, weights)


def recentre_limit(pi, profiles, weights, floor=0.0):
    """The largest lam at which every recentred entry is at least floor; inf if none decreases."""
    d = _deviations(profiles, weights)
    falling = d < 0
    room = np.asarray(pi, dtype=float) - floor
    return float(np.min(np.where(falling, room / np.where(falling, -d, 1.0), math.inf)))


def mean_flux(components):
    """sum_k w_k f_k over the components (weights normalized): a balanced flux of total 1, the
    generating reference of a composition or profile mixture, a generator at pi* by flux.generator."""
    total = sum(c.weight for c in components)
    return sum(c.weight / total * flux.flux_of(c.Q, c.pi) for c in components)


def column_digits(k, base=STATES):
    """Every column of k letters, as a (base^k, k) array in itertools.product order."""
    return (np.arange(base ** k)[:, None] // base ** np.arange(k - 1, -1, -1)) % base


def _letters(gaps):
    return np.vstack([np.eye(STATES), np.ones(STATES)]) if gaps else np.eye(STATES)


def columns(names, gaps=False):
    """Every column over names as leaf name -> tip vector, in column_digits order, for
    nqc.Objective.from_patterns; with gaps, the 21st letter is the masked residue. The vectors are
    shared rows of one array."""
    letters = list(_letters(gaps))
    return [dict(zip(names, (letters[d] for d in row))) for row in column_digits(len(names), len(letters))]


def exact_weights(components, names, sites=1e6):
    """sites times each column's probability under the mixture of the components, over the 20^k
    columns of columns(names); every tree has the leaves names and is rooted at its component's pi."""
    letters, digits = _letters(False), column_digits(len(names))
    tips = {name: letters[digits[:, k]] for k, name in enumerate(names)}
    total, wsum = np.zeros(len(digits)), sum(c.weight for c in components)
    for c in components:
        if sorted(leaf.name for leaf in c.tree.leaves()) != sorted(names):
            raise ValueError("every component's tree needs exactly the leaves names")
        site = likelihood.site_log_likelihoods(c.tree, None, c.Q, c.pi, c.rates, c.rate_weights, tips=tips)
        total += c.weight / wsum * np.exp(site)
    return sites * total


def mask(weights, prob):
    """Weights over the 21^k columns of columns(names, gaps=True) from weights over the 20^k columns
    of columns(names): each leaf's residue masked independently with probability prob[column, leaf].
    Exact, as a sum over the 2^k masked subsets of each column."""
    weights, prob = np.asarray(weights, dtype=float), np.asarray(prob, dtype=float)
    m, k = prob.shape
    digits = column_digits(k)
    if m != len(digits) or len(weights) != m:
        raise ValueError("weights and prob need one row per column of 20^k")
    place = (STATES + 1) ** np.arange(k - 1, -1, -1)
    out = np.zeros((STATES + 1) ** k)
    for subset in itertools.product((False, True), repeat=k):
        s = np.array(subset)
        factor = np.prod(np.where(s, prob, 1.0 - prob), axis=1)
        out += np.bincount(np.where(s, GAP, digits) @ place, weights * factor, minlength=len(out))
    return out


def residue_mask(k, p_by_state):
    """prob for mask: each residue masked with its own probability."""
    return np.asarray(p_by_state, dtype=float)[column_digits(k)]


def context_mask(k, p):
    """prob for mask: a residue matching no other leaf's residue in its column masked with
    probability p."""
    d = column_digits(k)
    alone = (d[:, :, None] == d[:, None, :]).sum(axis=2) == 1
    return np.where(alone, p, 0.0)


def rarity_mask(states, s, edges, p, u):
    """SIM-2's masking of simulated alignments (decision 045), True where masked: a residue is masked
    where u < min(1, s p[k]), k the bin edges[k] <= r < edges[k+1] (r = 1 in the last bin) of r, the
    share of the other residues in its column identical to it. states is (taxa, sites) in 0..19, without
    gaps; u holds one uniform draw per residue, and sharing it across strengths s nests the masks."""
    states, u = np.asarray(states), np.asarray(u, dtype=float)
    edges, p = np.asarray(edges, dtype=float), np.asarray(p, dtype=float)
    if states.ndim != 2 or states.shape[0] < 2 or np.any((states < 0) | (states >= STATES)):
        raise ValueError("states need two or more taxa and every residue in 0..19")
    if u.shape != states.shape:
        raise ValueError("u needs one draw per residue")
    if (len(edges) != len(p) + 1 or edges[0] != 0.0 or edges[-1] != 1.0 or np.any(np.diff(edges) <= 0)
            or np.any(p < 0)):
        raise ValueError("edges must rise from 0 to 1, one more than the non-negative p")
    taxa, sites = states.shape
    counts = np.zeros((STATES, sites), dtype=int)
    np.add.at(counts, (states, np.broadcast_to(np.arange(sites), states.shape)), 1)
    r = (counts[states, np.arange(sites)] - 1) / (taxa - 1)
    k = np.minimum(np.searchsorted(edges, r, side="right") - 1, len(p) - 1)
    return u < np.minimum(1.0, s * p[k])


def masked_fraction(weights, prob):
    """The expected fraction of residues mask removes."""
    weights = np.asarray(weights, dtype=float)
    return float(weights @ np.asarray(prob, dtype=float).mean(axis=1) / weights.sum())
