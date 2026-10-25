"""The NQC model as the S2 checks see it (decision 030): coordinates z -> Q at pi* by the log-ratio
chart -> IQ-TREE's rooted likelihood with pi* as the root distribution.

The objective builds the site patterns once, so a gradient costs one likelihood evaluation per
coordinate plus one, as in IQ-TREE.
"""
import math

import numpy as np

from . import builtin, chart, likelihood, optimize

LOWER, UPPER = math.log(1e-5), math.log(1e5)  # the default domain (decisions 007 and 027)
# decision 007's first default, the benchmark box, kept for the records made in it: the decision 017
# toy comparison and S3c's stopping study
BOX_LOWER, BOX_UPPER = math.log(1e-5), math.log(10.0)
GRADIENT_EPSILON = 1e-4  # ModelFactory's default, passed to the model fit (model/modelfactory.h:204-205)


def lg_seed(pi):
    """LG rebuilt at pi (D05) and its row-maximum references (decision 011)."""
    Q = chart.reversible_seed(builtin.exchangeabilities("LG")[0], pi)
    return Q, chart.references(Q)


class Objective:
    """f(z) = -lnL of the Q built from z at pi, on a fixed tree with fixed rate parameters; 1e30
    when the build fails, as the class's targetFunk returns. A partition rate (-p, -q) multiplies
    every rate category, which gives the likelihood of the partition tree IQ-TREE forms by scaling
    every branch by that rate."""

    def __init__(self, alignment, tree, pi, ref, shape=None, ncat=4, p_invar=0.0, build=chart.build, free=None,
                 part_rate=1.0):
        self.tree, self.pi, self.ref, self.build = tree, np.asarray(pi, dtype=float), ref, build
        patterns, self.counts, self.columns = likelihood.compress(alignment, tree)
        self.rates, self.weights = likelihood.rate_categories(shape, ncat, p_invar, free)
        self.rates = self.rates * part_rate
        self.invariant = ([likelihood.invariant_term(c, self.pi, p_invar) for c in patterns]
                          if p_invar else None)
        self.evaluations = 0

    def log_likelihood(self, Q):
        site = likelihood.site_log_likelihoods(self.tree, self.columns, Q, self.pi, self.rates,
                                               self.weights, self.invariant)
        return float(self.counts @ site)

    def __call__(self, z):
        self.evaluations += 1
        try:
            Q = self.build(z, self.pi, self.ref)[0]
        except chart.BuildError:
            return 1.0e30
        return -self.log_likelihood(Q)


class LinkedObjective:
    """The linked objective of decision 005, as PartitionModel::targetFunk forms it: each partition's
    Objective at the same z, summed from 0 in IQ-TREE's partition order (part_order)."""

    def __init__(self, parts):
        self.parts = list(parts)
        self.evaluations = 0

    def __call__(self, z):
        self.evaluations += 1
        total = 0.0
        for f in self.parts:
            total += f(z)
        return total


def normalize_part_rates(rates, sites):
    """PhyloSuperTreePlen::normalizePartRate's statements in order (tree/phylosupertreeplen.cpp): the
    rates a -q partition file gives as {x}, divided by their mean weighted by each partition's sites."""
    total, nsite = 0.0, 0
    for r, n in zip(rates, sites):
        total += r * n
        nsite += n
    total /= nsite
    total = 1.0 / total
    return [r * total for r in rates]


def scaled_steps(z, upper=UPPER):
    """The displacement decision 005's rule takes at each coordinate."""
    return np.array([optimize.scaled_step(x, upper) for x in z])


def legacy_steps(z):
    """The displacement IQ-TREE's legacy rule takes (utils/optimization.cpp:916-939)."""
    out = np.empty(len(z))
    for k, x in enumerate(z):
        h = optimize.ERROR_X * abs(x) or optimize.ERROR_X
        out[k] = (x + h) - x
    return out


def central(f, z, s):
    """Central first and second differences with h = s * max(1, |z_k|): the reference gradient and
    the curvature f''_kk of decision 030. optimize.central_derivative gives only the first."""
    x = np.array(z, dtype=float)
    f0 = f(x)
    g, c = np.empty(len(x)), np.empty(len(x))
    for k in range(len(x)):
        temp = x[k]
        h = s * max(1.0, abs(temp))
        x[k] = temp + h
        up = f(x)
        x[k] = temp - h
        down = f(x)
        x[k] = temp
        g[k] = (up - down) / (2 * h)
        c[k] = (up - 2 * f0 + down) / (h * h)
    return g, c


def fit(f, z0):
    """One model fit as ModelMarkov::optimizeParameters makes it, with the scaled step and the
    benchmark box: the oracle's side of the decision 017 toy comparison, calibrated in that box."""
    n = len(z0)
    lower, upper = np.full(n, BOX_LOWER), np.full(n, BOX_UPPER)
    return optimize.optimize_model(f, z0, lower, upper, GRADIENT_EPSILON,
                                   derivative=lambda g, x: optimize.scaled_derivative(g, x, upper))
