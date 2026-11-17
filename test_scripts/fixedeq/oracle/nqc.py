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
        self.tips = likelihood.tip_arrays(tree, self.columns)
        self.rates, self.weights = likelihood.rate_categories(shape, ncat, p_invar, free)
        self.rates = self.rates * part_rate
        self.invariant = ([likelihood.invariant_term(c, self.pi, p_invar) for c in patterns]
                          if p_invar else None)
        self.evaluations = 0

    @classmethod
    def from_patterns(cls, tree, columns, weights, pi, ref, shape=None, ncat=4, p_invar=0.0, build=chart.build,
                      free=None, part_rate=1.0):
        """Tip-vector columns (leaf name -> vector of len(pi)) with real weights in place of an
        alignment: observed counts, or exact pattern frequencies (S5, SIM-1). +I is refused, since
        its term reads characters."""
        if p_invar:
            raise ValueError("from_patterns refuses +I")
        names = {leaf.name for leaf in tree.leaves()}
        if any(set(col) != names or any(len(v) != len(pi) for v in col.values()) for col in columns):
            raise ValueError("every column needs one vector of len(pi) per leaf")
        self = cls.__new__(cls)
        self.tree, self.pi, self.ref, self.build = tree, np.asarray(pi, dtype=float), ref, build
        self.columns, self.counts = list(columns), np.asarray(weights, dtype=float)
        self.tips = likelihood.tip_arrays(tree, self.columns)
        self.rates, self.weights = likelihood.rate_categories(shape, ncat, 0.0, free)
        self.rates = self.rates * part_rate
        self.invariant = None
        self.evaluations = 0
        return self

    def log_likelihood(self, Q):
        site = likelihood.site_log_likelihoods(self.tree, self.columns, Q, self.pi, self.rates,
                                               self.weights, self.invariant, tips=self.tips)
        return float(self.counts @ site)

    def log_likelihood_gradient(self, Q):
        """lnL and G = dlnL/dQ (likelihood.log_likelihood_gradient); +I is refused."""
        if self.invariant is not None:
            raise ValueError("the analytic gradient does not cover +I")
        return likelihood.log_likelihood_gradient(self.tree, self.columns, Q, self.pi, self.rates, self.weights,
                                                  self.counts, tips=self.tips)

    def gradient(self, z):
        """f(z) and its analytic gradient in z; a failed build raises chart.BuildError."""
        Q = self.build(z, self.pi, self.ref)[0]
        lnL, G = self.log_likelihood_gradient(Q)
        return -lnL, -z_gradient(z, self.pi, self.ref, G)

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

    @property
    def pi(self):
        return self.parts[0].pi

    @property
    def ref(self):
        return self.parts[0].ref

    def __call__(self, z):
        self.evaluations += 1
        total = 0.0
        for f in self.parts:
            total += f(z)
        return total

    def log_likelihood(self, Q):
        total = 0.0
        for f in self.parts:
            total += f.log_likelihood(Q)
        return total

    def log_likelihood_gradient(self, Q):
        total, G = 0.0, 0.0
        for f in self.parts:
            lnL, g = f.log_likelihood_gradient(Q)
            total += lnL
            G = G + g
        return total, G

    def gradient(self, z):
        Q = self.parts[0].build(z, self.pi, self.ref)[0]
        lnL, G = self.log_likelihood_gradient(Q)
        return -lnL, -z_gradient(z, self.pi, self.ref, G)


def z_gradient(z, pi, ref, G):
    """<G, dQ/dz_k> for every coordinate k, as the adjoint of chart.derivative's linear map (tested
    against it): with H_ij = G_ij - G_ii, a_i = sum_j H_ij K_ij / pi_i and c solving
    (I - K + 1 nu^T) c = a, the pairing is sum_ij beta_ij dK_ij with beta_ij = H_ij nu_i / pi_i + nu_i c_j,
    and dK_ij / dz_ik = K_ij (delta_jk - K_ik)."""
    n = len(ref)
    pi, G = np.asarray(pi, dtype=float), np.asarray(G, dtype=float)
    K = chart.jump_matrix(z, ref)
    nu = chart.generator_from_jump(K, pi)[1]
    H = G - np.diag(G)[:, None]
    np.fill_diagonal(H, 0.0)
    a = (H * K).sum(axis=1) / pi
    c = np.linalg.solve(np.eye(n) - K + np.outer(np.ones(n), nu), a)
    beta = H * (nu / pi)[:, None] + np.outer(nu, c)
    row = (beta * K).sum(axis=1)
    return np.concatenate([K[i, free] * (beta[i, free] - row[i])
                           for i, free in enumerate(chart.free_destinations(ref))])


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


def export_gradient(f, z, eta=chart.FIRST_ORDER_ETA):
    """The export diagnostic's central differences (decision 009), probe for probe as the C++ takes
    them: h = eta max(1, |z_k|), probes at z_k + h and z_k - h, g_k = (f+ - f-) / (h+ + h-) with h+-
    the displacements the probes make. Returns f at z, g and the displacement sums h+ + h-."""
    x = np.array(z, dtype=float)
    f0 = f(x)
    g, span = np.empty(len(x)), np.empty(len(x))
    for k in range(len(x)):
        temp = x[k]
        h = eta * max(1.0, abs(temp))
        x[k] = temp + h
        up_step = x[k] - temp
        up = f(x)
        x[k] = temp - h
        down_step = temp - x[k]
        down = f(x)
        x[k] = temp
        span[k] = up_step + down_step
        g[k] = (up - down) / span[k]
    return f0, g, span


def fit(f, z0):
    """One model fit as ModelMarkov::optimizeParameters makes it, with the scaled step and the
    benchmark box: the oracle's side of the decision 017 toy comparison, calibrated in that box."""
    n = len(z0)
    lower, upper = np.full(n, BOX_LOWER), np.full(n, BOX_UPPER)
    return optimize.optimize_model(f, z0, lower, upper, GRADIENT_EPSILON,
                                   derivative=lambda g, x: optimize.scaled_derivative(g, x, upper))
