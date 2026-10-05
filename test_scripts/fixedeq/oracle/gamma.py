"""Discrete gamma rates as IQ-TREE computes them (model/rategamma.cpp at 63c330d9).

The four numerical routines are ports of PAML's, which IQ-TREE uses (rategamma.cpp:302-475);
SciPy's own functions would not reproduce their stated accuracies (1e-8 for the incomplete
gamma, a relative change of 0.5e-6 for the chi-square quantile).
"""
from math import exp, log, sqrt

import numpy as np

MIN_GAMMA_SHAPE = 0.02  # utils/tools.h:605


def ln_gamma(alpha):
    """cmpLnGamma (305-322): Stirling's series, accurate to 10 decimal places."""
    x, f = alpha, 0.0
    if x < 7:
        f, z = 1.0, x - 1
        while True:
            z += 1
            if not z < 7:
                break
            f *= z
        x, f = z, -log(f)
    z = 1 / (x * x)
    return (f + (x - 0.5) * log(x) - x + .918938533204673
            + (((-.000595238095238 * z + .000793650793651) * z - .002777777777778) * z
               + .083333333333333) / x)


def incomplete_gamma(x, alpha, ln_gamma_alpha):
    """cmpIncompleteGamma (325-385), AS 32: series expansion or continued fraction."""
    p, g = alpha, ln_gamma_alpha
    accurate, overflow = 1e-8, 1e30
    if x == 0:
        return 0.0
    if x < 0 or p <= 0:
        return -1.0
    factor = exp(p * log(x) - x - g)
    if not (x > 1 and x >= p):
        gin, term, rn = 1.0, 1.0, p
        while True:
            rn += 1
            term *= x / rn
            gin += term
            if not term > accurate:
                break
        return gin * factor / p
    a = 1 - p
    b = a + x + 1
    term = 0.0
    pn = [1.0, x, x + 1, x * b, 0.0, 0.0]
    gin = pn[2] / pn[3]
    while True:
        a += 1
        b += 2
        term += 1
        an = a * term
        pn[4] = b * pn[2] - an * pn[0]
        pn[5] = b * pn[3] - an * pn[1]
        if pn[5] != 0:
            rn = pn[4] / pn[5]
            dif = abs(gin - rn)
            if dif <= accurate and dif <= accurate * rn:
                return 1 - factor * gin
            gin = rn
        pn[0:4] = pn[2:6]
        if not abs(pn[4]) < overflow:
            pn[0:4] = [v / overflow for v in pn[0:4]]


def point_normal(prob):
    """cmpPointNormal (391-414), AS 111."""
    a0, a1, a2, a3 = -.322232431088, -1, -.342242088547, -.0204231210245
    a4, b0, b1 = -.453642210148e-4, .0993484626060, .588581570495
    b2, b3, b4 = .531103462366, .103537752850, .0038560700634
    p = prob
    p1 = p if p < 0.5 else 1 - p
    if p1 < 1e-20:
        return -9999.0
    y = sqrt(log(1 / (p1 * p1)))
    z = y + ((((y * a4 + a3) * y + a2) * y + a1) * y + a0) / ((((y * b4 + b3) * y + b2) * y + b1) * y + b0)
    return -z if p < 0.5 else z


def point_chi2(prob, v):
    """cmpPointChi2 (420-475), AS 91."""
    e, aa, p = .5e-6, .6931471805, prob
    if p < .000002 or p > .999998 or v <= 0:
        return -1.0
    g = ln_gamma(v / 2)
    xx = v / 2
    c = xx - 1
    if v < -1.24 * log(p):
        ch = (p * xx * exp(g + xx * aa)) ** (1 / xx)
        if ch - e < 0:
            return ch
    elif v <= .32:
        ch, a = 0.4, log(1 - p)
        while True:
            q = ch
            p1 = 1 + ch * (4.67 + ch)
            p2 = ch * (6.73 + ch * (6.66 + ch))
            t = -0.5 + (4.67 + 2 * ch) / p1 - (6.73 + ch * (13.32 + 3 * ch)) / p2
            ch -= (1 - exp(a + g + .5 * ch + c * aa) * p2 / p1) / t
            if abs(q / ch - 1) - .01 <= 0:
                break
    else:
        x = point_normal(p)
        p1 = 0.222222 / v
        ch = v * (x * sqrt(p1) + 1 - p1) ** 3.0
        if ch > 2.2 * v + 6:
            ch = -2 * (log(1 - p) - c * log(.5 * ch) + g)
    while True:
        q = ch
        p1 = .5 * ch
        t = incomplete_gamma(p1, xx, g)
        if t < 0:
            return -1.0
        p2 = p - t
        t = p2 * exp(xx * aa + g + p1 - c * log(ch))
        b = t / ch
        a = 0.5 * t - b * c
        s1 = (210 + a * (140 + a * (105 + a * (84 + a * (70 + 60 * a))))) / 420
        s2 = (420 + a * (735 + a * (966 + a * (1141 + 1278 * a)))) / 2520
        s3 = (210 + a * (462 + a * (707 + 932 * a))) / 2520
        s4 = (252 + a * (672 + 1182 * a) + c * (294 + a * (889 + 1740 * a))) / 5040
        s5 = (84 + 264 * a + c * (175 + 606 * a)) / 2520
        s6 = (120 + c * (346 + 127 * c)) / 5040
        ch += t * (1 + 0.5 * t * s1 - b * c * (s1 - b * (s2 - b * (s3 - b * (s4 - b * (s5 - b * s6))))))
        if not abs(q / ch - 1) > e:
            return ch


def mean_rates(shape, ncat):
    """computeRatesMean (158-171): the mean rate of each of ncat equal-probability categories."""
    lnga1 = ln_gamma(shape + 1)
    cut = [point_chi2((i + 1.0) / ncat, 2.0 * shape) / (2.0 * shape) for i in range(ncat - 1)]
    freq = [incomplete_gamma(f * shape, shape + 1, lnga1) for f in cut]
    rates = [0.0] * ncat
    rates[0] = freq[0] * ncat
    rates[ncat - 1] = (1 - freq[ncat - 2]) * ncat
    for i in range(1, ncat - 1):
        rates[i] = (freq[i] - freq[i - 1]) * ncat
    return rates


def _compute_rates(previous, shape, ncat):
    """computeRates (98-150) with gamma_median false: new rates rescaled to the previous sum."""
    if ncat == 1:
        return [1.0]
    cur = 0.0
    for r in previous:
        cur += r
    rates = mean_rates(shape, ncat)
    new = 0.0
    for r in rates:
        new += r
    return [r * (cur / new) for r in rates] if new != cur else rates


def category_rates(shape, ncat, p_invar=0.0):
    """Rates and weights for +G (RateGamma) or +I+G (RateGammaInvar) with a fixed shape.

    The constructors set every rate to 1, or to 1/(1-p) under +I+G, before computeRates
    (rategamma.cpp:31, 74-82; rategammainvar.cpp:24-35); the weights are 1/ncat and
    (1-p)/ncat (rategamma.h:120, rategammainvar.h:69)."""
    shape = max(MIN_GAMMA_SHAPE, abs(shape))
    rates = _compute_rates([1.0] * ncat, shape, ncat)
    if not p_invar:
        return np.array(rates), np.full(ncat, 1.0 / ncat)
    rates = _compute_rates([1.0 / (1.0 - p_invar)] * ncat, shape, ncat)
    return np.array(rates), np.full(ncat, (1.0 - p_invar) / ncat)


def invariant_rates(p_invar):
    """+I alone: one category with rate 1/(1-p) and weight 1-p (rateinvar.h:77, 84)."""
    return np.array([1.0 / (1.0 - p_invar)]), np.array([1.0 - p_invar])
