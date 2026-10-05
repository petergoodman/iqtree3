"""IQ-TREE's BFGS path, ported from utils/optimization.cpp at 63c330d9 (decision 017), and the
scaled derivative step of decision 005.

The algorithm, constants and control flow follow the C++, including its quirks; arrays are
0-indexed here where the C++ is 1-indexed. Floating-point identity with the C++ is not expected,
because NumPy and the C++ loops sum in different orders. Arithmetic runs on NumPy floats so that
division by zero gives IEEE infinities, as in the C++, instead of raising.
"""
from dataclasses import dataclass

import numpy as np

ERROR_X = 1.0e-4                    # utils/optimization.cpp:23
ALF, TOLX_LNSRCH = 1.0e-4, 1.0e-7   # 639-640
MAX_ITER = 3                        # 723
ITMAX, EPS, STPMX = 200, 3.0e-8, 100.0  # 781-785
TOLX_DFPMIN = 4 * EPS
MIN_RATE, TOL_RATE, MAX_RATE = 1e-4, 1e-4, 100.0  # model/modelmarkov.h:30-32


@dataclass
class Result:
    x: np.ndarray
    f: float
    iterations: int
    reason: str            # "TOLX", "gtol" or "ITMAX"; the C++ does not report it
    failed_searches: int   # line searches that set check = 1; the C++ never reads the flag
    runs: int = 1          # dfpmin calls made by minimize_multi_dimen


@dataclass
class ModelFit:
    x: np.ndarray
    score: float
    gtol: float
    result: Result


def fix_bound(x, lower, upper):
    """Clamp x in place (149-156)."""
    np.clip(x, lower, upper, out=x)


def legacy_derivative(f, x):
    """Forward difference with step 1e-4|x|, or 1e-4 when that is 0, and no bound handling
    (916-939). Returns f(x) and the gradient; x is restored."""
    fx = np.float64(f(x))
    n = len(x)
    h, df = np.empty(n), np.empty(n)
    for d in range(n):
        temp = x[d]
        h[d] = ERROR_X * abs(temp)
        if h[d] == 0.0:
            h[d] = ERROR_X
        x[d] = temp + h[d]
        h[d] = x[d] - temp
        df[d] = f(x)
        x[d] = temp
    return fx, (df - fx) / h


def scaled_step(x, upper, eta=ERROR_X, scale=1.0):
    """Decision 005: h = eta * max(scale, |x|), exactly representable, backward when x + h would
    leave the upper bound."""
    h = eta * max(scale, abs(x))
    if x + h > upper:
        h = -h
    return (x + h) - x


def scaled_derivative(f, x, upper, eta=ERROR_X, scale=1.0):
    """The legacy loop with decision 005's signed step."""
    fx = np.float64(f(x))
    upper = np.broadcast_to(np.asarray(upper, dtype=float), x.shape)
    g = np.empty(len(x))
    for d in range(len(x)):
        temp = x[d]
        h = scaled_step(temp, upper[d], eta, scale)
        x[d] = temp + h
        g[d] = (f(x) - fx) / h
        x[d] = temp
    return fx, g


def central_derivative(f, x, steps):
    """Central differences with h = s * max(1, |x_k|) for each s in steps; validation only."""
    out = np.empty((len(steps), len(x)))
    for a, s in enumerate(steps):
        for d in range(len(x)):
            temp = x[d]
            h = s * max(1.0, abs(temp))
            x[d] = temp + h
            up = f(x)
            x[d] = temp - h
            down = f(x)
            x[d] = temp
            out[a, d] = (up - down) / (2 * h)
    return out


def lnsrch(f, xold, fold, g, p, stpmax, lower, upper):
    """Backtracking line search (645-718). Returns (x, f, check). On failure x is xold and f is
    the value at the last, rejected trial (687-690)."""
    with np.errstate(divide="ignore", invalid="ignore"):
        p = np.array(p, dtype=float)
        total = np.sqrt(p @ p)
        if total > stpmax:
            p *= stpmax / total
        slope = np.float64(g @ p)
        test = np.max(np.abs(p) / np.maximum(np.abs(xold), 1.0))
        alamin = np.float64(TOLX_LNSRCH) / test
        alam, alam2, f2, fold2 = np.float64(1.0), 0.0, 0.0, 0.0
        first = True
        while True:
            x = xold + alam * p
            fix_bound(x, lower, upper)
            fv = np.float64(f(x))
            if alam < alamin:
                return np.array(xold, dtype=float), fv, 1
            if fv <= fold + ALF * alam * slope:
                return x, fv, 0
            if first:
                tmplam = -slope / (2.0 * (fv - fold - slope))
            else:
                rhs1 = fv - fold - alam * slope
                rhs2 = f2 - fold2 - alam2 * slope
                a = (rhs1 / (alam * alam) - rhs2 / (alam2 * alam2)) / (alam - alam2)
                b = (-alam2 * rhs1 / (alam * alam) + alam * rhs2 / (alam2 * alam2)) / (alam - alam2)
                if a == 0.0:
                    tmplam = -slope / (2.0 * b)
                else:
                    disc = b * b - 3.0 * a * slope
                    if disc < 0.0:
                        tmplam = 0.5 * alam
                    elif b <= 0.0:
                        tmplam = (-b + np.sqrt(disc)) / (3.0 * a)
                    else:
                        tmplam = -slope / (b + np.sqrt(disc))
                if tmplam > 0.5 * alam:
                    tmplam = 0.5 * alam
            alam2, f2, fold2 = alam, fv, fold
            alam = max(tmplam, 0.1 * alam)
            first = False


def update_inverse_hessian(H, dg, xi):
    """The BFGS update of dfpmin (863-889). It runs when fac^2 > EPS*sumdg*sumxi (875), so also
    for negative curvature, where Numerical Recipes tests fac > sqrt(EPS*sumdg*sumxi)."""
    with np.errstate(divide="ignore", invalid="ignore"):
        hdg = H @ dg
        fac, fae = dg @ xi, dg @ hdg
        if not fac * fac > EPS * (dg @ dg) * (xi @ xi):
            return H, False
        fac, fad = 1.0 / fac, 1.0 / fae
        u = fac * xi - fad * hdg
        return H + fac * np.outer(xi, xi) - fad * np.outer(hdg, hdg) + fae * np.outer(u, u), True


def dfpmin(f, p, lower, upper, gtol, derivative=legacy_derivative, hessian=None, itmax=ITMAX):
    """BFGS on the inverse Hessian (793-900). `derivative(f, x)` returns (f(x), gradient)."""
    p = np.array(p, dtype=float)
    n = len(p)
    fp, g = derivative(f, p)
    xi = -g
    H = np.eye(n) if hessian is None else np.array(hessian, dtype=float).reshape(n, n)
    stpmax = STPMX * max(np.sqrt(p @ p), float(n))
    failed = 0
    fret = fp
    for its in range(1, itmax + 1):
        pnew, fret, check = lnsrch(f, p, fp, g, xi, stpmax, lower, upper)
        failed += check
        fp = fret
        xi = pnew - p
        p = pnew
        if np.max(np.abs(xi) / np.maximum(np.abs(p), 1.0)) < TOLX_DFPMIN:
            return Result(p, fret, its, "TOLX", failed)
        dg = g.copy()
        _, g = derivative(f, p)
        den = max(abs(fret), 1.0)
        if np.max(np.abs(g) * np.maximum(np.abs(p), 1.0) / den) < gtol:
            return Result(p, fret, its, "gtol", failed)
        H, _ = update_inverse_hessian(H, g - dg, xi)
        xi = -(H @ g)
    return Result(p, fret, itmax, "ITMAX", failed)


def restart_parameters(guess, lower, upper, bound_check, iteration, random_double):
    """726-748: up to MAX_ITER restarts when a checked coordinate is within 1e-4 of a bound; every
    coordinate is then redrawn. IQ-TREE draws with rand(), so a caller must supply the draws."""
    if iteration > MAX_ITER:
        return False
    near = (np.abs(guess - lower) < 1e-4) | (np.abs(guess - upper) < 1e-4)
    if not np.any(np.asarray(bound_check, dtype=bool) & near):
        return False
    if random_double is None:
        raise ValueError("a restart needs random_double; IQ-TREE's rand() draws cannot be reproduced")
    for i in range(len(guess)):
        guess[i] = random_double() * (upper[i] - lower[i]) / 3 + lower[i]
    return True


def minimize_multi_dimen(f, guess, lower, upper, bound_check, gtol, hessian=None,
                         derivative=legacy_derivative, random_double=None):
    """750-778: dfpmin, repeated after each restart; the best result is returned."""
    guess = np.array(guess, dtype=float)
    minf, minx, count = 1e12, None, 0
    total_failed = 0
    while True:
        r = dfpmin(f, guess, lower, upper, gtol, derivative, hessian)
        guess, fret = r.x, r.f
        total_failed += r.failed_searches
        if fret < minf:
            minf, minx = fret, guess.copy()
        count += 1
        if not restart_parameters(guess, lower, upper, bound_check, count, random_double):
            break
    if count > 1:
        if minx is None:
            raise ValueError("every run returned at least 1e12; the C++ would read uninitialised memory")
        guess, fret = minx, minf
    return Result(guess, fret, r.iterations, r.reason, total_failed, count)


def optimize_model(f, x0, lower, upper, gradient_epsilon, derivative=legacy_derivative):
    """The call as ModelMarkov::optimizeParameters makes it (model/modelmarkov.cpp:1199-1235):
    gtol = max(gradient_epsilon, TOL_RATE), no bound checks, so no restarts, and the score is
    recomputed when the returned variables differ from the model's last evaluated state."""
    last = {}

    def target(x):
        last["x"] = np.array(x, dtype=float)
        return f(x)

    gtol = max(gradient_epsilon, TOL_RATE)
    r = minimize_multi_dimen(target, x0, lower, upper, np.zeros(len(x0), bool), gtol, derivative=derivative)
    score = -r.f
    if "x" not in last or not np.array_equal(last["x"], r.x) or score == -1.0e30:
        score = -np.float64(target(r.x))
    return ModelFit(r.x, score, gtol, r)
