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


# ---------------------------------------------------------------- a reference climber, not IQ-TREE's

@dataclass
class Climb:
    x: np.ndarray
    f: float
    evaluations: int
    phases: list           # (derivative, iterations, SciPy's message) per L-BFGS-B pass
    gradient: np.ndarray   # central differences at x
    pg_max: float          # the largest projected-gradient component of that gradient
    at_lower: np.ndarray   # indices exactly at a bound
    at_upper: np.ndarray
    analytic_pg: float = None  # the analytic pass's projected gradient at x (S5), else None
    gradients: int = 0         # analytic gradient calls


def projected_gradient(x, g, lower, upper):
    """x - clip(x - g, lower, upper): g in the interior, 0 where a bound blocks descent."""
    return x - np.clip(x - g, lower, upper)


def probes(x, steps):
    """One point per coordinate, x with that coordinate moved by its step."""
    pts = np.repeat(np.asarray(x, dtype=float)[None, :], len(x), axis=0)
    pts[np.arange(len(x)), np.arange(len(x))] += steps
    return list(pts)


def forward_gradient(evaluate, x, upper, fx):
    """scaled_derivative's gradient with every probe handed to evaluate at once."""
    h = np.array([scaled_step(v, u) for v, u in zip(x, upper)])
    return (np.array(evaluate(probes(x, h)), dtype=float) - fx) / h


def central_gradient(evaluate, x, step):
    """central_derivative's gradient for one step, with every probe handed to evaluate at once."""
    h = step * np.maximum(1.0, np.abs(x))
    values = np.array(evaluate(probes(x, h) + probes(x, -h)), dtype=float)
    return (values[:len(x)] - values[len(x):]) / (2 * h)


def lbfgsb_fit(f, x0, lower, upper, pgtol=1e-3, max_iter=2000, central_step=ERROR_X, evaluate=None,
               gradient=None, maxcor=20, max_seconds=None, gradient_map=None):
    """A test reference for where a fit's peak lies (PLAN.md risks 22 and 24), not IQ-TREE's method:
    SciPy's L-BFGS-B minimizing f in the box, stopped by its projected-gradient test at pgtol. The
    first pass uses decision 005's forward differences; when the central-difference gradient at its
    end still fails pgtol, a second pass uses central differences, whose error does not grow with
    the curvature as the forward step's does. Central probes may leave the box by one step, which
    the NQC chart allows (constraint 1). evaluate(points) returns f at each point, in order; a
    process pool's map spreads a gradient's probes, with the same values as the serial default.

    With gradient (x -> f, g analytic), one pass replaces both (S5, decision 042): L-BFGS-B with
    that gradient and maxcor memory, stopped by the analytic projected gradient at pgtol (ftol 0),
    by max_iter iterations or by max_seconds, then _newton_finish while still above pgtol;
    gradient_map(points) returns gradient at each point, in order (serially by default). The
    central-difference gradient is still taken at x."""
    from scipy.optimize import Bounds, minimize

    if gradient is not None:
        return _analytic_lbfgsb(f, gradient, x0, lower, upper, pgtol, max_iter, central_step, evaluate, maxcor,
                                max_seconds, gradient_map or (lambda pts: [gradient(p) for p in pts]))

    n = len(x0)
    lower = np.broadcast_to(np.asarray(lower, dtype=float), (n,)).copy()
    upper = np.broadcast_to(np.asarray(upper, dtype=float), (n,)).copy()
    evaluate = evaluate or (lambda pts: [f(p) for p in pts])
    count = [0]

    def counted(pts):
        count[0] += len(pts)
        return evaluate(pts)

    def forward(x):
        x = np.array(x, dtype=float)
        fx = np.float64(counted([x])[0])
        return fx, forward_gradient(counted, x, upper, fx)

    def central(x):
        x = np.array(x, dtype=float)
        return np.float64(counted([x])[0]), central_gradient(counted, x, central_step)

    x = np.clip(np.asarray(x0, dtype=float), lower, upper)
    phases, g = [], None
    for name, fun in (("forward", forward), ("central", central)):
        if g is not None and np.max(np.abs(projected_gradient(x, g, lower, upper))) <= pgtol:
            break
        r = minimize(fun, x, jac=True, method="L-BFGS-B", bounds=Bounds(lower, upper),
                     options={"maxiter": max_iter, "maxfun": 4 * max_iter, "gtol": pgtol, "ftol": 1e-15,
                              "maxcor": 20})
        x = np.clip(r.x, lower, upper)
        phases.append((name, int(r.nit), str(r.message)))
        _, g = central(x)
    fx = float(counted([x])[0])
    pg = projected_gradient(x, g, lower, upper)
    return Climb(x, fx, count[0], phases, g, float(np.max(np.abs(pg))), np.flatnonzero(x == lower),
                 np.flatnonzero(x == upper))


def _analytic_lbfgsb(f, gradient, x0, lower, upper, pgtol, max_iter, central_step, evaluate, maxcor, max_seconds,
                     gradient_map):
    import time

    from scipy.optimize import Bounds, minimize

    n = len(x0)
    lower = np.broadcast_to(np.asarray(lower, dtype=float), (n,)).copy()
    upper = np.broadcast_to(np.asarray(upper, dtype=float), (n,)).copy()
    evaluate = evaluate or (lambda pts: [f(p) for p in pts])
    calls = [0]
    deadline = None if max_seconds is None else time.monotonic() + max_seconds

    def fun(x):
        calls[0] += 1
        fx, g = gradient(np.array(x, dtype=float))
        return np.float64(fx), np.asarray(g, dtype=float)

    def stop(intermediate_result):
        if deadline is not None and time.monotonic() > deadline:
            raise StopIteration

    x = np.clip(np.asarray(x0, dtype=float), lower, upper)
    r = minimize(fun, x, jac=True, method="L-BFGS-B", bounds=Bounds(lower, upper), callback=stop,
                 options={"maxiter": max_iter, "maxfun": 4 * max_iter, "gtol": pgtol, "ftol": 0.0,
                          "maxcor": maxcor})
    x = np.clip(r.x, lower, upper)
    fx, ga = fun(x)
    phases, probes = [("analytic", int(r.nit), str(r.message))], 0
    if np.max(np.abs(projected_gradient(x, ga, lower, upper))) > pgtol:
        x, fx, ga, steps, ending, probes = _newton_finish(fun, gradient_map, x, fx, ga, lower, upper, pgtol, deadline)
        phases.append(("newton", steps, ending))
    g = central_gradient(evaluate, x, central_step)
    pg = projected_gradient(x, g, lower, upper)
    return Climb(x, float(fx), 2 * n + calls[0] + probes, phases, g,
                 float(np.max(np.abs(pg))), np.flatnonzero(x == lower), np.flatnonzero(x == upper),
                 float(np.max(np.abs(projected_gradient(x, ga, lower, upper)))), calls[0] + probes)


NEWTON_STEPS = 20   # steps of the analytic pass's projected Newton finish


def _newton_direction(H, g):
    """-H^-1 g by Cholesky, H + tau I with tau raised by 10 from 1e-12 max|diag H| while H is not
    positive definite, as flux.py's finish."""
    from scipy.linalg import LinAlgError, cho_factor, cho_solve

    if not len(g):
        return np.zeros(0)
    tau, scale = 0.0, max(1.0, float(np.max(np.abs(np.diag(H)))))
    while True:
        try:
            return -cho_solve(cho_factor(H + tau * np.eye(len(H))), g)
        except LinAlgError:
            tau = 1e-12 * scale if tau == 0.0 else 10 * tau


def _newton_finish(fun, gradient_map, x, fx, g, lower, upper, pgtol, deadline):
    """The analytic pass's finish (S5, decision 042): projected Newton in the box, for the tail where
    L-BFGS-B's steps gain at most a unit of f's spacing and its zero-reduction test ends the pass above
    pgtol (explore/z_climb_diagnosis.py). Each step takes the Hessian by central differences of the
    analytic gradient (step HESSIAN_STEP max(1, |x|)) on the free coordinates, those not at a bound
    the gradient pushes against; holds for that step a free coordinate at a bound that the Newton step
    would push out; cuts the step at the first bound it reaches (coordinates that reach one within
    ALPHA_MIN of the step are first put on it); and accepts as flux.py's finish: Armijo above
    ROUNDING |f|, below it no rise above that, with the decrement at least halving from one full step
    to the next. Ends at pgtol, a stall, the deadline or NEWTON_STEPS steps, keeping its entry point
    if lower when pgtol is not met. Returns x, f, g, steps, ending and the Hessian's gradients taken."""
    import time

    from .flux import ALPHA_MIN, ARMIJO, HESSIAN_STEP, ROUNDING

    entry, probes, last_pure, ending, step = (x, fx, g), 0, None, f"{NEWTON_STEPS} steps", 0
    while step < NEWTON_STEPS:
        if np.max(np.abs(projected_gradient(x, g, lower, upper))) <= pgtol:
            break
        if deadline is not None and time.monotonic() > deadline:
            ending = "time cap"
            break
        step += 1
        free = np.flatnonzero(~(((x <= lower) & (g > 0)) | ((x >= upper) & (g < 0))))
        h = HESSIAN_STEP * np.maximum(1.0, np.abs(x[free]))
        pts = []
        for k, hk in zip(free, h):
            for s in (hk, -hk):
                v = x.copy()
                v[k] += s
                pts.append(v)
        G = np.array([gk for _, gk in gradient_map(pts)])
        probes += len(pts)
        H = ((G[0::2] - G[1::2]) / (2 * h)[:, None])[:, free]
        H = (H + H.T) / 2
        keep = np.ones(len(free), dtype=bool)
        while True:
            sub = free[keep]
            d_sub = _newton_direction(H[np.ix_(keep, keep)], g[sub])
            out = ((x[sub] <= lower[sub]) & (d_sub < 0)) | ((x[sub] >= upper[sub]) & (d_sub > 0))
            if not out.any():
                break
            keep[np.flatnonzero(keep)[out]] = False
        d = np.zeros_like(x)
        d[sub] = d_sub
        ratio = np.full(len(x), np.inf)
        ratio[d < 0] = (lower - x)[d < 0] / d[d < 0]
        ratio[d > 0] = (upper - x)[d > 0] / d[d > 0]
        if np.min(ratio) < ALPHA_MIN:   # within rounding of a bound: put those coordinates on it
            x = np.where(ratio < ALPHA_MIN, np.where(d < 0, lower, upper), x)
            fx, g = fun(x)
            last_pure = None
            continue
        slope = float(g @ d)
        decrement = -slope / 2
        rounding = ROUNDING * max(1.0, abs(fx))
        pure = decrement <= rounding
        if decrement <= 0 or (pure and last_pure is not None and decrement > 0.5 * last_pure):
            ending = "stalled"
            break
        alpha = min(1.0, float(np.min(ratio)))
        while alpha >= ALPHA_MIN:
            xn = np.clip(x + alpha * d, lower, upper)
            fn, gn = fun(xn)
            if fn <= fx + (rounding if pure else ARMIJO * alpha * slope):
                break
            alpha *= 0.5
        else:
            ending = "no acceptable step"
            break
        last_pure = decrement if pure and alpha == 1.0 else None
        x, fx, g = xn, fn, gn
    if np.max(np.abs(projected_gradient(x, g, lower, upper))) <= pgtol:
        return x, fx, g, step, "reached pgtol", probes
    if entry[1] < fx:
        x, fx, g = entry
        ending += "; kept its entry point"
    return x, fx, g, step, ending, probes


# ---------------------------------------------------------------- SIM-2's stage B

LOW_FLUX_SHARE = 0.4   # the two lowest quintiles


@dataclass
class StageB:
    x: np.ndarray
    f: float
    f_start: float
    gain: float            # f_start - f
    iterations: int
    gradients: int         # calls of the gradient function
    seconds: float
    stop: str              # "stalled", "time cap", or SciPy's message
    message: str           # SciPy's message
    pg_max: float          # the largest analytic projected-gradient component at x
    at_lower: np.ndarray
    at_upper: np.ndarray
    window: dict           # the stall rule's three movements over the last window, and its span
    trace: list            # (iteration, seconds, f, C) per iterate, the start first


def low_flux(flux, share=LOW_FLUX_SHARE):
    """Indices of the rates in the lowest share of flux: 152 of 380 at the default."""
    return np.argsort(np.asarray(flux, dtype=float), kind="stable")[:int(round(share * len(flux)))]


def _window_moves(history):
    """The stall rule's movements from the oldest kept iterate to the newest: -f's rise, C's change and
    the largest log-rate change among the newest iterate's low-flux rates."""
    f_old, c_old, lr_old, _ = history[0]
    f_new, c_new, lr_new, fl_new = history[-1]
    low = low_flux(fl_new)
    return {"span": len(history) - 1, "lnl_rise": float(f_old - f_new), "c_move": float(abs(c_new - c_old)),
            "low_flux_move": float(np.max(np.abs(lr_new[low] - lr_old[low]))) if len(low) else 0.0}


def stage_b_climb(gradient, x0, lower, upper, watch, max_seconds, window=100, gain=0.1, c_move=1e-4,
                  low_move=0.01, pgtol=1e-6, maxcor=20, max_iter=10 ** 6, clock=None):
    """SIM-2's stage B (decision 039; CODE_PLAN.md 3.8): SciPy's L-BFGS-B on gradient (x -> f, analytic
    g) from x0 clipped into the box, ftol 0, memory maxcor. It ends by SciPy's own tests (pgtol, a failed
    line search, max_iter), by the time cap, checked after each iteration, or by the stall rule: between
    the current iterate and the one window iterations earlier, -f rises less than gain, C moves less than
    c_move, and no rate among the current iterate's low-flux rates moves by more than low_move in log.
    watch(x) returns C, the log rates and their fluxes. Returns the lowest iterate, never above the start;
    no Newton finish, no central differences."""
    import time
    from collections import deque

    from scipy.optimize import Bounds, minimize

    clock = clock or time.monotonic
    start = clock()
    n = len(x0)
    lower = np.broadcast_to(np.asarray(lower, dtype=float), (n,)).copy()
    upper = np.broadcast_to(np.asarray(upper, dtype=float), (n,)).copy()
    calls, last = [0], {}

    def fun(x):
        calls[0] += 1
        fx, g = gradient(np.array(x, dtype=float))
        last.update(x=np.array(x, dtype=float), f=np.float64(fx), g=np.asarray(g, dtype=float))
        return last["f"], last["g"]

    def observe(x, fx):
        c, lr, fl = watch(x)
        return float(fx), float(c), np.asarray(lr, dtype=float), np.asarray(fl, dtype=float)

    x = np.clip(np.asarray(x0, dtype=float), lower, upper)
    f_start, g = fun(x)
    best = {"x": x, "f": f_start, "g": g}
    history = deque([observe(x, f_start)], maxlen=window + 1)
    trace = [(0, clock() - start, float(f_start), history[0][1])]
    state = {"stop": None}

    def callback(intermediate_result):
        xk, fk = np.array(intermediate_result.x, dtype=float), np.float64(intermediate_result.fun)
        now = clock() - start
        history.append(observe(xk, fk))
        trace.append((len(trace), now, float(fk), history[-1][1]))
        if fk < best["f"]:
            best.update(x=xk, f=fk, g=last["g"] if np.array_equal(last.get("x"), xk) else None)
        if len(history) == window + 1:
            m = _window_moves(history)
            if m["lnl_rise"] < gain and m["c_move"] < c_move and m["low_flux_move"] <= low_move:
                state["stop"] = "stalled"
                raise StopIteration
        if now > max_seconds:
            state["stop"] = "time cap"
            raise StopIteration

    r = minimize(fun, x, jac=True, method="L-BFGS-B", bounds=Bounds(lower, upper), callback=callback,
                 options={"maxiter": max_iter, "maxfun": 4 * max_iter, "gtol": pgtol, "ftol": 0.0,
                          "maxcor": maxcor})
    if best["g"] is None:
        best["f"], best["g"] = fun(best["x"])
    xb = best["x"]
    return StageB(xb, float(best["f"]), float(f_start), float(f_start - best["f"]), len(trace) - 1, calls[0],
                  clock() - start, state["stop"] or str(r.message), str(r.message),
                  float(np.max(np.abs(projected_gradient(xb, best["g"], lower, upper)))),
                  np.flatnonzero(xb == lower), np.flatnonzero(xb == upper), _window_moves(history), trace)
