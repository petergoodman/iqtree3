"""The direct balanced-flux reference of S5 (D10; decision 042; CODE_PLAN.md 3.1).

Fluxes f_ij = pi_i q_ij on the n(n - 1) ordered pairs, in t3._edges order. Balance and unit mean
rate are Cf = b with C = [D; 1^T] and b = e_n; the polytope does not depend on pi, which enters only
through q_ij = f_ij / pi_i. A fit maximizes an objective of f (a log-likelihood through Q, or the
concave validation objective) by a log-barrier Newton method in a re-centred scaled basis
f = f_c (1 + N_c u), N_c spanning ker(C diag(f_c)), restarted at its endpoint until a restart gains
less than GAIN_STOP. The matched domain is NQC's box written exactly as the linear inequalities
e^lo f_i,r(i) <= f_ij <= e^hi f_i,r(i); the open domain is f > 0.

The fitting path uses only t3._edges, cases.generator_from_flux, the objective, NumPy and SciPy's
linear algebra. The references and the chart enter only through the matched domain, the start's
inward move and the stopping test's projected z-gradient.
"""
import math
import time
from collections import OrderedDict
from dataclasses import dataclass, field

import numpy as np
from scipy.linalg import LinAlgError, cho_factor, cho_solve, cholesky, lstsq, null_space, solve_triangular

from . import cases, chart, nqc, t3
from .optimize import projected_gradient

MU_STOP = 1e-10        # the declared stop (decision 042): barrier parameter at most this
PG_STOP = 1e-6         # and the analytic projected z-gradient at most this
GAIN_STOP = 1e-4       # a restart gaining less ends the fit
MAX_RESTARTS = 10
MAX_ITER = 50_000      # Newton steps per fit, over its restarts
MAX_SECONDS = 4 * 3600.0
INWARD = 1e-9          # a start's coordinates within this of a bound are moved this far inside
ACTIVE = 1e-6          # normalized slack at or below which a constraint counts in m (A5 iii)
MU0, MU_DECAY = 0.1, 0.1   # the barrier parameter of each run, from MU0 down to MU_STOP
TAPIA = 0.6            # the finish takes a row whose slack fell to at most this over the last level
NEAR = 1e-2            # and lies within this normalized slack (about this in z) of its bound
RELEASE = 1e-7         # the finish releases a face whose multiplier is below -RELEASE
FINISH_STEPS = 300     # Newton steps of the active-set finish
CLOSED = 1e-12         # a working row with a larger slack is still being closed
CENTRED = 0.1          # a subproblem ends when its Newton decrement lambda^2 / 2 <= CENTRED * mu
INNER_MAX = 100        # Newton steps per barrier parameter
TO_BOUNDARY = 0.995    # a step keeps this fraction of every slack
ARMIJO, ALPHA_MIN = 1e-4, 1e-12
ROUNDING = 1e-12       # a decrement below this times |phi| is under the objective's rounding
HESSIAN_STEP = 1e-6    # step in u of the Hessian's differences (forward; central in the finish)
RECENTRE = 0.5         # the basis is re-centred once a flux has moved this fraction from the centre
KKT_THRESHOLDS = (0.0, 1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.0, math.inf)


def geometry(n):
    i, j, D = t3._edges(n)
    b = np.zeros(n)
    b[-1] = 1.0
    return i, j, np.vstack([D, np.ones(len(i))]), b


def edge_index(n):
    """index[i, j] of the ordered pair (i, j) in t3._edges order, -1 on the diagonal."""
    i, j, _ = t3._edges(n)
    idx = np.full((n, n), -1)
    idx[i, j] = np.arange(len(i))
    return idx


def flux_of(Q, pi):
    Q, pi = np.asarray(Q, dtype=float), np.asarray(pi, dtype=float)
    i, j, _ = t3._edges(len(pi))
    return pi[i] * Q[i, j]


def generator(f, pi):
    """Q = f / pi off the diagonal, the diagonal from row sums, never renormalized. A non-finite or
    non-positive flux is refused before anything else."""
    f, pi = np.asarray(f, dtype=float), np.asarray(pi, dtype=float)
    if not np.all(np.isfinite(f)) or np.any(f <= 0):
        raise chart.BuildError("fluxes must be finite and strictly positive")
    n = len(pi)
    i, j, _ = t3._edges(n)
    F = np.zeros((n, n))
    F[i, j] = f
    return cases.generator_from_flux(F, pi)


def residuals(f, pi):
    d = chart.diagnostics(generator(f, pi), pi)
    return {k: d[k] for k in ("stationarity", "row_sums", "mean_rate_error")}


def canonical(n):
    """f0 uniform and Z an orthonormal basis of ker C: the canonical chart f = f0 + Zy."""
    return np.full(n * (n - 1), 1.0 / (n * (n - 1))), null_space(geometry(n)[2])


def z_of_flux(f, ref):
    """z_ij = log(f_ij / f_i,r(i)), in the chart's coordinate order (K_ij / K_ir = f_ij / f_ir)."""
    idx, f = edge_index(len(ref)), np.asarray(f, dtype=float)
    return np.concatenate([np.log(f[idx[i, free]] / f[idx[i, ref[i]]])
                           for i, free in enumerate(chart.free_destinations(ref))])


@dataclass
class Domain:
    """Linear inequalities A f >= 0. Matched: rows f_ij - e^lo f_ir, then rows e^hi f_ir - f_ij, in
    the chart's coordinate order, coef the bound's coefficient and ref_edge the reference flux of
    each row. Open: A = I."""
    kind: str
    ref: np.ndarray
    lower: float
    upper: float
    A: np.ndarray
    coef: np.ndarray
    ref_edge: np.ndarray

    def scale(self, f):
        """Each row's unit: coef f_i,r(i) (matched); 1 (open), so a flux's slack is itself."""
        return self.coef * f[self.ref_edge] if self.kind == "matched" else np.ones(len(self.A))

    def slack(self, f):
        return (self.A @ f) / self.scale(f)

    def contains(self, f):
        return bool(np.all(np.isfinite(f)) and np.all(f > 0) and np.all(self.A @ f > 0))


def matched(ref, lower=nqc.LOWER, upper=nqc.UPPER):
    ref = np.asarray(ref)
    n = len(ref)
    idx = edge_index(n)
    E = n * (n - 1)
    lo_rows, hi_rows, refs = [], [], []
    for i, free in enumerate(chart.free_destinations(ref)):
        r = idx[i, ref[i]]
        for j in free:
            row = np.zeros(E)
            row[idx[i, j]], row[r] = 1.0, -math.exp(lower)
            lo_rows.append(row)
            row = np.zeros(E)
            row[idx[i, j]], row[r] = -1.0, math.exp(upper)
            hi_rows.append(row)
            refs.append(r)
    m = len(refs)
    return Domain("matched", ref, lower, upper, np.array(lo_rows + hi_rows),
                  np.r_[np.full(m, math.exp(lower)), np.full(m, math.exp(upper))], np.array(refs + refs))


def open_domain(ref):
    E = len(ref) * (len(ref) - 1)
    return Domain("open", np.asarray(ref), -math.inf, math.inf, np.eye(E), np.ones(E), np.arange(E))


def inward(Q, pi, domain):
    """A strictly feasible start: coordinates within INWARD of a bound, or beyond it, set INWARD
    inside and the flux rebuilt by the chart; otherwise Q's own flux. Returns the flux and the
    number of coordinates moved."""
    return inward_flux(flux_of(Q, pi), domain, pi)


def inward_flux(f, domain, pi=None):
    """inward for a flux; the chart's flux nu_i K_ij does not depend on pi (uniform by default)."""
    if domain.kind == "open":
        return f, 0
    z = z_of_flux(f, domain.ref)
    lo, hi = domain.lower + INWARD, domain.upper - INWARD
    moved = (z < lo) | (z > hi)
    if not moved.any():
        return f, 0
    pi = np.full(len(domain.ref), 1.0 / len(domain.ref)) if pi is None else pi
    return flux_of(chart.build(np.clip(z, lo, hi), pi, domain.ref)[0], pi), int(moved.sum())


class Basis:
    """f = f_c (1 + N u), N an orthonormal basis of ker(C diag(f_c)); the domain's rows in u,
    G u >= h, each divided by its unit at the centre (the open domain's by f_c)."""

    def __init__(self, fc, domain):
        n = len(domain.ref)
        self.fc = np.array(fc, dtype=float)
        self.N = null_space(geometry(n)[2] * self.fc[None, :])
        if self.N.shape[1] != n * (n - 2):
            raise ValueError(f"ker(C diag(f_c)) has dimension {self.N.shape[1]}, not {n * (n - 2)}")
        unit = domain.scale(self.fc) if domain.kind == "matched" else self.fc
        self.G = (domain.A * self.fc[None, :]) @ self.N / unit[:, None]
        self.h = -(domain.A @ self.fc) / unit

    def flux(self, u):
        return self.fc * (1.0 + self.N @ u)

    def coordinates(self, f):
        return self.N.T @ (np.asarray(f, dtype=float) / self.fc - 1.0)

    def chain(self, gf):
        """d/du from d/df."""
        return self.N.T @ (self.fc * gf)

    def slack(self, u):
        return self.G @ u - self.h


class LikelihoodFlux:
    """f -> (lnL, dlnL/df) for an objective with log_likelihood_gradient(Q) (nqc.Objective,
    nqc.LinkedObjective): dlnL/df_ij = (G_ij - G_ii) / pi_i."""

    def __init__(self, objective, pi):
        self.objective, self.pi = objective, np.asarray(pi, dtype=float)
        self.i, self.j, _ = t3._edges(len(self.pi))

    def __call__(self, f):
        lnL, G = self.objective.log_likelihood_gradient(generator(f, self.pi))
        return lnL, (G[self.i, self.j] - G[self.i, self.i]) / self.pi[self.i]


class Concave:
    """V1 and V1b: sum_e N_e log f_e."""

    def __init__(self, N):
        self.N = np.asarray(N, dtype=float)

    def __call__(self, f):
        return float(self.N @ np.log(f)), self.N / f


def z_gradient(f, gf, ref):
    """d(objective)/dz at a flux, by nqc.z_gradient with G_ij = pi_i df_ij / dq_ij (pi cancels,
    since the chart's flux nu_i K_ij does not depend on it; uniform here)."""
    n = len(ref)
    pi = np.full(n, 1.0 / n)
    i, j, _ = t3._edges(n)
    G = np.zeros((n, n))
    G[i, j] = pi[i] * gf
    return nqc.z_gradient(z_of_flux(f, ref), pi, ref, G)


def projected_z_gradient(f, gf, domain):
    """max |projected gradient| of -objective in z, against the box (none in the open domain)."""
    z = z_of_flux(f, domain.ref)
    pg = projected_gradient(z, -z_gradient(f, gf, domain.ref), domain.lower, domain.upper)
    return float(np.max(np.abs(pg)))


def _multipliers(f, gf, domain, threshold, divisor=None):
    """eta and nu >= 0 by bounded least squares on d objective / df = C^T eta - A_act^T nu, the
    active rows those with normalized slack at most threshold, each equation divided by divisor."""
    from scipy.optimize import lsq_linear

    C = geometry(len(domain.ref))[2]
    act = np.flatnonzero(domain.slack(f) <= threshold)
    M = np.hstack([C.T, -domain.A[act].T])
    lb = np.r_[np.full(len(C), -np.inf), np.zeros(len(act))]
    w = np.ones(len(gf)) if divisor is None else divisor
    x = lsq_linear(M / w[:, None], gf / w, bounds=(lb, np.full(len(lb), np.inf)), method="bvls").x
    return act, x[:len(C)], x[len(C):], M @ x


def kkt(f, gf, domain):
    """The KKT residual at f, recorded per fit: r = d objective / df - C^T eta + A_act^T nu with the
    multipliers of _multipliers, for the threshold in KKT_THRESHOLDS giving the smallest ||r||_inf.
    A record, not a test: the "all rows" fit drops complementarity (risk 27's probe: 7.6e-7 at an
    NQC endpoint 1.48 below aa_example's peak); certificate is the test."""
    best = None
    for t in KKT_THRESHOLDS:
        act, eta, nu, fitted = _multipliers(f, gf, domain, t)
        r = float(np.max(np.abs(gf - fitted)))
        if best is None or r < best["residual"]:
            best = {"residual": r, "threshold": t if math.isfinite(t) else "all rows", "active": int(len(act)),
                    "complementarity": float(nu @ (domain.A[act] @ f)), "gradient_norm": float(np.max(np.abs(gf)))}
    return best


def certificate(f, gf, domain, threshold=1e-12):
    """Decision 039's flux-space first-order test (line in 021): max |d objective / df - C^T eta +
    A_act^T nu|, nu >= 0, the active rows fixed as those with normalized slack <= threshold. Unlike a
    z-gradient, it does not shrink with a flux, so it fails where a near-zero flux should grow
    (PLAN.md risk 27, retired)."""
    act, _, nu, fitted = _multipliers(f, gf, domain, threshold)
    return {"residual": float(np.max(np.abs(gf - fitted))), "active": int(len(act)),
            "complementarity": float(nu @ (domain.A[act] @ f))}


def duality_gap(f, N, domain, mu=None):
    """V1 and V1b's certificate (decision 021): D(eta, nu) - F(f) for F = sum N_e log f_e, with
    D = b^T eta + sum_e N_e (log(N_e / w_e) - 1) and w = C^T eta - A^T nu > 0, nu >= 0. By weak duality
    D bounds max F over the domain for any such pair, so D - F(f) bounds f's optimality gap, which
    near the optimum is about nu . A f + sum_e N_e rho_e^2 / 2, rho = w / gf - 1. Pairs tried: the
    barrier's multipliers nu = mu / (A f) with eta by least squares (given mu), and bounded least
    squares at each threshold in KKT_THRESHOLDS, both with equation e divided by sqrt(N_e) / f_e so
    that they minimize the second term; the smallest bound is kept. Returns inf when no pair gives
    w > 0."""
    from scipy.optimize import lsq_linear

    N, f = np.asarray(N, dtype=float), np.asarray(f, dtype=float)
    F = float(N @ np.log(f))
    C, b = geometry(len(domain.ref))[2:]
    gf, divisor = N / f, np.sqrt(N) / f
    pairs = [_multipliers(f, gf, domain, t, divisor)[1:] for t in KKT_THRESHOLDS]
    if mu and np.all(domain.A @ f > 0):
        nu = mu / (domain.A @ f)
        eta = lsq_linear(C.T / divisor[:, None], (gf + domain.A.T @ nu) / divisor).x
        pairs.append((eta, nu, C.T @ eta - domain.A.T @ nu))
    best = math.inf
    for eta, _, w in pairs:
        if np.all(w > 0):
            best = min(best, float(b @ eta + N @ (np.log(N / w) - 1.0)) - F)
    return best


def dual_newton(N, n, tol=1e-12, max_iter=200):
    """V1's open optimum by the technical plan's normalized dual (section 9.2, d = 0): f_e =
    N_e / (C^T lambda)_e, lambda minimizing b^T lambda - sum N_e log (C^T lambda)_e by damped Newton."""
    N = np.asarray(N, dtype=float)
    C, b = geometry(n)[2:]
    lam = np.r_[np.zeros(n - 1), N.sum()]

    def dual(v):
        den = C.T @ v
        return b @ v - N @ np.log(den) if np.all(den > 0) else math.inf

    for _ in range(max_iter):
        den = C.T @ lam
        f = N / den
        g = b - C @ f
        if np.max(np.abs(g)) < tol:
            return f
        step = np.linalg.solve((C * (N / den ** 2)) @ C.T, -g)
        alpha = 1.0
        while dual(lam + alpha * step) > dual(lam) + 1e-4 * alpha * (g @ step) + 1e-13:
            alpha *= 0.5
            if alpha < 1e-14:
                raise RuntimeError("dual Newton: line search stalled")
        lam += alpha * step
    raise RuntimeError("dual Newton: no convergence")


def gradient_check(fun, x, g, directions, h=1e-4):
    """Decision 021's S5 gradient row along each direction d: |g.d - D(h)| <= 2 [|D(2h) - D(h)| / 3
    + 2e-12 |F| / h], D(h) the central difference of fun along d. Returns (error, line) pairs."""
    F = fun(x)
    out = []
    for d in directions:
        d1 = (fun(x + h * d) - fun(x - h * d)) / (2 * h)
        d2 = (fun(x + 2 * h * d) - fun(x - 2 * h * d)) / (4 * h)
        out.append((float(abs(g @ d - d1)), float(2 * (abs(d2 - d1) / 3 + 2e-12 * abs(F) / h))))
    return out


@dataclass
class FluxFit:
    status: str              # converged, unconverged or failed
    message: str
    f: np.ndarray            # the last accepted point
    value: float             # the objective there
    start_value: float
    runs: list = field(default_factory=list)
    seconds: float = 0.0
    mu: float = math.nan
    pg: float = math.nan     # analytic projected z-gradient at f
    m_active: int = 0
    kkt: dict = None

    @property
    def iterations(self):
        return sum(r["iterations"] for r in self.runs)

    def totals(self):
        keys = ("iterations", "evaluations", "hessian_gradients", "backtracks", "infeasible", "stalls", "damped",
                "recentres", "pure")
        return {k: sum(r[k] for r in self.runs) for k in keys}


def difference_hessian(gradients, basis, u, g0, central=False):
    """The Hessian of -objective in u, symmetrized, from the analytic gradient: forward differences
    against g0 (its gradient at u), or central differences, which the active-set finish takes since
    forward ones are indefinite along directions the data barely constrain and stalled V3
    (explore/v3_diagnosis.py). Each probe moves HESSIAN_STEP, shrunk by 10 until its fluxes stay
    positive (a probe may leave the box by up to the step). gradients(list of f) returns
    [(value, d value / df)] in order. Returns H and the number of gradients taken."""
    n = len(u)
    probes, steps = [], []
    for k in range(n):
        step = HESSIAN_STEP
        while True:
            v = u.copy()
            v[k] += step
            fv = basis.flux(v)
            if central:
                w = u.copy()
                w[k] -= step
                fw = basis.flux(w)
            if np.all(fv > 0) and (not central or np.all(fw > 0)):
                break
            step *= 0.1
            if step < 1e-15:
                raise RuntimeError("no Hessian probe keeps the fluxes positive")
        probes += [fv, fw] if central else [fv]
        steps.append(v[k] - w[k] if central else v[k] - u[k])
    g = [-basis.chain(np.asarray(gf, dtype=float)) for _, gf in gradients(probes)]
    H = np.empty((n, n))
    for k in range(n):
        H[:, k] = (g[2 * k] - g[2 * k + 1]) / steps[k] if central else (g[k] - g0) / steps[k]
    return (H + H.T) / 2, len(probes)


def fit(objective, f_start, domain, objective_map=None, max_iter=MAX_ITER, max_seconds=MAX_SECONDS, progress=None):
    """Maximize objective(f) -> (value, d value / df) over the domain from a strictly feasible
    f_start. Each run is a log-barrier Newton method from u = 0 at its centre, barrier parameter
    from MU0 down by MU_DECAY; the declared stop is the barrier parameter at most MU_STOP with the
    analytic projected z-gradient at most PG_STOP, checked after every Newton step. Converged: that
    stop, on a restart that gained less than GAIN_STOP over the previous run's endpoint (each restart
    records its inward move's lnL change). Any other ending is unconverged (caps, the active-set finish
    without the stop, the open domain's barrier floor, restarts exhausted) or failed (an exception).
    objective_map(list of f) -> results spreads the Hessian's probes, serially by default."""
    t0 = time.monotonic()
    deadline = t0 + max_seconds
    f = np.array(f_start, dtype=float)
    if not domain.contains(f):
        raise ValueError("the start is not strictly inside the domain")
    out = FluxFit("unconverged", f"{MAX_RESTARTS} restarts each gained at least {GAIN_STOP}", f, math.nan, math.nan)
    prev_f = prev_value = None
    for restart in range(MAX_RESTARTS + 1):
        moved = 0
        if restart:
            f, moved = inward_flux(f, domain)   # an active-set finish ends exactly on its faces
        if not domain.contains(f):
            out.status, out.message = "failed", "a restart centre is not strictly inside the domain"
            break
        try:
            basis = Basis(f, domain)
        except Exception as e:  # noqa: BLE001  a failed fit is recorded, not raised
            out.status, out.message = "failed", f"{type(e).__name__}: {e}"
            break
        run = _run(objective, basis, domain, objective_map, max_iter - out.iterations, deadline, progress)
        run["moved"] = moved
        if restart:   # the inward move's lnL change, and the gain over the previous endpoint, not the moved centre
            run["inward_change"] = run["value_start"] - prev_value
            run["gain"] = run["value_end"] - prev_value
        out.runs.append(run)
        if restart == 0:
            out.start_value = run["value_start"]
        if run["ending"] != "converged":
            f = run["f"]
            out.status = "failed" if run["ending"] == "failed" else "unconverged"
            out.message = f"{run['ending']}: {run['message']}"
            break
        if restart > 0 and run["gain"] < -GAIN_STOP:   # the restart ended lower: keep the higher endpoint
            f = prev_f
            out.status, out.message = "unconverged", f"restart {restart} ended {-run['gain']:.3g} lower"
            break
        f = prev_f = run["f"]
        prev_value = run["value_end"]
        if restart > 0 and run["gain"] < GAIN_STOP:
            out.status, out.message = "converged", f"declared stop; restart {restart} gained {run['gain']:.3g}"
            break
    out.f = f
    out.seconds = time.monotonic() - t0
    out.mu = out.runs[-1]["mu"] if out.runs else math.nan
    try:
        value, gf = objective(f)
        out.value = float(value)
        out.pg = projected_z_gradient(f, gf, domain)
        out.kkt = kkt(f, gf, domain)
        out.m_active = int(np.sum(domain.slack(f) <= ACTIVE))
    except Exception as e:  # noqa: BLE001
        out.status, out.message = "failed", f"{out.message}; endpoint not evaluable: {type(e).__name__}: {e}"
    return out


def _run(objective, basis, domain, objective_map, max_iter, deadline, progress=None):
    rec = {"evaluations": 0, "hessian_gradients": 0, "backtracks": 0, "infeasible": 0, "stalls": 0, "damped": 0,
           "recentres": 0, "pure": 0, "iterations": 0, "mu": MU0, "pg": math.nan, "ending": None, "message": "", "value_start": math.nan,
           "value_end": math.nan, "gain": math.nan}
    t0 = time.monotonic()
    cache = OrderedDict()

    def evaluate(u, face=None):
        """(F, dF/du, f, value, gf) with F = -value; None outside the domain. With face, the rows
        it marks are held at zero slack by the caller and only the others must be positive."""
        key = u.tobytes()
        if key in cache:
            return cache[key]
        res = None
        s = basis.slack(u)
        if np.all(s > 0) if face is None else np.all(s[~face] > 0):
            f = basis.flux(u)
            if np.all(np.isfinite(f)) and np.all(f > 0):
                rec["evaluations"] += 1
                value, gf = objective(f)
                if math.isfinite(value) and np.all(np.isfinite(gf)):
                    res = (-float(value), -basis.chain(gf), f, float(value), np.asarray(gf, dtype=float))
        cache[key] = res
        if len(cache) > 4:
            cache.popitem(last=False)
        return res

    last_hessian = [None, None]

    def hessian(u, g0, central=False):
        """difference_hessian, kept for reuse while u and the kind of difference are unchanged."""
        key = (u.tobytes(), central)
        if last_hessian[0] == key:
            return last_hessian[1]
        H, count = difference_hessian(
            lambda fs: objective_map(fs) if objective_map else [objective(p) for p in fs], basis, u, g0, central)
        rec["hessian_gradients"] += count
        last_hessian[:] = [key, H]
        return H

    def direction(HF, gF, s, mu):
        """Newton direction of the barrier subproblem, (HF + mu G^T S^-2 G) d = -(gF - mu G^T S^-1 1),
        solved as the least-squares problem [L^T; sqrt(mu) S^-1 G] d = [-L^-1 gF; sqrt(mu) 1] with
        HF + tau I = L L^T, whose condition is the square root of the normal equations' (the barrier
        terms of binding rows reach 1e10 and more). tau is raised by 10 from 1e-12 max|diag HF| only
        when HF is not positive definite (the likelihood is not concave)."""
        tau, scale = 0.0, max(1.0, float(np.max(np.abs(np.diag(HF)))))
        while True:
            try:
                L = cholesky(HF + tau * np.eye(len(HF)), lower=True)
                break
            except LinAlgError:
                tau = 1e-12 * scale if tau == 0.0 else 10 * tau
        if tau:
            rec["damped"] += 1
        root = math.sqrt(mu)
        A = np.vstack([L.T, (root / s)[:, None] * basis.G])
        c = np.r_[-solve_triangular(L, gF, lower=True), np.full(len(s), root)]
        return lstsq(A, c)[0]

    def capped():
        if rec["iterations"] >= max_iter:
            rec["ending"], rec["message"] = "iteration cap", f"{MAX_ITER} Newton steps per fit"
        elif time.monotonic() > deadline:
            rec["ending"], rec["message"] = "time cap", f"{MAX_SECONDS:.0f} s per fit"
        return rec["ending"] is not None

    def recentre(f):
        """The finish's re-centring at f, as the barrier's: a new basis, u = 0, caches dropped."""
        nonlocal basis
        basis = Basis(f, domain)
        cache.clear()
        last_hessian[:] = [None, None]
        rec["recentres"] += 1
        return np.zeros(basis.N.shape[1])

    def finish(u, cur, level_slack):
        """The active-set finish (decision 042, Peter 2026-11-08), a primal active-set Newton method
        entered when the barrier subproblem at MU_STOP ends without the declared stop (never in the
        open domain). The working set starts as the rows within ACTIVE of their bound, or within NEAR
        of it with a slack that fell to at most TAPIA of its value one barrier level earlier (binding
        faces fall about tenfold, faces with zero multipliers about threefold). Each step solves the
        Newton problem, its Hessian by central differences, with the working rows as equalities,
        closing their slacks in the same step (the closing's effect on the other coordinates is part
        of the step); the rows that block a step join the set, with no step when one is at its bound;
        at a stationary point the rows whose multiplier is below -RELEASE leave it. The basis is
        re-centred as in the barrier. No barrier: mu is recorded as 0, and the declared stop is tested
        after every step. Without that stop, the finish never ends below the point it started from."""
        rec["mu"], rec["finish_steps"], rec["released"] = 0.0, 0, 0
        slack = level_slack[-1]
        falling = (slack <= TAPIA * level_slack[-2]) & (slack <= NEAR) if len(level_slack) > 1 else False
        W = (slack <= ACTIVE) | falling
        cache.clear()
        cur = evaluate(u, W)
        entry = cur
        last_pure = None
        for _ in range(FINISH_STEPS):
            if declared_stop(cur, 0.0):
                rec["ending"] = "converged"
                return u, cur
            if capped():
                break
            s, g = basis.slack(u), cur[1]
            HF = hessian(u, g, central=True)
            GW = basis.G[W]
            dp = lstsq(GW, -s[W])[0] if W.any() else np.zeros(len(u))
            P = null_space(GW) if W.any() else np.eye(len(u))
            d = dp + P @ newton_step(P.T @ HF @ P, P.T @ (g + HF @ dp))
            rec["iterations"] += 1
            rec["finish_steps"] += 1
            rounding = ROUNDING * max(1.0, abs(cur[0]))
            closing = W.any() and float(np.max(s[W])) > CLOSED
            decrement = float(-(g @ d) / 2)
            Gd = basis.G @ d
            block = ~W & (Gd < 0)
            ratio = np.full(len(s), math.inf)
            ratio[block] = -s[block] / Gd[block]
            a_block = float(np.min(ratio))
            blockers = block & (ratio <= a_block * (1 + 1e-9))   # exact ties join together
            if a_block < ALPHA_MIN:   # a blocking row is at its bound: it joins with no step
                W = W | blockers
                last_pure = None
                cache.clear()
                continue
            if not closing:
                pure = decrement <= rounding
                stalled = pure and last_pure is not None and decrement > 0.5 * last_pure
                if pure:
                    last_pure = decrement
                    rec["pure"] += 1
                if stalled or decrement <= 1e-6 * rounding:
                    lam = lstsq(GW.T, HF @ d + g)[0] if W.any() else np.zeros(0)
                    out = np.flatnonzero(W)[lam < -RELEASE]
                    if not len(out):
                        break                     # stationary with every multiplier kept
                    W[out] = False
                    rec["released"] += len(out)
                    pin = W.copy()
                    pin[out] = True       # released rows moved INWARD off their faces, the others kept
                    u = u + lstsq(basis.G[pin], np.where(W[pin], 0.0, INWARD - s[pin]))[0]
                    last_pure = None
                    cache.clear()
                    cur = evaluate(u, W)
                    if cur is None:
                        raise RuntimeError("the active-set finish left the domain at a release")
                    continue
            alpha = min(1.0, a_block)
            while alpha >= ALPHA_MIN:
                Wn = W | blockers if alpha >= a_block else W.copy()   # the rows that block the step join
                nxt = evaluate(u + alpha * d, Wn)
                if nxt is None:
                    rec["infeasible"] += 1
                elif closing or nxt[0] <= (cur[0] + rounding if decrement <= rounding
                                           else cur[0] + ARMIJO * alpha * (g @ d)):
                    break
                alpha *= 0.5
                rec["backtracks"] += 1
            else:
                rec["stalls"] += 1
                break
            if Wn.sum() > W.sum() or alpha < 1.0:
                last_pure = None   # a cut step's decrement does not bound the next one's
            if Wn.sum() > W.sum():
                cache.clear()
            u, cur, W = u + alpha * d, nxt, Wn
            if float(np.max(np.abs(basis.N @ u))) > RECENTRE:
                u = recentre(basis.flux(u))
                cur = evaluate(u, W)
            if progress:
                progress({"step": rec["iterations"], "mu": 0.0, "decrement": decrement, "alpha": alpha,
                          "pure": decrement <= rounding, "value": cur[3], "recentres": rec["recentres"],
                          "pg": rec["pg"], "faces": int(W.sum()), "closing": closing})
        if declared_stop(cur, 0.0):
            rec["ending"] = "converged"
            return u, cur
        if rec["ending"] is None:
            rec["ending"], rec["message"] = "finish", f"the active-set finish ended without the declared stop " \
                                                       f"(projected z-gradient {rec['pg']:.3g}, {int(W.sum())} faces)"
        if entry[3] > cur[3]:
            rec["message"] += f"; kept its entry point, {entry[3] - cur[3]:.3g} higher"
            u = recentre(entry[2])
            cur = evaluate(u)
        return u, cur

    def newton_step(H, g):
        tau, scale = 0.0, max(1.0, float(np.max(np.abs(np.diag(H))))) if len(H) else 1.0
        while True:
            try:
                c = cho_factor(H + tau * np.eye(len(H)))
                break
            except LinAlgError:
                tau = 1e-12 * scale if tau == 0.0 else 10 * tau
        if tau:
            rec["damped"] += 1
        return -cho_solve(c, g)

    def declared_stop(cur, mu):
        if mu > MU_STOP * (1 + 1e-9):   # mu is 0.1 times 0.1 nine times, 1.0000000000000006e-10 at MU_STOP
            return False
        rec["pg"] = projected_z_gradient(cur[2], cur[4], domain)
        return rec["pg"] <= PG_STOP

    u = np.zeros(basis.N.shape[1])
    try:
        cur = evaluate(u)
        if cur is None:
            raise ValueError("the run's centre is outside the domain")
        rec["value_start"] = cur[3]
        mu = MU0
        level_slack = []
        while rec["ending"] is None:
            rec["mu"] = mu
            last_pure = None
            for _ in range(INNER_MAX):
                if rec["iterations"] >= max_iter:
                    rec["ending"], rec["message"] = "iteration cap", f"{MAX_ITER} Newton steps per fit"
                    break
                if time.monotonic() > deadline:
                    rec["ending"], rec["message"] = "time cap", f"{MAX_SECONDS:.0f} s per fit"
                    break
                s = basis.slack(u)
                g = cur[1] - mu * (basis.G.T @ (1.0 / s))
                d = direction(hessian(u, cur[1]), cur[1], s, mu)
                rec["iterations"] += 1
                decrement = float(-(g @ d) / 2)
                if decrement <= CENTRED * mu:
                    break
                Gd = basis.G @ d
                shrink = Gd < 0
                alpha = min(1.0, TO_BOUNDARY * float(np.min(-s[shrink] / Gd[shrink]))) if shrink.any() else 1.0
                phi0 = cur[0] - mu * float(np.sum(np.log(s)))
                # Below phi's rounding the value cannot judge a step and the gradient still can: the
                # Newton step is taken without the Armijo test, until it stops halving the decrement.
                pure = decrement <= ROUNDING * max(1.0, abs(phi0))
                if pure:
                    if last_pure is not None and decrement > 0.5 * last_pure:
                        rec["stalls"] += 1
                        break
                    last_pure = decrement
                    rec["pure"] += 1
                while alpha >= ALPHA_MIN:
                    v = u + alpha * d
                    nxt = evaluate(v)
                    if nxt is None:
                        rec["infeasible"] += 1
                    else:
                        phi = nxt[0] - mu * float(np.sum(np.log(basis.slack(v))))
                        if (phi <= phi0 + ROUNDING * max(1.0, abs(phi0)) if pure
                                else phi <= phi0 + ARMIJO * alpha * (g @ d)):
                            break
                    alpha *= 0.5
                    rec["backtracks"] += 1
                else:
                    rec["stalls"] += 1      # no decrease above rounding: the subproblem ends here
                    break
                if progress:
                    progress({"step": rec["iterations"], "mu": mu, "decrement": decrement, "alpha": alpha,
                              "pure": pure, "value": nxt[3], "recentres": rec["recentres"], "pg": rec["pg"]})
                u, cur = v, nxt
                if float(np.max(np.abs(basis.N @ u))) > RECENTRE:
                    basis = Basis(basis.flux(u), domain)
                    cache.clear()
                    last_hessian[:] = [None, None]
                    u = np.zeros(basis.N.shape[1])
                    cur = evaluate(u)
                    rec["recentres"] += 1
                if declared_stop(cur, mu):
                    rec["ending"] = "converged"
                    break
            if rec["ending"] is None and declared_stop(cur, mu):
                rec["ending"] = "converged"
            if rec["ending"] is None:
                level_slack.append(domain.slack(cur[2]))
                if mu <= MU_STOP * (1 + 1e-9):
                    if domain.kind == "open":   # its faces f = 0 lie outside it
                        rec["ending"], rec["message"] = "barrier floor", "the open domain has no active-set finish"
                    else:
                        u, cur = finish(u, cur, level_slack)
                else:
                    mu *= MU_DECAY
        if rec["ending"] == "converged":
            rec["message"] = f"mu {rec['mu']:.3g} <= {MU_STOP}, projected z-gradient {rec['pg']:.3g} <= {PG_STOP}"
    except Exception as e:  # noqa: BLE001  a failed fit is recorded, not raised
        rec["ending"], rec["message"] = "failed", f"{type(e).__name__}: {e}"
        cur = cache.get(u.tobytes()) or (None, None, basis.flux(u), math.nan, None)   # no second evaluation
    rec["f"], rec["value_end"] = cur[2], cur[3]
    rec["gain"] = rec["value_end"] - rec["value_start"]
    rec["seconds"] = time.monotonic() - t0
    return rec
