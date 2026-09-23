"""
Compare positive-ratio (C2) versus log-ratio (logit) coordinates for the
jump-chain chart of the fixed-pi nonreversible family, under a faithful
Python port of IQ-TREE 3's optimizer path (utils/optimization.cpp at commit
63c330d9): dfpmin (BFGS on the inverse Hessian), lnsrch with fixBound
clamping, derivativeFunk with h = 1e-4*|x| (1e-4 if x == 0), convergence
tests with max(|x|,1) scaling, gtol = 1e-4 (TOL_RATE).

Experiments
  E1  gradient accuracy of IQ-TREE's forward-difference rule on both charts,
      against a central-difference reference with a step-size plateau check,
      at n = 20 states on a rooted tree likelihood.
  E2  full BFGS runs from the same seed on identical physical feasible sets
      (ratio box [1e-4, 100] with pin 10  <->  log box [ln 1e-5, ln 10] with
      pin 0), n states small enough to afford several runs, comparing final
      log-likelihood, iterations, likelihood evaluations, bound hits, and
      the third variant "log coordinates + scale-aware step".

Usage: python chart_fd_bfgs_compare.py [--n 8] [--seed 20260922] [--ntips 8]
                                        [--nsites 2000]
Outputs JSON to stdout with provenance.
"""
import argparse, json, sys, time, subprocess, datetime
import numpy as np
from scipy.linalg import expm

ERROR_X = 1e-4
MIN_RATE, MAX_RATE = 1e-4, 100.0
PIN = 10.0                       # AA_SCALE used by rescaleRates
ITMAX, EPS, TOLX_DFP, STPMX = 200, 3.0e-8, 4*3.0e-8, 100.0
ALF, TOLX_LN = 1.0e-4, 1.0e-7

# ---------------------------------------------------------------- charts
def stationary(P):
    n = len(P)
    A = np.vstack((np.ones(n), (P - np.eye(n)).T))
    b = np.zeros(n+1); b[0] = 1.0
    nu, *_ = np.linalg.lstsq(A, b, rcond=None)      # column-pivoted QR in IQ-TREE; same solution
    return nu

def build_from_P(P, pi):
    nu = stationary(P)
    lam = nu / pi
    Q = lam[:, None] * P
    np.fill_diagonal(Q, -lam)
    return Q

def offdiag_idx(n):
    i, j = np.where(~np.eye(n, dtype=bool))
    return i, j

class RatioChart:
    """C2: raw positive rates, one pinned per row at PIN, row-normalized."""
    def __init__(self, n, pi, ref):
        self.n, self.pi, self.ref = n, pi, ref
        self.free = [(i, j) for i in range(n) for j in range(n) if j != i and j != ref[i]]
    def Q(self, x):
        n = self.n
        W = np.zeros((n, n))
        for k, (i, j) in enumerate(self.free):
            W[i, j] = x[k]
        W[np.arange(n), self.ref] = PIN
        P = W / W.sum(1, keepdims=True)
        return build_from_P(P, self.pi)
    def from_Q(self, Q):
        n = self.n
        P = Q / (-np.diag(Q))[:, None]; np.fill_diagonal(P, 0)
        W = P * (PIN / P[np.arange(n), self.ref])[:, None]
        return np.array([W[i, j] for (i, j) in self.free])
    def bounds(self):
        m = len(self.free)
        return np.full(m, MIN_RATE), np.full(m, MAX_RATE)

class LogChart:
    """Logits z, reference logit 0, row softmax."""
    def __init__(self, n, pi, ref):
        self.n, self.pi, self.ref = n, pi, ref
        self.free = [(i, j) for i in range(n) for j in range(n) if j != i and j != ref[i]]
    def Q(self, z):
        n = self.n
        Z = np.full((n, n), -np.inf)
        for k, (i, j) in enumerate(self.free):
            Z[i, j] = z[k]
        Z[np.arange(n), self.ref] = 0.0
        Z = Z - Z.max(1, keepdims=True)
        E = np.exp(Z); np.fill_diagonal(E, 0)
        P = E / E.sum(1, keepdims=True)
        return build_from_P(P, self.pi)
    def from_Q(self, Q):
        n = self.n
        P = Q / (-np.diag(Q))[:, None]
        return np.array([np.log(P[i, j] / P[i, self.ref[i]]) for (i, j) in self.free])
    def bounds(self):
        # identical physical box to the ratio chart: ratio/PIN in [1e-5, 10]
        m = len(self.free)
        return np.full(m, np.log(MIN_RATE/PIN)), np.full(m, np.log(MAX_RATE/PIN))

# ---------------------------------------------------------------- likelihood
class RootedTreeLik:
    """Rooted tree as parent pointers; root distribution = pi; pattern counts."""
    def __init__(self, n, ntips, rng):
        # random rooted binary tree by coalescing
        self.n = n
        nodes = list(range(ntips)); nxt = ntips
        parent = {}; blen = {}
        while len(nodes) > 1:
            a, b = rng.choice(len(nodes), 2, replace=False)
            a, b = nodes[a], nodes[b]
            parent[a] = nxt; parent[b] = nxt
            blen[a] = rng.uniform(0.05, 0.8); blen[b] = rng.uniform(0.05, 0.8)
            nodes = [x for x in nodes if x not in (a, b)] + [nxt]; nxt += 1
        self.root = nodes[0]; self.parent = parent; self.blen = blen
        self.ntips = ntips; self.nnodes = nxt
        self.children = {}
        for c, p in parent.items(): self.children.setdefault(p, []).append(c)
        self.postorder = self._postorder(self.root)
    def _postorder(self, v):
        out = []
        for c in self.children.get(v, []): out += self._postorder(c)
        return out + [v]
    def simulate(self, Q, pi, nsites, rng):
        n = self.n
        states = {self.root: rng.choice(n, size=nsites, p=pi)}
        for v in reversed(self.postorder):
            for c in self.children.get(v, []):
                P = expm(self.blen[c] * Q)
                s = states[v]
                states[c] = np.array([rng.choice(n, p=P[a]) for a in s])
        tips = np.stack([states[t] for t in range(self.ntips)], 1)   # nsites x ntips
        pats, counts = np.unique(tips, axis=0, return_counts=True)
        self.pats, self.counts = pats, counts
        self.ncalls = 0
    def loglik(self, Q, pi):
        self.ncalls += 1
        n = self.n; npat = len(self.pats)
        part = {}
        # eigen-based transition matrices, as IQ-TREE's nonreversible path does
        ev, V = np.linalg.eig(Q); Vinv = np.linalg.inv(V)
        def trans(t):
            return np.real(V @ (np.exp(ev*t)[:, None] * Vinv))
        for t in range(self.ntips):
            L = np.zeros((npat, n)); L[np.arange(npat), self.pats[:, t]] = 1.0
            part[t] = L
        for v in self.postorder:
            if v < self.ntips: continue
            L = np.ones((npat, n))
            for c in self.children[v]:
                P = trans(self.blen[c])          # P_ab = Pr(b at child | a at parent)
                L *= part[c] @ P.T
            part[v] = L
        site = part[self.root] @ pi
        return float(np.sum(self.counts * np.log(site)))

# ---------------------------------------------------------------- IQ-TREE optimizer port
class IQOpt:
    def __init__(self, f, lower, upper, fd_rule='iqtree', track=None):
        self.f, self.lower, self.upper, self.fd_rule = f, lower, upper, fd_rule
        self.nfev = 0
        self.track = track
    def targetFunk(self, x):
        self.nfev += 1
        return self.f(x)
    def derivativeFunk(self, x):
        n = len(x); h = np.empty(n); df = np.empty(n)
        fx = self.targetFunk(x)
        for d in range(n):
            t = x[d]
            if self.fd_rule == 'iqtree':
                hd = ERROR_X * abs(t)
                if hd == 0.0: hd = ERROR_X
            else:                       # scale-aware: h = eta*max(1,|x|)
                hd = ERROR_X * max(1.0, abs(t))
            x[d] = t + hd; hd = x[d] - t
            df[d] = self.targetFunk(x); x[d] = t; h[d] = hd
        return fx, (df - fx) / h
    def fixBound(self, x):
        np.clip(x, self.lower, self.upper, out=x)
    def lnsrch(self, xold, fold, g, p, stpmax):
        n = len(xold)
        s = np.sqrt(np.sum(p*p))
        if s > stpmax: p *= stpmax / s
        slope = float(g @ p)
        test = np.max(np.abs(p) / np.maximum(np.abs(xold), 1.0))
        alamin = TOLX_LN / test
        alam = 1.0; first = True; alam2 = f2 = fold2 = 0.0
        while True:
            x = xold + alam * p
            self.fixBound(x)
            f = self.targetFunk(x)
            if alam < alamin:
                return xold.copy(), fold, 1
            elif f <= fold + ALF * alam * slope:
                return x, f, 0
            else:
                if first:
                    tmplam = -slope / (2.0 * (f - fold - slope))
                else:
                    rhs1 = f - fold - alam * slope
                    rhs2 = f2 - fold2 - alam2 * slope
                    a = (rhs1/(alam*alam) - rhs2/(alam2*alam2)) / (alam - alam2)
                    b = (-alam2*rhs1/(alam*alam) + alam*rhs2/(alam2*alam2)) / (alam - alam2)
                    if a == 0.0: tmplam = -slope / (2.0*b)
                    else:
                        disc = b*b - 3.0*a*slope
                        if disc < 0.0: tmplam = 0.5*alam
                        elif b <= 0.0: tmplam = (-b + np.sqrt(disc)) / (3.0*a)
                        else: tmplam = -slope / (b + np.sqrt(disc))
                    if tmplam > 0.5*alam: tmplam = 0.5*alam
            alam2 = alam; f2 = f; fold2 = fold
            alam = max(tmplam, 0.1*alam); first = False
    def dfpmin(self, p, gtol):
        n = len(p)
        fp, g = self.derivativeFunk(p)
        hessin = np.eye(n)
        xi = -g
        stpmax = STPMX * max(np.sqrt(np.sum(p*p)), float(n))
        its_done = 0; reason = 'ITMAX'
        for its in range(1, ITMAX+1):
            its_done = its
            pnew, fret, check = self.lnsrch(p, fp, g, xi, stpmax)
            fp = fret
            xi = pnew - p; p = pnew
            if self.track is not None: self.track.append(-fp)
            test = np.max(np.abs(xi) / np.maximum(np.abs(p), 1.0))
            if test < TOLX_DFP: reason = 'TOLX'; break
            dg = g.copy()
            _, g = self.derivativeFunk(p)
            den = max(abs(fret), 1.0)
            test = np.max(np.abs(g) * np.maximum(np.abs(p), 1.0) / den)
            if test < gtol: reason = 'gtol'; break
            dg = g - dg
            hdg = hessin @ dg
            fac = float(dg @ xi); fae = float(dg @ hdg)
            sumdg = float(dg @ dg); sumxi = float(xi @ xi)
            if fac*fac > EPS*sumdg*sumxi:
                fac = 1.0/fac; fad = 1.0/fae
                u = fac*xi - fad*hdg
                hessin += fac*np.outer(xi, xi) - fad*np.outer(hdg, hdg) + fae*np.outer(u, u)
            xi = -(hessin @ g)
        return p, fp, its_done, reason

# ---------------------------------------------------------------- experiments
def make_problem(n, ntips, nsites, rng):
    pi = rng.dirichlet(np.full(n, 3.0))
    # true nonreversible Q stationary at pi: balanced flux from symmetric part plus directed cycles
    F = np.zeros((n, n))
    for a in range(n):
        for b in range(a+1, n):
            F[a, b] = F[b, a] = rng.lognormal(0, 1.0)
    for _ in range(3*n):
        cyc = rng.choice(n, size=3, replace=False); w = rng.lognormal(0, 1.0)
        for a, b in zip(cyc, np.roll(cyc, -1)): F[a, b] += w
    F /= F.sum()
    Qtrue = F / pi[:, None]; np.fill_diagonal(Qtrue, 0); np.fill_diagonal(Qtrue, -Qtrue.sum(1))
    tree = RootedTreeLik(n, ntips, rng)
    tree.simulate(Qtrue, pi, nsites, rng)
    # seed: reversible member of the fiber (R = symmetric part of F / (pi_i pi_j)), seed A analogue
    R = (F + F.T) / 2 / np.outer(pi, pi)
    Qseed = R * pi[None, :]; np.fill_diagonal(Qseed, 0); np.fill_diagonal(Qseed, -Qseed.sum(1))
    Qseed /= -(pi @ np.diag(Qseed))
    ref = np.argmax(np.where(np.eye(n, dtype=bool), -np.inf, Qseed), axis=1)
    return pi, Qtrue, Qseed, ref, tree

def central_ref(fun, x, d, steps):
    """central differences at several relative/absolute steps; return the plateau value"""
    vals = []
    for h in steps:
        xp = x.copy(); xm = x.copy(); xp[d] += h; xm[d] -= h
        vals.append((fun(xp) - fun(xm)) / (2*h))
    return vals

def experiment_fd(n, ntips, nsites, rng):
    pi, Qtrue, Qseed, ref, tree = make_problem(n, ntips, nsites, rng)
    ratio, logc = RatioChart(n, pi, ref), LogChart(n, pi, ref)
    # evaluate at a generic interior point: the true Q perturbed
    x_r = ratio.from_Q(Qtrue); z_l = logc.from_Q(Qtrue)
    fr = lambda x: -tree.loglik(ratio.Q(x), pi)
    fl = lambda z: -tree.loglik(logc.Q(z), pi)
    out = {'n': n, 'npatterns': int(len(tree.pats)), 'logL_at_true': -fr(x_r),
           'ratio_coord_range': [float(x_r.min()), float(x_r.max())],
           'log_coord_range': [float(z_l.min()), float(z_l.max())]}
    # IQ-TREE FD rule on both charts, every coordinate, against central plateau reference
    o_r = IQOpt(fr, *ratio.bounds()); o_l = IQOpt(fl, *logc.bounds())
    _, g_r = o_r.derivativeFunk(x_r.copy()); _, g_l = o_l.derivativeFunk(z_l.copy())
    # reference gradient via chain: dL/dz_k = dL/dx_k * dx_k/dz_k where x_k = PIN*exp(z_k) -> both
    # charts have the same physical gradient content; compute central refs on each chart directly
    rel_err_r, rel_err_l, abs_err_l, zsmall = [], [], [], []
    for k in range(len(x_r)):
        ref_r = central_ref(fr, x_r, k, [1e-5*abs(x_r[k]), 1e-4*abs(x_r[k])])
        ref_l = central_ref(fl, z_l, k, [1e-5, 1e-4])
        rr = ref_r[1]; rl = ref_l[1]
        rel_err_r.append(abs(g_r[k]-rr)/max(abs(rr), 1e-12))
        rel_err_l.append(abs(g_l[k]-rl)/max(abs(rl), 1e-12))
        abs_err_l.append(abs(g_l[k]-rl))
        zsmall.append(abs(z_l[k]))
    rel_err_r, rel_err_l, zsmall = map(np.array, (rel_err_r, rel_err_l, zsmall))
    out['ratio_chart_rel_err_median'] = float(np.median(rel_err_r))
    out['ratio_chart_rel_err_max'] = float(rel_err_r.max())
    out['log_chart_rel_err_median'] = float(np.median(rel_err_l))
    out['log_chart_rel_err_max'] = float(rel_err_l.max())
    out['log_chart_rel_err_for_|z|<0.01'] = [float(v) for v in rel_err_l[zsmall < 1e-2]]
    out['log_chart_|z|<0.01'] = [float(v) for v in z_l[zsmall < 1e-2]]
    # controlled sweep: push one log coordinate to exactly 1e-2, 1e-4, 1e-6, 1e-8 and 0
    k = int(np.argmin(zsmall)); sweep = {}
    for zval in [1e-2, 1e-4, 1e-6, 1e-8, 0.0]:
        z = z_l.copy(); z[k] = zval
        _, g = IQOpt(fl, *logc.bounds()).derivativeFunk(z.copy())
        rl = central_ref(fl, z, k, [1e-4])[0]
        sweep[str(zval)] = {'iqtree_fd': float(g[k]), 'central_ref': float(rl),
                            'rel_err': float(abs(g[k]-rl)/max(abs(rl), 1e-12))}
        # the same physical point in ratio coordinates
        x = ratio.from_Q(logc.Q(z))
        _, gr = IQOpt(fr, *ratio.bounds()).derivativeFunk(x.copy())
        rr = central_ref(fr, x, k, [1e-4*abs(x[k])])[0]
        sweep[str(zval)]['ratio_chart_same_point_rel_err'] = float(abs(gr[k]-rr)/max(abs(rr), 1e-12))
        sweep[str(zval)]['ratio_coordinate_value'] = float(x[k])
    out['log_coordinate_near_zero_sweep'] = sweep
    return out

def experiment_bfgs(n, ntips, nsites, rng, nrep):
    reps = []
    for r in range(nrep):
        pi, Qtrue, Qseed, ref, tree = make_problem(n, ntips, nsites, rng)
        ratio, logc = RatioChart(n, pi, ref), LogChart(n, pi, ref)
        x0 = ratio.from_Q(Qseed); z0 = logc.from_Q(Qseed)
        assert np.allclose(ratio.Q(x0), logc.Q(z0), atol=1e-12)
        gtol = 1e-4
        rec = {'rep': r, 'logL_seed': tree.loglik(Qseed, pi), 'logL_true': tree.loglik(Qtrue, pi)}
        for name, chart, xs, rule in [('ratio_iqtree_rule', ratio, x0, 'iqtree'),
                                      ('log_iqtree_rule', logc, z0, 'iqtree'),
                                      ('log_scale_aware_rule', logc, z0, 'scaled')]:
            f = lambda v, chart=chart: -tree.loglik(chart.Q(v), pi)
            lo, hi = chart.bounds()
            opt = IQOpt(f, lo, hi, fd_rule=rule)
            t0 = time.time()
            xf, ff, its, reason = opt.dfpmin(xs.copy(), gtol)
            at_bound = int(np.sum((np.abs(xf-lo) < 1e-4*np.maximum(1, np.abs(lo))) | (np.abs(xf-hi) < 1e-4*np.maximum(1, np.abs(hi)))))
            Qf = chart.Q(xf)
            rec[name] = {'logL_final': -ff, 'iterations': its, 'stop': reason, 'nfev': opt.nfev,
                         'coords_at_bound': at_bound,
                         'stationarity_resid': float(np.max(np.abs(pi @ Qf))),
                         'mean_rate_minus_1': float(-(pi @ np.diag(Qf)) - 1),
                         'max_abs_Q_diff_from_true': float(np.max(np.abs(Qf - Qtrue))),
                         'seconds': round(time.time()-t0, 1)}
        reps.append(rec)
        print(json.dumps(rec), file=sys.stderr)
    return reps

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--n', type=int, default=8)
    ap.add_argument('--ntips', type=int, default=8)
    ap.add_argument('--nsites', type=int, default=2000)
    ap.add_argument('--seed', type=int, default=20260922)
    ap.add_argument('--nrep', type=int, default=3)
    ap.add_argument('--skip-bfgs', action='store_true')
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    res = {'script': __file__, 'date': datetime.datetime.now().isoformat(), 'seed': a.seed,
           'numpy': np.__version__, 'args': vars(a)}
    res['E1_fd_accuracy_n20'] = experiment_fd(20, a.ntips, a.nsites, rng)
    if not a.skip_bfgs:
        res['E2_bfgs'] = experiment_bfgs(a.n, a.ntips, a.nsites, rng, a.nrep)
    print(json.dumps(res, indent=1))
