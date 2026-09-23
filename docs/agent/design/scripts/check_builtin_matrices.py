#!/usr/bin/env python3
"""
check_builtin_matrices.py

Reads the LG and NQ.PFAM definitions out of IQ-TREE 3's model/modelprotein.cpp
(path given as argv[1], default ./iqtree3/model/modelprotein.cpp) and reports
the quantities that matter for the constrained (fixed-pi) chart design:
within-row dynamic range of the seed Q, stationarity of the shipped pi,
mean rate, C2 (row-shape) coordinates, gauge-fixed C1 tilt magnitudes.
Writes check_builtin_matrices_results.json.  No randomness.
"""
import sys, re, json
import numpy as np

src = sys.argv[1] if len(sys.argv) > 1 else "iqtree3/model/modelprotein.cpp"
text = open(src).read()

def read_block(name):
    m = re.search(r"^model %s=\n(.*?);" % re.escape(name), text, re.S | re.M)
    if not m: raise SystemExit("model %s not found" % name)
    rows = [r for r in m.group(1).strip().split("\n") if r.strip() and not r.strip().startswith("[")]
    return [list(map(float, r.split())) for r in rows]

def stationary(Q):
    n = Q.shape[0]; A = np.vstack([np.ones(n), Q.T]); b = np.zeros(n + 1); b[0] = 1
    return np.linalg.lstsq(A, b, rcond=None)[0]

def generator_from_off(off):
    Q = off.copy(); np.fill_diagonal(Q, 0); np.fill_diagonal(Q, -Q.sum(1)); return Q

R = {}
# ---- LG: lower-triangular exchangeabilities then frequencies
lg = read_block("LG"); n = 20
Rm = np.zeros((n, n))
for i in range(1, n):
    for j in range(i):
        Rm[i, j] = Rm[j, i] = lg[i - 1][j]
pi_lg = np.array(lg[n - 1]); pi_lg /= pi_lg.sum()
Q_lg = generator_from_off(Rm * pi_lg[None, :]); Q_lg /= -(pi_lg @ np.diag(Q_lg))
off = ~np.eye(n, dtype=bool)
row_range = [float(Q_lg[i, off[i]].max() / Q_lg[i, off[i]].min()) for i in range(n)]
R["LG"] = {"stationarity_residual_at_shipped_pi": float(np.max(np.abs(pi_lg @ Q_lg))),
           "exchangeability_min_max": [float(Rm[off].min()), float(Rm[off].max())],
           "Q_offdiag_min_max": [float(Q_lg[off].min()), float(Q_lg[off].max())],
           "within_row_dynamic_range_max": float(max(row_range)),
           "within_row_dynamic_range_median": float(np.median(row_range)),
           "global_dynamic_range": float(Q_lg[off].max() / Q_lg[off].min())}
# IQ-TREE's seed handling: rescale so max rate = 10 (AA_SCALE), then floor 1e-4 applies as a bound
scaled = Q_lg[off] * (10.0 / Q_lg[off].max())
R["LG"]["fraction_of_seed_rates_below_1e-4_after_AA_SCALE"] = float(np.mean(scaled < 1e-4))

# ---- NQ.PFAM: full Q rows then frequencies
nq = read_block("NQ.PFAM")
Q_nq = np.array(nq[:n]); pi_nq_file = np.array(nq[n]); pi_nq_file /= pi_nq_file.sum()
pi_nq = stationary(Q_nq)
R["NQ.PFAM"] = {"row_sums_max_abs": float(np.max(np.abs(Q_nq.sum(1)))),
                "mean_rate_at_solved_pi": float(-(pi_nq @ np.diag(Q_nq))),
                "shipped_pi_vs_solved_pi_max_abs": float(np.max(np.abs(pi_nq - pi_nq_file))),
                "Q_offdiag_min_max": [float(Q_nq[off].min()), float(Q_nq[off].max())],
                "within_row_dynamic_range_max": float(max(Q_nq[i, off[i]].max() / Q_nq[i, off[i]].min() for i in range(n))),
                "global_dynamic_range": float(Q_nq[off].max() / Q_nq[off].min())}
# C2 coordinates of NQ.PFAM: row shapes, reference = row argmax -> ratios in (0,1]
P = Q_nq / (-np.diag(Q_nq))[:, None]; np.fill_diagonal(P, 0)
ratios = P / P.max(1, keepdims=True)
R["NQ.PFAM"]["c2_min_ratio_to_row_max"] = float(ratios[off].min())
R["NQ.PFAM"]["c2_fraction_ratios_below_1e-4"] = float(np.mean(ratios[off] < 1e-4))
# retarget NQ.PFAM to LG's pi by C2 (row shapes kept): stationarity check
nu = stationary(P - np.eye(n)); lam = nu / pi_lg; Q_re = lam[:, None] * P; np.fill_diagonal(Q_re, -lam)
R["NQ.PFAM"]["c2_retarget_to_LG_pi_residual"] = float(np.max(np.abs(pi_lg @ Q_re)))
R["NQ.PFAM"]["c2_retarget_mean_rate"] = float(-(pi_lg @ np.diag(Q_re)))
R["NQ.PFAM"]["c2_retarget_max_rel_change_vs_original"] = float(np.max(np.abs(Q_re[off] - Q_nq[off]) / Q_nq[off]))
# C1 tilt coordinates of NQ.PFAM at its own pi: raw h and gauge-fixed h' (star at state 0)
M = Q_nq / pi_nq[None, :]; np.fill_diagonal(M, 1)
h = 0.5 * np.log(M / M.T); np.fill_diagonal(h, 0)
c = h[0].copy(); hg = h + c[:, None] - c[None, :]
R["NQ.PFAM"]["c1_raw_h_max_abs"] = float(np.max(np.abs(h)))
R["NQ.PFAM"]["c1_gauge_fixed_h_max_abs"] = float(np.max(np.abs(hg)))
R["NQ.PFAM"]["c1_r_min_max"] = [float(np.sqrt(M * M.T)[off].min()), float(np.sqrt(M * M.T)[off].max())]
# flux antisymmetry summary: fraction of pairs with |F_ij - F_ji|/(F_ij+F_ji) > 0.2
F = pi_nq[:, None] * Q_nq; np.fill_diagonal(F, 0)
asym = np.abs(F - F.T) / (F + F.T + 1e-300)
R["NQ.PFAM"]["fraction_pairs_flux_asymmetry_gt_0.2"] = float(np.mean(asym[np.triu_indices(n, 1)] > 0.2))

json.dump(R, open(__file__.replace(".py", "_results.json"), "w"), indent=1)
print(json.dumps(R, indent=1))
