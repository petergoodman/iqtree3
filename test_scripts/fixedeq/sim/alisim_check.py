"""Confirmatory: SIM-1's AliSim exact check (CODE_PLAN.md 3.8, "AliSim exact check"; decision 021's SIM-1
AliSim rows).

  python test_scripts/fixedeq/sim/alisim_check.py run --frozen ~/iqtree3-baseline/iqtree3 \
      --gate ~/iqtree3-binaries/25d78d6e/iqtree3 --git git.exe [--out DIR] [--length 10000000] [--allow-dirty]

Two taxa, four configurations: M1 (NQ.bird at pi*, written by differential.write_matrix) on T_joint
(a:0.2,b:0.6), with and without +G4{0.5}; M2 (the .NQC.qmat of an NQC Level 1 fit on aa_example at pi*
by the gate binary) on T_joint; M1 on T_root (a:0,b:0.6), whose tip a copies the root. Each is simulated
by the frozen binary and then by the gate binary with the same command and seed, each in a new
directory. Asserted: Pearson X^2 of the 400 tip-pair counts against the oracle's exact law (the -m FILE
convention of the same file, site_log_likelihoods on the 400 columns), p >= 1e-3, per T_joint
configuration; T_root's tip a against pi* over 20 cells, p >= 1e-3; only the 20 residues; the gate's
output byte-identical to the frozen binary's; M2's fit exits 0 with its export. Reported before any
simulation: the noncentrality and power of named alternatives at the run's length. A failing
configuration is rerun with +F{pi*} and the same seed as a diagnosis only; its row stays failed.
"""
import argparse
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
from scipy import special, stats
from scipy.linalg import expm

HERE = Path(__file__).resolve().parent
for _p in ("", "regression", "probes", "nqc", "differential"):
    sys.path.insert(0, str(HERE.parent / _p))
sys.path.insert(0, str(HERE))
from oracle import likelihood as lk, readers, simulate  # noqa: E402
from regress import REPO, read, sha256, write  # noqa: E402
from probes import AA, Runner, baseline_runs, tree_text  # noqa: E402
from nqcchecks import check  # noqa: E402
from differential import write_matrix  # noqa: E402
import common  # noqa: E402

SCRIPT = Path(__file__).resolve().relative_to(REPO).as_posix()
LENGTH = 10_000_000
SEED = 20261117
ALPHA = 0.5
P_MIN = 1e-3            # 021: p >= 1e-3 per configuration
POOL = 5.0              # cells expected below this are pooled
T_JOINT, T_ROOT = "(a:0.2,b:0.6);", "(a:0,b:0.6);"
LENGTHS = (0.2, 0.6)
CONFIGS = {"m1": ("sim_m1.qmat", "", T_JOINT), "m1g4": ("sim_m1.qmat", f"+G4{{{ALPHA}}}", T_JOINT),
           "m2": ("sim_m2.qmat", "", T_JOINT), "m1root": ("sim_m1.qmat", "", T_ROOT)}
# noncentralities at N = 1e6 measured while planning (docs/agent/archive/SIM1_plan_2026-11-16.md), shown beside
PLANNING_1E6 = {"start at tip a": 651, "start at tip b": 5440, "tips swapped": 2380, "reversed process": 2520,
                "reversible part": 631, "median gamma": 502, "continuous gamma": 1120, "scale +3%": 390,
                "scale +1%": 44}
CODE = np.full(256, -1, dtype=np.int64)
for _k, _c in enumerate(lk.AMINO):
    CODE[ord(_c)] = _k


# ---------------------------------------------------------------- laws and statistics

def direct_law(Q, pi, ta, tb, rates=(1.0,), weights=(1.0,)):
    """J_xy = sum_k w_k [P(r_k t_a)^T diag(pi) P(r_k t_b)]_xy: the root above both tips, drawn from pi."""
    J = np.zeros((len(pi), len(pi)))
    for r, w in zip(rates, weights):
        J += w * expm(Q * r * ta).T @ (np.asarray(pi)[:, None] * expm(Q * r * tb))
    return J


def exact_law(tree, Q, pi, rates=(1.0,), weights=(1.0,)):
    """J[x_a, x_b] through the oracle's likelihood on the 400 columns of simulate.columns(("a", "b"))."""
    site = lk.site_log_likelihoods(readers.parse_newick(tree), simulate.columns(("a", "b")), Q, pi, rates, weights)
    return np.exp(site).reshape(20, 20)


def reversed_process(Q, pi):
    """The time reversal q~_ij = pi_j q_ji / pi_i, stationary at pi."""
    return np.asarray(Q).T * np.asarray(pi)[None, :] / np.asarray(pi)[:, None]


def median_categories(shape, ncat=4):
    """Equal-weight categories at the gamma's interval medians, scaled to mean 1 (an alternative only)."""
    r = stats.gamma.ppf((2 * np.arange(ncat) + 1) / (2 * ncat), a=shape, scale=1 / shape)
    return r / r.mean(), np.full(ncat, 1 / ncat)


def continuous_gamma(shape, nodes=64):
    """The continuous gamma of mean 1 by generalized Gauss-Laguerre quadrature: rates x_i / shape."""
    x, w = special.roots_genlaguerre(nodes, shape - 1)
    return x / shape, w / w.sum()


def pearson(observed, expected, pool=POOL):
    """Pearson X^2. Cells expected below pool are merged into one, and that cell, if still below pool,
    into the smallest other. No parameter is estimated, so df = cells - 1."""
    o, e = np.ravel(observed).astype(float), np.ravel(expected).astype(float)
    small = e < pool
    if small.any():
        o, e = np.r_[o[~small], o[small].sum()], np.r_[e[~small], e[small].sum()]
        if e[-1] < pool and len(e) > 1:
            k = int(np.argmin(e[:-1]))
            o[k], e[k] = o[k] + o[-1], e[k] + e[-1]
            o, e = o[:-1], e[:-1]
    x2 = float(((o - e) ** 2 / e).sum())
    df = len(e) - 1
    return {"x2": x2, "df": df, "p": float(stats.chi2.sf(x2, df)), "cells": len(e),
            "pooled": int(small.sum()), "smallest_expected": float(e.min())}


def codes(seq):
    """Residue indices 0..19 in IQ-TREE's order; -1 for any other character."""
    return CODE[np.frombuffer(seq.encode("ascii"), dtype=np.uint8)]


def residue_counts(seq):
    c = codes(seq)
    return np.bincount(c[c >= 0], minlength=20), int((c < 0).sum())


def pair_counts(a, b):
    ca, cb = codes(a), codes(b)
    ok = (ca >= 0) & (cb >= 0)
    return np.bincount(20 * ca[ok] + cb[ok], minlength=400).reshape(20, 20)


def power_rows(Q, pi, n):
    """Noncentrality n sum (p - p0)^2 / p0 and power against X^2_399's 1e-3 critical value, for the named
    alternatives, from the direct law."""
    ta, tb = LENGTHS
    rates, weights = lk.rate_categories(shape=ALPHA)
    null, null_g = direct_law(Q, pi, ta, tb), direct_law(Q, pi, ta, tb, rates, weights)
    alternatives = {
        "start at tip a": (null, direct_law(Q, pi, 0.0, ta + tb)),
        "start at tip b": (null, direct_law(Q, pi, ta + tb, 0.0)),
        "tips swapped": (null, null.T),
        "reversed process": (null, direct_law(reversed_process(Q, pi), pi, ta, tb)),
        "reversible part": (null, direct_law(simulate.dial(Q, pi, 0.0), pi, ta, tb)),
        **{f"scale {s:+.0%}": (null, direct_law(Q, pi, ta * (1 + s), tb * (1 + s))) for s in (0.01, -0.01, 0.03, -0.03)},
        "median gamma": (null_g, direct_law(Q, pi, ta, tb, *median_categories(ALPHA))),
        "continuous gamma": (null_g, direct_law(Q, pi, ta, tb, *continuous_gamma(ALPHA))),
        **{f"shape {s:g}": (null_g, direct_law(Q, pi, ta, tb, *lk.rate_categories(shape=s))) for s in (0.45, 0.55)}}
    critical = float(stats.chi2.isf(P_MIN, 399))
    rows = []
    for name, (p0, p1) in alternatives.items():
        lam = float(n * ((p1 - p0) ** 2 / p0).sum())
        rows.append(check(f"power: {name}", None, length=n, noncentrality=lam, critical=critical,
                          power=float(stats.ncx2.sf(critical, 399, lam)), noncentrality_1e6=lam * 1e6 / n,
                          planning_1e6=PLANNING_1E6.get(name)))
    return rows


# ---------------------------------------------------------------- running

def alisim(binary, rundir, name, matrix, rate, tree, seed, length):
    """One AliSim run in the new directory rundir; returns (exit code, output path, printed seed, argv)."""
    rundir.mkdir()
    shutil.copyfile(matrix, rundir / matrix.name)
    write(rundir / "tree.nwk", tree + "\n")
    argv = [str(binary), "--alisim", name, "-pre", name, "--seqtype", "AA", "-m", f"{matrix.name}{rate}",
            "-t", "tree.nwk", "--length", str(length), "--seed", str(seed), "-af", "phy", "-nt", "1",
            "--simulation-thresh", "0"]
    t0 = time.perf_counter()
    with open(rundir / "stdout.txt", "w", encoding="utf-8") as fh:
        code = subprocess.run(argv, cwd=rundir, stdout=fh, stderr=subprocess.STDOUT).returncode
    printed = re.search(r"^Seed:\s+(\d+)", read(rundir / "stdout.txt"), re.M)
    print(f"  {rundir.name}: exit {code}, {time.perf_counter() - t0:.1f} s", flush=True)
    return code, rundir / f"{name}.phy", int(printed.group(1)) if printed else None, argv


def analyse(name, path, law, pi_star):
    """Rows for one frozen-binary output: its residues, and its X^2 against the exact law."""
    aln = readers.read_phylip(path)
    if sorted(aln) != ["a", "b"]:
        return [check(f"{name}: the output holds taxa a and b", False, names=sorted(aln))], None
    rows, other = [], sum(residue_counts(s)[1] for s in aln.values())
    rows.append(check(f"{name}: only the 20 residues appear", other == 0, other_characters=other))
    if CONFIGS[name][2] == T_ROOT:
        counts = residue_counts(aln["a"])[0]
        test = pearson(counts, counts.sum() * pi_star / pi_star.sum())
        rows.append(check(f"{name}: tip a (the root) against pi*, X^2 over 20 cells, p >= {P_MIN:g}",
                          test["p"] >= P_MIN, **test))
    else:
        counts = pair_counts(aln["a"], aln["b"])
        test = pearson(counts, counts.sum() * law)
        rows.append(check(f"{name}: joint law, X^2 over 400 cells, p >= {P_MIN:g}", test["p"] >= P_MIN, **test))
    return rows, test


def cmd_run(a):
    frozen, gate = (Path(p).expanduser().resolve() for p in (a.frozen, a.gate))
    seeds = {name: SEED + k for k, name in enumerate(CONFIGS)}
    prov = common.provenance(a, SCRIPT, "alisim-check", "confirmatory",
                             frozen={"path": str(frozen), "sha256": sha256(frozen)},
                             gate={"path": str(gate), "sha256": sha256(gate)},
                             seeds=seeds, length=a.length, alpha=ALPHA, p_min=P_MIN, pool=POOL,
                             familywise_false_failure=1 - (1 - P_MIN) ** 4)
    out = common.new_out(a.out, prov)
    text = common.target_text()
    pi_star = common.target_values()
    prov["target"] = {"text": text, "sha256": sha256(common.TARGET)}

    m1 = out / "sim_m1.qmat"
    Q1 = simulate.truths()[1]["NQ.bird"]
    write_matrix(m1, Q1, pi_star)
    Qm1, pim1 = lk.file_model(readers.read_qmat(read(m1))[0])
    rows = power_rows(Qm1, pim1, a.length)

    fit = Runner(gate, out, a.timeout).run("m2-fit", ["-s", AA, "-m", f"NQC+F{{{text}}}+G4{{{ALPHA}}}",
                                                      "-te", "input.tree", "-blfix"],
                                           files={"input.tree": tree_text(baseline_runs(), "run02")})
    export = fit.file(".NQC.qmat")
    fitted = fit.code == 0 and export.exists()
    rows.append(check("M2: NQC fit on aa_example at pi* exits 0 with its export", fitted, exit_code=fit.code))
    matrices = {"sim_m1.qmat": m1}
    if fitted:
        matrices["sim_m2.qmat"] = out / "sim_m2.qmat"
        shutil.copyfile(export, matrices["sim_m2.qmat"])
        pi_row = readers.read_qmat(read(export))[1]
        rows.append(check("M2: the export's pi row against pi*", None, bitwise=bool(np.array_equal(pi_row, pi_star)),
                          max_abs=float(np.max(np.abs(pi_row - pi_star)))))

    for name, (fname, rate, tree) in CONFIGS.items():
        if fname not in matrices:
            rows.append(check(f"{name}: simulated", False, reason="M2's fit failed"))
            continue
        Q, pi = lk.file_model(readers.read_qmat(read(matrices[fname]))[0])
        rates, weights = lk.rate_categories(shape=ALPHA) if rate else (np.ones(1), np.ones(1))
        law = exact_law(tree, Q, pi, rates, weights) if tree == T_JOINT else None
        if law is not None:
            ta, tb = LENGTHS
            rows.append(check(f"{name}: exact law, site_log_likelihoods against the direct product", None,
                              max_relative=float(np.max(np.abs(law - direct_law(Q, pi, ta, tb, rates, weights)) / law)),
                              total=float(law.sum())))
        runs = {}
        for which, binary in (("frozen", frozen), ("gate", gate)):
            runs[which] = alisim(binary, out / f"{which}-{name}", name, matrices[fname], rate, tree, seeds[name],
                                 a.length)
            code, path, printed, argv = runs[which]
            rows.append(check(f"{name}: {which} AliSim exits 0 with its output", code == 0 and path.exists(),
                              exit_code=code, seed=seeds[name], printed_seed=printed, argv=argv))
        code, path = runs["frozen"][:2]
        if not (code == 0 and path.exists()):
            continue
        found, test = analyse(name, path, law, pi_star)
        rows += found
        gpath = runs["gate"][1]
        if gpath.exists():
            gate_other = sum(residue_counts(s)[1] for s in readers.read_phylip(gpath).values())
            rows.append(check(f"{name}: gate output byte-identical to the frozen binary's",
                              sha256(gpath) == sha256(path), frozen=sha256(path), gate=sha256(gpath),
                              gate_other_characters=gate_other))
        if test is not None and test["p"] < P_MIN:
            code, fpath, printed, argv = alisim(frozen, out / f"frozen-{name}-fallback", name, matrices[fname],
                                                f"{rate}+F{{{text}}}", tree, seeds[name], a.length)
            diag = analyse(name, fpath, law, pi_star)[1] if code == 0 and fpath.exists() else None
            rows.append(check(f"{name}: fallback with +F{{pi*}}, same seed (diagnosis only)", None, exit_code=code,
                              test=diag))
    common.finish(out, prov, "SIM-1 AliSim exact check", rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="step", required=True)
    r = sub.add_parser("run")
    r.add_argument("--frozen", required=True, help="the frozen unmodified binary")
    r.add_argument("--gate", required=True, help="the gate binary")
    r.add_argument("--out", help="new directory outside the repository (default under ~/iqtree3-runs/sim1/)")
    r.add_argument("--length", type=int, default=LENGTH)
    r.add_argument("--timeout", type=float, default=3600)
    r.add_argument("--git", default="git", help="git executable (git.exe under WSL)")
    r.add_argument("--allow-dirty", action="store_true")
    a = ap.parse_args()
    cmd_run(a)


if __name__ == "__main__":
    main()
