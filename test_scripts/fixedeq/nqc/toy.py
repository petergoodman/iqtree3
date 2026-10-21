"""Decision 030's toy problem: 8 taxa and 3000 sites, simulated without rate variation by the frozen
binary's AliSim from NQ.pfam's jump chain transferred to decision 023's target, on a rooted tree
drawn from a fixed seed. Used by explore/toy_calibration.py and nqc/nqcchecks.py.
"""
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE.parent / "differential"))
sys.path.insert(0, str(HERE.parent / "regression"))
from oracle import builtin, chart  # noqa: E402
from differential import write_matrix  # noqa: E402
from regress import sha256, write  # noqa: E402

SEED = 20261020
TAXA, SITES = 8, 3000
LENGTHS = (0.05, 0.3)


def tree_text(seed=SEED):
    """A rooted binary tree on T1..T8, built by joining two random subtrees at a time; every branch
    length is drawn from LENGTHS and written with 4 decimals, which both IQ-TREE and the oracle read."""
    rng = np.random.default_rng(seed)
    nodes = [f"T{i + 1}:{rng.uniform(*LENGTHS):.4f}" for i in range(TAXA)]
    while len(nodes) > 2:
        a, b = sorted(rng.choice(len(nodes), size=2, replace=False), reverse=True)
        joined = f"({nodes[b]},{nodes[a]})"
        del nodes[a], nodes[b]
        nodes.append(f"{joined}:{rng.uniform(*LENGTHS):.4f}")
    return f"({nodes[0]},{nodes[1]});\n"


def truth(pi):
    """NQ.pfam's jump chain at pi; its stored diagonal is not exactly minus the row sum, so it is
    rebuilt first, as make_fixtures.py does."""
    Q0 = builtin.full_matrix("NQ.PFAM")[0].copy()
    np.fill_diagonal(Q0, 0.0)
    np.fill_diagonal(Q0, -Q0.sum(axis=1))
    return chart.transfer_seed(Q0, pi)


def make(out, alisim, pi):
    """Write toy.nwk, toy.qmat and the simulated toy.phy into the existing directory out."""
    out = Path(out)
    write(out / "toy.nwk", tree_text())
    write_matrix(out / "toy.qmat", truth(pi), pi)
    argv = [str(alisim), "--alisim", "toy", "--seqtype", "AA", "-m", "toy.qmat", "-t", "toy.nwk",
            "--length", str(SITES), "--seed", str(SEED), "-af", "phy", "-nt", "1"]
    with open(out / "alisim.log", "w", encoding="utf-8") as fh:
        code = subprocess.run(argv, cwd=out, stdout=fh, stderr=subprocess.STDOUT).returncode
    if code != 0 or not (out / "toy.phy").exists():
        sys.exit(f"AliSim failed with exit {code}; see {out / 'alisim.log'}")
    return {"argv": argv, "sha256": {n: sha256(out / n) for n in ("toy.nwk", "toy.qmat", "toy.phy")}}
