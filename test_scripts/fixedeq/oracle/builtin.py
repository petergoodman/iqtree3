"""IQ-TREE's built-in protein matrices, read from model/modelprotein.cpp as check_builtin_matrices.py does."""
import re
from pathlib import Path

import numpy as np

SOURCE = Path(__file__).resolve().parents[3] / "model" / "modelprotein.cpp"
STATES = "ARNDCQEGHILKMFPSTWYV"


def _numbers(name, path):
    text = Path(path).read_text(encoding="utf-8")
    m = re.search(r"^model %s=\s*\n(.*?);" % re.escape(name), text, re.S | re.M)
    if not m:
        raise KeyError(f"model {name} not found in {path}")
    return np.array([float(t) for t in re.sub(r"\[.*?\]", " ", m.group(1), flags=re.S).split()])


def exchangeabilities(name, path=SOURCE):
    """R (symmetric, zero diagonal) and the frequencies of a model stored as a lower triangle."""
    v = _numbers(name, path)
    n = 20
    if v.size != n * (n - 1) // 2 + n:
        raise ValueError(f"{name}: {v.size} numbers, expected {n * (n - 1) // 2 + n}")
    R = np.zeros((n, n))
    R[np.tril_indices(n, -1)] = v[:-n]
    R = R + R.T
    return R, v[-n:] / v[-n:].sum()


def full_matrix(name="NQ.PFAM", path=SOURCE):
    """Q (20 rows) and the frequency row of a model stored as a full matrix."""
    v = _numbers(name, path)
    if v.size != 420:
        raise ValueError(f"{name}: {v.size} numbers, expected 420")
    return v[:400].reshape(20, 20), v[400:] / v[400:].sum()
