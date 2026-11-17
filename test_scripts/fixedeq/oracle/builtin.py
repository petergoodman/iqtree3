"""IQ-TREE's built-in protein matrices, read from model/modelprotein.cpp as check_builtin_matrices.py does."""
import re
from pathlib import Path

import numpy as np

SOURCE = Path(__file__).resolve().parents[3] / "model" / "modelprotein.cpp"
MIXTURE_SOURCE = SOURCE.parent / "modelmixture.cpp"
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


def frequency_profiles(name="C20", path=MIXTURE_SOURCE):
    """The profiles (one row each) and weights of a fixed-weight profile mixture
    `model NAME = ...FMIX{P1:1:w1,...};`, in FMIX order, each profile from its `frequency P = ...;`
    line; profiles and weights normalized to sum 1."""
    text = Path(path).read_text(encoding="utf-8")
    m = re.search(r"^model %s = .*?FMIX\{([^}]*)\}" % re.escape(name), text, re.M)
    if not m:
        raise KeyError(f"model {name} with FMIX not found in {path}")
    profiles, weights = [], []
    for entry in m.group(1).split(","):
        p, rate, w = entry.split(":")
        if float(rate) != 1.0:
            raise ValueError(f"{name}: component {p} has rate {rate}, not 1")
        f = re.search(r"^frequency %s =([^;]*);" % re.escape(p), text, re.M)
        if not f:
            raise KeyError(f"frequency {p} not found in {path}")
        v = np.array([float(t) for t in f.group(1).split()])
        if v.size != 20:
            raise ValueError(f"frequency {p}: {v.size} numbers, expected 20")
        profiles.append(v / v.sum())
        weights.append(float(w))
    return np.array(profiles), np.array(weights) / sum(weights)
