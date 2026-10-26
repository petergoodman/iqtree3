"""Readers for the alignment, tree and partition files the oracle uses, and for NQC's export."""
import math
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np


def read_phylip(path):
    """Interleaved or sequential PHYLIP: names in the first block only, whitespace ignored."""
    lines = [line for line in Path(path).read_text(encoding="utf-8").splitlines() if line.strip()]
    ntax, nsite = (int(t) for t in lines[0].split()[:2])
    names, seqs = [], []
    for k, line in enumerate(lines[1:]):
        if k < ntax:
            name, *rest = line.split()
            names.append(name)
            seqs.append("".join(rest))
        else:
            seqs[k % ntax] += "".join(line.split())
    if len(names) != ntax or any(len(s) != nsite for s in seqs):
        raise ValueError(f"{path}: expected {ntax} sequences of {nsite} characters")
    return dict(zip(names, seqs))


def read_fasta(path):
    seqs, name = {}, None
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = ""
        elif line.strip():
            seqs[name] += "".join(line.split())
    return seqs


def read_charsets(path):
    """`charset NAME = a-b c ...;` lines, as 0-based site index arrays in file order."""
    text = Path(path).read_text(encoding="utf-8")
    sets = {}
    for name, spec in re.findall(r"charset\s+(\S+)\s*=\s*([^;]+);", text, re.I):
        sites = []
        for token in spec.split():
            if "\\" in token:
                raise ValueError(f"charset {name}: step ranges are not supported")
            lo, _, hi = token.partition("-")
            sites.extend(range(int(lo) - 1, int(hi or lo)))
        sets[name] = np.array(sites)
    return sets


@dataclass
class Node:
    name: str = None
    length: float = None
    children: list = field(default_factory=list)

    def preorder(self):
        yield self
        for c in self.children:
            yield from c.preorder()

    def leaves(self):
        return [n for n in self.preorder() if not n.children]


def parse_newick(text):
    s = "".join(text.split())
    if not s.endswith(";"):
        raise ValueError("a Newick tree ends with ';'")
    pos = 0
    label = re.compile(r"[^:,();]*")
    number = re.compile(r"[-+0-9.eE]+")

    def node():
        nonlocal pos
        nd = Node()
        if s[pos] == "(":
            pos += 1
            nd.children.append(node())
            while s[pos] == ",":
                pos += 1
                nd.children.append(node())
            if s[pos] != ")":
                raise ValueError(f"Newick: expected ')' at {pos}")
            pos += 1
        m = label.match(s, pos)
        nd.name, pos = m.group(0) or None, m.end()
        if s[pos] == ":":
            m = number.match(s, pos + 1)
            nd.length, pos = float(m.group(0)), m.end()
        return nd

    root = node()
    if s[pos:] != ";":
        raise ValueError(f"Newick: unexpected text at {pos}")
    return root


def read_newick(path):
    return parse_newick(Path(path).read_text(encoding="utf-8"))


def read_qmat(text):
    """An exported <prefix>.NQC.qmat (decision 010): n rows of the full Q, then the frequency row,
    whitespace-separated; returns (Q, pi) with all n(n + 1) numbers."""
    values = np.array([float(v) for v in text.split()])
    n = int(round((math.sqrt(1 + 4 * values.size) - 1) / 2))
    if n < 2 or n * (n + 1) != values.size:
        raise ValueError(f".qmat: {values.size} numbers do not make n rows of n and a row of n")
    return values[:n * n].reshape(n, n), values[n * n:]


def read_info(text):
    """An exported <prefix>.NQC.info: one `key: value` per line, lines starting with # ignored;
    returns the values as strings. A repeated key or a line without ': ' is an error."""
    out = {}
    for line in text.splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        key, sep, value = line.partition(": ")
        if not sep or not key or key in out:
            raise ValueError(f".info: bad or repeated line {line!r}")
        out[key] = value
    return out
