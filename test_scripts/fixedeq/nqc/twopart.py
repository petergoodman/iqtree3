"""S3b's two-partition data (CODE_PLAN.md section 4, S3b): one seeded rooted 12-taxon tree;
partition A simulated by the frozen binary's AliSim from NQ.pfam's jump chain transferred to
decision 023's target, partition B from the same chain at that target cubed and renormalized; 600
sites each, with a two-charset NEXUS file. Used by nqc/nqcchecks.py.

No taxon is missing and no two sequences are identical in either partition, since -p removes
identical sequences; -p's partition trees then equal the shared tree, and the oracle needs no
pruning.
"""
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE.parent / "regression"))
from oracle import readers  # noqa: E402
from oracle.likelihood import UNKNOWN  # noqa: E402
from regress import sha256, write  # noqa: E402
import toy  # noqa: E402

SEED_TREE, SEED_A, SEED_B = 20261023, 20261024, 20261025
TAXA, SITES = 12, 600
CHARSETS = ("partA", "partB")


def cubed(pi):
    """Decision 023's target cubed and renormalized, partition B's composition."""
    p = np.asarray(pi, dtype=float) ** 3
    return p / p.sum()


def nexus(order=CHARSETS):
    """The charsets, in the given order, over the concatenated alignment."""
    start = {"partA": 1, "partB": SITES + 1}
    lines = [f"charset {name} = {start[name]}-{start[name] + SITES - 1};" for name in order]
    return "#nexus\nbegin sets;\n" + "\n".join(lines) + "\nend;\n"


def problems(aln):
    """Missing or all-unknown sequences, and pairs of identical sequences."""
    out = [f"{n} has no known character" for n, s in aln.items() if all(c in UNKNOWN for c in s.upper())]
    names = sorted(aln)
    out += [f"{a} and {b} are identical" for i, a in enumerate(names) for b in names[i + 1:] if aln[a] == aln[b]]
    return out


def make(out, alisim, pi):
    """Write partA.*, partB.*, and the concatenated twopart.phy, twopart.nex and twopart.nwk into the
    existing directory out; stop if a partition lacks a taxon or holds identical sequences."""
    out = Path(out)
    tree = toy.tree_text(SEED_TREE, TAXA)
    made = {}
    for name, p, seed in (("partA", np.asarray(pi, dtype=float), SEED_A), ("partB", cubed(pi), SEED_B)):
        made[name] = toy.simulate(out, alisim, name, tree, toy.truth(p), p, SITES, seed)
    parts = [readers.read_phylip(out / f"{name}.phy") for name in CHARSETS]
    taxa = {f"T{i + 1}" for i in range(TAXA)}
    found = [f"{name}: taxa {sorted(taxa ^ set(a))} missing or extra" for name, a in zip(CHARSETS, parts)
             if set(a) != taxa]
    # distinct within each partition implies distinct in the whole
    for name, a in zip(CHARSETS, parts):
        found += [f"{name}: {p}" for p in problems(a)]
    if found:
        sys.exit("two-partition data unusable: " + "; ".join(found))
    whole = {n: parts[0][n] + parts[1][n] for n in sorted(taxa, key=lambda n: int(n[1:]))}
    write(out / "twopart.phy", f"{TAXA} {2 * SITES}\n" + "".join(f"{n} {s}\n" for n, s in whole.items()))
    write(out / "twopart.nex", nexus())
    write(out / "twopart.nwk", tree)
    made["sha256"] = {n: sha256(out / n) for n in ("twopart.phy", "twopart.nex", "twopart.nwk")}
    made["seeds"] = {"tree": SEED_TREE, "partA": SEED_A, "partB": SEED_B}
    return made
