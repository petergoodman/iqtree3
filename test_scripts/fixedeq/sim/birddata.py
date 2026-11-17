"""SIM-1's bird data and scaffold (CODE_PLAN.md 3.8, "Bird data and scaffold"): nQMaker's Dryad deposit
(doi:10.5061/dryad.3tx95x6hx, CC0), raw and read-only outside the repository; derived data only from here.

  python test_scripts/fixedeq/sim/birddata.py fetch --register FILE [FILE ...] --git git.exe [--raw DIR]
  python test_scripts/fixedeq/sim/birddata.py inventory --git git.exe [--raw DIR] [--out DIR]
  python test_scripts/fixedeq/sim/birddata.py scaffold --binary FROZEN --git git.exe [--threads 6] [--raw DIR] [--out DIR]

fetch: Dryad refuses scripted downloads (API 401 without a token, file stream 403, 2026-11-17), so the files
are downloaded in a browser and registered: each is copied into the raw directory (default
~/iqtree3-data/raw/nqmaker-dryad/), checked against the SHA-256 Dryad publishes (API version 167087, read
2026-11-17), made read-only and recorded in fetch.json there. A registered file is checked, never overwritten.

The deposit's bird set (05_clades/bird): one concatenated protein alignment (alignment.nex inside
alignment.nex.zip; 52 taxa, 4,519,041 columns, X missing, - gap) with 1,000 training and 6,295 test loci as
charsets (train.nex, test.nex); its 52 taxa are 48 birds and the outgroups HUMAN, ALLIG, ANOCA and CHEMY.

inventory: every locus's length, taxa with a residue, gap fraction and residue counts, against QMaker's
"52 taxa, mean 545 columns" (reported); writes the derived table outside the repository and the committed
scaffold/compositions.tsv (locus number, length, 20 counts) for the screen's S2 and SIM-2.
scaffold: nQMaker's bird ASTRAL_NEW topology, rooted on HUMAN (the outgroups make Palaeognathae the wrong
root here; Peter to confirm), lengths dropped; one -p fit of NQ.bird+G4 under -te (each locus its own +G4
shape and partition rate, branch lengths estimated) on a seeded random 100 loci by the frozen binary; writes
scaffold/tree.nwk (the fitted rooted tree, lengths at partition rate 1) and scaffold/loci.tsv (locus,
length, shape, rate), with provenance header lines, for sim/pipeline.py.
"""
import argparse
import datetime
import io
import json
import os
import re
import shutil
import stat
import sys
import zipfile
from collections import defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
for _p in ("", "regression", "probes", "nqc"):
    sys.path.insert(0, str(HERE.parent / _p))
sys.path.insert(0, str(HERE))
from oracle import readers  # noqa: E402
from oracle.likelihood import AMINO  # noqa: E402
from regress import REPO, read, sha256, write  # noqa: E402
from probes import Runner  # noqa: E402
from nqcchecks import check, q_route  # noqa: E402
import common  # noqa: E402

SCRIPT = Path(__file__).resolve().relative_to(REPO).as_posix()
RAW = Path("~/iqtree3-data/raw/nqmaker-dryad")
SCAFFOLD = HERE / "scaffold"
DOI = "10.5061/dryad.3tx95x6hx"
# Dryad's published digests (https://datadryad.org/api/v2/versions/167087/files), the three files SIM-1 uses
DRYAD_SHA256 = {
    "05_clades.zip": "c34baac62d3c720943fb667fe712658256644953bd68fec054c3cfba1de9d03f",
    "ASTRAL_trees.zip": "a5d6cbee380777d8b432c1de82acc6e0a9e2d56f4bd6572c053a1311f265f86e",
    "README.txt": "cf66fedfd2a4e550e5eae737d0a0e1186a368c6d187d7b43e6855a63307ca7c2",
}
BIRD = "05_clades/bird/"
TREE = "ASTRAL trees/bird/ASTRAL_NEW.treefile"
ROOT = "HUMAN"
OUTGROUPS = ("HUMAN", "ALLIG", "ANOCA", "CHEMY")
SEED = 20261117
SCAFFOLD_LOCI = 100
QMAKER = {"taxa": 52, "mean_columns": 545}     # reported by QMaker for the bird set (not verified here)
CODE = np.full(256, -1, dtype=np.int64)
for _k, _c in enumerate(AMINO):
    CODE[ord(_c)] = CODE[ord(_c.lower())] = _k
CHARSET = re.compile(r"^\s*charset\s+(\S+)\s*=\s*(\d+)-(\d+)\s*;\s*$", re.I | re.M)


# ---------------------------------------------------------------- fetch

def register(files, raw):
    """Copy each named file into raw unless it is there, check every file against Dryad's digest, make
    it read-only; returns the record and the problems found."""
    raw.mkdir(parents=True, exist_ok=True)
    record, problems = {}, []
    for src in map(Path, files):
        name = src.name
        if name not in DRYAD_SHA256:
            problems.append(f"{name}: not one of {sorted(DRYAD_SHA256)}")
            continue
        dest = raw / name
        copied = not dest.exists()
        if copied:
            shutil.copyfile(src, dest)
        digest = sha256(dest)
        if digest != DRYAD_SHA256[name]:
            problems.append(f"{name}: SHA-256 {digest} is not Dryad's {DRYAD_SHA256[name]}")
        os.chmod(dest, stat.S_IRUSR | stat.S_IRGRP | stat.S_IROTH)
        record[name] = {"source": str(src.resolve()) if copied else "registered earlier (see the first record)",
                        "sha256": digest, "bytes": dest.stat().st_size, "matches_dryad": digest == DRYAD_SHA256[name]}
    present = {p.name for p in raw.iterdir()}
    problems += [f"{name}: not registered yet" for name in sorted(set(DRYAD_SHA256) - present)]
    return record, problems


def cmd_fetch(a):
    prov = common.provenance(a, SCRIPT, "fetch", "raw data registration", doi=DOI)
    raw = Path(a.raw).expanduser().resolve()
    if raw.is_relative_to(REPO):
        sys.exit("raw data stay outside the repository")
    record, problems = register(a.register, raw)
    log = raw / "fetch.json"
    previous = json.loads(log.read_text(encoding="utf-8")) if log.exists() else {"registrations": []}
    previous["registrations"].append({"provenance": prov, "files": record, "problems": problems,
                                      "finished_utc": datetime.datetime.now(datetime.timezone.utc).isoformat()})
    if log.exists():
        os.chmod(log, stat.S_IRUSR | stat.S_IWUSR | stat.S_IRGRP | stat.S_IROTH)
    write(log, json.dumps(previous, indent=1))
    for p in problems:
        print(f"  {p}")
    print(f"registered {sorted(record)} in {raw}; {len(problems)} problems")
    sys.exit(1 if problems else 0)


# ---------------------------------------------------------------- reading the deposit

def verified(raw):
    bad = [n for n, h in DRYAD_SHA256.items() if not (raw / n).exists() or sha256(raw / n) != h]
    if bad:
        sys.exit(f"{raw}: {bad} missing or not matching Dryad's SHA-256; run fetch first")


def charsets(raw):
    """[(locus, set, start, end)], 1-based inclusive ranges, from train.nex then test.nex."""
    out = []
    with zipfile.ZipFile(raw / "05_clades.zip") as z:
        for kind in ("train", "test"):
            text = z.read(BIRD + f"{kind}.nex").decode("utf-8")
            found = CHARSET.findall(text)
            if len(found) != text.lower().count("charset"):
                raise ValueError(f"{kind}.nex: a charset that is not one plain range")
            out += [(name, kind, int(lo), int(hi)) for name, lo, hi in found]
    return out


def alignment(raw):
    """Taxon names and the matrix (taxa x columns, uint8) of the bird alignment, streamed from the nested zip."""
    with zipfile.ZipFile(raw / "05_clades.zip") as z:
        inner = zipfile.ZipFile(io.BytesIO(z.read(BIRD + "alignment.nex.zip")))
        with inner.open("alignment.nex") as fh:
            text = io.TextIOWrapper(fh, encoding="ascii")
            dims, names, rows, in_matrix = {}, [], [], False
            for line in text:
                s = line.strip()
                if not in_matrix:
                    for key in ("ntax", "nchar"):
                        m = re.search(rf"{key}\s*=\s*(\d+)", s, re.I)
                        if m:
                            dims[key] = int(m.group(1))
                    in_matrix = s.lower() == "matrix"
                    continue
                if s.startswith(";"):
                    break
                if s:
                    name, seq = s.split(None, 1)
                    seq = seq.replace(" ", "")
                    names.append(name)
                    rows.append(np.frombuffer(seq.encode("ascii"), dtype=np.uint8))
    if len(names) != dims.get("ntax") or any(len(r) != dims.get("nchar") for r in rows):
        raise ValueError(f"alignment.nex: expected {dims} as one sequential line per taxon")
    return names, np.vstack(rows)


def astral_tree(raw):
    with zipfile.ZipFile(raw / "ASTRAL_trees.zip") as z:
        return z.read(TREE).decode("utf-8").strip()


def root_on(text, leaf):
    """The topology of a Newick tree rooted on the pendant edge of leaf, without lengths or labels."""
    tree = readers.parse_newick(text)
    adj, by_name = defaultdict(list), {}
    for nd in tree.preorder():
        for ch in nd.children:
            adj[id(nd)].append(ch)
            adj[id(ch)].append(nd)
        if not nd.children:
            by_name[nd.name] = nd
    target = by_name[leaf]
    (parent,) = adj[id(target)]

    def side(nd, came):
        nbrs = [x for x in adj[id(nd)] if x is not came]
        if len(nbrs) == 1:         # a degree-2 node (the old root) is passed through, not kept unary
            return side(nbrs[0], nd)
        return nd.name if not nbrs else "(" + ",".join(side(x, nd) for x in nbrs) + ")"

    return f"({leaf},{side(parent, target)});"


# ---------------------------------------------------------------- inventory

def inventory_table(names, matrix, sets):
    rows = []
    for name, kind, lo, hi in sets:
        block = CODE[matrix[:, lo - 1:hi]]
        counts = np.bincount(block[block >= 0], minlength=20)
        rows.append({"locus": name, "number": name.split(".")[0], "set": kind, "start": lo, "end": hi,
                     "length": hi - lo + 1, "taxa": int(np.count_nonzero((block >= 0).any(axis=1))),
                     "residue_fraction": float(np.count_nonzero(block >= 0) / block.size), "counts": counts.tolist()})
    return rows


def header(prov, content, **extra):
    keys = ("script", "step", "git_commit", "git_clean", "started_utc", "doi")
    lines = [f"# content: {content}"] + [f"# {k}: {prov[k]}" for k in keys]
    lines += [f"# raw {n} sha256: {h}" for n, h in DRYAD_SHA256.items()]
    lines += [f"# {k}: {v}" for k, v in extra.items()]
    return "\n".join(lines) + "\n"


def refuse_existing(*paths):
    found = [str(p) for p in paths if p.exists()]
    if found:
        sys.exit(f"{found} exist; derived data are never regenerated in place (remove them deliberately first)")


def cmd_inventory(a):
    refuse_existing(SCAFFOLD / "compositions.tsv")
    prov = common.provenance(a, SCRIPT, "inventory", "derived data", doi=DOI)
    raw = Path(a.raw).expanduser().resolve()
    verified(raw)
    out = common.new_out(a.out, prov)
    names, matrix = alignment(raw)
    sets = charsets(raw)
    numbers = [name.split(".")[0] for name, _, _, _ in sets]
    if len(set(numbers)) != len(numbers):
        sys.exit("locus numbers repeat across train.nex and test.nex")
    rows = inventory_table(names, matrix, sets)
    lengths = {k: np.array([r["length"] for r in rows if r["set"] == k]) for k in ("train", "test")}
    summary = {"taxa": len(names), "columns": int(matrix.shape[1]), "outgroups_present": [t for t in OUTGROUPS if t in names],
               "loci": {k: len(v) for k, v in lengths.items()},
               "mean_columns": {k: float(v.mean()) for k, v in lengths.items()},
               "mean_columns_all": float(np.concatenate(list(lengths.values())).mean()),
               "loci_with_all_taxa": sum(r["taxa"] == len(names) for r in rows), "qmaker_reported": QMAKER,
               "covered_columns": int(sum(r["length"] for r in rows))}
    write(out / "inventory.json", json.dumps({"provenance": prov, "summary": summary, "loci": rows}, indent=0))
    SCAFFOLD.mkdir(exist_ok=True)
    write(SCAFFOLD / "compositions.tsv",
          header(prov, "every bird locus's length and residue counts (S2 of the screen, SIM-2's compositions)")
          + "locus\tlength\t" + "\t".join(AMINO) + "\n"
          + "".join(f"{r['number']}\t{r['length']}\t" + "\t".join(map(str, r["counts"])) + "\n" for r in rows))
    common.finish(out, prov, "SIM-1 bird inventory", [check("inventory", None, **summary)])


# ---------------------------------------------------------------- scaffold

def cmd_scaffold(a):
    refuse_existing(SCAFFOLD / "tree.nwk", SCAFFOLD / "loci.tsv")
    prov = common.provenance(a, SCRIPT, "scaffold", "derived data", doi=DOI, seed=SEED, loci=SCAFFOLD_LOCI,
                             topology=TREE, root=ROOT, threads=a.threads)
    raw = Path(a.raw).expanduser().resolve()
    verified(raw)
    out = common.new_out(a.out, prov)
    binary = Path(a.binary).expanduser().resolve()
    prov["binary"] = {"path": str(binary), "sha256": sha256(binary)}
    names, matrix = alignment(raw)
    sets = charsets(raw)
    pick = sorted(np.random.default_rng(SEED).choice(len(sets), SCAFFOLD_LOCI, replace=False).tolist())
    chosen = [sets[i] for i in pick]
    blocks = [matrix[:, lo - 1:hi] for _, _, lo, hi in chosen]
    sub = np.hstack(blocks)
    phy = f"{len(names)} {sub.shape[1]}\n" + "".join(f"{n} {row.tobytes().decode('ascii')}\n" for n, row in zip(names, sub))
    sets_text, start = [], 1
    for (name, _, lo, hi) in chosen:
        sets_text.append(f"charset L{name.split('.')[0]} = {start}-{start + hi - lo};")
        start += hi - lo + 1
    parts = "#nexus\nbegin sets;\n" + "\n".join(sets_text) + "\ncharpartition scaffold = " + ", ".join(
        f"NQ.bird+G4: L{name.split('.')[0]}" for name, _, _, _ in chosen) + ";\nend;\n"
    rooted = root_on(astral_tree(raw), ROOT)
    write(out / "sub.phy", phy)
    run = Runner(binary, out, a.timeout).run("fit", ["-s", out / "sub.phy", "-p", "parts.nex", "-te", "rooted.nwk",
                                                     "-keep-ident"], {"parts.nex": parts, "rooted.nwk": rooted + "\n"},
                                            threads=a.threads)
    rows = [check("scaffold fit exits 0", run.code == 0, exit_code=run.code, wall_s=run.wall)]
    labels = [f"L{n.split('.')[0]}" for n, _, _, _ in chosen]
    lengths = [hi - lo + 1 for _, _, lo, hi in chosen]
    # shapes and partition rates from the checkpoint (best_model.nex holds each partition's tree length in braces)
    route = q_route(run, labels) if run.code == 0 else None
    rows.append(check("the fit's checkpoint holds the tree, every shape and the partition rates", route is not None))
    if route is not None:
        # IQ-TREE writes a rooted tree with a stem, the branch to its virtual root; the root is drawn from pi
        # solved from NQ.bird's Q, stationary, so the stem changes no likelihood and its length is not
        # identified. It is dropped (recorded), since AliSim would simulate it as a branch.
        fitted = route["tree"].strip()
        stem = readers.parse_newick(fitted).length
        tree = re.sub(r"\)\s*:\s*[-+0-9.eE]+\s*;\s*$", ");", fitted)
        root = readers.parse_newick(tree)
        rows.append(check("the fitted tree's root stem, dropped (report)", None, stem=stem))
        shapes = [float(v) for v in route["shapes"]]
        rates = [float(v) for v in route["rates"]]
        mean = float(np.dot(rates, lengths) / np.sum(lengths))
        rows += [check("the fitted tree is rooted, without a stem, on all 52 taxa",
                       len(root.children) == 2 and not root.length and len(root.leaves()) == len(names)),
                 check("one shape and one partition rate per scaffold locus",
                       len(shapes) == len(rates) == SCAFFOLD_LOCI, shapes=len(shapes), rates=len(rates)),
                 check("the partition rates' site-weighted mean is 1 (IQ-TREE's -p normalization)",
                       abs(mean - 1.0) <= 1e-6, mean=mean)]
        if all(r["passed"] is not False for r in rows):     # reported rows (None) do not block
            SCAFFOLD.mkdir(exist_ok=True)
            extra = {"seed": SEED, "loci": SCAFFOLD_LOCI, "topology": TREE, "root": ROOT,
                     "binary_sha256": prov["binary"]["sha256"]}
            write(SCAFFOLD / "tree.nwk", header(prov, f"bird {TREE} rooted on {ROOT}, lengths from the scaffold fit",
                                                **extra) + tree + "\n")
            write(SCAFFOLD / "loci.tsv", header(prov, "the scaffold fit's loci: length, +G4 shape, partition rate "
                                                      "(checkpoint, 10 significant digits)", **extra)
                  + "locus\tlength\tshape\trate\n"
                  + "".join(f"{loc}\t{n}\t{s!r}\t{r!r}\n" for loc, n, s, r in zip(labels, lengths, shapes, rates)))
            rows.append(check("scaffold written", None, tree_length=sum(n.length or 0 for n in root.preorder()),
                              shapes=shapes, rates=rates))
    common.finish(out, prov, "SIM-1 bird scaffold", rows, {"loci": [c[0] for c in chosen], "rooted_input": rooted})


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="step", required=True)
    f = sub.add_parser("fetch")
    f.add_argument("--register", nargs="+", required=True, help="the files downloaded in a browser")
    i = sub.add_parser("inventory")
    s = sub.add_parser("scaffold")
    s.add_argument("--binary", required=True, help="the frozen unmodified binary")
    s.add_argument("--threads", type=int, default=6)
    s.add_argument("--timeout", type=float, default=4 * 3600)
    for p in (f, i, s):
        p.add_argument("--raw", default=str(RAW))
        p.add_argument("--git", default="git")
        p.add_argument("--allow-dirty", action="store_true")
    for p in (i, s):
        p.add_argument("--out", help="new directory outside the repository (default under ~/iqtree3-runs/sim1/)")
    a = ap.parse_args()
    {"fetch": cmd_fetch, "inventory": cmd_inventory, "scaffold": cmd_scaffold}[a.step](a)


if __name__ == "__main__":
    main()
