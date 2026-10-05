#!/usr/bin/env python3
"""S0 runtime probes (a), (b), (c), (e), (f) and (g) on an IQ-TREE binary (docs/agent/CODE_PLAN.md
section 4). NONREV stands in for NQC, which does not exist yet.

  python probes.py --binary ~/iqtree3-baseline/iqtree3 --out DIR [--only a,b,...] [--git git.exe]

Each probe's outcome as predicted from reading the source is written below before any run. A probe
whose observations disagree is marked "differs"; that is a finding to report, not something to
tune away. Exit status 1 if any probe differs or a run fails unexpectedly.
"""
import argparse
import datetime
import gzip
import json
import platform
import re
import subprocess
import sys
import time
import traceback
from collections import defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE.parent / "regression"))
from oracle import readers, target  # noqa: E402
from oracle.readers import Node  # noqa: E402
from regress import FREQ, REPO, capture, compare_run, extract_run, read, sha256, split_sections, write  # noqa: E402

SCRIPT = Path(__file__).resolve().relative_to(REPO).as_posix()
AA = REPO / "example" / "aa_example.phy"
TURTLE = REPO / "test_scripts" / "test_data" / "turtle_aa.fasta"
NEX = REPO / "test_scripts" / "test_data" / "turtle_aa.nex"
BASELINE = REPO / "test_scripts" / "fixedeq" / "regression" / "baseline" / "baseline.json"
PREFIX = "probe"
SEED = 20261004
FROZEN_TOL = 1e-10    # branch lengths are printed with 10 decimals
SPLIT_TOL = 1e-9      # root-adjacent lengths against the midpoint prediction
CKP_DIGITS = 10       # utils/checkpoint.h:22
BRACKET_ERROR = "Close bracket not found in +F{0.08e"

PREDICTIONS = {
    "a": "The named route (+F<name> from --mdef) delivers the vector exactly as the literal +F{...} does, "
         "under -m and under --model-joint with -p and -S (modelfactory.cpp:599-604; 284-305): a1 matches "
         "baseline run 3 in every compared item, a2 = a3 and a4 = a5 item by item, and every checkpoint "
         "state_freq equals the vector at 10 significant digits.",
    "b": "The +F split is not brace-aware (modelfactory.cpp:258-267, 295-303): b1 and b3 exit with code 2 and "
         f"'{BRACKET_ERROR}'; b2 (8e-02) and b4 (the e+ vector through --mdef) match baseline run 3.",
    "c": "Level 1 with -te, -blfix and braced rate parameters freezes every nuisance for one alignment (c1), "
         "under -S (c2) and under -q (c4; c5 with partition rates fixed through the charpartition): topology, "
         "branch lengths, gamma shape, p_invar and partition rates unchanged. Under -p (c3) the partition rates "
         "are re-optimized (partitionmodelplen.cpp:150-157), so Level 1 is not frozen there.",
    "e": "-S differs from -p beyond tree linkage: partition rates only under -p; no padding of missing taxa in "
         "the pooled frequencies under -S (partitionmodel.cpp:132-133); -S roots each partition tree separately; "
         "the free-parameter count adds partition rates only under -p; -S writes one tree per partition. The "
         "+I+G restart (modelfactory.cpp:1379-1538) and identical-sequence handling also differ but are not "
         "exercised by these runs.",
    "f": "Under -te the root never moves to another branch (no tree search; --root-find off). A rooted input "
         "keeps its root; an unrooted input is rooted once at the midpoint of the longest path "
         "(phylotree.cpp:5906-5960), per partition under -S and on the shared tree under -p and -q. At Level 1 "
         "the root-adjacent lengths are unchanged (or the midpoint split); at Level 2 the root may slide along "
         "its own branch (phylotree.cpp:2639-2640).",
    "g": "g0, fitted on a fixed tree (-te), has no final model optimization (phyloanalysis.cpp:3880), so its "
         "checkpoint holds the final model: every checkpoint rate equals the report's PAML block within the "
         "block's rounding (5e-7) and the checkpoint gamma shape equals the report's within 5e-5. On g0: gB "
         "(checkpoint rates, 10 digits, -m FILE) and gC (GTR20{...}) agree within 1e-9 in lnL; gA, gA2 and gA3 "
         "(the report's block, 6 decimal places) differ from gB by more than the spread of gD (gB with rates "
         "perturbed within their 10th digit); |gB - g0's printed lnL| <= 1e-4. gS, fitted with a tree search, "
         "runs the final model optimization, which no model checkpoint save follows (phyloanalysis.cpp:3895-3899; "
         "iqtree.cpp:134-150): its checkpoint rates differ from its report block, and gS_B, its checkpoint model "
         "on its final tree, lies within 0.01 of the final step's printed initial lnL and at least 0.1 below its "
         "printed final lnL. (The gS prediction was written after the 2026-10-05 dry run observed this.)",
}


# ---------------------------------------------------------------- helpers (unit-tested)

def read_checkpoint(path):
    """A checkpoint as a flat 'Struct!key' -> value map, following Checkpoint::dump (utils/checkpoint.cpp:128-144)."""
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        lines = fh.read().splitlines()
    flat, struct = {}, None
    for line in lines:
        if line.startswith("---") or not line.strip():
            continue
        if line.startswith(" "):
            key, _, value = line[1:].partition(": ")
            flat[f"{struct}!{key}"] = value
        elif line.endswith(":"):
            struct = line[:-1]
        else:
            key, _, value = line.partition(": ")
            flat[key] = value
    return flat


def ckp_values(value):
    return [float(v) for v in value.split(",")]


def leafset(node):
    return frozenset(n.name for n in node.leaves())


def canonical(side, all_leaves):
    """The side of a split that does not hold the alphabetically first taxon."""
    return all_leaves - side if min(all_leaves) in side else side


def edge_map(tree):
    """Unrooted split -> branch length; a rooted tree's two root edges are one split, summed."""
    all_leaves, out = leafset(tree), {}
    for node in tree.preorder():
        if node is not tree:
            s = canonical(leafset(node), all_leaves)
            out[s] = out.get(s, 0.0) + node.length
    return out


def root_info(tree):
    """The root split (canonical side) and the two root-adjacent lengths, canonical side first; None if unrooted."""
    if len(tree.children) != 2:
        return None
    a, b = tree.children
    side = canonical(leafset(a), leafset(tree))
    return side, ((a.length, b.length) if side == leafset(a) else (b.length, a.length))


def unroot(tree):
    """Join the two root edges into one, under the internal root child."""
    a, b = tree.children
    if not a.children:
        a, b = b, a
    return Node(children=list(a.children) + [Node(b.name, a.length + b.length, b.children)])


def to_newick(tree):
    def rec(n):
        s = "(" + ",".join(rec(c) for c in n.children) + ")" if n.children else ""
        s += n.name or ""
        return s + (f":{n.length:.10f}" if n is not tree and n.length is not None else "")
    return rec(tree) + ";"


def midpoint(tree, start):
    """Where convertToRooted puts the root of an unrooted tree (phylotree.cpp:5918-5934; node.cpp:163-170):
    u farthest from the leaf `start`, v farthest from u, the root on the path u-v at half its length, at
    distance (half - curlen) from the u end of its edge. Returns root_info's (split, lengths)."""
    nodes = list(tree.preorder())
    idx = {id(n): i for i, n in enumerate(nodes)}
    adj = defaultdict(list)
    for n in nodes:
        for c in n.children:
            adj[idx[id(n)]].append((idx[id(c)], c.length))
            adj[idx[id(c)]].append((idx[id(n)], c.length))
    leaf = {i: n.name for i, n in enumerate(nodes) if not n.children}

    def sweep(src):
        dist, parent, stack = {src: 0.0}, {src: None}, [src]
        while stack:
            x = stack.pop()
            for y, w in adj[x]:
                if y not in dist:
                    dist[y], parent[y] = dist[x] + w, x
                    stack.append(y)
        far = max((i for i in leaf), key=lambda i: dist[i])
        return far, dist, parent

    s = next(i for i, name in leaf.items() if name == start)
    u, _, _ = sweep(s)
    v, dist, parent = sweep(u)
    path = [v]
    while parent[path[-1]] is not None:
        path.append(parent[path[-1]])
    path.reverse()  # u ... v
    half, curlen, k = dist[v] / 2, 0.0, 0
    weight = {(x, y): w for x in adj for y, w in adj[x]}
    while path[k] != v and curlen + weight[(path[k], path[k + 1])] < half:
        curlen += weight[(path[k], path[k + 1])]
        k += 1
    x, y = path[k], path[k + 1]
    node_len = half - curlen
    dad_len = weight[(x, y)] - node_len
    # leaves on x's side of the edge x-y
    seen, stack = {x, y}, [x]
    side = set()
    while stack:
        z = stack.pop()
        if z in leaf:
            side.add(leaf[z])
        for w_, _ in adj[z]:
            if w_ not in seen:
                seen.add(w_)
                stack.append(w_)
    all_leaves = frozenset(leaf.values())
    side = frozenset(side)
    canon = canonical(side, all_leaves)
    return canon, ((node_len, dad_len) if canon == side else (dad_len, node_len))


def paml_rows(rates, n=20):
    """IQ-TREE's internal pair order to PAML's lower-triangle rows (model/modelprotein.cpp:1281-1302)."""
    return [[rates[col * (2 * n - col - 1) // 2 + (row - col - 1)] for col in range(row)] for row in range(1, n)]


def mdef(name, vector):
    return f"#nexus\nbegin models;\nfrequency {name} = {vector};\nend;\n"


def with_charpartition(nex_text, parts):
    """turtle_aa.nex plus a charpartition of (model, charset, rate or None) entries (nclextra/msetsblock.cpp:184-242)."""
    entries = ", ".join(f"{m}: {cs}" + (f"{{{r}}}" if r is not None else "") for m, cs, r in parts)
    return nex_text.replace("end;", f"charpartition probe = {entries};\nend;", 1)


def partition_table(report_text):
    """The partition rows of the SUBSTITUTION PROCESS section: (third column's title, [(model, value, parameters)])."""
    section = split_sections(report_text).get("SUBSTITUTION PROCESS", "")
    lines = section.splitlines()
    head = next((i for i, l in enumerate(lines) if re.match(r"^\s+ID\s+Model\s+", l)), None)
    if head is None:
        return None, []
    title = lines[head].split()[2]
    rows = []
    for line in lines[head + 1:]:
        m = re.match(r"^\s+\d+\s+(\S+)\s+(\S+)\s+(.*)$", line)
        if not m:
            break
        rows.append((m.group(1), m.group(2), m.group(3).strip()))
    return title, rows


def turtle_target():
    """Decision 023's test target at 17 significant digits."""
    pi = target.pooled_frequencies(readers.read_fasta(TURTLE), readers.read_charsets(NEX), "p")
    return ",".join(f"{v:.17g}" for v in pi)


# ---------------------------------------------------------------- running

class Run:
    def __init__(self, rundir, name, args, code, wall):
        self.dir, self.name, self.args, self.code, self.wall = rundir, name, args, code, wall

    def file(self, suffix):
        return self.dir / f"{PREFIX}{suffix}"

    def stdout(self):
        return read(self.dir / "stdout.txt")

    def report(self):
        p = self.file(".iqtree")
        return read(p) if p.exists() else ""

    def trees(self):
        p = self.file(".treefile")
        return [l for l in read(p).splitlines() if l.strip()] if p.exists() else []

    def checkpoint(self):
        p = self.file(".ckp.gz")
        return read_checkpoint(p) if p.exists() else {}

    def extract(self):
        return extract_run(self.dir, PREFIX, self.code)

    def record(self):
        return {"name": self.name, "args": self.args, "exit_code": self.code, "wall_s": self.wall}


class Runner:
    def __init__(self, binary, out, timeout):
        self.binary, self.out, self.timeout, self.runs = binary, out, timeout, []

    def run(self, name, args, files=None):
        rundir = self.out / name
        rundir.mkdir()
        for fname, text in (files or {}).items():
            write(rundir / fname, text)
        argv = [str(self.binary), *[str(a) for a in args], "-seed", "1", "-T", "1", "--prefix", PREFIX]
        t0 = time.perf_counter()
        with open(rundir / "stdout.txt", "w", encoding="utf-8") as fh:
            try:
                code = subprocess.run(argv, cwd=rundir, stdout=fh, stderr=subprocess.STDOUT,
                                      timeout=self.timeout).returncode
            except subprocess.TimeoutExpired:
                code = "timeout"
        r = Run(rundir, name, [str(a) for a in args], code, round(time.perf_counter() - t0, 1))
        self.runs.append(r.record())
        print(f"  {name}: exit {code}, {r.wall} s", flush=True)
        return r


def baseline_runs():
    return json.loads(read(BASELINE))["runs"]


def tree_text(base, run):
    return base[run]["items"]["(treefile)"]["text"]


def unrooted_text(text):
    return "\n".join(to_newick(unroot(readers.parse_newick(l))) for l in text.splitlines() if l.strip()) + "\n"


def exact(rec):
    return {"items": {k: {"mode": "exact", "text": v} for k, v in rec["items"].items()}}


def freq_problems(run, vector_text):
    """Checkpoint state_freq entries that differ from the vector at 10 significant digits, or None when
    the run wrote no state_freq."""
    want = [float(f"{float(v):.{CKP_DIGITS}g}") for v in re.split(r"[,\s/]+", vector_text.strip()) if v]
    found = {k: ckp_values(v) for k, v in run.checkpoint().items() if k.endswith("state_freq")}
    return [k for k, v in found.items() if v != want] if found else None


def lnl_evidence(run):
    text = run.stdout()
    first = re.search(r"Initial log-likelihood: (\S+)", text)
    optimal = re.findall(r"Optimal log-likelihood: (\S+)", text)
    final = re.search(r"^Log-likelihood of the tree: (\S+)", run.report(), re.M)
    return {"initial": first.group(1) if first else None, "optimal": optimal[-1] if optimal else None,
            "final": final.group(1) if final else None}


def tree_comparison(in_line, out_line, start=None):
    tin, tout = readers.parse_newick(in_line), readers.parse_newick(out_line)
    ein, eout = edge_map(tin), edge_map(tout)
    common = [s for s in ein if s in eout]
    ratios = [eout[s] / ein[s] for s in common if ein[s] > 1e-3]
    rin, rout = root_info(tin), root_info(tout)
    expected = rin if rin else (midpoint(tin, start) if start else None)
    return {
        "same_topology": set(ein) == set(eout),
        "max_abs_change": max(abs(eout[s] - ein[s]) for s in common) if common else None,
        "ratio_min": min(ratios) if ratios else None, "ratio_max": max(ratios) if ratios else None,
        "input_rooted": rin is not None,
        "root_expected": [sorted(expected[0]), list(expected[1])] if expected else None,
        "root_out": [sorted(rout[0]), list(rout[1])] if rout else None,
        "root_split_kept": bool(expected and rout and rout[0] == expected[0]),
        "root_lengths_kept": bool(expected and rout and all(abs(a - b) <= SPLIT_TOL for a, b in zip(rout[1], expected[1]))),
    }


def first_taxa(path, charsets=None):
    """The first taxon IQ-TREE's midpoint sweep starts from: of the alignment, or of each partition with its
    all-unknown taxa removed (alignment/superalignment.cpp:495-499)."""
    aln = readers.read_phylip(path) if path.suffix == ".phy" else readers.read_fasta(path)
    if charsets is None:
        return [next(iter(aln))]
    out = []
    for sites in charsets.values():
        part = {n: "".join(s[i] for i in sites) for n, s in aln.items()}
        out.append(next(iter(target.remove_gappy(part))))
    return out


def nuisances(run):
    """Rate parameters from the checkpoint (10 digits), else from the report's partition table or its
    single-model lines (rategamma.cpp:248-252, rateinvar.cpp:120-122)."""
    ck, report = run.checkpoint(), run.report()
    title, rows = partition_table(report)
    gamma = [float(v) for k, v in ck.items() if k.endswith("gamma_shape")]
    pinv = [float(v) for k, v in ck.items() if k.endswith("p_invar")]
    source = "checkpoint"
    if not gamma:
        gamma = [float(g) for r in rows for g in re.findall(r"G4\{([^}]*)\}", r[2])]
        gamma += [float(v) for v in re.findall(r"^Gamma shape alpha: (\S+)", report, re.M)]
        source = "report"
    if not pinv:
        pinv = [float(v) for v in re.findall(r"^Proportion of invariable sites: (\S+)", report, re.M)]
    return {"gamma_shape": gamma, "p_invar": pinv, "rate_source": source,
            "part_rates_ckp": ck.get("PartitionModelPlen!part_rates"),
            "table_column": title, "table_values": [r[1] for r in rows]}


# ---------------------------------------------------------------- probes

def probe_a(r, ctx):
    base = ctx["base"]
    f = {"freq.nex": mdef("TESTFREQ", FREQ)}
    a1 = r.run("a1", ["-s", AA, "-m", "NONREV+FTESTFREQ", "--mdef", "freq.nex"], f)
    runs, comparisons = [a1], {"a1 against baseline run03": compare_run(base["run03"], a1.extract())}
    for mode, run, lit_name, named_name in (("-p", "run08", "a2", "a3"), ("-S", "run09", "a4", "a5")):
        tree = {"input.tree": tree_text(base, run)}
        common = ["-s", TURTLE, mode, NEX, "-te", "input.tree", "-blfix", "--model-joint"]
        lit = r.run(lit_name, common + ["NONREV+F{" + FREQ + "}"], tree)
        named = r.run(named_name, common + ["NONREV+FTESTFREQ", "--mdef", "freq.nex"], {**tree, **f})
        comparisons[f"{lit_name} against {named_name}"] = compare_run(exact(lit.extract()), named.extract())
        runs += [lit, named]
    obs = {"exit_codes": {x.name: x.code for x in runs}, "comparisons": comparisons,
           "state_freq_problems": {x.name: freq_problems(x, FREQ) for x in runs}}
    ok = (all(x.code == 0 for x in runs) and not any(comparisons.values())
          and all(v == [] for v in obs["state_freq_problems"].values()))
    return obs, ok


def probe_b(r, ctx):
    base = ctx["base"]
    head, rest = FREQ.split(",", 1)
    bad = "0.08e+00," + rest
    b1 = r.run("b1", ["-s", AA, "-m", "NONREV+F{" + bad + "}"])
    b2 = r.run("b2", ["-s", AA, "-m", "NONREV+F{8e-02," + rest + "}"])
    b3 = r.run("b3", ["-s", TURTLE, "-p", NEX, "--model-joint", "NONREV+F{" + bad + "}"])
    b4 = r.run("b4", ["-s", AA, "-m", "NONREV+FTESTFREQ", "--mdef", "freq.nex"], {"freq.nex": mdef("TESTFREQ", bad)})
    obs = {}
    for x in (b1, b3):
        line = next((l for l in x.stdout().splitlines() if "ERROR" in l), None)
        obs[x.name] = {"exit_code": x.code, "error_line": line, "expected_message": BRACKET_ERROR in x.stdout()}
    for x in (b2, b4):
        obs[x.name] = {"exit_code": x.code, "vs_run03": compare_run(base["run03"], x.extract())}
    ok = (all(obs[n]["exit_code"] == 2 and obs[n]["expected_message"] for n in ("b1", "b3"))
          and all(obs[n]["exit_code"] == 0 and not obs[n]["vs_run03"] for n in ("b2", "b4")))
    return obs, ok


def joint_run(r, name, mode, rate, tree, blfix, partfile=None):
    """A --model-joint NONREV run on turtle. The rate model is given by -m, or, if that is rejected, by a
    charpartition (the plan's pre-declared fallback); `partfile` (name, text) replaces the partition file
    from the start. Returns the run and the attempts made."""
    extra = ["-blfix"] if blfix else []
    files = {"input.tree": tree}

    def args(nex):
        return ["-s", TURTLE, mode, nex, "--model-joint", "NONREV", "-te", "input.tree"] + extra

    if partfile:
        return r.run(name, args(partfile[0]), {**files, partfile[0]: partfile[1]}), ["charpartition"]
    first = r.run(name, args(NEX) + ["-m", rate], files)
    if first.code == 0:
        return first, ["-m"]
    cp = with_charpartition(read(NEX), [(rate, cs, None) for cs in readers.read_charsets(NEX)])
    return (r.run(name + "_charpartition", args("parts.nex"), {**files, "parts.nex": cp}),
            ["-m (failed)", "charpartition"])


def judge_frozen(run, in_text, gamma=None, pinv=None, rates=None):
    """Whether every nuisance of a Level 1 run kept its input value."""
    outs = run.trees()
    ins = [l for l in in_text.splitlines() if l.strip()]
    trees = [tree_comparison(i, o) for i, o in zip(ins, outs)] if len(outs) == len(ins) else []
    nz = nuisances(run)
    checks = {
        "exit_ok": run.code == 0,
        "trees": bool(trees) and all(t["same_topology"] and t["max_abs_change"] <= FROZEN_TOL for t in trees),
    }
    if gamma is not None:
        checks["gamma"] = bool(nz["gamma_shape"]) and all(v == gamma for v in nz["gamma_shape"])
    if pinv is not None:
        checks["p_invar"] = bool(nz["p_invar"]) and all(v == pinv for v in nz["p_invar"])
    if rates is not None:
        got = [float(v) for v in nz["table_values"]]
        checks["part_rates"] = len(got) == len(rates) and all(abs(a - b) <= 5e-5 + 1e-12 for a, b in zip(got, rates))
    return {"checks": checks, "trees": trees, "nuisances": nz, "lnl": lnl_evidence(run)}, all(checks.values())


def probe_c(r, ctx):
    base = ctx["base"]
    t02, t08, t09 = (tree_text(base, k) for k in ("run02", "run08", "run09"))
    sizes = [len(s) for s in readers.read_charsets(NEX).values()]
    _, rows = partition_table(base["run08"]["items"]["SUBSTITUTION PROCESS"]["text"])
    speeds = [float(v) for _, v, _ in rows]
    mean = sum(s * n for s, n in zip(speeds, sizes)) / sum(sizes)
    normalized = [float(f"{s / mean:.4f}") for s in speeds]
    rates_nex = with_charpartition(read(NEX), [("LG+G4{0.5}", cs, f"{s:.4f}")
                                               for cs, s in zip(readers.read_charsets(NEX), speeds)])
    c1 = r.run("c1", ["-s", AA, "-m", "NONREV+I{0.2}+G4{0.5}", "-te", "input.tree", "-blfix"], {"input.tree": t02})
    c2, at2 = joint_run(r, "c2", "-S", "LG+G4{0.5}", t09, True)
    c3, at3 = joint_run(r, "c3", "-p", "LG+G4{0.5}", t08, True)
    c4, at4 = joint_run(r, "c4", "-q", "LG+G4{0.5}", t08, True)
    c5, at5 = joint_run(r, "c5", "-q", None, t08, True, partfile=("rates.nex", rates_nex))
    ctx.update(c1=c1, c2=c2, c3=c3, c4=c4)
    obs, oks = {"attempts": {"c2": at2, "c3": at3, "c4": at4, "c5": at5},
                "c5_input_speeds": speeds, "c5_expected_speeds": normalized}, []
    for run, tree, kw in ((c1, t02, dict(gamma=0.5, pinv=0.2)), (c2, t09, dict(gamma=0.5)),
                          (c4, t08, dict(gamma=0.5, rates=[1.0, 1.0, 1.0])),
                          (c5, t08, dict(gamma=0.5, rates=normalized))):
        obs[run.name], ok = judge_frozen(run, tree, **kw)
        oks.append(ok)
    obs[c3.name], frozen3 = judge_frozen(c3, t08, gamma=0.5)
    rates3 = c3.checkpoint().get("PartitionModelPlen!part_rates")
    moved3 = rates3 is not None and any(abs(v - 1.0) > 1e-6 for v in ckp_values(rates3))
    obs["c3_part_rates_moved"] = moved3
    return obs, all(oks) and c3.code == 0 and moved3 and not frozen3


def probe_f(r, ctx):
    base = ctx["base"]
    t02, t08, t09 = (tree_text(base, k) for k in ("run02", "run08", "run09"))
    u02, u08, u09 = unrooted_text(t02), unrooted_text(t08), unrooted_text(t09)
    aa_start = first_taxa(AA)
    sup_start = first_taxa(TURTLE)
    part_start = first_taxa(TURTLE, readers.read_charsets(NEX))
    runs = {
        "f1=c1": (ctx.get("c1"), t02, aa_start, 1),
        "f2": (r.run("f2", ["-s", AA, "-m", "NONREV+G4", "-te", "input.tree"], {"input.tree": t02}), t02, aa_start, 2),
        "f3": (r.run("f3", ["-s", AA, "-m", "NONREV+I{0.2}+G4{0.5}", "-te", "input.tree", "-blfix"],
                     {"input.tree": u02}), u02, aa_start, 1),
        "f4": (r.run("f4", ["-s", AA, "-m", "NONREV+G4", "-te", "input.tree"], {"input.tree": u02}), u02, aa_start, 2),
        "f5": (joint_run(r, "f5", "-p", "LG+G4", t08, False)[0], t08, sup_start, 2),
        "f6": (joint_run(r, "f6", "-p", "LG+G4", u08, False)[0], u08, sup_start, 2),
        "f-c3": (ctx.get("c3"), t08, sup_start, 1),
        "f-c4": (ctx.get("c4"), t08, sup_start, 1),
        "f7": (joint_run(r, "f7", "-q", "LG+G4{0.5}", u08, True)[0], u08, sup_start, 1),
        "f-c2": (ctx.get("c2"), t09, part_start, 1),
        "f8": (joint_run(r, "f8", "-S", "LG+G4", t09, False)[0], t09, part_start, 2),
        "f9": (joint_run(r, "f9", "-S", "LG+G4{0.5}", u09, True)[0], u09, part_start, 1),
        "f10": (joint_run(r, "f10", "-S", "LG+G4", u09, False)[0], u09, part_start, 2),
    }
    obs, oks = {}, []
    for label, (run, in_text, starts, level) in runs.items():
        if run is None:
            obs[label] = "not run (needs probe c)"
            continue
        ins = [l for l in in_text.splitlines() if l.strip()]
        outs = run.trees()
        if run.code != 0 or len(outs) != len(ins):
            obs[label] = {"run": run.name, "exit_code": run.code, "trees_out": len(outs)}
            oks.append(False)
            continue
        starts = starts if len(starts) == len(ins) else starts[:1] * len(ins)
        comps = [tree_comparison(i, o, s) for i, o, s in zip(ins, outs, starts)]
        kept = all(c["root_split_kept"] for c in comps)
        lengths = all(c["root_lengths_kept"] for c in comps)
        ok = kept and (lengths if level == 1 and label != "f-c3" else True)
        obs[label] = {"run": run.name, "level": level, "root_split_kept": kept, "root_lengths_kept": lengths,
                      "trees": comps}
        oks.append(ok)
    return obs, all(oks)


def probe_e(r, ctx):
    raw = ctx["baseline_raw"]
    obs = {"inputs": {}}
    files = {}
    for run in ("run08", "run09"):
        d = raw / run
        for suffix in (".iqtree", ".treefile", ".ckp.gz", ".log"):
            p = d / f"{run}{suffix}"
            files[(run, suffix)] = p
            obs["inputs"][str(p)] = sha256(p) if p.exists() else None
    if not all(p.exists() for p in files.values()):
        return {**obs, "error": "baseline raw outputs not found"}, False
    rep = {k: read(files[(k, ".iqtree")]) for k in ("run08", "run09")}
    ck = {k: read_checkpoint(files[(k, ".ckp.gz")]) for k in ("run08", "run09")}
    trees = {k: [l for l in read(files[(k, ".treefile")]).splitlines() if l.strip()] for k in ("run08", "run09")}
    cols = {k: partition_table(rep[k])[0] for k in rep}
    items = {}
    items["partition rates only under -p"] = ("PartitionModelPlen!part_rates" in ck["run08"]
                                              and not any(k.endswith("part_rates") for k in ck["run09"])
                                              and cols == {"run08": "Speed", "run09": "TreeLen"})
    aln, cs = readers.read_fasta(TURTLE), readers.read_charsets(NEX)
    pad = float(np.max(np.abs(target.pooled_frequencies(aln, cs, "p") - target.pooled_frequencies(aln, cs, "S"))))
    obs["pooled_p_minus_S_max_abs"] = pad
    items["no padding under -S (pooled vectors differ)"] = pad > 0
    infos = {k: [root_info(readers.parse_newick(t)) for t in trees[k]] for k in trees}
    roots = {k: [sorted(i[0]) if i else None for i in v] for k, v in infos.items()}
    obs["root_splits"] = roots
    items["one rooted tree under -p, one rooted tree per partition under -S"] = (
        len(trees["run08"]) == 1 and len(trees["run09"]) == len(cs) and all(roots["run09"]) and all(roots["run08"]))
    df = {}
    for k in ("run08", "run09"):
        m = re.search(r"^Number of free parameters \(#branches \+ #model parameters\): (\d+)", rep[k], re.M)
        edges = sum(len(edge_map(readers.parse_newick(t))) + 1 for t in trees[k])  # a rooted tree's root split counts twice
        df[k] = {"free_parameters": int(m.group(1)) if m else None, "edges": edges,
                 "residual": (int(m.group(1)) - 379 - edges) if m else None}
    obs["free_parameters"] = df
    items["partition rates counted only under -p"] = df["run09"]["residual"] == 0 and (df["run08"]["residual"] or 0) > 0
    obs["items"] = items
    obs["not_exercised"] = ["+I+G restart under -S", "identical sequences kept under -S",
                            "'logl got worse' abort skipped under -S"]
    return obs, all(items.values())


FINAL_OPT = "Performs final model parameters optimization"


def fitted_gtr20(run):
    """A GTR20+G4 fit's model as the checkpoint and the report hold it, and the log-likelihoods it printed."""
    ck, report, out = run.checkpoint(), run.report(), run.stdout()
    rates = ck.get("ModelProtein!rates", "").split(", ")
    alpha = next((v for k, v in ck.items() if k.endswith("gamma_shape")), None)
    lines = split_sections(report).get("SUBSTITUTION PROCESS", "").splitlines()
    at = next((i for i, l in enumerate(lines) if l.startswith("Substitution parameters (lower-diagonal)")), None)
    block = [l for l in lines[at + 1:] if l.strip()][:20] if at is not None else []
    if run.code != 0 or len(rates) != 190 or alpha is None or len(block) != 20:
        return None
    printed = [float(t) for row in block[:19] for t in row.split()]
    from_ckp = [float(v) for row in paml_rows(rates) for v in row]
    report_alpha = re.search(r"^Gamma shape alpha: (\S+)", report, re.M).group(1)
    final_opt = out.find(FINAL_OPT)
    initial = re.search(r"^1\. Initial log-likelihood: (\S+)", out[final_opt:], re.M) if final_opt >= 0 else None
    return {"rates": rates, "alpha": alpha, "block": block, "report_alpha": report_alpha,
            "final_printed": re.search(r"^Log-likelihood of the tree: (\S+)", report, re.M).group(1),
            "final_opt_ran": final_opt >= 0, "final_opt_initial": initial.group(1) if initial else None,
            "candidate_00": ck.get("CandidateSet!00", "").split(" ")[0] or None,
            "rates_off_block": sum(abs(a - b) > 5e-7 + 1e-12 for a, b in zip(from_ckp, printed)),
            "alpha_off_report": abs(float(alpha) - float(report_alpha)) > 5e-5 + 1e-12}


def probe_g(r, ctx):
    pi = turtle_target()
    model = "GTR20+F{" + pi + "}+G4"
    g0 = r.run("g0", ["-s", AA, "-m", model, "-te", "input.tree"], {"input.tree": tree_text(ctx["base"], "run05")})
    gs = r.run("gS", ["-s", AA, "-m", model])
    fits = {"g0": fitted_gtr20(g0), "gS": fitted_gtr20(gs)}
    obs = {"target": pi, "fit_exit_codes": {"g0": g0.code, "gS": gs.code},
           "fits": {k: (None if f is None else {x: f[x] for x in f if x not in ("rates", "block")})
                    for k, f in fits.items()}}
    if None in fits.values():
        obs["error"] = "a fit gave no checkpoint rates, gamma shape or report block"
        return obs, False
    pif = "+F{" + pi + "}"
    pi_list = pi.split(",")

    def paml(values):
        return "\n".join(" ".join(row) for row in paml_rows(values)) + "\n" + " ".join(pi_list) + "\n"

    def show(name, fit_run, model_text, files):
        tree = {"input.tree": "\n".join(fit_run.trees()) + "\n"}
        run = r.run(name, ["-s", AA, "-m", model_text, "-te", "input.tree", "--show-lh"], {**tree, **files})
        m = re.search(r"^1\. Initial log-likelihood: (\S+)", run.stdout(), re.M)
        return run, (m.group(1) if m else None)

    f0 = fits["g0"]
    g4, g4_report = "+G4{" + f0["alpha"] + "}", "+G4{" + f0["report_alpha"] + "}"
    r6 = {"r6.paml": "\n".join(f0["block"]) + "\n"}
    flat = [float(v) for row in paml_rows(f0["rates"]) for v in row]
    scaled = ",".join(repr(v / flat[-1]) for v in flat[:-1])
    results = {
        "gA": show("gA", g0, "r6.paml" + g4, r6),
        "gA2": show("gA2", g0, "r6.paml" + pif + g4, r6),
        "gA3": show("gA3", g0, "r6.paml" + pif + g4_report, r6),
        "gB": show("gB", g0, "r10.paml" + pif + g4, {"r10.paml": paml(f0["rates"])}),
        "gC": show("gC", g0, "GTR20{" + scaled + "}" + pif + g4, {}),
    }
    rng = np.random.default_rng(SEED)
    draws = []
    for k in range(10):
        u = rng.uniform(-5e-10, 5e-10, len(f0["rates"]))
        values = [repr(float(v) * (1.0 + float(e))) for v, e in zip(f0["rates"], u)]
        draws.append(show(f"gD{k}", g0, "r10.paml" + pif + g4, {"r10.paml": paml(values)})[1])
    fs = fits["gS"]
    stale = show("gS_B", gs, "r10.paml" + pif + "+G4{" + fs["alpha"] + "}", {"r10.paml": paml(fs["rates"])})
    lnl = {k: v[1] for k, v in results.items()}
    obs.update(lnl=lnl, gD=draws, gS_B=stale[1], exit_codes={k: v[0].code for k, v in [*results.items(), ("gS_B", stale)]},
               state_freq_problems={k: freq_problems(v[0], pi) for k, v in results.items() if k != "gA"})
    if None in lnl.values() or None in draws or stale[1] is None or fs["final_opt_initial"] is None:
        return obs, False
    gb = float(lnl["gB"])
    spread = max(map(float, draws + [lnl["gB"]])) - min(map(float, draws + [lnl["gB"]]))
    obs.update(diff_from_gB={k: float(v) - gb for k, v in lnl.items()}, gD_spread=spread,
               gB_minus_g0_printed=gb - float(f0["final_printed"]),
               gS_B_minus_final_opt_initial=float(stale[1]) - float(fs["final_opt_initial"]),
               gS_final_printed_minus_gS_B=float(fs["final_printed"]) - float(stale[1]))
    g0_ok = (not f0["final_opt_ran"] and f0["rates_off_block"] == 0 and not f0["alpha_off_report"]
             and abs(float(lnl["gC"]) - gb) <= 1e-9 and abs(gb - float(f0["final_printed"])) <= 1e-4
             and all(abs(float(lnl[k]) - gb) > spread for k in ("gA", "gA2", "gA3"))
             and not any(obs["state_freq_problems"].values()))
    gs_ok = (fs["final_opt_ran"] and fs["rates_off_block"] > 0
             and abs(obs["gS_B_minus_final_opt_initial"]) <= 0.01 and obs["gS_final_printed_minus_gS_B"] >= 0.1)
    obs["verdict_parts"] = {"g0 (fixed tree)": g0_ok, "gS (tree search)": gs_ok}
    return obs, g0_ok and gs_ok


PROBES = {"a": ("Does a target vector named in an --mdef file reach the model as the literal +F{...} does?", probe_a),
          "b": ("Does a number such as 1e+00 inside +F{...} break parsing?", probe_b),
          "c": ("Do existing flags give a Level 1 fit with every nuisance frozen?", probe_c),
          "f": ("Does the root move under -te at Level 1 and Level 2, from rooted and unrooted input?", probe_f),
          "e": ("Does -S behave as -p apart from tree linkage?", probe_e),
          "g": ("How precisely does a fitted GTR20+F{pi*} carry into a new run without source changes?", probe_g)}


def provenance(a, binary, out):
    import numpy
    return {
        "script": SCRIPT, "command": sys.argv,
        "started_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
        "git_commit": capture(a.git, "rev-parse", "HEAD", cwd=REPO),
        "git_clean": capture(a.git, "status", "--porcelain", cwd=REPO) == "",
        "binary": {"path": str(binary), "sha256": sha256(binary)},
        "out": str(out), "python": sys.version, "numpy": numpy.__version__, "platform": platform.platform(),
        "inputs": {str(p.relative_to(REPO)): sha256(p) for p in (AA, TURTLE, NEX, BASELINE)},
        "settings": {"seed": SEED, "frozen_tol": FROZEN_TOL, "split_tol": SPLIT_TOL},
    }


def markdown(prov, results):
    lines = ["# S0 runtime probes", "",
             f"Generated by `{SCRIPT}` on {prov['started_utc']} at commit `{prov['git_commit']}` "
             f"(clean: {prov['git_clean']}); binary SHA-256 `{prov['binary']['sha256']}`.", ""]
    for key, res in results.items():
        lines += [f"## Probe ({key}): {res['question']}", "", f"Prediction: {res['prediction']}", "",
                  f"Verdict: **{res['verdict']}**", "", "```json",
                  json.dumps(res["observations"], indent=1, default=str)[:20000], "```", ""]
    return "\n".join(lines) + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--binary", required=True)
    ap.add_argument("--out", required=True, help="new directory outside the repository")
    ap.add_argument("--only", help="comma-separated probes (default: all; f reuses c's runs)")
    ap.add_argument("--baseline-raw", default="~/iqtree3-runs/baseline/rep1", help="raw outputs of baseline rep1, for probe e")
    ap.add_argument("--git", default="git", help="git executable (git.exe under WSL)")
    ap.add_argument("--allow-dirty", action="store_true")
    ap.add_argument("--timeout", type=float, default=7200, help="seconds per run")
    a = ap.parse_args()
    out = Path(a.out).expanduser().resolve()
    if out.exists() or out.is_relative_to(REPO):
        sys.exit(f"{out}: needs a new directory outside the repository")
    binary = Path(a.binary).expanduser().resolve()
    prov = provenance(a, binary, out)
    if not prov["git_clean"] and not a.allow_dirty:
        sys.exit(f"the repository is not clean ({prov['git_commit']}); commit first or pass --allow-dirty")
    keys = a.only.split(",") if a.only else list(PROBES)
    if set(keys) - set(PROBES):
        sys.exit(f"unknown probes: {sorted(set(keys) - set(PROBES))}")
    out.mkdir(parents=True)
    runner = Runner(binary, out, a.timeout)
    ctx = {"base": baseline_runs(), "baseline_raw": Path(a.baseline_raw).expanduser()}
    results = {}
    for key in keys:
        question, fn = PROBES[key]
        print(f"probe ({key}): {question}", flush=True)
        try:
            obs, ok = fn(runner, ctx)
            verdict = "as predicted" if ok else "differs"
        except Exception:  # a crash in a probe is recorded, and the other probes still run
            obs, verdict = {"exception": traceback.format_exc()}, "error"
        results[key] = {"question": question, "prediction": PREDICTIONS[key], "verdict": verdict, "observations": obs}
        print(f"probe ({key}): {verdict}", flush=True)
    prov["finished_utc"] = datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")
    write(out / "results.json", json.dumps({"provenance": prov, "runs": runner.runs, "probes": results},
                                           indent=1, default=str))
    write(out / "results.md", markdown(prov, results))
    print(f"wrote {out / 'results.md'}")
    sys.exit(0 if all(r["verdict"] == "as predicted" for r in results.values()) else 1)


if __name__ == "__main__":
    main()
