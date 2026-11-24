"""SIM-1's simulation pipeline (CODE_PLAN.md 3.8, "Pipeline"; decisions 024, 040, 043, 044). One replicate
per invocation, with paths, binaries and threads as arguments, so it runs unchanged on a cluster.

  pipeline.py simulate --scaffold DIR --truth NAME --genes G --replicate R --alisim BIN --git git.exe [--seed BASE] [--out DIR]
  pipeline.py fit --data DIR --binary BIN --threads T [--level2] [--relaxed] [--gtr20] [--rounds R] --git git.exe [--out DIR]
  pipeline.py checks --scaffold DIR --alisim BIN --binary BIN --truth NAME --git git.exe [--out DIR]
  pipeline.py metrics --fit FITDIR [--heldout SIMDIR --binary BIN] --git git.exe [--out DIR]
  pipeline.py pilot --scaffold DIR --alisim BIN --binary BIN [--threads 6] [--cells main1,ram,taxa26,genes100,variants,main] --git git.exe [--out DIR]

A scaffold directory holds tree.nwk (rooted, no stem, lengths at partition rate 1) and loci.tsv (a header,
then locus, length, +G4 shape and partition rate per line, tab-separated); sim/birddata.py scaffold writes
the bird one. Truths at pi* (sim/target.txt): NQ.bird, Q.bird, LG, or dial:LAMBDA (simulate.dial of NQ.bird).

simulate draws each gene's (length, shape, rate rho_k) with replacement from the scaffold's loci, seeded by
SeedSequence(base, replicate); simulates gene k by AliSim on T rho_k (-nt 1, --simulation-thresh 0, its own
directory and -pre, seed from SeedSequence(base, replicate, k) in [1, 2^31 - 1025]); and writes the
concatenation in the tree's leaf order, charsets, the -q file with true shapes and rates rho_k / m at 17
digits (m = sum rho_k n_k / sum n_k, which IQ-TREE's site-weighted renormalization leaves unchanged to
rounding: asserted), the fit tree T m, and a manifest with every file's SHA-256.
fit runs NQC Level 1 (-q, --model-joint NQC+F{pi*}, -te T m, -blfix, -keep-ident, default domain), with
options for -p Level 2, relaxed settings (-me 1e-6, 300-round caps), GTR20+F{pi*} (on the tree unrooted)
and a cap on the linked rounds.
checks asserts decision 021's two pipeline rows for one truth: the same seed gives byte-identical genes,
also with AliSim confined to 2 CPUs (taskset), and IQ-TREE's --show-lh at the true Q under -q equals the
oracle's linked objective built from the trees AliSim used (T rho_k at partition rate 1).
metrics reads an NQC fit's export: C, per-rate log error (RMS, maximum, twofold counts), net-flux error, and
with a held-out replicate the --show-lh totals (-q, true shapes and rates) of the fitted and the true Q.
pilot (exploratory; CODE_PLAN.md 3.8, "Pilot"): one fit at a time, caps 3 h (relaxed 4 h), fits at a cap
counted, results saved after each, a fit's record before its metrics: main1 (30 genes, replicate 1 of
Q.bird and NQ.bird), ram (IQ-TREE's -memcheck estimate at 300 and 1000 genes), main (replicates 2 to 5),
taxa26 (a seeded prune to 26 taxa keeping one on each side of the root), genes100, variants (on NQ.bird's
replicate 1: -p Level 2, relaxed, GTR20+F{pi*}), default1 (NQ.bird's replicate 1 at IQ-TREE's default
rounds, beside a capped main1); held-out lnL on 10 genes; each fit's own peak memory and CPU (os.wait4); the
SDs over main replicates with their upper 80% bounds. --rounds caps the linked rounds of every fit but
default1 and relaxed; --replicates and --variants trim the main cell and the variants (Peter, 2026-11-17:
--rounds 20 --replicates 3 --variants gtr20, to finish within a day).
"""
import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import time
import traceback
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
for _p in ("", "regression", "probes", "nqc", "differential", "explore", "profile"):
    sys.path.insert(0, str(HERE.parent / _p))
sys.path.insert(0, str(HERE))
from oracle import builtin, chart, flux, likelihood as lk, nqc, readers, simulate  # noqa: E402
from regress import REPO, read, sha256, write  # noqa: E402
from probes import PREFIX, Run, Runner, paml_rows, read_checkpoint, unroot, with_charpartition  # noqa: E402
from nqcchecks import check, gamma_shapes  # noqa: E402
from differential import write_matrix  # noqa: E402
from convergence import (LOOP_DEFAULT, ME_DEFAULT, ME_RELAXED, NPARAM_DEFAULT, ended_by, run_state,  # noqa: E402
                         trajectory)
from profiling import scale_tree  # noqa: E402
import twopart  # noqa: E402
import common  # noqa: E402

SCRIPT = Path(__file__).resolve().relative_to(REPO).as_posix()
SEED = 20261117
SEED_MAX = 2 ** 31 - 1025          # AliSim reads a signed 32-bit seed and adds up to twice the CPU count
ROUNDS = 300
TOL_RATES = 1e-12                   # the -q rates against IQ-TREE's renormalization of them
TOL_LNL = 1e-8                      # 021's wiring row
MIN_LNL = 5e4                       # so that --show-lh's 3-decimal total is within the line
TOTAL = re.compile(r"^Initial log-likelihood: (\S+)", re.M)   # -q's total, 3 decimals (nqcchecks.PARTITION_INITIAL)
TREE_RATE1 = "tree_rate1.nwk"       # the tree AliSim used at partition rate 1 (the oracle's, stage B's)
DATA_FILES = ("truth.qmat", "alignment.phy", "sets.nex", "parts.nex", "tree.nwk", TREE_RATE1)


# ---------------------------------------------------------------- inputs

def _lines(path):
    """A scaffold file's lines, '#' provenance lines and blank lines left out."""
    return [line for line in read(path).splitlines() if line.strip() and not line.startswith("#")]


def read_scaffold(path):
    path = Path(path)
    tree = "".join(_lines(path / "tree.nwk")).strip()
    check_rooted(tree, path / "tree.nwk")
    rows = [line.split("\t") for line in _lines(path / "loci.tsv")[1:]]
    if not rows or any(len(r) != 4 for r in rows):
        raise ValueError(f"{path / 'loci.tsv'}: expected lines of locus, length, shape, rate")
    loci = [{"locus": r[0], "length": int(r[1]), "shape": float(r[2]), "rate": float(r[3])} for r in rows]
    return tree, loci


def check_rooted(tree, where):
    root = readers.parse_newick(tree)
    if len(root.children) != 2:
        raise ValueError(f"{where}: not rooted (the root has {len(root.children)} children)")
    if root.length:
        raise ValueError(f"{where}: a root stem, which AliSim would simulate as a branch")


def truth(name, pi):
    """The truth's Q (unit mean rate) and its stationary distribution: NQ.bird, Q.bird, LG and dial:LAMBDA at
    pi; shift:MODEL:TV (SIM-2's cell B), R0 = dial:0 transported (decision 044) to pi shifted toward MODEL's
    composition to total variation TV, or to that composition itself with TV = own."""
    t = simulate.truths()[1]
    if name in t:
        return t[name], pi
    if name.startswith("dial:"):
        return simulate.dial(t["NQ.bird"], pi, float(name[5:])), pi
    if name.startswith("shift:"):
        _, model, tv = name.split(":")
        toward = np.asarray(builtin.exchangeabilities(model)[1], dtype=float)
        toward = toward / toward.sum()
        pi_d = toward if tv == "own" else simulate.shift_to_tv(pi, toward, float(tv))[1]
        return simulate.t3_transport(simulate.dial(t["NQ.bird"], pi, 0.0), pi, pi_d), pi_d
    raise ValueError(f"unknown truth {name}")


def at_target(name):
    """Whether the truth generating the data lies at pi* (rates are then compared with it)."""
    return not name.startswith(("shift:", "profiles"))


def profile_components(pi):
    """SIM-2's cell C (CODE_PLAN 3.8): R0 = dial:0 transported to each C20 profile recentred at pi with half
    the positivity limit, weighted as C20; and the reference, the mean-flux generator at pi (C 0)."""
    profiles, weights = builtin.frequency_profiles("C20")
    lam = simulate.recentre_limit(pi, profiles, weights) / 2
    r0 = simulate.dial(simulate.truths()[1]["NQ.bird"], pi, 0.0)
    comps = [simulate.Component(w, None, simulate.t3_transport(r0, pi, p), p)
             for w, p in zip(weights, simulate.recentre(pi, profiles, weights, lam))]
    return comps, flux.generator(simulate.mean_flux(comps), pi)


def gene_seed(base, replicate, gene, *more):
    """AliSim's seed for one gene (or one of its components, more = (c,)), in [1, SEED_MAX]."""
    v = int(np.random.SeedSequence([base, replicate, gene, *more]).generate_state(1, np.uint64)[0])
    return 1 + v % SEED_MAX


COUNTS_KEY = 1000   # SeedSequence pads short entropy with zeros, so a key of 0 would repeat the gene's own seed


def component_counts(base, replicate, gene, length, weights):
    """Cell C: a gene's sites per component, multinomial(length, weights), seeded by (base, replicate, gene,
    COUNTS_KEY); the components' AliSim seeds take the keys 1 to K."""
    rng = np.random.default_rng(np.random.SeedSequence([base, replicate, gene, COUNTS_KEY]))
    return [int(c) for c in rng.multinomial(length, np.asarray(weights, dtype=float))]


def draw(loci, genes, base, replicate):
    rng = np.random.default_rng(np.random.SeedSequence([base, replicate]))
    return [dict(loci[i], index=int(i)) for i in rng.integers(len(loci), size=genes)]


def tree_factor(genes):
    """m = sum rho_k n_k / sum n_k, and each gene's -q rate rho_k / m."""
    n = np.array([g["length"] for g in genes], dtype=float)
    rho = np.array([g["rate"] for g in genes], dtype=float)
    m = float(rho @ n / n.sum())
    return m, [float(r) / m for r in rho]


def leaf_order(tree):
    return [leaf.name for leaf in readers.parse_newick(tree).leaves()]


def concatenate(alignments, order):
    """Sequential PHYLIP of the genes' alignments, in the tree's leaf order, and their charsets."""
    length = sum(len(next(iter(a.values()))) for a in alignments)
    text = f"{len(order)} {length}\n" + "".join(f"{t} {''.join(a[t] for a in alignments)}\n" for t in order)
    sets, start = [], 1
    for k, a in enumerate(alignments, start=1):
        n = len(next(iter(a.values())))
        sets.append(f"charset g{k:04d} = {start}-{start + n - 1};")
        start += n
    return text, "#nexus\nbegin sets;\n" + "\n".join(sets) + "\nend;\n"


def q_file(sets_text, genes, rates, model="LG"):
    return with_charpartition(sets_text, [(f"{model}+G4{{{g['shape']:.17g}}}", f"g{k:04d}", f"{r:.17g}")
                                          for k, (g, r) in enumerate(zip(genes, rates), start=1)])


# ---------------------------------------------------------------- AliSim

def alisim(binary, rundir, name, matrix, shape, tree, length, seed, prefix=()):
    """One gene in the new directory rundir; returns (code, output path, printed seed, argv)."""
    rundir.mkdir()
    if any(rundir.glob("*.ckp.gz")):
        raise RuntimeError(f"{rundir}: a stale checkpoint")
    shutil.copyfile(matrix, rundir / matrix.name)
    write(rundir / "tree.nwk", tree + "\n")
    argv = [*prefix, str(binary), "--alisim", name, "-pre", name, "--seqtype", "AA", "-m",
            f"{matrix.name}+G4{{{shape:.17g}}}", "-t", "tree.nwk", "--length", str(length), "--seed", str(seed),
            "-af", "phy", "-nt", "1", "--simulation-thresh", "0"]
    with open(rundir / "stdout.txt", "w", encoding="utf-8") as fh:
        code = subprocess.run(argv, cwd=rundir, stdout=fh, stderr=subprocess.STDOUT).returncode
    printed = re.search(r"^Seed:\s+(\d+)", read(rundir / "stdout.txt"), re.M)
    return code, rundir / f"{name}.phy", int(printed.group(1)) if printed else None, argv


def simulate_genes(out, alisim_bin, tree, genes, Q, pi, base, replicate, prefix=()):
    """Every gene's alignment, with the rows of its checks."""
    write_matrix(out / "truth.qmat", Q, pi)
    order, rows, alns = leaf_order(tree), [], []
    for k, g in enumerate(genes, start=1):
        seed = gene_seed(base, replicate, k)
        code, path, printed, argv = alisim(alisim_bin, out / f"gene{k:04d}", f"g{k:04d}", out / "truth.qmat",
                                           g["shape"], scale_tree(tree, g["rate"]), g["length"], seed, prefix)
        aln = readers.read_phylip(path) if code == 0 and path.exists() else {}
        g.update({"seed": seed, "printed_seed": printed, "sha256": sha256(path) if aln else None,
                  "identical": [p for p in twopart.problems(aln) if "identical" in p] if aln else None})
        ok = code == 0 and bool(aln) and sorted(aln) == sorted(order) and printed == seed
        if not ok:
            rows.append(check(f"gene {k}: AliSim exits 0 with every taxon and its seed", False, exit_code=code,
                              printed_seed=printed, seed=seed, argv=argv))
        alns.append(aln)
    return alns, rows


def write_phylip(path, aln, order):
    write(path, f"{len(order)} {len(aln[order[0]])}\n" + "".join(f"{t} {aln[t]}\n" for t in order))


def simulate_profile_genes(out, alisim_bin, tree, genes, pi, base, replicate, prefix=()):
    """Cell C's genes: each gene's sites split among the C20 components (component_counts), each part simulated
    by AliSim from its component's matrix with seed gene_seed(base, replicate, gene, c), the parts joined in
    component order; truth.qmat holds the reference, the mean-flux generator at pi."""
    comps, reference = profile_components(pi)
    write_matrix(out / "truth.qmat", reference, pi)
    for c, comp in enumerate(comps, start=1):
        write_matrix(out / f"component{c:02d}.qmat", comp.Q, comp.pi)
    order, rows, alns = leaf_order(tree), [], []
    for k, g in enumerate(genes, start=1):
        gdir = out / f"gene{k:04d}"
        gdir.mkdir()
        counts = component_counts(base, replicate, k, g["length"], [c.weight for c in comps])
        parts, seeds, ok = [], [], True
        for c, n in enumerate(counts, start=1):
            if not n:
                continue
            seed = gene_seed(base, replicate, k, c)
            code, path, printed, argv = alisim(alisim_bin, gdir / f"c{c:02d}", f"g{k:04d}c{c:02d}",
                                               out / f"component{c:02d}.qmat", g["shape"], scale_tree(tree, g["rate"]),
                                               n, seed, prefix)
            part = readers.read_phylip(path) if code == 0 and path.exists() else {}
            if not (code == 0 and part and sorted(part) == sorted(order) and printed == seed):
                ok = False
                rows.append(check(f"gene {k}, component {c}: AliSim exits 0 with every taxon and its seed", False,
                                  exit_code=code, printed_seed=printed, seed=seed, argv=argv))
            parts.append(part)
            seeds.append(seed)
        aln = {t: "".join(p[t] for p in parts) for t in order} if ok else {}
        if aln:
            write_phylip(gdir / f"g{k:04d}.phy", aln, order)
        g.update({"counts": counts, "seeds": seeds, "sha256": sha256(gdir / f"g{k:04d}.phy") if aln else None,
                  "identical": [p for p in twopart.problems(aln) if "identical" in p] if aln else None})
        alns.append(aln)
    return alns, rows


# ---------------------------------------------------------------- steps

def simulate_into(out, alisim_bin, tree, loci, truth_name, genes_n, base, replicate, prov=None):
    """One replicate's data in the existing directory out: data.json (manifest) and the rows of its checks.
    profiles:C20 simulates cell C's mixture (simulate_profile_genes); every other name is a truth (truth())."""
    pi = common.target_values()
    genes = draw(loci, genes_n, base, replicate)
    if truth_name.startswith("profiles"):
        alns, rows = simulate_profile_genes(out, alisim_bin, tree, genes, pi, base, replicate)
    else:
        alns, rows = simulate_genes(out, alisim_bin, tree, genes, *truth(truth_name, pi), base, replicate)
    m, rates = tree_factor(genes)
    renorm = nqc.normalize_part_rates(rates, [g["length"] for g in genes])
    worst = float(np.max(np.abs(np.array(renorm) - rates) / np.array(rates)))
    rows.append(check("the -q rates are unchanged by IQ-TREE's site-weighted renormalization", worst <= TOL_RATES,
                      max_relative=worst, m=m))
    if all(alns):
        text, sets = concatenate(alns, leaf_order(tree))
        write(out / "alignment.phy", text)
        write(out / "sets.nex", sets)
        write(out / "parts.nex", q_file(sets, genes, rates))
        write(out / "tree.nwk", scale_tree(tree, m) + "\n")
    write(out / TREE_RATE1, tree + "\n")
    identical = sum(len(g["identical"] or []) for g in genes)
    rows.append(check("identical sequences within genes (counted)", None, pairs=identical))
    manifest = {"truth": truth_name, "replicate": replicate, "seed": base, "genes": genes, "m": m, "rates": rates,
                "target": common.target_text(), "taxa": len(leaf_order(tree)), "at_target": at_target(truth_name),
                "files": {f: sha256(out / f) for f in DATA_FILES if (out / f).exists()}}
    write(out / "data.json", json.dumps({"provenance": prov, **manifest}, indent=1, default=common._json))
    return manifest, rows


STATES = "ARNDCQEGHILKMFPSTWYV"     # IQ-TREE's amino-acid order
MASK_KEY = 99                       # cell D's uniform draws: SeedSequence([base, replicate, gene, MASK_KEY])


def read_curve(path):
    """CLOAK's binned removal curve (sim/cloak.py curve): the bin edges and each bin's removal probability."""
    rows = [line.split("\t") for line in _lines(path)]
    head, body = rows[0], rows[1:]
    lo, hi, p = (head.index(c) for c in ("bin_lo", "bin_hi", "p"))
    edges = [float(body[0][lo])] + [float(r[hi]) for r in body]
    probs = np.array([float(r[p]) for r in body])
    if not np.all(np.isfinite(probs)):
        raise ValueError(f"{path}: a bin without a finite removal probability (no residues in it)")
    return np.array(edges), probs


def derived_dataset(src, dst, alns, order, extra):
    """A dataset made from src's genes (masked or pruned): the genes' alignments, their concatenation, src's
    partition files, truth and rate-1 tree unless replaced, and a manifest recording extra."""
    src, dst = Path(src), Path(dst)
    man = json.loads(read(src / "data.json"))
    for k, aln in enumerate(alns, start=1):
        (dst / f"gene{k:04d}").mkdir()
        write_phylip(dst / f"gene{k:04d}" / f"g{k:04d}.phy", aln, order)
    text, sets = concatenate(alns, order)
    write(dst / "alignment.phy", text)
    for f in DATA_FILES:
        if not (dst / f).exists() and f != "alignment.phy" and (src / f).exists():
            shutil.copyfile(src / f, dst / f)
    man.update({"source": str(src), "taxa": len(order), **extra,
                "files": {f: sha256(dst / f) for f in DATA_FILES if (dst / f).exists()}})
    write(dst / "data.json", json.dumps(man, indent=1, default=common._json))
    return man, sets == read(src / "sets.nex")


def mask_dataset(src, dst, s, edges, p, curve_sha256=None):
    """Cell D (decision 045): src's genes masked by simulate.rarity_mask at strength s, with one uniform draw per
    residue from SeedSequence([base, replicate, gene, MASK_KEY]) in the tree's leaf order, so strengths nest.
    Records the masked fraction and the total variation between the masked and the unmasked composition."""
    src = Path(src)
    man = json.loads(read(src / "data.json"))
    order = leaf_order(read(src / TREE_RATE1).strip())
    alns, kept, before, after = [], 0, np.zeros(20), np.zeros(20)
    for k, aln in enumerate(gene_alignments_of(src, len(man["genes"])), start=1):
        states = np.array([[STATES.index(c) for c in aln[t]] for t in order])
        u = np.random.default_rng(np.random.SeedSequence([man["seed"], man["replicate"], k, MASK_KEY])).uniform(
            size=states.shape)
        masked = simulate.rarity_mask(states, s, edges, p, u)
        before += np.bincount(states.ravel(), minlength=20)
        after += np.bincount(states[~masked], minlength=20)
        kept += int((~masked).sum())
        alns.append({t: "".join("-" if m else c for c, m in zip(aln[t], row)) for t, row in zip(order, masked)})
    total = int(before.sum())
    return derived_dataset(src, dst, alns, order, {"masking": {
        "s": s, "curve_sha256": curve_sha256, "masked_fraction": 1 - kept / total,
        "composition_tv": simulate.total_variation(after / after.sum(), before / before.sum())}})


def prune_dataset(src, dst, keep):
    """Cell E: src's genes restricted to the taxa in keep, on src's rate-1 tree pruned (prune), the fit tree that
    pruned tree times m; returns the manifest and whether the charsets are unchanged."""
    src, dst = Path(src), Path(dst)
    man = json.loads(read(src / "data.json"))
    pruned = prune(read(src / TREE_RATE1).strip(), set(keep))
    write(dst / TREE_RATE1, pruned + "\n")
    write(dst / "tree.nwk", scale_tree(pruned, man["m"]) + "\n")
    order = leaf_order(pruned)
    alns = [{t: a[t] for t in order} for a in gene_alignments_of(src, len(man["genes"]))]
    return derived_dataset(src, dst, alns, order, {"pruned": sorted(keep)})


def gene_alignments_of(data, genes):
    return [readers.read_phylip(Path(data) / f"gene{k:04d}" / f"g{k:04d}.phy") for k in range(1, genes + 1)]


def cmd_simulate(a):
    prov = common.provenance(a, SCRIPT, f"simulate-{a.truth.replace(':', '')}-g{a.genes}-r{a.replicate}",
                             "exploratory", truth=a.truth, genes=a.genes, replicate=a.replicate, seed=a.seed)
    out = common.new_out(a.out, prov)
    alisim_bin = Path(a.alisim).expanduser().resolve()
    prov["alisim"] = {"path": str(alisim_bin), "sha256": sha256(alisim_bin)}
    tree, loci = read_scaffold(a.scaffold)
    prov["scaffold"] = scaffold_record(a.scaffold)
    manifest, rows = simulate_into(out, alisim_bin, tree, loci, a.truth, a.genes, a.seed, a.replicate, prov)
    common.finish(out, prov, f"SIM-1 pipeline simulate ({a.truth}, {a.genes} genes, replicate {a.replicate})", rows,
                  manifest)


def scaffold_record(path):
    return {"path": str(Path(path).resolve()), "sha256": {f: sha256(Path(path) / f) for f in ("tree.nwk", "loci.tsv")}}


def fit_args(data, model, level2=False, relaxed=False, rounds=None, tree=None):
    """IQ-TREE's arguments for one fit; rounds caps the linked rounds (-nparam rounds + 1:
    PartitionModelPlen::optimizeParameters runs rounds 1 to num_param_iterations - 1, default 100); tree
    replaces the data's rooted tree.nwk."""
    args = _fit_args(data, model, level2, relaxed, tree or data / "tree.nwk")
    return args + ["-nparam", str(rounds + 1)] if rounds else args


def _fit_args(data, model, level2, relaxed, tree):
    if level2:
        args = ["-s", data / "alignment.phy", "-p", "parts.nex", "--model-joint", model, "-te", tree, "-keep-ident"]
    else:
        args = ["-s", data / "alignment.phy", "-q", "parts.nex", "--model-joint", model, "-te", tree, "-blfix",
                "-keep-ident"]
    if relaxed:
        args += ["-me", f"{ME_RELAXED:g}", "--loop-model", str(ROUNDS), "-nparam", str(ROUNDS + 2)]
    return args


def measured_run(binary, out, name, args, files, threads, timeout):
    """probes.Runner.run's command in the new directory out/name (-seed 1 -T threads --prefix PREFIX), with
    this child's own resource use from os.wait4 (RUSAGE_CHILDREN would hold the largest child so far); at
    the timeout IQ-TREE is killed and the code is "timeout". Returns the probes.Run and the usage."""
    rundir = out / name
    rundir.mkdir()
    for fname, content in (files or {}).items():
        write(rundir / fname, content)
    argv = [str(binary), *[str(x) for x in args], "-seed", "1", "-T", str(threads), "--prefix", PREFIX]
    t0 = time.monotonic()
    with open(rundir / "stdout.txt", "w", encoding="utf-8") as fh:
        proc = subprocess.Popen(argv, cwd=rundir, stdout=fh, stderr=subprocess.STDOUT)
        code = None
        while code is None:
            pid, status, usage = os.wait4(proc.pid, os.WNOHANG)
            if pid:
                code = os.waitstatus_to_exitcode(status)
            elif time.monotonic() - t0 > timeout:
                proc.kill()
                _, status, usage = os.wait4(proc.pid, 0)
                code = "timeout"
            else:
                time.sleep(0.5)
    wall = round(time.monotonic() - t0, 1)
    print(f"  {name}: exit {code}, {wall} s", flush=True)
    run = Run(rundir, name, [str(x) for x in args], code, wall, threads)
    return run, {"wall_s": wall, "cpu_s": round(usage.ru_utime + usage.ru_stime, 1),
                 "peak_rss_mb": round(usage.ru_maxrss / 1024, 1)}


def unrooted(tree_text):
    """The rooted tree with its two root edges joined (probes.unroot), at 17 digits: on a rooted -te tree
    IQ-TREE's GTR20 segfaults under --model-joint (maintainer question 17), and a reversible likelihood
    does not depend on the root."""
    return newick17(unroot(readers.parse_newick(tree_text)))


def round_ending(run, relaxed=False, rounds=None):
    """The rule that ended a -q fit's linked rounds, from its printed rounds and the cap its -nparam set
    (convergence.run_state tests IQ-TREE's default cap)."""
    vals, done = trajectory(run, "-q")
    nparam = rounds + 1 if rounds else (ROUNDS + 2 if relaxed else NPARAM_DEFAULT)
    return ended_by("-q", done, list(np.diff(vals)), ME_RELAXED if relaxed else ME_DEFAULT, LOOP_DEFAULT, nparam)


def fit_step(level2, relaxed, gtr20):
    kind = {"counted": "gtr20f"}.get(gtr20, "gtr20") if gtr20 else "nqc"
    return f"fit-{kind}-{'p2' if level2 else 'q1'}{'-relaxed' if relaxed else ''}"


def fit_into(out, data, binary, threads, timeout, level2=False, relaxed=False, gtr20=False, rounds=None,
             init_model=None, extra=()):
    """One fit of a replicate's data in the existing directory out (its run in out/fit); the fit's record
    (a timeout is counted as a fit at its cap). gtr20: True or "pi" fits GTR20+F{pi*}, "counted" GTR20+F (SIM-2's
    cell D), both on the tree unrooted; init_model = (file name, text) starts NQC from that matrix (--init-model:
    a PAML file is rebuilt at pi*, a full -m FILE matrix keeps its jump chain)."""
    data = Path(data)
    manifest = json.loads(read(data / "data.json"))
    bad = [f for f, h in manifest["files"].items() if sha256(data / f) != h]
    text = manifest["target"]
    model = ("GTR20+F" if gtr20 == "counted" else f"GTR20+F{{{text}}}") if gtr20 else f"NQC+F{{{text}}}"
    # Level 2 (-p) estimates the shapes and partition rates, so its charpartition fixes neither
    parts = (with_charpartition(read(data / "sets.nex"), [("LG+G4", f"g{k:04d}", None)
                                                          for k in range(1, len(manifest["genes"]) + 1)])
             if level2 else read(data / "parts.nex"))
    files = {"parts.nex": parts}
    if gtr20:
        files[UNROOTED] = unrooted(read(data / "tree.nwk")) + "\n"
    args = fit_args(data, model, level2, relaxed, rounds, UNROOTED if gtr20 else None)
    if init_model is not None:
        files[init_model[0]] = init_model[1]
        args += ["--init-model", init_model[0]]
    run, usage = measured_run(binary, out, "fit", args + list(extra), files, threads, timeout)
    state = run_state(run, "-q", not gtr20, "relaxed" if relaxed else "default", ROUNDS) if run.code == 0 else {}
    if run.code == 0:
        state["ended_by"] = round_ending(run, relaxed, rounds)
    ram = RAM_NOTE.search(run.stdout())
    info = readers.read_info(read(run.file(".NQC.info"))) if run.file(".NQC.info").exists() else {}
    return {"step": fit_step(level2, relaxed, gtr20), "rounds_cap": rounds, "data_match": not bad, "mismatched": bad,
            "exit_code": run.code,
            "at_cap": run.code == "timeout", **usage, "ram_note_mb": int(ram.group(1)) if ram else None,
            "threads": threads, "columns": sum(g["length"] for g in manifest["genes"]), "taxa": manifest.get("taxa"),
            "evaluations_key": int(info["evaluations"]) if "evaluations" in info else None,
            "diagnostic_evaluations": int(info["diagnostic_evaluations"]) if "diagnostic_evaluations" in info else None,
            "stop": info.get("stop"), "start_score": float(info["start_score"]) if "start_score" in info else None,
            "init_model": init_model[0] if init_model else None,
            **{k: state.get(k) for k in ("final_lnl", "initial_lnl", "rounds", "ended_by", "evaluations", "fits",
                                         "switches", "falls", "assertion_or_abort", "bounds", "outside_box")}}


def gtr20_state(ck):
    """A linked GTR20 fit's state from its checkpoint dict (probes.read_checkpoint): the exchangeabilities R,
    symmetric, from the 190 rates IQ-TREE stores row by row over the upper triangle (also kept as their 10-digit
    strings, for a PAML carry), the linked model's frequencies and each partition's. The linked model is named
    GTR20 for GTR20+F{pi*} and GTR20+F for counted frequencies, which --model-joint pools (SIM-2's smoke run)."""
    def linked(k, key):
        return k.endswith(f"!ModelProtein!{key}") and k.split("!")[1] in ("GTR20", "GTR20+F")

    rates = [v for k, v in ck.items() if linked(k, "rates")]
    freqs = [v for k, v in ck.items() if linked(k, "state_freq")]
    if len(rates) != 1 or len(freqs) != 1:
        raise ValueError(f"expected one linked GTR20 rates and frequencies, found {len(rates)} and {len(freqs)}")
    text = [v.strip() for v in rates[0].split(",")]
    if len(text) != 190:
        raise ValueError(f"{len(text)} GTR20 rates, expected 190")
    R = np.zeros((20, 20))
    R[np.triu_indices(20, 1)] = [float(v) for v in text]
    parts = {k.split("!")[1]: np.array([float(x) for x in v.split(",")]) for k, v in ck.items()
             if k.endswith("!ModelProtein!state_freq") and not linked(k, "state_freq")}
    return {"R": R + R.T, "rates_text": text, "freq": np.array([float(x) for x in freqs[0].split(",")]),
            "part_freqs": parts}


def gtr20_paml(state, target_text):
    """decision 025's carry as nqcchecks.case_nesting writes it: the 19 lower-triangle rows of the checkpoint's
    10-digit rates in PAML order, then pi*'s 20 values; --init-model rebuilds it at the target."""
    return "\n".join(" ".join(row) for row in paml_rows(state["rates_text"])) + "\n" + target_text.replace(",", " ") \
        + "\n"


MAX_RATE = 100.0          # ModelMarkov's upper bound on a free rate (model/modelmarkov.h, MAX_RATE)
FLOOR_HELD = 1.0001e-4    # a free rate at ModelMarkov's lower bound, MIN_RATE = 1e-4


def gtr20_seed(p, target_text):
    """Decision 046's --init-model file for GTR20, written as gtr20_paml writes a carry: every exchangeability 10
    and Tyr-Val, IQ-TREE's last rate (the one GTR20 keeps fixed), p; then pi*. Its largest entry is 10, so
    rescaleRates scales by exactly 1 and Tyr-Val stays exactly p."""
    return gtr20_paml({"rates_text": ["10"] * 189 + [f"{p:.10g}"]}, target_text)


def rate_bounds(state):
    """A GTR20 fit's rates against IQ-TREE's box (decision 046): the fixed last rate's checkpoint text; the
    largest free rate and its headroom (MAX_RATE over it); the free rates above isUnstableParameters' 0.99 x
    MAX_RATE; those at the floor."""
    free = np.array([float(v) for v in state["rates_text"][:-1]])
    largest = float(free.max())
    return {"fixed_text": state["rates_text"][-1], "largest_free": largest, "headroom": MAX_RATE / largest,
            "above_upper": int(np.sum(free > 0.99 * MAX_RATE)), "at_floor": int(np.sum(free <= FLOOR_HELD))}


def level2_state(fitdir, taxa):
    """A -p fit's state for the oracle (SIM-2's cell G): the tree from the checkpoint (PhyloTree!newick, 10
    significant digits where .treefile has 10 decimals; taxon IDs index taxa, the alignment's rows), its paths
    checked against .treefile's; each partition's rate (PartitionModelPlen!part_rates) and gamma shape, in
    partition order g0001, g0002, ..."""
    fitdir = Path(fitdir)
    ck = read_checkpoint(fitdir / f"{PREFIX}.ckp.gz")
    tree = readers.parse_newick(ck["PhyloTree!newick"])
    for leaf in tree.leaves():
        leaf.name = taxa[int(leaf.name)]
    printed = readers.parse_newick(read(fitdir / f"{PREFIX}.treefile").strip())
    a, b = leaf_depths(tree), leaf_depths(printed)
    if sorted(a) != sorted(b) or any(abs(a[t] - b[t]) > 1e-9 * max(1.0, b[t]) + 1e-9 for t in a):
        raise ValueError(f"{fitdir}: the checkpoint's tree does not match .treefile")
    rates = [float(v) for v in ck["PartitionModelPlen!part_rates"].split(",")]
    names = [f"g{k:04d}" for k in range(1, len(rates) + 1)]
    shapes = gamma_shapes(ck, names)
    if None in shapes:
        raise ValueError(f"{fitdir}: a partition's gamma shape is missing")
    return {"tree": newick17(tree), "rates": rates, "shapes": [float(s) for s in shapes]}


def leaf_depths(tree):
    """Each leaf's path length from the root."""
    out = {}

    def rec(node, depth):
        if not node.children:
            out[node.name] = depth
        for c in node.children:
            rec(c, depth + (c.length or 0.0))

    rec(tree, 0.0)
    return out


def cmd_fit(a):
    data = Path(a.data).expanduser().resolve()
    step = fit_step(a.level2, a.relaxed, a.gtr20)
    prov = common.provenance(a, SCRIPT, step, "exploratory", data=str(data), threads=a.threads,
                             level2=a.level2, relaxed=a.relaxed, gtr20=a.gtr20, rounds_cap=a.rounds)
    out = common.new_out(a.out, prov)
    binary = Path(a.binary).expanduser().resolve()
    prov["binary"] = {"path": str(binary), "sha256": sha256(binary)}
    rec = fit_into(out, data, binary, a.threads, a.timeout, a.level2, a.relaxed, a.gtr20, a.rounds)
    rows = [check("data match their manifest", rec["data_match"], mismatched=rec["mismatched"]),
            check("the fit exits 0", rec["exit_code"] == 0, **{k: v for k, v in rec.items() if k != "mismatched"})]
    common.finish(out, prov, f"SIM-1 pipeline {step}", rows)


def heldout_lnl(binary, rundir_parent, name, heldout, qmat_text, timeout):
    """IQ-TREE's --show-lh total under -q on a held-out replicate at a fixed matrix (true shapes and rates)."""
    heldout = Path(heldout)
    man = json.loads(read(heldout / "data.json"))
    files = {"parts.nex": q_file(read(heldout / "sets.nex"), man["genes"], man["rates"], model="model.qmat"),
             "model.qmat": qmat_text}
    run = Runner(binary, rundir_parent, timeout).run(name, ["-s", heldout / "alignment.phy", "-q", "parts.nex", "-te",
                                                            heldout / "tree.nwk", "--show-lh", "-keep-ident"], files)
    printed = TOTAL.findall(run.stdout())
    return float(printed[-1]) if run.code == 0 and printed else None


def metrics_of(fitdir, data, heldout=None, binary=None, timeout=3600, rundir=None):
    """matrix_metrics of an NQC fit's export in fitdir/fit, its held-out runs in rundir (default fitdir)."""
    qmat = Path(fitdir) / "fit" / f"{PREFIX}.NQC.qmat"
    if not qmat.exists():
        return {"fitted": False}
    return matrix_metrics(readers.read_qmat(read(qmat))[0], data, heldout, binary, timeout,
                          Path(fitdir) if rundir is None else Path(rundir), qmat_text=read(qmat))


def matrix_metrics(Qh, data, heldout=None, binary=None, timeout=3600, rundir=None, name="fit", qmat_text=None,
                   pi_fit=None):
    """A fitted matrix against its cell's reference (data's truth.qmat, at its own pi): C and its error, and
    with a truth at pi* (data.json's at_target, true when absent) rate_metrics; with a held-out replicate, the
    --show-lh totals (-q, true shapes and rates, 3 decimals) of the fitted and the reference matrix, its runs in
    rundir as heldout-NAME and heldout-truth. pi_fit is the fitted matrix's stationary distribution, pi* by
    default; qmat_text, the file scored, write_matrix's of (Qh, pi_fit) by default."""
    data = Path(data)
    manifest = json.loads(read(data / "data.json"))
    pi = np.array([float(v) for v in manifest["target"].split(",")])
    pi_fit = pi if pi_fit is None else np.asarray(pi_fit, dtype=float)
    Qt, pi_t = readers.read_qmat(read(data / "truth.qmat"))
    out = {"fitted": True, "C_fit": simulate.circulation(Qh, pi_fit), "C_truth": simulate.circulation(Qt, pi_t)}
    out["C_error"] = out["C_fit"] - out["C_truth"]
    if manifest.get("at_target", True):
        out.update(rate_metrics(Qh, Qt, pi, pi_fit))
    if heldout is not None:
        if qmat_text is None:
            write_matrix(Path(rundir) / f"{name}.qmat", Qh, pi_fit)
            qmat_text = read(Path(rundir) / f"{name}.qmat")
        lf = heldout_lnl(binary, rundir, f"heldout-{name}", heldout, qmat_text, timeout)
        lt = heldout_lnl(binary, rundir, "heldout-truth", heldout, read(data / "truth.qmat"), timeout) \
            if not (Path(rundir) / "heldout-truth").exists() else _heldout_total(Path(rundir) / "heldout-truth")
        out.update({"heldout_lnl_fit": lf, "heldout_lnl_truth": lt,
                    "heldout_difference": None if lf is None or lt is None else lf - lt})
    return out


def _heldout_total(rundir):
    printed = TOTAL.findall(read(rundir / "stdout.txt"))
    return float(printed[-1]) if printed else None


def rate_metrics(Qh, Qt, pi, pi_fit=None):
    """Per-rate log error of Qh against the truth Qt at pi: RMS, maximum, counts off twofold and tenfold, RMS by
    quintile of true flux and flux-weighted RMS; net-flux error; and risk 28's split, the pair-sum flux's log
    error (RMS) and the direction share's error (mean absolute) by quintile of true pair-sum flux (archive
    SIM2_plan, A4). Qh's fluxes are taken at its own stationary distribution pi_fit (pi by default), so a
    reversible fit at other frequencies has none net. Quintiles run from the highest true flux to the lowest."""
    pi = np.asarray(pi, dtype=float)
    pi_fit = pi if pi_fit is None else np.asarray(pi_fit, dtype=float)
    off = ~np.eye(len(pi), dtype=bool)
    F, Fh = pi[:, None] * Qt, pi_fit[:, None] * Qh
    e, f = np.log(Qh[off] / Qt[off]), F[off]
    iu = np.triu_indices(len(pi), 1)
    S, Sh = (F + F.T)[iu], (Fh + Fh.T)[iu]
    share_error = np.abs((Fh / (Fh + Fh.T))[iu] - (F / (F + F.T))[iu])
    pair_error = np.log(Sh / S)
    net = 2 * (simulate.flux_parts(Qh, pi_fit)[1] - simulate.flux_parts(Qt, pi)[1])

    def rms(v):
        return float(np.sqrt(np.mean(v ** 2)))

    by_flux = np.array_split(np.argsort(-f, kind="stable"), 5)
    by_pair = np.array_split(np.argsort(-S, kind="stable"), 5)
    return {"log_error_rms": rms(e), "log_error_max": float(np.max(np.abs(e))),
            "twofold": int(np.sum(np.abs(e) > np.log(2))), "tenfold": int(np.sum(np.abs(e) > np.log(10))),
            "log_error_rms_by_flux_quintile": [rms(e[k]) for k in by_flux],
            "log_error_rms_flux_weighted": float(np.sqrt(f @ e ** 2 / f.sum())),
            "pair_sum_log_error_by_quintile": [rms(pair_error[k]) for k in by_pair],
            "direction_share_error_by_quintile": [float(np.mean(share_error[k])) for k in by_pair],
            "net_flux_error_max": float(np.max(np.abs(net))), "net_flux_error_rms": rms(net)}


def cmd_metrics(a):
    fitdir = Path(a.fit).expanduser().resolve()
    data = Path(json.loads(read(fitdir / "results.json"))["provenance"]["data"])
    if a.heldout and not a.binary:
        sys.exit("--heldout needs --binary")
    prov = common.provenance(a, SCRIPT, "metrics", "exploratory", fit=str(fitdir), data=str(data), heldout=a.heldout)
    out = common.new_out(a.out, prov)
    binary = Path(a.binary).expanduser().resolve() if a.binary else None
    m = metrics_of(fitdir, data, Path(a.heldout).expanduser().resolve() if a.heldout else None, binary, a.timeout,
                   rundir=out)
    common.finish(out, prov, "SIM-1 pipeline metrics", [check("metrics", None, **m)])


SEEDS_USED = re.compile(r"^- Random seeds used:((?: \d+)+)", re.M)


def seed_list(binary, rundir, matrix, tree, prefix=()):
    """AliSim's -vvv list of per-thread seeds (2 x the CPU cores it counts), from a short run."""
    rundir.mkdir()
    shutil.copyfile(matrix, rundir / matrix.name)
    write(rundir / "tree.nwk", tree + "\n")
    argv = [*prefix, str(binary), "--alisim", "probe", "-pre", "probe", "--seqtype", "AA", "-m", matrix.name, "-t",
            "tree.nwk", "--length", "10", "--seed", "1", "-af", "phy", "-nt", "1", "--simulation-thresh", "0", "-vvv"]
    with open(rundir / "stdout.txt", "w", encoding="utf-8") as fh:
        subprocess.run(argv, cwd=rundir, stdout=fh, stderr=subprocess.STDOUT)
    m = SEEDS_USED.search(read(rundir / "stdout.txt"))
    return [int(v) for v in m.group(1).split()] if m else None


def wiring_objective(alns, tree, genes, pi):
    """The oracle's linked objective from the trees AliSim used: gene k on T rho_k at partition rate 1,
    rooted at the -m FILE convention's pi, with its true shape."""
    ref = nqc.lg_seed(pi)[1]
    return nqc.LinkedObjective([nqc.Objective(a, readers.parse_newick(scale_tree(tree, g["rate"])), pi, ref,
                                              shape=g["shape"]) for a, g in zip(alns, genes)])


def two_cpus(affinity=None):
    """taskset's argument for AliSim's two-CPU check: the two lowest CPUs this process may use (a cluster job's
    CPUs need not include 0 and 1)."""
    cpus = sorted(os.sched_getaffinity(0) if affinity is None else affinity)
    if len(cpus) < 2:
        raise ValueError("the two-CPU check needs at least two CPUs")
    return f"{cpus[0]},{cpus[1]}"


def cmd_checks(a):
    prov = common.provenance(a, SCRIPT, f"checks-{a.truth.replace(':', '')}", "confirmatory", truth=a.truth,
                             genes=a.genes, seed=a.seed)
    out = common.new_out(a.out, prov)
    alisim_bin, binary = (Path(p).expanduser().resolve() for p in (a.alisim, a.binary))
    prov["alisim"] = {"path": str(alisim_bin), "sha256": sha256(alisim_bin)}
    prov["binary"] = {"path": str(binary), "sha256": sha256(binary)}
    tree, loci = read_scaffold(a.scaffold)
    prov["scaffold"] = scaffold_record(a.scaffold)
    rows = checks_rows(out, alisim_bin, binary, tree, loci, a.truth, a.genes, a.seed, a.timeout, prov)
    common.finish(out, prov, f"SIM-1 pipeline checks ({a.truth})", rows)


def checks_rows(out, alisim_bin, binary, tree, loci, truth_name, genes_n, seed, timeout, prov=None):
    """Decision 021's pipeline rows and the wiring row for one truth, run in the existing directory out."""
    pi = common.target_values()
    Q, pi_t = truth(truth_name, pi)
    taskset = ("taskset", "-c", two_cpus())
    first = out / "first"
    first.mkdir()
    manifest, rows = simulate_into(first, alisim_bin, tree, loci, truth_name, genes_n, seed, 0, prov)
    hashes = {"first": [g["sha256"] for g in manifest["genes"]]}
    alns = {}
    for name, prefix in (("again", ()), ("two-cpus", taskset)):
        d = out / name
        d.mkdir()
        genes = draw(loci, genes_n, seed, 0)
        alns[name], found = simulate_genes(d, alisim_bin, tree, genes, Q, pi_t, seed, 0, prefix)
        rows += found
        hashes[name] = [g["sha256"] for g in genes]
    made = None not in hashes["first"]
    rows.append(check("same seed: every gene byte-identical", made and hashes["first"] == hashes["again"],
                      first=hashes["first"], again=hashes["again"]))
    seeds = {name: seed_list(alisim_bin, out / f"seeds-{name}", first / "truth.qmat", tree, prefix)
             for name, prefix in (("all", ()), ("two-cpus", taskset))}
    rows.append(check("taskset confines AliSim to 2 CPUs (its -vvv seed list has 2 x 2 entries, more without it)",
                      seeds["two-cpus"] is not None and len(seeds["two-cpus"]) == 4 and seeds["all"] is not None
                      and len(seeds["all"]) > 4, **seeds))
    rows.append(check("AliSim on 2 CPUs against all: every gene byte-identical",
                      made and hashes["first"] == hashes["two-cpus"], two_cpus=hashes["two-cpus"]))
    wiring = f"wiring: IQ-TREE --show-lh under -q against the oracle, <= 1e-8 |lnL|, |lnL| >= {MIN_LNL:g}"
    if not (first / "parts.nex").exists():
        rows.append(check(wiring, False, reason="no data: a gene failed"))
    else:
        # the files a fit reads (simulate_into's), with the true matrix in place of LG
        files = {"parts.nex": read(first / "parts.nex").replace("LG+G4{", "truth.qmat+G4{"),
                 "truth.qmat": read(first / "truth.qmat")}
        run = Runner(binary, out, timeout).run("wiring", ["-s", first / "alignment.phy", "-q", "parts.nex", "-te",
                                                            first / "tree.nwk", "--show-lh", "-keep-ident"], files)
        printed = TOTAL.findall(run.stdout())
        Qf, pif = lk.file_model(readers.read_qmat(files["truth.qmat"])[0])
        gene_alns = [readers.read_phylip(first / f"gene{k:04d}" / f"g{k:04d}.phy")
                     for k in range(1, len(manifest["genes"]) + 1)]
        oracle = wiring_objective(gene_alns, tree, manifest["genes"], pif).log_likelihood(Qf)
        iq = float(printed[-1]) if printed else None
        rows.append(check(wiring, run.code == 0 and iq is not None and abs(oracle) >= MIN_LNL
                          and abs(iq - oracle) <= TOL_LNL * abs(oracle),
                          iqtree=iq, oracle=oracle, difference=None if iq is None else iq - oracle, exit_code=run.code,
                          printed=printed))
    return rows


# ---------------------------------------------------------------- pilot (exploratory)

PILOT_TRUTHS = ("Q.bird", "NQ.bird")
PILOT_REPLICATES = 5
PILOT_GENES = 30
HELDOUT_GENES = 10
HELDOUT_REPLICATE = 1000        # held-out data apart from every fitted replicate
PRUNED_TAXA = 26
CAP, CAP_RELAXED = 3 * 3600.0, 4 * 3600.0
RAM_GENES = (300, 1000)
CELLS = ("main1", "ram", "main", "taxa26", "genes100", "variants", "default1")
UNROOTED = "tree-unrooted.nwk"     # GTR20's -te tree, written into its run directory
RAM_NOTE = re.compile(r"^NOTE: (\d+) MB RAM", re.M)
RAM_ERROR = re.compile(r"^ERROR: Your RAM is below minimum requirement of (\S+) GB", re.M)
RAM_ALLOC = re.compile(r"^ERROR: Not enough memory, allocation of (\d+) bytes failed", re.M)


def newick17(tree):
    def rec(n):
        s = "(" + ",".join(rec(c) for c in n.children) + ")" if n.children else n.name
        return s + (f":{n.length:.17g}" if n is not tree and n.length is not None else "")
    return rec(tree) + ";"


def prune(tree_text, keep):
    """The rooted tree restricted to the leaves in keep, each degree-2 node merged into its child with the
    lengths summed, so every kept path keeps its length; no stem."""
    def rec(nd):
        if not nd.children:
            return readers.Node(nd.name, nd.length) if nd.name in keep else None
        kids = [k for k in (rec(c) for c in nd.children) if k is not None]
        if not kids:
            return None
        if len(kids) == 1:
            kids[0].length = (kids[0].length or 0.0) + (nd.length or 0.0)
            return kids[0]
        return readers.Node(None, nd.length, kids)
    new = rec(readers.parse_newick(tree_text))
    new.length = None
    return newick17(new)


def pruned_taxa(tree_text, k, seed):
    """k leaves: one drawn from each side of the root, the rest drawn from all others, seeded."""
    root = readers.parse_newick(tree_text)
    rng = np.random.default_rng(seed)
    sides = [[leaf.name for leaf in child.leaves()] for child in root.children]
    keep = {str(rng.choice(side)) for side in sides}
    rest = sorted(set(leaf_order(tree_text)) - keep)
    keep |= set(rng.choice(rest, k - len(keep), replace=False).tolist())
    return keep


def ram_probe(binary, out, name, data, threads, timeout=3600):
    """The Level 1 fit's command with -memcheck, which stops after IQ-TREE's memory estimate: its RAM note
    in MB, the GB of its "RAM below minimum" error when the estimate exceeds the machine, or the bytes of an
    allocation that failed before the estimate (the pilot's 1000 genes)."""
    manifest = json.loads(read(Path(data) / "data.json"))
    args = [*fit_args(Path(data), f"NQC+F{{{manifest['target']}}}"), "-memcheck"]
    run, usage = measured_run(binary, out, name, args, {"parts.nex": read(Path(data) / "parts.nex")}, threads, timeout)
    text = run.stdout()
    note, error, alloc = RAM_NOTE.search(text), RAM_ERROR.search(text), RAM_ALLOC.search(text)
    return {"exit_code": run.code, "ram_note_mb": int(note.group(1)) if note else None,
            "ram_below_minimum_gb": float(error.group(1)) if error else None,
            "allocation_failed_bytes": int(alloc.group(1)) if alloc else None, **usage}


def sd_bounds(values):
    """The sample SD and its one-sided upper 80% bound (chi-square), for replicate counts n = (SD / MCSE)^2."""
    v = np.array([x for x in values if x is not None], dtype=float)
    if len(v) < 2:
        return {"n": len(v), "sd": None, "sd_upper80": None}
    from scipy import stats
    sd = float(v.std(ddof=1))
    return {"n": len(v), "mean": float(v.mean()), "sd": sd,
            "sd_upper80": sd * float(np.sqrt((len(v) - 1) / stats.chi2.ppf(0.2, len(v) - 1)))}


def cmd_pilot(a):
    cells = a.cells.split(",")
    if any(c not in CELLS for c in cells):
        sys.exit(f"--cells takes a comma-separated subset of {', '.join(CELLS)}")
    prov = common.provenance(a, SCRIPT, "pilot", "exploratory", cells=cells, threads=a.threads, seed=a.seed,
                             rounds_cap=a.rounds, replicates=a.replicates, variants=a.variants,
                             caps={"default": CAP, "relaxed": CAP_RELAXED}, scaffold=scaffold_record(a.scaffold))
    out = common.new_out(a.out, prov)
    alisim_bin, binary = (Path(p).expanduser().resolve() for p in (a.alisim, a.binary))
    prov["alisim"] = {"path": str(alisim_bin), "sha256": sha256(alisim_bin)}
    prov["binary"] = {"path": str(binary), "sha256": sha256(binary)}
    tree, loci = read_scaffold(a.scaffold)
    records, rows, held = {}, [], {}
    title = f"SIM-1 pilot ({', '.join(cells)})"

    def heldout(name, tree_text, truth_name):
        """The held-out replicate for a tree and truth, simulated once; None (cached) if any gene failed."""
        if (name, truth_name) not in held:
            d = out / f"heldout-{name}-{truth_name}"
            d.mkdir()
            _, found = simulate_into(d, alisim_bin, tree_text, loci, truth_name, HELDOUT_GENES, a.seed,
                                     HELDOUT_REPLICATE)
            failed = [r["check"] for r in found if r["passed"] is False]
            if failed:
                save(f"heldout-{name}-{truth_name}", {"simulate_failed": failed})
            held[(name, truth_name)] = None if failed else d
        return held[(name, truth_name)]

    def measure(name, d, data, tree_text, truth_name, held_name):
        """The fit's metrics, apart from its record: a failure here keeps the fit's record."""
        try:
            return metrics_of(d, data, heldout(held_name, tree_text, truth_name), binary)
        except Exception:  # noqa: BLE001
            return {"metrics_error": traceback.format_exc()[-2000:]}

    def save(name, rec):
        records[name] = rec
        rows.append(check(name, None, **{k: v for k, v in rec.items() if k not in ("genes",)}))
        common.save(out, prov, title + ", in progress", rows, {"records": records})
        print(f"{name}: exit {rec.get('exit_code')}, {rec.get('wall_s')} s, ended by {rec.get('ended_by')}, "
              f"C {rec.get('C_fit')}", flush=True)

    def cell(name, tree_text, truth_name, genes, replicate, held_name="full", rounds=a.rounds):
        d = out / name
        (d / "data").mkdir(parents=True)
        try:
            manifest, found = simulate_into(d / "data", alisim_bin, tree_text, loci, truth_name, genes, a.seed, replicate)
            if any(r["passed"] is False for r in found):
                save(name, {"simulate_failed": [r["check"] for r in found if r["passed"] is False]})
                return None
            rec = fit_into(d, d / "data", binary, a.threads, CAP, rounds=rounds)
        except Exception:  # noqa: BLE001  a failed cell is recorded; the pilot goes on
            save(name, {"error": traceback.format_exc()[-2000:]})
            return None
        rec.update({"truth": truth_name, "genes": genes, "replicate": replicate})
        save(name, rec)                 # the fit's record first, so a metrics failure cannot lose it
        rec.update(measure(name, d, d / "data", tree_text, truth_name, held_name))
        save(name, rec)
        return d / "data"

    main_data = {}
    for c in cells:
        if c in ("main1", "main"):
            reps = [1] if c == "main1" else range(2 if "main1" in cells else 1, a.replicates + 1)
            for r in reps:
                for t in PILOT_TRUTHS:
                    main_data[(t, r)] = cell(f"main-{t}-r{r}", tree, t, PILOT_GENES, r)
        elif c == "default1":       # NQ.bird's replicate 1 at IQ-TREE's default rounds, beside a capped main1
            cell("default-NQ.bird-r1", tree, "NQ.bird", PILOT_GENES, 1, rounds=None)
        elif c == "ram":
            for g in RAM_GENES:
                d = out / f"ram-{g}"
                (d / "data").mkdir(parents=True)
                try:
                    _, found = simulate_into(d / "data", alisim_bin, tree, loci, "NQ.bird", g, a.seed, 400 + g)
                    failed = [r["check"] for r in found if r["passed"] is False]
                    rec = ({"simulate_failed": failed} if failed
                           else ram_probe(binary, d, "probe", d / "data", a.threads))
                    save(f"ram-{g}", {"genes_count": g, **rec})
                except Exception:  # noqa: BLE001
                    save(f"ram-{g}", {"error": traceback.format_exc()[-2000:]})
        elif c == "taxa26":
            try:
                t26 = prune(tree, pruned_taxa(tree, PRUNED_TAXA, a.seed))
            except Exception:  # noqa: BLE001
                save("taxa26-NQ.bird", {"error": traceback.format_exc()[-2000:]})
                continue
            write(out / "tree26.nwk", t26 + "\n")
            cell("taxa26-NQ.bird", t26, "NQ.bird", PILOT_GENES, 301, held_name="taxa26")
        elif c == "genes100":
            cell("genes100-NQ.bird", tree, "NQ.bird", 100, 201)
        elif c == "variants":
            data = main_data.get(("NQ.bird", 1))
            if data is None:
                save("variants", {"skipped": "main replicate 1 of NQ.bird was not simulated in this run"})
                continue
            for kind, opts, cap in (("p2", {"level2": True}, CAP), ("relaxed", {"relaxed": True}, CAP_RELAXED),
                                    ("gtr20", {"gtr20": True}, CAP)):
                if kind not in a.variants.split(","):
                    continue
                d = out / f"variant-{kind}-NQ.bird-r1"
                d.mkdir()
                try:
                    rec = fit_into(d, data, binary, a.threads, cap, rounds=None if kind == "relaxed" else a.rounds,
                                   **opts)
                except Exception:  # noqa: BLE001
                    save(f"variant-{kind}", {"error": traceback.format_exc()[-2000:]})
                    continue
                save(f"variant-{kind}", rec)
                if kind != "gtr20":
                    rec.update(measure(f"variant-{kind}", d, data, tree, "NQ.bird", "full"))
                    save(f"variant-{kind}", rec)
    summary = {}
    for t in PILOT_TRUTHS:
        reps = [records[k] for k in records if k.startswith(f"main-{t}-r")]
        summary[t] = {m: sd_bounds([r.get(m) for r in reps])
                      for m in ("C_fit", "C_error", "log_error_rms", "heldout_difference", "wall_s", "evaluations")}
    rows.append(check("pilot summary: SDs and their upper 80% bounds over main replicates", None, **summary))
    common.finish(out, prov, title, rows, {"records": records, "summary": summary})


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="step", required=True)
    s = sub.add_parser("simulate")
    s.add_argument("--scaffold", required=True)
    s.add_argument("--truth", required=True)
    s.add_argument("--genes", type=int, required=True)
    s.add_argument("--replicate", type=int, required=True)
    s.add_argument("--alisim", required=True)
    s.add_argument("--seed", type=int, default=SEED)
    f = sub.add_parser("fit")
    f.add_argument("--data", required=True)
    f.add_argument("--binary", required=True)
    f.add_argument("--threads", type=int, default=1)
    f.add_argument("--level2", action="store_true")
    f.add_argument("--relaxed", action="store_true")
    f.add_argument("--gtr20", action="store_true")
    f.add_argument("--rounds", type=int, help="cap the linked rounds")
    f.add_argument("--timeout", type=float, default=4 * 3600)
    c = sub.add_parser("checks")
    c.add_argument("--scaffold", required=True)
    c.add_argument("--truth", required=True)
    c.add_argument("--alisim", required=True)
    c.add_argument("--binary", required=True)
    c.add_argument("--genes", type=int, default=20)    # about 10,800 columns: |lnL| well above 021's 5e4
    c.add_argument("--seed", type=int, default=SEED)
    c.add_argument("--timeout", type=float, default=3600)
    m = sub.add_parser("metrics")
    m.add_argument("--fit", required=True, help="a fit step's output directory")
    m.add_argument("--heldout", help="a simulate step's output directory of held-out genes")
    m.add_argument("--binary", help="needed with --heldout")
    m.add_argument("--timeout", type=float, default=3600)
    p_ = sub.add_parser("pilot")
    p_.add_argument("--scaffold", required=True)
    p_.add_argument("--alisim", required=True)
    p_.add_argument("--binary", required=True)
    p_.add_argument("--threads", type=int, default=6)
    p_.add_argument("--cells", default=",".join(CELLS))
    p_.add_argument("--rounds", type=int, help="cap the linked rounds of every fit but default1 and relaxed")
    p_.add_argument("--replicates", type=int, default=PILOT_REPLICATES, help="main cell replicates per truth")
    p_.add_argument("--variants", default="p2,relaxed,gtr20")
    p_.add_argument("--seed", type=int, default=SEED)
    for p in (s, f, c, m, p_):
        p.add_argument("--out", help="new directory outside the repository (default under ~/iqtree3-runs/sim1/)")
        p.add_argument("--git", default="git", help="git executable (git.exe under WSL)")
        p.add_argument("--allow-dirty", action="store_true")
    a = ap.parse_args()
    {"simulate": cmd_simulate, "fit": cmd_fit, "checks": cmd_checks, "metrics": cmd_metrics,
     "pilot": cmd_pilot}[a.step](a)


if __name__ == "__main__":
    main()
