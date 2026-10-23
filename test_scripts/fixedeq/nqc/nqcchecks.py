#!/usr/bin/env python3
"""S2, S3a and S3b checks of the compiled NQC model (docs/agent/CODE_PLAN.md section 4; decision 021).

  python nqcchecks.py --binary ~/iqtree3-build/iqtree3 --alisim ~/iqtree3-baseline/iqtree3 \
      --calibration CALIBRATION/results.json --older-run S2_GATE/nqc/fit-level1 --out DIR \
      [--cases root,seed,...] [--git git.exe]

Cases:
  root      a rooted input tree keeps its root and an unrooted one is rooted at the midpoint
  seed      the seed Q and the first log-likelihood against the oracle (decision 021), +R4 included
  gradient  the first gradient against the oracle, by the same rule and against a central
            reference, at the LG seed and at a start near zero (decision 021)
  target    targets refused for their sum and for an entry below the minimum
  refusals  every combination decision 008 refuses
  fits      Level 1 and Level 2 fits on aa_example, +I+G4 and free +R4 included, and the domain
            option
  nesting   NQC seeded at a fitted GTR20+F{pi*}+G4 in a widened domain (decision 025)
  skewed    a Level 1 fit at the test target cubed and renormalized, reported (decision 023)
  restart   --undo with another target or domain refused, with the same command accepted, and
            from S2's checkpoint (needs --older-run)
  risk22    Level 1 and Level 2 +G4 at --nqc-domain 1e-8,1e7 beside the default box (PLAN.md
            risk 22), reported
  toy       the compiled optimizer against the oracle's port (decisions 017, 021), using the
            line measured by explore/toy_calibration.py (needs --alisim and --calibration)
  short     a Level 1 and a +I+G4 Level 2 fit on 300 sites of the toy data, one round each, for
            the sanitizer build (CODE_PLAN.md section 3.6; needs --alisim)

S3b's linked cases (CODE_PLAN.md sections 3.5 and 4, S3b), on nqc/twopart.py's data (needs
--alisim) unless they name turtle:
  linked-seed      -lnL at the first linked gradient against the oracle's linked objective, under -S
                   and -p with +G4{0.5} and the given +R4; the printed steps against decision 005's
                   rule, bitwise; the gradient lines; and a -S start near zero
  linked-fits      Level 1 under -S, -q (decision 024, nuisances from the -p Level 2 +G4 fit) and
                   --link-model (-S); Level 2 under -S and -p with +G4, +I+G4 and free +R4; at every
                   numbered checkpoint pi*, coordinates, references and residuals in every struct
  linked-refusals  --link-model with another target and with the target at 8 decimals, refused at
                   the first linked fit; a charset without a model and no -m; a charset with +H4
  permuted         twopart -S Level 1 with the charsets reversed, bitwise; turtle at the seed under -p
                   and -S with the charsets (and -S trees) reversed, by decision 038's lines
  threads          turtle at the seed, -p at -T 4 and -S at -T 3, 10 runs each against -T 1 (decision
                   038); 5 -p Level 2 fits at -T 4 without rate variation, the setup of NONREV's
                   measured spread (explore/thread_counts.py "pL2"), reported
  linked-restart   -p and -S Level 2 +G4 resumed from numbered checkpoint #1; --undo of the -S Level
                   1 fit with the same command, another target and another domain
  risk22-turtle    turtle Level 1 -S (+G4{0.5}) and Level 2 +G4 under -p and -S, in the default box
                   and at --nqc-domain 1e-8,1e7, reported (PLAN.md risk 22)
  warm-starts      two Level 1 continuations, in the box and at 1e-8,1e7, from the final state of
                   each default-box fit of risk22, nesting and risk22-turtle, reported
  linked-short     a Level 1 -S and a Level 2 -p +G4 fit on 150 sites of each twopart partition, for
                   the sanitizer build

Every run has its own directory under DIR and uses -seed 1 and -T 1 unless a case says otherwise.
Final states are read from the class's checkpoint keys; the seed Q and the first gradient from its
verbose output (decision 021). How often something is printed is reported, not asserted. Exit
status 1 if an asserted check fails; checks marked "report" are not asserted.
"""
import argparse
import datetime
import json
import math
import platform
import re
import sys
import traceback
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE.parent / "probes"))
sys.path.insert(0, str(HERE.parent / "differential"))
from oracle import chart, likelihood as lk, nqc, optimize, readers, target  # noqa: E402
from probes import (AA, NEX, PREFIX, TURTLE, Runner, baseline_runs, first_taxa, fitted_gtr20, mdef,  # noqa: E402
                    paml_rows, read_checkpoint, tree_comparison, tree_text, turtle_target, unrooted_text,
                    with_charpartition)
from differential import datasets, significant_digits, write_matrix  # noqa: E402
from regress import REPO, SCALARS, capture, read, sha256, write  # noqa: E402
import toy  # noqa: E402
import twopart  # noqa: E402

SCRIPT = Path(__file__).resolve().relative_to(REPO).as_posix()
DNA = REPO / "example" / "example.phy"
TOL_LNL = 1e-8          # decision 021
RESIDUAL = 1e-10        # decision 026
TAU = 1e-12             # decision 030: relative disagreement allowed between two likelihoods
CENTRAL_STEPS = (1e-2, 1e-3, 1e-4)
ALPHA, PINV = 0.5, 0.2
NEAR_ZERO = 1e-9
BRANCHES_AA = 34        # free branch lengths of the rooted aa_example tree (baseline runs 2 and 4)
NAMED = "NQCTARGET"
# dyadic FreeRate values with sum(w) = 1 and sum(w r) = 1, so RateFree's normalization is the identity
R4_GIVEN = (0.25, 0.125, 0.25, 0.625, 0.25, 1.125, 0.25, 2.125)
R4 = "+R4{" + ",".join(f"{v:g}" for v in R4_GIVEN) + "}"
WIDE = "1e-8,1e7"       # the nesting test's aa_example domain, reused for risk 22 (CODE_PLAN.md section 4)

INITIAL = re.compile(r"^1\. Initial log-likelihood: (\S+)", re.M)
SEED_Q = "NQC seed Q after rescaling (17 significant digits):"
GRADIENT_HEAD = re.compile(r"NQC first gradient at -lnL (\S+) ", re.M)
GRADIENT = re.compile(r"^NQC gradient (\d+) (\S+) (\S+) (\S+)$", re.M)
FIT = re.compile(r"NQC fit (\d+): start lnL (\S+), optimizer lnL (\S+), re-scored lnL (\S+)$", re.M)
RESIDUALS = re.compile(r"NQC residuals: \|pi\*Q\| (\S+), row sums (\S+), \|mean rate - 1\| (\S+), nu (\S+), "
                       r"minimum rate (\S+)$", re.M)
BOUNDS = re.compile(r"NQC coordinates at the lower bound (\d+), at the upper bound (\d+), outside the domain (\d+)$",
                    re.M)
SWITCH = "INFO: Switch to scaling-squaring"
FALL = re.compile(r"WARNING: NQC fit (\d+): the re-scored lnL is below the start by (\S+)", re.M)
BUG_ABORT = "BUG: Tree logl gets worse after model optimization"
SUMMARY = "NQC summary of the final model (17 significant digits):"
RESTORE_REFUSAL = "the checkpoint was written with a different"


# ---------------------------------------------------------------- pass lines (decisions 026 and 030)

def q_line(kappa):
    return max(1e-12, 1e-14 * kappa)


def same_rule_bound(lnl, h):
    return 2 * TAU * abs(lnl) / np.abs(h)


def step_bound(lnl, h, curvature):
    return 2 * (np.abs(h) / 2 * np.abs(curvature) + 2 * TAU * abs(lnl) / np.abs(h))


# ---------------------------------------------------------------- reading a run

def seed_q(text):
    lines = text.splitlines()
    i = next(k for k, line in enumerate(lines) if line.endswith(SEED_Q))
    return np.array([[float(v) for v in row.split(",")] for row in lines[i + 1:i + 21]])


def first_gradient(text):
    rows = GRADIENT.findall(text)
    head = GRADIENT_HEAD.search(text)
    if not rows or not head:
        return None
    a = np.array([[float(v) for v in r[1:]] for r in rows])
    return float(head.group(1)), a[:, 0], a[:, 1], a[:, 2]


def nqc_structs(flat):
    """The class's keys in a flat checkpoint, one dict per ModelNonrevFixedEq struct, keyed by the
    prefix before the struct's name: "" for one alignment; for partitioned data, the linked model's
    name (NQC) or a partition's, under the PartitionModel struct."""
    structs = {}
    for k, v in flat.items():
        if "ModelNonrevFixedEq!" in k:
            prefix, key = k.split("ModelNonrevFixedEq!", 1)
            structs.setdefault(prefix, {})[key] = v
    return structs


def nqc_checkpoint(run):
    return nqc_structs(run.checkpoint())


def single_checkpoint(run):
    """The keys of a single-alignment run's one NQC struct, or None if it has not exactly one."""
    structs = nqc_checkpoint(run)
    return next(iter(structs.values())) if len(structs) == 1 else None


def linked_struct(structs):
    """The linked model's struct, which PartitionModel saves under the model's name, or None."""
    return next((s for prefix, s in structs.items() if prefix.endswith("!NQC!")), None)


def final_struct(run):
    """A run's final NQC state: the one struct of a single alignment, else the linked struct."""
    return single_checkpoint(run) or linked_struct(nqc_checkpoint(run))


def summary(text):
    """The end-of-run summary's lines, or None if it was not printed."""
    lines = text.splitlines()
    i = next((k for k, line in enumerate(lines) if line == SUMMARY), None)
    if i is None:
        return None
    block = []
    for line in lines[i + 1:]:
        if not line.startswith("NQC "):
            break
        block.append(line)
    return block


def floats(text):
    return np.array([float(v) for v in text.split(",")])


def check(name, passed, **values):
    """passed is True or False for an asserted check, None for a reported one."""
    return {"check": name, "passed": passed, **values}


def failed_run(run):
    return [check(f"{run.name} completes", False, exit_code=run.code, tail=run.stdout()[-1500:])]


# ---------------------------------------------------------------- the oracle's side

class Context:
    def __init__(self, out):
        self.target_text = turtle_target()
        self.pi = target.parse_target(self.target_text)[0]
        self.seed, self.ref = nqc.lg_seed(self.pi)
        self.z = chart.inverse(self.seed, self.ref)
        self.kappa = chart.nu_condition(chart.jump_of(self.seed))
        base = baseline_runs()
        self.t02 = tree_text(base, "run02")
        self.u02 = unrooted_text(self.t02)
        self.aa = readers.read_phylip(AA)
        self.tree02 = readers.parse_newick(self.t02)
        data_dir = out / "data"
        data_dir.mkdir()
        self.data = datasets(data_dir)
        self.model = f"NQC+F{{{self.target_text}}}"

    def named_files(self):
        return {"target.nex": mdef(NAMED, self.target_text)}


def forward(f, x, h):
    """Forward differences at the steps the compiled class took."""
    fx = f(x)
    y = np.array(x, dtype=float)
    g = np.empty(len(x))
    for k in range(len(x)):
        y[k] = x[k] + h[k]
        g[k] = (f(y) - fx) / h[k]
        y[k] = x[k]
    return fx, g


def gradient_checks(label, run, objective):
    """Decision 030's two gradient lines at the point the class printed."""
    got = first_gradient(run.stdout())
    if got is None:
        return [check(f"{label}: gradient printed", False)], None, None, None, None
    f_iq, x, h, g = got
    f_or, g_or = forward(objective, x, h)
    lnl = -f_or
    out = [check(f"{label}: -lnL at the point, against the oracle (decision 021)",
                 abs(f_iq - f_or) <= TOL_LNL * abs(f_or), iqtree=f_iq, oracle=f_or)]
    h_rule = nqc.scaled_steps(x)
    out.append(check(f"{label}: steps against decision 005's rule (report)", None,
                     max_relative_difference=float(np.max(np.abs(h - h_rule) / np.abs(h_rule)))))
    ratio = np.abs(g - g_or) / same_rule_bound(lnl, h)
    out.append(check(f"{label}: same rule, |g - g_oracle| / line", bool(np.all(ratio <= 1)),
                     worst=float(ratio.max()), worst_index=int(ratio.argmax())))
    refs = {s: nqc.central(objective, x, s) for s in CENTRAL_STEPS}
    g_ref, curvature = refs[1e-3]
    ratio = np.abs(g - g_ref) / step_bound(lnl, h, curvature)
    spread = max(float(np.max(np.abs(refs[s][0] - g_ref) / np.maximum(1.0, np.abs(g_ref)))) for s in (1e-2, 1e-4))
    out.append(check(f"{label}: against the central reference, |g - g_ref| / line", bool(np.all(ratio <= 1)),
                     worst=float(ratio.max()), worst_index=int(ratio.argmax()),
                     central_spread_relative=spread))
    return out, x, g_ref, curvature, lnl


# ---------------------------------------------------------------- cases

def case_root(r, c):
    out = []
    start = first_taxa(AA)[0]
    for label, text in (("rooted", c.t02), ("unrooted", c.u02)):
        run = r.run(f"root-{label}", ["-s", AA, "-m", c.model + f"+G4{{{ALPHA}}}", "-te", "input.tree", "-blfix",
                                      "-nparam", "2"], {"input.tree": text})
        if run.code != 0 or not run.trees():
            out += failed_run(run)
            continue
        comp = tree_comparison(text.strip(), run.trees()[0], start)
        out.append(check(f"{label} input: root split as NONREV places it", comp["root_split_kept"],
                         root_expected=comp["root_expected"], root_out=comp["root_out"]))
        out.append(check(f"{label} input: root-adjacent lengths", comp["root_lengths_kept"]))
    return out


def case_seed(r, c):
    out, printed = [], {}
    cases = [("aa-none", "aa_example", "", None, 0.0, None),
             ("aa-G4", "aa_example", f"+G4{{{ALPHA}}}", ALPHA, 0.0, None),
             ("aa-I", "aa_example", f"+I{{{PINV}}}", None, PINV, None),
             ("aa-IG4", "aa_example", f"+I{{{PINV}}}+G4{{{ALPHA}}}", ALPHA, PINV, None),
             ("aa-R4", "aa_example", R4, None, 0.0, R4_GIVEN),
             ("turtle1-G4", "turtle_part1", f"+G4{{{ALPHA}}}", ALPHA, 0.0, None)]
    Q_or = chart.build(c.z, c.pi, c.ref)[0]
    for name, dname, suffix, shape, p_invar, free in cases:
        aln_path, tree_path = c.data[dname]
        run = r.run(f"seed-{name}", ["-s", aln_path, "-m", c.model + suffix, "-te", tree_path, "--show-lh",
                                     "-nparam", "2"])
        m = INITIAL.search(run.stdout())
        if run.code != 0 or not m:
            out += failed_run(run)
            continue
        printed[name] = m.group(1)
        aln = readers.read_phylip(aln_path) if aln_path.suffix == ".phy" else readers.read_fasta(aln_path)
        oracle = lk.alignment_log_likelihood(aln, readers.read_newick(tree_path), Q_or, c.pi, shape=shape,
                                             p_invar=p_invar, free=free)
        iq = float(m.group(1))
        out.append(check(f"{name}: first log-likelihood against the oracle (decision 021)",
                         significant_digits(m.group(1)) >= 10 and abs(iq - oracle) <= TOL_LNL * abs(iq),
                         iqtree=m.group(1), oracle=oracle, relative=abs(iq - oracle) / abs(iq)))
        Q_iq = seed_q(run.stdout())
        rel = float(np.max(np.abs(Q_iq - Q_or) / np.abs(Q_or)))
        out.append(check(f"{name}: seed Q against the oracle (decision 026)", rel <= q_line(c.kappa),
                         max_relative=rel, line=q_line(c.kappa), kappa=c.kappa))
    aln_path, tree_path = c.data["aa_example"]
    run = r.run("seed-aa-G4-named", ["-s", aln_path, "-m", f"NQC+F{NAMED}+G4{{{ALPHA}}}", "-te", tree_path,
                                     "--show-lh", "-nparam", "2", "--mdef", "target.nex"], c.named_files())
    m = INITIAL.search(run.stdout())
    if run.code != 0 or not m:
        out += failed_run(run)
    else:
        out.append(check("+F<name> through --mdef gives the same first log-likelihood as +F{...}",
                         m.group(1) == printed.get("aa-G4"), named=m.group(1), literal=printed.get("aa-G4")))
    return out


def level1_args(model):
    return ["-s", AA, "-m", model + f"+G4{{{ALPHA}}}", "-te", "input.tree", "-blfix", "-vv"]


def level2_args(model, suffix):
    return ["-s", AA, "-m", model + suffix, "-te", "input.tree", "-v"]


def level1(r, c, ctx):
    """The Level 1 fit that the gradient, fits, restart and risk22 cases share."""
    if "level1" not in ctx:
        ctx["level1"] = r.run("fit-level1", level1_args(c.model), {"input.tree": c.t02})
    return ctx["level1"]


def case_gradient(r, c, ctx):
    out = []
    objective = nqc.Objective(c.aa, c.tree02, c.pi, c.ref, shape=ALPHA)
    run = level1(r, c, ctx)
    if run.code != 0:
        return failed_run(run)
    checks, x, *_ = gradient_checks("seed", run, objective)
    out += checks
    if x is not None:
        out.append(check("seed: printed coordinates against the oracle's LG seed (report)", None,
                         max_abs_difference=float(np.max(np.abs(x - c.z)))))

    z = c.z.copy()
    z[::18] = NEAR_ZERO   # the first free coordinate of each row
    text_dir = r.out / "nearzero-input"
    text_dir.mkdir()
    write_matrix(text_dir / "nearzero.qmat", chart.build(z, c.pi, c.ref)[0], c.pi)
    run = r.run("gradient-nearzero", ["-s", AA, "-m", c.model + f"+G4{{{ALPHA}}}", "-te", "input.tree", "-blfix",
                                      "-nparam", "3", "-vv", "--init-model", "nearzero.qmat"],
                {"input.tree": c.t02, "nearzero.qmat": read(text_dir / "nearzero.qmat")})
    if run.code != 0:
        return out + failed_run(run)
    checks, x, g_ref, curvature, lnl = gradient_checks("near zero", run, objective)
    out += checks
    if x is None:
        return out
    near = np.abs(x) < 1e-8
    out.append(check("near zero: coordinates within 1e-8 of zero", int(near.sum()) == 20, count=int(near.sum())))
    # The control D02 rests on: relative error at the near-zero coordinates, scaled step against
    # IQ-TREE's legacy step. Decision 030's step-size line is not used here, because its rounding
    # term grows as 1/h and so forgives the legacy step's tiny steps.
    _, g_legacy = optimize.legacy_derivative(objective, np.array(x))
    got = first_gradient(run.stdout())
    scale = np.maximum(1.0, np.abs(g_ref[near]))
    out.append(check("near zero: relative error against the central reference, scaled and legacy steps (report)",
                     None, scaled=float(np.max(np.abs(got[3][near] - g_ref[near]) / scale)),
                     legacy=float(np.max(np.abs(g_legacy[near] - g_ref[near]) / scale)),
                     legacy_steps=[float(v) for v in nqc.legacy_steps(x[near])[:3]]))
    return out


def case_target(r, c):
    out = []
    five = ",".join(f"{v:.5f}" for v in c.pi)
    low = c.pi.copy()
    low[0] += low[17] - 5e-5
    low[17] = 5e-5
    for name, text, why in (("sum", five, "a sum more than 1e-6 from 1"),
                            ("minimum", ",".join(f"{v:.17g}" for v in low), "an entry below --min-freq")):
        run = r.run(f"target-{name}", ["-s", AA, "-m", f"NQC+F{{{text}}}", "-te", "input.tree", "-nparam", "2"],
                    {"input.tree": c.t02})
        out.append(check(f"refused: {why}", run.code not in (0, "timeout") and "NQC: invalid target" in run.stdout(),
                         exit_code=run.code))
    out.append(check("accepted through +F{...} and through --mdef: see the seed case", None))
    return out


def case_refusals(r, c):
    t = c.target_text
    two = c.t02.strip() + "\n" + c.t02.strip() + "\n"
    named = c.named_files()
    fmix = {"target.nex": f"#nexus\nbegin models;\nfrequency {NAMED} = {t};\nfrequency {NAMED}2 = {t};\nend;\n"}
    mixdef = {"target.nex": f"#nexus\nbegin models;\nfrequency {NAMED} = {t};\n"
                            f"model NQCMIX = MIX{{NQC+F{NAMED},LG}};\nend;\n"}
    aa = ["-s", AA, "-te", "input.tree", "-nparam", "2"]
    tree = {"input.tree": c.t02}
    runs = [
        ("dna", ["-s", DNA, "-m", c.model, "-nparam", "2"], {}),
        ("freq-F", aa + ["-m", "NQC+F"], tree),
        ("freq-FQ", aa + ["-m", "NQC+FQ"], tree),
        ("freq-FO", aa + ["-m", "NQC+FO"], tree),
        ("freq-FU", aa + ["-m", "NQC+FU"], tree),
        ("freq-none", aa + ["-m", "NQC"], tree),
        ("braces", aa + ["-m", f"NQC{{1,2}}+F{{{t}}}"], tree),
        ("optfromgiven", aa + ["-m", c.model, "-optfromgiven"], tree),
        ("jointopt", aa + ["-m", c.model, "-jointopt"], tree),
        # --model-joint and --link-model were refused until S3b; the linked cases now fit them
        ("tree-mixture", ["-s", AA, "-m", c.model + "+T", "-te", "two.trees", "-nparam", "2"], {"two.trees": two}),
        ("ghost", aa + ["-m", c.model + "+H4"], tree),
        ("site-freq", aa + ["-m", c.model, "-ft", "input.tree"], tree),
        # IQ-TREE itself refuses the flag unless GTR or GTR20 is in the model
        ("link-exchange", aa + ["-m", f"MIX{{NQC+F{NAMED},GTR20}}", "--mdef", "target.nex", "--link-exchange-rates"],
         {**named, **tree}),
        ("eigen", aa + ["-m", c.model, "--eigen"], tree),
        ("modelfinder", ["-s", AA, "-m", "MFP", "-mset", f"NQC+F{NAMED}", "--mdef", "target.nex", "-te", "input.tree"],
         {**named, **tree}),
        # S3a: ModelFinder with no -m at all, which S2's test of the -m string missed
        ("modelfinder-no-m", ["-s", AA, "-mset", f"NQC+F{NAMED}", "--mdef", "target.nex", "-te", "input.tree"],
         {**named, **tree}),
        # IQ-TREE reads the site file before it builds the model, so the file must be valid
        ("site-freq-fs", aa + ["-m", c.model, "-fs", "site.freq"], {**tree, "site.freq": f"1 {t.replace(',', ' ')}\n"}),
        ("hmm", ["-s", AA, "-m", c.model + "+T", "-te", "two.trees", "-hmmster", "-nparam", "2"], {"two.trees": two}),
        ("alisim", ["--alisim", "sim", "-m", c.model, "-t", "input.tree", "--seqtype", "AA"], tree),
        ("mix", aa + ["-m", f"MIX{{NQC+F{NAMED},LG}}", "--mdef", "target.nex"], {**named, **tree}),
        ("fmix", aa + ["-m", f"NQC+FMIX{{{NAMED},{NAMED}2}}", "--mdef", "target.nex"], {**fmix, **tree}),
        ("mdef-mix", aa + ["-m", "NQCMIX", "--mdef", "target.nex"], {**mixdef, **tree}),
        ("unlinked-p", ["-s", TURTLE, "-p", NEX, "-m", c.model, "-nparam", "2"], {}),
        ("unlinked-q", ["-s", TURTLE, "-q", NEX, "-m", c.model, "-nparam", "2"], {}),
        ("unlinked-S", ["-s", TURTLE, "-S", NEX, "-m", c.model, "-nparam", "2"], {}),
    ]
    out = []
    for name, args, files in runs:
        run = r.run(f"refuse-{name}", args, files)
        message = nqc_error(run)
        out.append(check(f"refused: {name}", run.code not in (0, "timeout") and message is not None,
                         exit_code=run.code, message=message))
    return out


def nqc_error(run):
    """The first error line that names NQC, or None."""
    return next((l for l in run.stdout().splitlines() if "NQC" in l and "ERROR" in l), None)


START_LNL = re.compile(r"Init pinv, alpha: \S+ \S+ / Estimate: \S+ \S+ / LogL: (\S+)")
BEST_LNL = re.compile(r"Optimal pinv,alpha: \S+ \S+ / LogL: (\S+)")


def restore_check(label, text):
    """The +I+G path keeps the best of its starts through the model's checkpoint (decision 033): the
    model it ends with must score as the best start did, at the 3 decimals IQ-TREE prints. For one
    alignment, which takes the path once; block_restore_check covers -S."""
    starts, best = [float(v) for v in START_LNL.findall(text)], BEST_LNL.search(text)
    ok = bool(starts and best and abs(float(best.group(1)) - max(starts)) <= 5e-4)
    return check(f"{label}: the restored model scores as the best +I+G start", ok,
                 best_start=max(starts) if starts else None, final=best.group(1) if best else None)


def block_restore_check(label, text):
    """restore_check for -S, which keeps IQ-TREE's +I+G multi-start on (only -p, -q and -blfix turn
    it off), so each partition takes the path in every round of the first model optimization: every
    "Optimal pinv,alpha" line must score as the best start printed since the one before it."""
    blocks, starts = [], []
    for line in text.splitlines():
        m = START_LNL.search(line)
        if m:
            starts.append(float(m.group(1)))
            continue
        m = BEST_LNL.search(line)
        if m:
            blocks.append((max(starts) if starts else None, float(m.group(1))))
            starts = []
    bad = [i for i, (s, b) in enumerate(blocks) if s is None or abs(b - s) > 5e-4]
    gaps = [abs(b - s) for s, b in blocks if s is not None]
    return check(f"{label}: in each block of +I+G starts, the restored model scores as the best start", bool(blocks)
                 and not bad, blocks=len(blocks), failing=bad, largest_gap=max(gaps, default=None))


KEYS = ("target", "residuals", "coordinates", "domain", "bounds")
RECORD_KEYS = ("fits", "accepted_fit", "evaluations", "start_score", "optimizer_score", "accepted_score", "bounds")


def key_checks(label, ck, pi):
    """Asserted on the class's checkpoint keys (decisions 003, 021): pi* bitwise, the three residuals
    of decision 021's line, and the bound counts against a recount from the coordinates and domain.
    pi* after every earlier fit is the class's own guard, observed as exit 0."""
    missing = [k for k in KEYS if ck is None or k not in ck]
    if missing:
        return [check(f"{label}: one NQC struct in the checkpoint, with its keys", False, missing=missing)]
    res = floats(ck["residuals"])  # fixedeq::Diagnostics order
    z, (lo, hi) = floats(ck["coordinates"]), floats(ck["domain"])
    recount = [int(np.sum(z == lo)), int(np.sum(z == hi)), int(np.sum((z < lo) | (z > hi)))]
    keys = [int(v) for v in ck["bounds"].split(",")]
    return [
        check(f"{label}: target key bitwise equal to the parsed target", bool(np.array_equal(floats(ck["target"]), pi))),
        check(f"{label}: residual keys at most {RESIDUAL:g} (decision 021)", bool(np.all(res[:3] <= RESIDUAL)),
              stationarity_row_sums_mean_rate=res[:3].tolist(), minimum_rate=res[3], minimum_flux=res[4], nu=res[5],
              exit_rate_mismatch=res[6]),
        check(f"{label}: bound keys equal a recount from the key coordinates", keys == recount, keys=keys,
              recount=recount),
    ]


def monitoring(label, run, ck):
    """CODE_PLAN.md section 4's monitoring, reported: the fit record, falls below the start, the
    optimizer's score against the re-scored one, IQ-TREE's hard stops and scaling-squaring."""
    text = run.stdout()
    fits = FIT.findall(text)  # printed under -v only
    gaps = [abs(float(f[2]) - float(f[3])) for f in fits]
    return [
        check(f"{label}: fit record keys (report)", None, **{k: (ck or {}).get(k) for k in RECORD_KEYS}),
        check(f"{label}: fits ending below their start (report; decision 009)", None, records_printed=len(fits),
              falls=FALL.findall(text)),
        check(f"{label}: optimizer score against the re-scored point, largest gap over printed records (report)",
              None, largest=max(gaps) if gaps else None),
        check(f"{label}: summary printed (report)", None, lines=summary(text)),
        check(f"{label}: scaling-squaring switches (report)", None, count=text.count(SWITCH)),
        check(f"{label}: IQ-TREE's ASSERT or 'Tree logl gets worse' abort (report)", None,
              assertion="Assertion" in text or "ASSERT" in text, bug_abort=BUG_ABORT in text, exit_code=run.code),
    ]


def free_parameters_check(label, run, expected):
    m = re.search(SCALARS["free_parameters"], run.report(), re.M)
    return check(f"{label}: free parameters", bool(m) and int(m.group(1)) == expected,
                 reported=m.group(1) if m else None, expected=expected)


def fit_checks(label, run, c, expected_parameters, level):
    """A single-alignment fit's checks."""
    if run.code != 0:
        return failed_run(run) + monitoring(label, run, None)
    text, ck = run.stdout(), single_checkpoint(run)
    out = key_checks(label, ck, c.pi)
    out.append(free_parameters_check(label, run, expected_parameters))
    comp = tree_comparison(c.t02.strip(), run.trees()[0]) if run.trees() else None
    out.append(check(f"{label}: root split unchanged", bool(comp and comp["root_split_kept"])))
    if level == 1:
        out.append(check(f"{label}: root-adjacent lengths unchanged", bool(comp and comp["root_lengths_kept"])))
    if "+I+G" in label:
        out.append(restore_check(label, text))
    return out + monitoring(label, run, ck)


def slug(label):
    return label.lower().replace(" ", "-").replace("+", "")


def case_fits(r, c, ctx):
    out = []
    run1 = level1(r, c, ctx)
    out += fit_checks("Level 1", run1, c, 360, 1)
    ck = single_checkpoint(run1) if run1.code == 0 else None
    if not (ck and "accepted_score" in ck):
        out.append(check("Level 1: accepted score against the oracle's re-evaluation (decisions 021, 030)", False,
                         reason="no accepted_score key"))
    else:
        z = floats(ck["coordinates"])
        ref_ok = [int(v) for v in ck["references"].split(",")] == [int(v) for v in c.ref]
        out.append(check("Level 1: references in the checkpoint equal the oracle's", ref_ok))
        oracle = -nqc.Objective(c.aa, c.tree02, c.pi, c.ref, shape=ALPHA)(z)
        accepted = float(ck["accepted_score"])
        out.append(check("Level 1: accepted score against the oracle's re-evaluation (decisions 021, 030)",
                         abs(accepted - oracle) <= TOL_LNL * abs(oracle), accepted=accepted, oracle=oracle))
        printed = re.search(SCALARS["log_likelihood"], run1.report(), re.M)
        out.append(check("Level 1: printed final score against the accepted score (report)", None,
                         printed=printed.group(1) if printed else None,
                         gap=abs(float(printed.group(1)) - accepted) if printed else None))
    # free rate parameters: +G4 1, +I+G4 2, +R4 6 (CODE_PLAN.md section 3.5)
    for label, suffix, extra in (("Level 2 +G4", "+G4", 1), ("Level 2 +I+G4", "+I+G4", 2), ("Level 2 +R4", "+R4", 6)):
        if label not in ctx:  # risk22 may have run it
            ctx[label] = r.run("fit-" + slug(label), level2_args(c.model, suffix), {"input.tree": c.t02})
        out += fit_checks(label, ctx[label], c, 360 + extra + BRANCHES_AA, 2)
    run = r.run("fit-domain", ["-s", AA, "-m", c.model + f"+G4{{{ALPHA}}}", "-te", "input.tree", "-blfix",
                               "-nparam", "3", "--nqc-domain", "1e-6,100"], {"input.tree": c.t02})
    if run.code != 0:
        out += failed_run(run)
    else:
        domain = floats((single_checkpoint(run) or {}).get("domain", "nan,nan"))
        want = np.array([math.log(1e-6), math.log(100.0)])
        out.append(check("--nqc-domain 1e-6,100 sets the domain in the checkpoint", bool(np.array_equal(domain, want)),
                         domain=domain.tolist()))
    return out


def widened_domain(z):
    """CODE_PLAN.md section 3.5: one decade beyond the start's smallest and largest jump-weight
    ratios, rounded outward to powers of ten, and never inside the default box."""
    ratios = np.exp(z)
    lo = min(10.0 ** math.floor(math.log10(ratios.min() / 10)), 1e-5)
    hi = max(10.0 ** math.ceil(math.log10(ratios.max() * 10)), 10.0)
    return f"{lo:g},{hi:g}", float(ratios.min()), float(ratios.max())


def case_nesting(r, c, ctx):
    """Decision 025: NQC seeded at a fixed-tree GTR20+F{pi*}+G4 fit, carried from its checkpoint,
    starts at that fit's --show-lh evaluation and ends no lower, in a domain widened to hold the
    start; the same NQC fit in the default box is reported, and kept for the warm-starts case."""
    pif = f"+F{{{c.target_text}}}"
    gtr = r.run("nesting-gtr20", ["-s", AA, "-m", "GTR20" + pif + "+G4", "-te", "input.tree"], {"input.tree": c.u02})
    fit = fitted_gtr20(gtr)
    if fit is None:
        return failed_run(gtr)
    files = {"gtr.treefile": "\n".join(gtr.trees()) + "\n",
             "start.paml": "\n".join(" ".join(row) for row in paml_rows(fit["rates"])) + "\n"
                           + c.target_text.replace(",", " ") + "\n"}
    g4 = "+G4{" + fit["alpha"] + "}"
    ctx["nesting"] = {"files": files, "g4": g4}
    ev = r.run("nesting-evaluation", ["-s", AA, "-m", "start.paml" + pif + g4, "-te", "gtr.treefile", "--show-lh"],
               files)
    m = INITIAL.search(ev.stdout())
    if ev.code != 0 or not m:
        return failed_run(ev)
    reference = float(m.group(1))
    R = np.zeros((20, 20))
    R[np.triu_indices(20, 1)] = [float(v) for v in fit["rates"]]  # IQ-TREE's pair order
    z = chart.inverse(chart.reversible_seed(R + R.T, c.pi), c.ref)
    domain, ratio_min, ratio_max = widened_domain(z)
    line = TOL_LNL * abs(reference)
    out = [check("carried start: jump-weight ratios and the widened domain (report)", None, ratio_min=ratio_min,
                 ratio_max=ratio_max, domain=domain, outside_default_box=int(np.sum((z < nqc.LOWER) | (z > nqc.UPPER))),
                 reference=m.group(1), alpha=fit["alpha"])]
    for label, extra in (("widened domain", ["--nqc-domain", domain]), ("default box", [])):
        run = r.run("nesting-nqc-" + label.replace(" ", "-"),
                    ["-s", AA, "-m", c.model + g4, "-te", "gtr.treefile", "-blfix", "-vv", "--init-model", "start.paml",
                     *extra], files)
        got = first_gradient(run.stdout()) if run.code == 0 else None
        ck = single_checkpoint(run) if run.code == 0 else None
        first = -got[0] if got else None
        final = float(ck["accepted_score"]) if ck and "accepted_score" in ck else None
        values = {"exit_code": run.code, "first_minus_reference": None if first is None else first - reference,
                  "final_minus_reference": None if final is None else final - reference, "line": line,
                  "bounds": (ck or {}).get("bounds"), "fits": (ck or {}).get("fits"), "falls": FALL.findall(run.stdout())}
        if label == "default box":
            ctx["nesting"]["default box"] = run
            out.append(check(f"{label}: first and final log-likelihood against the GTR20 evaluation (report)", None,
                             **values))
        else:
            out.append(check(f"{label}: first log-likelihood equals the GTR20 evaluation within 1e-8 |lnL|",
                             first is not None and abs(first - reference) <= line, **values))
            out.append(check(f"{label}: final log-likelihood not below the GTR20 evaluation by more than 1e-8 |lnL|",
                             final is not None and final >= reference - line))
        out += monitoring(label, run, ck)
    return out


def skewed_target(c):
    """Decision 023's reported case: the test target cubed and renormalized, at 17 digits."""
    p = c.pi ** 3
    return ",".join(f"{v:.17g}" for v in p / p.sum())


def case_skewed(r, c):
    text = skewed_target(c)
    pi = target.parse_target(text)[0]
    z = chart.inverse(*nqc.lg_seed(pi))
    run = r.run("skewed", level1_args(f"NQC+F{{{text}}}")[:-1] + ["-v"], {"input.tree": c.t02})
    ck = single_checkpoint(run) if run.code == 0 else None
    fits = FIT.findall(run.stdout())
    res = floats(ck["residuals"]).tolist() if ck and "residuals" in ck else None
    return [check("target cubed and renormalized, aa_example Level 1, default box (report; decision 023)", None,
                  exit_code=run.code, target_min=float(pi.min()), start_outside_default_box=int(
                      np.sum((z < nqc.LOWER) | (z > nqc.UPPER))), first_start=fits[0][1] if fits else None,
                  accepted_score=(ck or {}).get("accepted_score"), bounds=(ck or {}).get("bounds"), residuals=res,
                  fits=(ck or {}).get("fits"), falls=FALL.findall(run.stdout()), error=nqc_error(run))
            ] + monitoring("skewed", run, ck)


def staged(rundir, inputs=("input.tree",)):
    """A finished run's inputs and checkpoint, to rerun in a new directory."""
    return {**{name: read(rundir / name) for name in inputs},
            f"{PREFIX}.ckp.gz": (rundir / f"{PREFIX}.ckp.gz").read_bytes()}


def case_restart(r, c, ctx, older):
    """CODE_PLAN.md section 4, S3a: NQC's restore refuses a checkpoint written with another target or
    domain, accepts the same command (exit 0 is all the keys can show, since IQ-TREE re-saves no
    model after this restore), and reads S2's checkpoint, which lacks S3a's keys."""
    base = level1(r, c, ctx)
    if base.code != 0:
        return failed_run(base)
    args = level1_args(c.model)
    out = []
    for name, a in (("target", level1_args(f"NQC+F{{{skewed_target(c)}}}")),
                    ("domain", args + ["--nqc-domain", "1e-6,100"])):
        run = r.run(f"restart-other-{name}", a + ["--undo"], staged(base.dir))
        message = nqc_error(run)
        # refused by the restore, not by an earlier NQC check
        out.append(check(f"--undo with another {name}: refused by NQC's restore",
                         run.code not in (0, "timeout") and message is not None and RESTORE_REFUSAL in message,
                         exit_code=run.code, message=message))
    run = r.run("restart-undo", args + ["--undo"], staged(base.dir))
    out.append(check("--undo with the same command: exit 0", run.code == 0, exit_code=run.code))
    out.append(check("--undo: summary against the original's (report)", None,
                     equal=summary(run.stdout()) == summary(base.stdout()), undo=summary(run.stdout())))
    run = r.run("restart-older", args + ["--undo"], staged(older))
    out.append(check("S2's checkpoint, --undo with this binary (report)", None, exit_code=run.code,
                     summary=summary(run.stdout()), error=nqc_error(run)))
    return out


def fit_state(run):
    ck = final_struct(run) if run.code == 0 else None
    printed = re.search(SCALARS["log_likelihood"], run.report(), re.M)
    return {"exit_code": run.code, "accepted_score": (ck or {}).get("accepted_score"),
            "printed_final": printed.group(1) if printed else None, "bounds": (ck or {}).get("bounds"),
            "fits": (ck or {}).get("fits"), "falls": len(FALL.findall(run.stdout()))}


def case_risk22(r, c, ctx):
    """PLAN.md risk 22's evidence on aa_example, reported: the default box against 1e-8,1e7."""
    tree = {"input.tree": c.t02}
    if "Level 2 +G4" not in ctx:
        ctx["Level 2 +G4"] = r.run("fit-level-2-g4", level2_args(c.model, "+G4"), tree)
    out = []
    for label, default, args in (("Level 1", level1(r, c, ctx), level1_args(c.model)),
                                 ("Level 2 +G4", ctx["Level 2 +G4"], level2_args(c.model, "+G4"))):
        wide = r.run("risk22-" + slug(label), args + ["--nqc-domain", WIDE], tree)
        out.append(check(f"{label}: default box against --nqc-domain {WIDE} (report)", None,
                         default=fit_state(default), widened=fit_state(wide)))
        out += monitoring(f"{label} at {WIDE}", wide, single_checkpoint(wide) if wide.code == 0 else None)
    return out


def case_toy(r, c, calibration_path, alisim):
    cal = json.loads(read(calibration_path))
    data = r.out / "toy-data"
    data.mkdir()
    made = toy.make(data, alisim, c.pi)
    out = [check("toy data identical to the calibration's", made["sha256"] == cal["provenance"]["toy"]["sha256"]),
           check("toy target identical to the calibration's", cal["target_text"] == c.target_text)]
    run = r.run("toy", ["-s", data / "toy.phy", "-m", c.model, "-te", data / "toy.nwk", "-blfix", "-nparam", "3", "-vv"])
    if run.code != 0:
        return out + failed_run(run)
    plain = cal["fits"]["plain"]
    m = INITIAL.search(run.stdout())
    iq_start = float(m.group(1)) if m else None
    out.append(check("start score against the oracle's (decision 021)",
                     iq_start is not None and abs(iq_start - plain["start_lnl"]) <= TOL_LNL * abs(plain["start_lnl"]),
                     iqtree=iq_start, oracle=plain["start_lnl"]))
    objective = nqc.Objective(readers.read_phylip(data / "toy.phy"), readers.read_newick(data / "toy.nwk"), c.pi, c.ref)
    got = first_gradient(run.stdout())
    if got is None:
        out.append(check("first gradient printed", False))
    else:
        f_iq, x, h, g = got
        f_or, g_or = forward(objective, x, h)
        ratio = np.abs(g - g_or) / same_rule_bound(-f_or, h)
        out.append(check("first gradient, same rule, |g - g_oracle| / line", bool(np.all(ratio <= 1)),
                         worst=float(ratio.max())))
    ck = single_checkpoint(run)
    if not (ck and "accepted_score" in ck):
        return out + [check("final score against the oracle's port (decision 030's calibrated line)", False,
                            reason="no single NQC struct with accepted_score")]
    accepted = float(ck["accepted_score"])
    out.append(check("final score against the oracle's port (decision 030's calibrated line)",
                     abs(accepted - plain["final_lnl"]) <= cal["line"], iqtree=accepted, oracle=plain["final_lnl"],
                     difference=abs(accepted - plain["final_lnl"]), line=cal["line"]))
    z = floats(ck["coordinates"])
    Q_iq = chart.build(z, c.pi, c.ref)[0]
    out.append(check("final Q against the oracle's port, largest relative difference (report)", None,
                     max_relative=float(np.max(np.abs(Q_iq - np.array(plain["Q"])) / np.abs(np.array(plain["Q"])))),
                     max_coordinate_difference=float(np.max(np.abs(z - np.array(plain["x"])))),
                     oracle_iterations=plain["iterations"], oracle_stop=plain["stop"],
                     oracle_evaluations=plain["evaluations"]))
    return out


def case_short(r, c, alisim):
    data = r.out / "short-data"
    data.mkdir()
    toy.make(data, alisim, c.pi)
    aln = readers.read_phylip(data / "toy.phy")
    write(data / "short.phy", f"{len(aln)} 300\n" + "".join(f"{n} {s[:300]}\n" for n, s in aln.items()))
    out = []
    for label, suffix, extra in (("Level 1", "", ["-blfix"]), ("Level 2 +I+G4", "+I+G4", [])):
        run = r.run("short-" + label.lower().replace(" ", "-").replace("+", ""),
                    ["-s", data / "short.phy", "-m", c.model + suffix, "-te", data / "toy.nwk", "-nparam", "3", "-v",
                     *extra])
        if run.code != 0:
            out += failed_run(run)
            continue
        if suffix:
            out.append(restore_check(label, run.stdout()))
        ck = single_checkpoint(run)
        out += key_checks(label, ck, c.pi) + monitoring(label, run, ck)
    return out


# ---------------------------------------------------------------- linked cases (S3b)

# every linked fit dumped and numbered: otherwise a checkpoint is dumped at most once a minute
# (CODE_PLAN.md section 3.5)
LINKED_CKP = ["-cptime", "0", "--all-checkpoint"]
SHARE_REFUSAL = "linked partitions must share one model"
CPU = re.compile(r"^Total CPU time used: (\S+) sec", re.M)
WALL = re.compile(r"^Total wall-clock time used: (\S+) sec", re.M)
SEED_RUNS = 10          # decision 038: runs of each multi-threaded case at the seed
LEVEL2_THREADED = 5     # CODE_PLAN.md section 3.5: Level 2 fits at -T 4, reported
TURTLE_SEED = {"-p": ("run08", 4), "-S": ("run09", 3)}   # decision 038's trees and thread counts
# CODE_PLAN.md section 3.5: 360, plus rate dimensions and branches per partition under -S, plus the
# super tree's branchNum and P - 1 partition rates under -p (none under -q)
TWOPART_COUNTS = {("-S", "seed"): 360, ("-p", "seed"): 384, ("-S", "Level 1"): 360, ("-q", "Level 1"): 383,
                  ("-S", "+G4"): 406, ("-S", "+I+G4"): 408, ("-S", "+R4"): 416,
                  ("-p", "+G4"): 386, ("-p", "+I+G4"): 388, ("-p", "+R4"): 396}
LINKED_CASES = ("linked-seed", "linked-fits", "linked-refusals", "permuted", "linked-restart", "linked-short")


def numbered_checkpoints(run):
    """The paths of the checkpoints --all-checkpoint numbered, in order (utils/checkpoint.cpp)."""
    found = []
    for p in run.dir.glob(f"{PREFIX}.*.ckp.gz"):
        n = p.name[len(PREFIX) + 1:-len(".ckp.gz")]
        if n.isdigit():
            found.append((int(n), p))
    return [p for _, p in sorted(found)]


def part_alignment(aln, sites):
    return {n: "".join(s[i] for i in sites) for n, s in aln.items()}


def tree_name(mode):
    return "partition.trees" if mode == "-S" else "input.tree"


def one_round(mode):
    """One round of model optimization: -S counts rounds with --loop-model, -p and -q with -nparam
    (CODE_PLAN.md section 5)."""
    return ["--loop-model", "1"] if mode == "-S" else ["-nparam", "2"]


def joint_args(aln, mode, nex, model, rate, trees, *extra):
    """A --model-joint run. A rate model comes in through -m, whose substitution model the joint
    model replaces (model/modelfactory.cpp, model_joint)."""
    return (["-s", aln, mode, nex, "--model-joint", model, "-te", trees] + (["-m", "LG" + rate] if rate else [])
            + list(extra))


def qmat_text(r, c, name, z):
    """The Q built from z at pi*, as the 17-digit -m FILE that --init-model reads."""
    d = r.out / f"{name}-input"
    d.mkdir()
    write_matrix(d / f"{name}.qmat", chart.build(z, c.pi, c.ref)[0], c.pi)
    return read(d / f"{name}.qmat")


def reversed_charsets(text):
    """A NEXUS sets block with its charset lines in reverse order."""
    lines = text.splitlines()
    at = [i for i, line in enumerate(lines) if line.strip().lower().startswith("charset")]
    for i, line in zip(at, [lines[i] for i in reversed(at)]):
        lines[i] = line
    return "\n".join(lines) + "\n"


class TwoPart:
    """nqc/twopart.py's data, made once per run of this driver."""

    def __init__(self, out, alisim, c):
        d = out / "twopart-data"
        d.mkdir()
        self.made = twopart.make(d, alisim, c.pi)
        self.phy, self.nex = d / "twopart.phy", d / "twopart.nex"
        self.tree = read(d / "twopart.nwk")
        self.aln = readers.read_phylip(self.phy)
        self.charsets = readers.read_charsets(self.nex)

    def files(self, mode):
        return {tree_name(mode): self.tree * (2 if mode == "-S" else 1)}

    def objective(self, c, shapes=None, free=None, rates=None, tree=None):
        """The oracle's linked objective, partitions in file order. IQ-TREE sums in part_order, which
        changes the sum by rounding only."""
        parsed = readers.parse_newick(tree or self.tree)
        parts = []
        for k, sites in enumerate(self.charsets.values()):
            shape = shapes[k] if isinstance(shapes, (list, tuple)) else shapes
            parts.append(nqc.Objective(part_alignment(self.aln, sites), parsed, c.pi, c.ref, shape=shape, free=free,
                                       part_rate=rates[k] if rates else 1.0))
        return nqc.LinkedObjective(parts)


def twopart_data(r, c, ctx, alisim):
    if "twopart" not in ctx:
        ctx["twopart"] = TwoPart(r.out, alisim, c)
    return ctx["twopart"]


def turtle_files(mode, reverse=False):
    lines = [l for l in tree_text(baseline_runs(), TURTLE_SEED[mode][0]).splitlines() if l.strip()]
    return {tree_name(mode): "\n".join(reversed(lines) if reverse else lines) + "\n"}


def gamma_shapes(ck, names):
    """Each partition's gamma shape as its checkpoint struct holds it (10 significant digits)."""
    return [next((v for k, v in ck.items() if k.endswith(f"!{n}!RateGamma!gamma_shape")), None) for n in names]


def q_route(run, names):
    """Decision 024's -q inputs from a -p fit: its tree file, and each partition's gamma shape and
    partition rate as its checkpoint holds them; None if the run left any of them out."""
    ck = run.checkpoint() if run.code == 0 else {}
    rates, shapes, trees = ck.get("PartitionModelPlen!part_rates"), gamma_shapes(ck, names), run.trees()
    if rates is None or None in shapes or not trees:
        return None
    return {"tree": trees[0] + "\n", "shapes": shapes, "rates": rates.split(", ")}


def q_charpartition(nex_text, names, route):
    return with_charpartition(nex_text, [(f"LG+G4{{{a}}}", n, x) for n, a, x in zip(names, route["shapes"],
                                                                                    route["rates"])])


def linked_invariants(label, run, c, partitions):
    """Asserted at every numbered checkpoint, or at the final one when the run wrote none
    (CODE_PLAN.md section 4, S3b): pi* bitwise in the linked struct and in every partition struct,
    coordinates bitwise equal across them, the oracle's references, and the residual keys within
    decision 021's line. A checkpoint holding no NQC struct is counted, not judged."""
    paths = numbered_checkpoints(run)
    cks = [read_checkpoint(p) for p in paths] or [run.checkpoint()]
    bad = {k: [] for k in ("structs", "target", "coordinates", "references", "residuals")}
    judged, without, fits, worst = 0, 0, [], 0.0
    ref = [int(v) for v in c.ref]
    for n, flat in enumerate(cks, 1):
        structs = nqc_structs(flat)
        if not structs:
            without += 1
            continue
        judged += 1
        linked = linked_struct(structs)
        if (linked is None or len(structs) != partitions + 1 or
                any(k not in s for s in structs.values() for k in ("target", "coordinates", "references", "residuals"))):
            bad["structs"].append(n)
            continue
        fits.append(linked.get("fits"))
        if not all(np.array_equal(floats(s["target"]), c.pi) for s in structs.values()):
            bad["target"].append(n)
        if len({s["coordinates"] for s in structs.values()}) != 1:   # 17 digits identify a double
            bad["coordinates"].append(n)
        if any([int(v) for v in s["references"].split(",")] != ref for s in structs.values()):
            bad["references"].append(n)
        res = max(float(np.max(floats(s["residuals"])[:3])) for s in structs.values())
        worst = max(worst, res)
        if res > RESIDUAL:
            bad["residuals"].append(n)
    where = f"{judged} numbered checkpoints" if paths else "the final checkpoint"
    ok = judged > 0 and not bad["structs"]
    return [
        check(f"{label}: the linked struct and {partitions} partition structs in {where}", ok,
              failing=bad["structs"], checkpoints_without_nqc=without),
        check(f"{label}: target bitwise equal to the parsed target in every struct", ok and not bad["target"],
              failing=bad["target"]),
        check(f"{label}: coordinates bitwise equal across the structs", ok and not bad["coordinates"],
              failing=bad["coordinates"]),
        check(f"{label}: references equal the oracle's in every struct", ok and not bad["references"],
              failing=bad["references"]),
        check(f"{label}: residual keys at most {RESIDUAL:g} in every struct (decision 021)",
              ok and not bad["residuals"], failing=bad["residuals"], worst=worst),
        check(f"{label}: linked fits key at each checkpoint (report)", None, fits=fits),
    ]


def linked_fit_checks(label, run, c, partitions, expected_parameters):
    """A linked fit's checks, with block_restore_check for -S +I+G runs."""
    if run.code != 0:
        return failed_run(run) + monitoring(label, run, None)
    final = linked_struct(nqc_checkpoint(run))
    out = linked_invariants(label, run, c, partitions)
    out.append(key_checks(f"{label}, final linked struct", final, c.pi)[-1])
    out.append(free_parameters_check(label, run, expected_parameters))
    if label.startswith("-S ") and "+I+G" in label:
        out.append(block_restore_check(label, run.stdout()))
    return out + monitoring(label, run, final)


def oracle_score_check(label, run, objective):
    """Decision 021's Level 1 line: the accepted score against the oracle's linked objective at the
    checkpointed coordinates."""
    name = f"{label}: accepted score against the oracle's linked re-evaluation (decisions 021, 030)"
    ck = linked_struct(nqc_checkpoint(run)) if run.code == 0 else None
    if not (ck and "accepted_score" in ck):
        return [check(name, False, reason="no linked accepted_score key")]
    oracle = -objective(floats(ck["coordinates"]))
    accepted = float(ck["accepted_score"])
    return [check(name, abs(accepted - oracle) <= TOL_LNL * abs(oracle), accepted=accepted, oracle=oracle)]


def linked_gradient_checks(label, run, objective):
    """Decision 021's lines at the printed point; the printed steps against decision 005's rule,
    bitwise, with the upper bound from the domain key; and the count of coordinates at which IQ-TREE's
    legacy step would differ (0 < |x| < 1), which must be at least 1 for the bitwise check to tell
    the two rules apart."""
    out, x, *_ = gradient_checks(label, run, objective)
    ck = linked_struct(nqc_checkpoint(run))
    if x is None or ck is None:
        return out + [check(f"{label}: first gradient and linked struct present", False)]
    h = first_gradient(run.stdout())[2]
    rule = nqc.scaled_steps(x, floats(ck["domain"])[1])
    out.append(check(f"{label}: printed steps equal decision 005's rule, bitwise", bool(np.array_equal(h, rule)),
                     differing=int(np.sum(h != rule))))
    differ = int(np.sum((np.abs(x) < 1) & (x != 0)))
    out.append(check(f"{label}: coordinates where the legacy step differs, at least 1", differ >= 1, count=differ))
    return out


def seed_args(tp, c, mode, rate):
    return joint_args(tp.phy, mode, tp.nex, c.model, rate, tree_name(mode), "-blfix", "-vv", *one_round(mode))


def case_linked_seed(r, c, ctx, alisim):
    """CODE_PLAN.md section 4, S3b, "Linked seed" and "Linked gradient": nothing moves before the
    first linked gradient, so its -lnL is the oracle's linked objective at the LG seed."""
    tp = twopart_data(r, c, ctx, alisim)
    out = []
    for mode in ("-S", "-p"):
        for name, rate, shape, free in (("+G4", f"+G4{{{ALPHA}}}", ALPHA, None), ("+R4", R4, None, R4_GIVEN)):
            label = f"{mode} {name} seed"
            run = r.run("linked-seed" + slug(f"{mode} {name}"), seed_args(tp, c, mode, rate), tp.files(mode))
            if run.code != 0:
                out += failed_run(run)
                continue
            out += linked_gradient_checks(label, run, tp.objective(c, shape, free))
            out.append(free_parameters_check(label, run, TWOPART_COUNTS[(mode, "seed")]))
    z = c.z.copy()
    z[::18] = NEAR_ZERO   # the first free coordinate of each row, as in the gradient case
    files = {**tp.files("-S"), "nearzero.qmat": qmat_text(r, c, "linked-nearzero", z)}
    run = r.run("linked-gradient-nearzero", seed_args(tp, c, "-S", f"+G4{{{ALPHA}}}") +
                ["--init-model", "nearzero.qmat"], files)
    if run.code != 0:
        return out + failed_run(run)
    out += linked_gradient_checks("-S near zero", run, tp.objective(c, ALPHA))
    got = first_gradient(run.stdout())
    near = int(np.sum(np.abs(got[1]) < 1e-8)) if got else None
    # as in the gradient case: shows --init-model reached the linked objects
    out.append(check("-S near zero: coordinates within 1e-8 of zero", near == 20, count=near))
    return out


def linked_level1(r, c, ctx, tp):
    """The -S Level 1 fit that the fits, permuted and restart cases share."""
    if "linked-level1" not in ctx:
        ctx["linked-level1"] = r.run("linked-fit-s-level1", joint_args(
            tp.phy, "-S", tp.nex, c.model, f"+G4{{{ALPHA}}}", tree_name("-S"), "-blfix", "-vv", *LINKED_CKP),
            tp.files("-S"))
    return ctx["linked-level1"]


def linked_level2(r, c, ctx, tp, mode, name):
    """A linked Level 2 fit. The +G4 runs print the first gradient, which the restart case reads
    from the same command."""
    key = ("linked-level2", mode, name)
    if key not in ctx:
        ctx[key] = r.run("linked-fit" + slug(f"{mode} level 2 {name}"), joint_args(
            tp.phy, mode, tp.nex, c.model, name, tree_name(mode), "-vv" if name == "+G4" else "-v", *LINKED_CKP),
            tp.files(mode))
    return ctx[key]


def case_linked_fits(r, c, ctx, alisim):
    tp = twopart_data(r, c, ctx, alisim)
    names = list(tp.charsets)
    out = []
    run = linked_level1(r, c, ctx, tp)
    out += linked_fit_checks("-S Level 1", run, c, 2, TWOPART_COUNTS[("-S", "Level 1")])
    if run.code == 0:
        out += oracle_score_check("-S Level 1", run, tp.objective(c, ALPHA))
    for mode in ("-S", "-p"):
        for name in ("+G4", "+I+G4", "+R4"):
            out += linked_fit_checks(f"{mode} Level 2 {name}", linked_level2(r, c, ctx, tp, mode, name), c, 2,
                                     TWOPART_COUNTS[(mode, name)])
    route = q_route(linked_level2(r, c, ctx, tp, "-p", "+G4"), names)
    if route is None:
        out.append(check("-q Level 1: the -p Level 2 +G4 fit's tree, shapes and partition rates", False))
    else:
        run = r.run("linked-fit-q-level1", ["-s", tp.phy, "-q", "parts.nex", "--model-joint", c.model, "-te",
                                            "input.tree", "-blfix", "-v", *LINKED_CKP],
                    {"parts.nex": q_charpartition(read(tp.nex), names, route), "input.tree": route["tree"]})
        out += linked_fit_checks("-q Level 1", run, c, 2, TWOPART_COUNTS[("-q", "Level 1")])
        if run.code == 0:
            rates = nqc.normalize_part_rates([float(x) for x in route["rates"]], [len(s) for s in tp.charsets.values()])
            out += oracle_score_check("-q Level 1", run, tp.objective(c, [float(a) for a in route["shapes"]],
                                                                      rates=rates, tree=route["tree"]))
            out.append(check("-q Level 1: inputs from the -p Level 2 +G4 fit (report)", None, shapes=route["shapes"],
                             rates=route["rates"], normalized=rates))
    parts = with_charpartition(read(tp.nex), [(c.model + f"+G4{{{ALPHA}}}", n, None) for n in names])
    run = r.run("linked-fit-link-model", ["-s", tp.phy, "-S", "parts.nex", "--link-model", "-te", tree_name("-S"),
                                          "-blfix", "-vv", *LINKED_CKP], {"parts.nex": parts, **tp.files("-S")})
    out += linked_fit_checks("--link-model -S Level 1", run, c, 2, TWOPART_COUNTS[("-S", "Level 1")])
    if run.code == 0:
        out += oracle_score_check("--link-model -S Level 1", run, tp.objective(c, ALPHA))
    a, b = final_struct(linked_level1(r, c, ctx, tp)) or {}, final_struct(run) or {}
    out.append(check("--link-model against --model-joint, -S Level 1 (report)", None,
                     coordinates_equal="coordinates" in a and a.get("coordinates") == b.get("coordinates"),
                     accepted=[a.get("accepted_score"), b.get("accepted_score")], fits=[a.get("fits"), b.get("fits")]))
    return out


def case_linked_refusals(r, c, ctx, alisim):
    tp = twopart_data(r, c, ctx, alisim)
    nex = read(tp.nex)
    cubed = ",".join(f"{v:.17g}" for v in twopart.cubed(c.pi))
    eight = ",".join(f"{v:.8f}" for v in c.pi)
    common = ["-s", tp.phy, "-S", "parts.nex", "--link-model", "-te", tree_name("-S"), "-blfix"]
    out = []
    for name, other in (("another target", cubed), ("the target written to 8 decimals", eight)):
        parts = with_charpartition(nex, [(c.model + f"+G4{{{ALPHA}}}", "partA", None),
                                         (f"NQC+F{{{other}}}+G4{{{ALPHA}}}", "partB", None)])
        run = r.run("linked-compat-" + ("other" if other is cubed else "eight"), common + ["--loop-model", "1"],
                    {"parts.nex": parts, **tp.files("-S")})
        message = nqc_error(run)
        # refused by startLinkedFit, not by an earlier NQC check
        out.append(check(f"--link-model with {name}: refused at the first linked fit",
                         run.code not in (0, "timeout") and message is not None and SHARE_REFUSAL in message,
                         exit_code=run.code, message=message))
    # each refusal must be the one meant: ModelFinder's, and the partition string's +H
    for name, why, entries in (("a charset without a model and no -m", "ModelFinder", [(c.model, "partA"), ("", "partB")]),
                               ("a charset NQC+F{...}+H4", "heterotachy", [(c.model + "+H4", "partA"),
                                                                          (c.model + f"+G4{{{ALPHA}}}", "partB")])):
        parts = with_charpartition(nex, [(m, n, None) for m, n in entries])
        run = r.run("linked-refuse-" + ("modelfinder" if why == "ModelFinder" else "ghost"), common,
                    {"parts.nex": parts, **tp.files("-S")})
        message = nqc_error(run)
        out.append(check(f"refused: --link-model with {name}, for {why}",
                         run.code not in (0, "timeout") and message is not None and why in message,
                         exit_code=run.code, message=message))
    return out


def turtle_seed_args(c, mode, nex=NEX):
    return joint_args(TURTLE, mode, nex, c.model, f"+G4{{{ALPHA}}}", tree_name(mode), "-blfix", "-vv",
                      *one_round(mode))


def turtle_seed(r, c, ctx, mode):
    """Turtle at the seed at -T 1, the reference of the threads and permuted cases."""
    key = ("turtle-seed", mode)
    if key not in ctx:
        ctx[key] = r.run("turtle-seed" + slug(mode), turtle_seed_args(c, mode), turtle_files(mode))
    return ctx[key]


def seed_against(label, reference, runs):
    """Decision 038's lines: -lnL within 1e-8 |lnL| and each gradient entry within 2e-12 |lnL| / |h_k|
    of the reference run's, at the same printed point."""
    ref = first_gradient(reference.stdout()) if reference.code == 0 else None
    if ref is None:
        return [check(f"{label}: reference run printed its first gradient", False, exit_code=reference.code)]
    f1, x1, h1, g1 = ref
    got = [first_gradient(run.stdout()) if run.code == 0 else None for run in runs]
    same_point = [g is not None and np.array_equal(g[1], x1) and np.array_equal(g[2], h1) for g in got]
    usable = [g for g, ok in zip(got, same_point) if ok]
    complete = len(usable) == len(runs)
    lnl_rel = [abs(g[0] - f1) / abs(f1) for g in usable]
    ratio = [float(np.max(np.abs(g[3] - g1) / same_rule_bound(f1, h1))) for g in usable]
    return [
        check(f"{label}: -lnL within 1e-8 |lnL| of the reference (decision 038)",
              complete and max(lnl_rel) <= TOL_LNL, runs=len(runs), at_the_same_point=len(usable),
              largest_relative=max(lnl_rel, default=None)),
        check(f"{label}: first gradient within 2e-12 |lnL| / |h| of the reference (decision 038)",
              complete and max(ratio) <= 1, worst=max(ratio, default=None)),
        check(f"{label}: runs bitwise equal to the reference (report)", None,
              lnl=sum(g[0] == f1 for g in usable), gradient=sum(bool(np.array_equal(g[3], g1)) for g in usable)),
    ]


def cpu_wall(run):
    cpu, wall = CPU.search(run.report()), WALL.search(run.report())
    return round(float(cpu.group(1)) / float(wall.group(1)), 2) if cpu and wall and float(wall.group(1)) > 0 else None


def case_permuted(r, c, ctx, alisim):
    """CODE_PLAN.md section 3.5, "Permuted partition order"."""
    tp = twopart_data(r, c, ctx, alisim)
    base = linked_level1(r, c, ctx, tp)
    run = r.run("permuted-twopart", joint_args(tp.phy, "-S", "reversed.nex", c.model, f"+G4{{{ALPHA}}}",
                                               tree_name("-S"), "-blfix", "-vv", *LINKED_CKP),
                {"reversed.nex": twopart.nexus(tuple(reversed(list(tp.charsets)))), **tp.files("-S")})
    a = first_gradient(base.stdout()) if base.code == 0 else None
    b = first_gradient(run.stdout()) if run.code == 0 else None
    za, zb = ((final_struct(x) or {}) if x.code == 0 else {} for x in (base, run))
    out = [check("twopart -S Level 1, charsets reversed: seed -lnL bitwise equal", bool(a and b and a[0] == b[0]),
                 exit_codes=[base.code, run.code], lnl=[a and a[0], b and b[0]]),
           check("twopart -S Level 1, charsets reversed: first gradient (point, steps, values) bitwise equal",
                 bool(a and b and all(np.array_equal(u, v) for u, v in zip(a[1:], b[1:])))),
           check("twopart -S Level 1, charsets reversed: final coordinates bitwise equal",
                 "coordinates" in za and za.get("coordinates") == zb.get("coordinates"))]
    for mode in ("-p", "-S"):
        ref = turtle_seed(r, c, ctx, mode)
        rev = r.run("permuted-turtle" + slug(mode), turtle_seed_args(c, mode, "reversed.nex"),
                    {"reversed.nex": reversed_charsets(read(NEX)), **turtle_files(mode, reverse=True)})
        out += seed_against(f"turtle {mode} at the seed, charsets reversed", ref, [rev])
    return out


def case_threads(r, c, ctx):
    """Decision 038 on turtle: the linked seed at more threads against -T 1, never bitwise; and
    Level 2 fits at -T 4, reported. Under -p IQ-TREE caps turtle at 3 threads, one per partition,
    and says so; nothing logs the allocation under -S, so the CPU-to-wall ratio is reported."""
    out = []
    for mode in ("-p", "-S"):
        threads = TURTLE_SEED[mode][1]
        ref = turtle_seed(r, c, ctx, mode)
        runs = [r.run(f"threads{slug(mode)}-t{threads}-{k}", turtle_seed_args(c, mode), turtle_files(mode), threads)
                for k in range(1, SEED_RUNS + 1)]
        out += seed_against(f"turtle {mode} at the seed, -T {threads} against -T 1", ref, runs)
        note = [l for l in runs[0].stdout().splitlines() if "thread" in l.lower()][:4]
        out.append(check(f"turtle {mode} -T {threads}: CPU-to-wall ratios and IQ-TREE's thread lines (report)", None,
                         ratios=[cpu_wall(run) for run in runs], reference_ratio=cpu_wall(ref), lines=note))
    args = joint_args(TURTLE, "-p", NEX, c.model, None, tree_name("-p"), "-v")
    fits = [r.run(f"threads-p-level2-t4-{k}", args, turtle_files("-p"), 4) for k in range(1, LEVEL2_THREADED + 1)]
    one = r.run("threads-p-level2-t1", args, turtle_files("-p"))
    states = [fit_state(run) for run in fits]
    finals = [float(s["printed_final"]) for s in states if s["printed_final"]]
    out.append(check(f"turtle -p Level 2, no rate variation, {LEVEL2_THREADED} fits at -T 4 beside -T 1 "
                     "(report; decision 038, NONREV's pL2 spread)", None, t4=states, t1=fit_state(one),
                     t4_spread=max(finals) - min(finals) if finals else None,
                     hard_stops=[("ASSERT" in run.stdout() or BUG_ABORT in run.stdout()) for run in fits + [one]],
                     switches=[run.stdout().count(SWITCH) for run in fits + [one]]))
    return out


def case_linked_restart(r, c, ctx, alisim):
    """CODE_PLAN.md section 3.5, "Restart", for linked runs: the -p and -S Level 2 +G4 runs resumed
    from numbered checkpoint #1, written after the first linked fit (under -p before any branch
    update; under -S after round 1's per-partition updates); no likelihood comparison, because the
    resumed round refits the rates before its first linked fit. The counter must continue: the
    first checkpoint whose coordinates differ from #1's, so written after the resumed run's first
    fit, holds #1's count plus 1, where a reset counter would hold 1. Then --undo of the finished -S
    Level 1 run."""
    tp = twopart_data(r, c, ctx, alisim)
    out = []
    for mode in ("-p", "-S"):
        base = linked_level2(r, c, ctx, tp, mode, "+G4")
        label = f"{mode} Level 2 +G4"
        paths = numbered_checkpoints(base) if base.code == 0 else []
        ck1 = linked_struct(nqc_structs(read_checkpoint(paths[0]))) if paths else None
        ok = bool(ck1) and paths[0].name == f"{PREFIX}.1.ckp.gz" and ck1.get("fits") == "1"
        out.append(check(f"{label}: checkpoint #1 holds the first linked fit", ok, exit_code=base.code,
                         numbered=len(paths), fits=(ck1 or {}).get("fits")))
        if not ok:   # the record's precondition; stop and ask rather than pick another checkpoint
            continue
        rerun = r.run(f"restart{slug(mode)}-level2", base.args, {**tp.files(mode),
                                                                f"{PREFIX}.ckp.gz": paths[0].read_bytes()})
        got = first_gradient(rerun.stdout()) if rerun.code == 0 else None
        fits, first_new = [], None
        for p in numbered_checkpoints(rerun):
            s = linked_struct(nqc_structs(read_checkpoint(p)))
            if s and "fits" in s:
                fits.append(int(s["fits"]))
                if first_new is None and s.get("coordinates") != ck1["coordinates"]:
                    first_new = int(s["fits"])
        saved = int(ck1["fits"])
        out += [check(f"{label}, resumed from #1: exit 0", rerun.code == 0, exit_code=rerun.code),
                check(f"{label}, resumed: first linked gradient at #1's coordinates, bitwise",
                      got is not None and np.array_equal(got[1], floats(ck1["coordinates"]))),
                check(f"{label}, resumed: fits continue from #1's count, never reset",
                      first_new == saved + 1 and min(fits) >= saved and fits == sorted(fits),
                      fits=fits, first_after_a_fit=first_new, saved=saved)]
        out += linked_invariants(f"{label}, resumed", rerun, c, 2)
    base = linked_level1(r, c, ctx, tp)
    if base.code != 0:
        return out + failed_run(base)
    files = staged(base.dir, (tree_name("-S"),))
    cubed = f"NQC+F{{{skewed_target(c)}}}"
    for name, args in (("target", [cubed if a == c.model else a for a in base.args]),
                       ("domain", base.args + ["--nqc-domain", "1e-6,100"])):
        run = r.run(f"restart-linked-other-{name}", args + ["--undo"], files)
        message = nqc_error(run)
        out.append(check(f"linked --undo with another {name}: refused by NQC's restore",
                         run.code not in (0, "timeout") and message is not None and RESTORE_REFUSAL in message,
                         exit_code=run.code, message=message))
    run = r.run("restart-linked-undo", base.args + ["--undo"], files)
    out.append(check("linked --undo with the same command: exit 0", run.code == 0, exit_code=run.code))
    out.append(check("linked --undo: summary against the original's (report)", None,
                     equal=summary(run.stdout()) == summary(base.stdout()), undo=summary(run.stdout())))
    return out


# PLAN.md risk 22 on turtle (interpretation recorded in CODE_PLAN.md section 3.5 at S3b's close)
TURTLE_FITS = {"turtle Level 1 -S": ("-S", f"+G4{{{ALPHA}}}", ["-blfix"]),
               "turtle Level 2 -p +G4": ("-p", "+G4", []),
               "turtle Level 2 -S +G4": ("-S", "+G4", [])}


def turtle_fit(r, c, ctx, label, wide=False):
    key = ("turtle-fit", label, wide)
    if key not in ctx:
        mode, rate, extra = TURTLE_FITS[label]
        args = joint_args(TURTLE, mode, NEX, c.model, rate, tree_name(mode), *extra, "-v")
        ctx[key] = r.run(("risk22-" if wide else "") + slug(label), args + (["--nqc-domain", WIDE] if wide else []),
                         turtle_files(mode))
    return ctx[key]


def case_risk22_turtle(r, c, ctx):
    out = []
    for label in TURTLE_FITS:
        box, wide = turtle_fit(r, c, ctx, label), turtle_fit(r, c, ctx, label, True)
        out.append(check(f"{label}: default box against --nqc-domain {WIDE} (report)", None,
                         default=fit_state(box), widened=fit_state(wide)))
        out += monitoring(label, box, final_struct(box) if box.code == 0 else None)
        out += monitoring(f"{label} at {WIDE}", wide, final_struct(wide) if wide.code == 0 else None)
    return out


def continuations(r, c, label, box, args, files):
    """Two Level 1 continuations from a default-box fit's final state, its coordinates passed as a
    17-digit Q by --init-model, one in the box and one at WIDE (CODE_PLAN.md section 4, S3b, "Risk 22
    warm starts"). The first lnL is -lnL at the first-gradient print."""
    ck = final_struct(box) if box.code == 0 else None
    if not (ck and "coordinates" in ck) or args is None:
        return [check(f"{label}: continuations (report)", None, reason="no final state to start from",
                      exit_code=box.code)]
    files = {**files, "warm.qmat": qmat_text(r, c, "warm-" + slug(label), floats(ck["coordinates"]))}
    printed = re.search(SCALARS["log_likelihood"], box.report(), re.M)
    values = {"box_accepted_score": ck.get("accepted_score"), "box_printed_final": printed and printed.group(1)}
    for name, extra in (("in the default box", []), (f"at {WIDE}", ["--nqc-domain", WIDE])):
        run = r.run(f"warm-{slug(label)}-{'wide' if extra else 'box'}", args + ["-vv", "--init-model", "warm.qmat",
                                                                               *extra], files)
        got = first_gradient(run.stdout()) if run.code == 0 else None
        end = final_struct(run) if run.code == 0 else {}
        values[name] = {"exit_code": run.code, "first_lnl": -got[0] if got else None,
                        "accepted_score": (end or {}).get("accepted_score"), "bounds": (end or {}).get("bounds"),
                        "fits": (end or {}).get("fits"), "falls": len(FALL.findall(run.stdout())),
                        "switches": run.stdout().count(SWITCH)}
    return [check(f"{label}: Level 1 continuations from the default-box fit (report; PLAN.md risk 22)", None,
                  **values)]


def case_warm_starts(r, c, ctx):
    tree = {"input.tree": c.t02}
    if "Level 2 +G4" not in ctx:
        ctx["Level 2 +G4"] = r.run("fit-level-2-g4", level2_args(c.model, "+G4"), tree)
    if "nesting" not in ctx or "default box" not in ctx["nesting"]:
        case_nesting(r, c, ctx)
    out = continuations(r, c, "aa_example Level 1", level1(r, c, ctx), level1_args(c.model)[:-1], tree)
    box = ctx["Level 2 +G4"]
    alpha = next((v for k, v in box.checkpoint().items() if k.endswith("gamma_shape")), None) if box.code == 0 else None
    args = ["-s", AA, "-m", c.model + f"+G4{{{alpha}}}", "-te", "input.tree", "-blfix"] if alpha else None
    out += continuations(r, c, "aa_example Level 2 +G4", box, args,
                         {"input.tree": box.trees()[0] + "\n" if box.trees() else ""})
    nest = ctx.get("nesting", {})
    if "default box" in nest:
        out += continuations(r, c, "aa_example nesting fit", nest["default box"],
                             ["-s", AA, "-m", c.model + nest["g4"], "-te", "gtr.treefile", "-blfix"], nest["files"])
    else:
        out.append(check("aa_example nesting fit: continuations (report)", None, reason="no nesting fit"))
    names = list(readers.read_charsets(NEX))
    out += continuations(r, c, "turtle Level 1 -S", turtle_fit(r, c, ctx, "turtle Level 1 -S"),
                         joint_args(TURTLE, "-S", NEX, c.model, f"+G4{{{ALPHA}}}", tree_name("-S"), "-blfix"),
                         turtle_files("-S"))
    box = turtle_fit(r, c, ctx, "turtle Level 2 -S +G4")
    shapes = gamma_shapes(box.checkpoint(), names) if box.code == 0 else [None]
    ok = None not in shapes and len(box.trees()) == len(names)
    out += continuations(r, c, "turtle Level 2 -S +G4", box,
                         ["-s", TURTLE, "-S", "parts.nex", "--model-joint", c.model, "-te", tree_name("-S"), "-blfix"]
                         if ok else None,
                         {"parts.nex": with_charpartition(read(NEX), [(f"LG+G4{{{a}}}", n, None)
                                                                      for n, a in zip(names, shapes)]),
                          tree_name("-S"): "\n".join(box.trees()) + "\n"} if ok else {})
    box = turtle_fit(r, c, ctx, "turtle Level 2 -p +G4")
    route = q_route(box, names)
    out += continuations(r, c, "turtle Level 2 -p +G4 (-q route)", box,
                         ["-s", TURTLE, "-q", "parts.nex", "--model-joint", c.model, "-te", "input.tree", "-blfix"]
                         if route else None,
                         {"parts.nex": q_charpartition(read(NEX), names, route), "input.tree": route["tree"]}
                         if route else {})
    return out


def case_linked_short(r, c, ctx, alisim):
    """Cut-down twopart runs for the sanitizer build (CODE_PLAN.md section 3.6): a Level 1 -S fit and
    a Level 2 -p +G4 fit on the first 150 sites of each partition."""
    tp = twopart_data(r, c, ctx, alisim)
    n = 150
    d = r.out / "linked-short-data"
    d.mkdir()
    cut = {name: s[:n] + s[twopart.SITES:twopart.SITES + n] for name, s in tp.aln.items()}
    write(d / "short.phy", f"{len(cut)} {2 * n}\n" + "".join(f"{k} {s}\n" for k, s in cut.items()))
    write(d / "short.nex", twopart.nexus(sites=n))
    out = [check("cut-down data: missing or identical sequences (report)", None,
                 problems=twopart.problems({k: s[:n] for k, s in cut.items()}) +
                 twopart.problems({k: s[n:] for k, s in cut.items()}))]
    # The -S run has no -v: at -v and above, IQ-TREE ends a partitioned run by printing NNI counts
    # through a cast to PhyloSuperTreePlen, which an -S tree is not, and AddressSanitizer stops it
    # there (an upstream finding, PLAN.md maintainer question 6). Its asserted checks read the
    # checkpoint, so only reported print lines are lost.
    for label, mode, rate, extra in (("Level 1 -S", "-S", f"+G4{{{ALPHA}}}", ["-blfix", "--loop-model", "1"]),
                                     ("Level 2 -p +G4", "-p", "+G4", ["-nparam", "3", "-v"])):
        run = r.run("linked-short" + slug(f"{mode} {label}"),
                    joint_args(d / "short.phy", mode, d / "short.nex", c.model, rate, tree_name(mode), *extra),
                    tp.files(mode))
        if run.code != 0:
            out += failed_run(run)
            continue
        out += linked_invariants(label, run, c, 2) + monitoring(label, run, linked_struct(nqc_checkpoint(run)))
    return out


# ---------------------------------------------------------------- driver

def provenance(a, binary, out):
    import scipy
    return {
        "script": SCRIPT, "command": sys.argv,
        "started_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
        "git_commit": capture(a.git, "rev-parse", "HEAD", cwd=REPO),
        "git_clean": capture(a.git, "status", "--porcelain", cwd=REPO) == "",
        "binary": {"path": str(binary), "sha256": sha256(binary)},
        "out": str(out), "python": sys.version, "numpy": np.__version__, "scipy": scipy.__version__,
        "platform": platform.platform(),
        "settings": {"tol_lnl": TOL_LNL, "residual": RESIDUAL, "tau": TAU, "central_steps": CENTRAL_STEPS,
                     "alpha": ALPHA, "p_invar": PINV, "near_zero": NEAR_ZERO, "r4_given": R4_GIVEN, "wide": WIDE,
                     "seed_runs": SEED_RUNS, "level2_threaded": LEVEL2_THREADED, "turtle_seed": TURTLE_SEED},
    }


def markdown(prov, results):
    lines = ["# S2, S3a and S3b checks of the compiled NQC model", "",
             f"Generated by `{SCRIPT}` on {prov['started_utc']} at commit `{prov['git_commit']}` "
             f"(clean: {prov['git_clean']}); binary SHA-256 `{prov['binary']['sha256']}`.", ""]
    for case, rows in results.items():
        lines += [f"## {case}", "", "| Check | Result | Values |", "|---|---|---|"]
        for row in rows:
            verdict = {True: "pass", False: "FAIL", None: "report"}[row["passed"]]
            values = {k: v for k, v in row.items() if k not in ("check", "passed")}
            lines.append(f"| {row['check']} | {verdict} | `{json.dumps(values, default=str)[:400]}` |")
        lines.append("")
    return "\n".join(lines) + "\n"


CASES = ("root", "seed", "gradient", "target", "refusals", "fits", "nesting", "skewed", "restart", "risk22", "toy",
         "short", "linked-seed", "linked-fits", "linked-refusals", "permuted", "threads", "linked-restart",
         "risk22-turtle", "warm-starts", "linked-short")
SANITIZER_ONLY = ("short", "linked-short")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--binary", required=True)
    ap.add_argument("--alisim", help="the frozen unmodified binary, which simulates the toy data")
    ap.add_argument("--calibration", help="results.json of explore/toy_calibration.py")
    ap.add_argument("--older-run", help="S2's gate run fit-level1, whose checkpoint lacks S3a's keys")
    ap.add_argument("--out", required=True, help="new directory outside the repository")
    ap.add_argument("--cases", help="comma-separated, from " + ",".join(CASES) +
                    " (default: all but short and linked-short)")
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
    keys = a.cases.split(",") if a.cases else [k for k in CASES if k not in SANITIZER_ONLY]
    if set(keys) - set(CASES):
        sys.exit(f"unknown cases: {sorted(set(keys) - set(CASES))}")
    if {"toy", "short", *LINKED_CASES} & set(keys) and not a.alisim or "toy" in keys and not a.calibration:
        sys.exit("the toy case needs --alisim and --calibration; the short case and the twopart cases ("
                 + ",".join(LINKED_CASES) + ") need --alisim")
    older = Path(a.older_run).expanduser().resolve() if a.older_run else None
    if "restart" in keys and not (older and (older / f"{PREFIX}.ckp.gz").exists()):
        sys.exit(f"the restart case needs --older-run, a run directory holding {PREFIX}.ckp.gz")
    out.mkdir(parents=True)
    alisim = Path(a.alisim).expanduser().resolve() if a.alisim else None
    calibration = Path(a.calibration).expanduser().resolve() if a.calibration else None
    if calibration:
        prov["calibration"] = {"path": str(calibration), "sha256": sha256(calibration)}
    if older:
        prov["older_run"] = {"path": str(older), "checkpoint_sha256": sha256(older / f"{PREFIX}.ckp.gz")}
    r = Runner(binary, out, a.timeout)
    c = Context(out)
    ctx, results = {}, {}
    for key in keys:
        print(f"case {key}", flush=True)
        try:
            if key == "root":
                results[key] = case_root(r, c)
            elif key == "seed":
                results[key] = case_seed(r, c)
            elif key == "gradient":
                results[key] = case_gradient(r, c, ctx)
            elif key == "target":
                results[key] = case_target(r, c)
            elif key == "refusals":
                results[key] = case_refusals(r, c)
            elif key == "fits":
                results[key] = case_fits(r, c, ctx)
            elif key == "nesting":
                results[key] = case_nesting(r, c, ctx)
            elif key == "skewed":
                results[key] = case_skewed(r, c)
            elif key == "restart":
                results[key] = case_restart(r, c, ctx, older)
            elif key == "risk22":
                results[key] = case_risk22(r, c, ctx)
            elif key == "toy":
                results[key] = case_toy(r, c, calibration, alisim)
            elif key == "short":
                results[key] = case_short(r, c, alisim)
            elif key == "linked-seed":
                results[key] = case_linked_seed(r, c, ctx, alisim)
            elif key == "linked-fits":
                results[key] = case_linked_fits(r, c, ctx, alisim)
            elif key == "linked-refusals":
                results[key] = case_linked_refusals(r, c, ctx, alisim)
            elif key == "permuted":
                results[key] = case_permuted(r, c, ctx, alisim)
            elif key == "threads":
                results[key] = case_threads(r, c, ctx)
            elif key == "linked-restart":
                results[key] = case_linked_restart(r, c, ctx, alisim)
            elif key == "risk22-turtle":
                results[key] = case_risk22_turtle(r, c, ctx)
            elif key == "warm-starts":
                results[key] = case_warm_starts(r, c, ctx)
            elif key == "linked-short":
                results[key] = case_linked_short(r, c, ctx, alisim)
        except Exception:  # a crash in one case is recorded, and the other cases still run
            results[key] = [check("case ran", False, exception=traceback.format_exc())]
        for row in results[key]:
            print(f"  {row['check']}: {({True: 'pass', False: 'FAIL', None: 'report'})[row['passed']]}", flush=True)
    if "twopart" in ctx:
        prov["twopart"] = {k: ctx["twopart"].made[k] for k in ("sha256", "seeds")}
    prov["finished_utc"] = datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")
    write(out / "results.json", json.dumps({"provenance": prov, "runs": r.runs, "results": results}, indent=1,
                                           default=str))
    write(out / "results.md", markdown(prov, results))
    failed = [row["check"] for rows in results.values() for row in rows if row["passed"] is False]
    print(f"{len(failed)} asserted checks failed; wrote {out / 'results.md'}")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
