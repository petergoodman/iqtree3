#!/usr/bin/env python3
"""S2 checks of the compiled NQC model (docs/agent/CODE_PLAN.md section 4; decision 030).

  python nqcchecks.py --binary ~/iqtree3-build/iqtree3 --alisim ~/iqtree3-baseline/iqtree3 \
      --calibration CALIBRATION/results.json --out DIR [--cases root,seed,...] [--git git.exe]

Cases:
  root      a rooted input tree keeps its root and an unrooted one is rooted at the midpoint
  seed      the seed Q and the first log-likelihood against the oracle (decisions 021, 026)
  gradient  the first gradient against the oracle, by the same rule and against a central
            reference, at the LG seed and at a start near zero (decision 030)
  target    targets refused for their sum and for an entry below the minimum
  refusals  every combination decisions 008, 032 and 034 refuse
  fits      Level 1 and Level 2 fits on aa_example, +I+G4 included, and the domain option
  toy       the compiled optimizer against the oracle's port (decisions 017, 030), using the
            line measured by explore/toy_calibration.py (needs --alisim and --calibration)
  short     a Level 1 and a +I+G4 Level 2 fit on 300 sites of the toy data, one round each, for
            the sanitizer build (CODE_PLAN.md section 3.6; needs --alisim)

Every run has its own directory under DIR and uses -seed 1 -T 1. Final states are read from the
class's checkpoint; the seed Q and the first gradient from its verbose output (decision 034).
Exit status 1 if an asserted check fails; checks marked "report" are not asserted.
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
from probes import (AA, NEX, TURTLE, Runner, baseline_runs, first_taxa, mdef, read_checkpoint,  # noqa: E402
                    tree_comparison, tree_text, turtle_target, unrooted_text)
from differential import datasets, significant_digits, write_matrix  # noqa: E402
from regress import REPO, SCALARS, capture, read, sha256, write  # noqa: E402
import toy  # noqa: E402

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


def nqc_checkpoint(run):
    """The class's checkpoint keys, without their struct prefix."""
    return {k.split("ModelNonrevFixedEq!", 1)[1]: v for k, v in run.checkpoint().items()
            if "ModelNonrevFixedEq!" in k}


def floats(text):
    return np.array([float(v) for v in text.split(",")])


def last(regex, text):
    found = regex.findall(text)
    return found[-1] if found else None


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
    cases = [("aa-none", "aa_example", "", None, 0.0), ("aa-G4", "aa_example", f"+G4{{{ALPHA}}}", ALPHA, 0.0),
             ("aa-I", "aa_example", f"+I{{{PINV}}}", None, PINV),
             ("aa-IG4", "aa_example", f"+I{{{PINV}}}+G4{{{ALPHA}}}", ALPHA, PINV),
             ("turtle1-G4", "turtle_part1", f"+G4{{{ALPHA}}}", ALPHA, 0.0)]
    Q_or = chart.build(c.z, c.pi, c.ref)[0]
    for name, dname, suffix, shape, p_invar in cases:
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
                                             p_invar=p_invar)
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


def case_gradient(r, c, ctx):
    out = []
    objective = nqc.Objective(c.aa, c.tree02, c.pi, c.ref, shape=ALPHA)
    run = ctx.get("level1") or r.run("fit-level1", ["-s", AA, "-m", c.model + f"+G4{{{ALPHA}}}", "-te", "input.tree",
                                                     "-blfix", "-vv"], {"input.tree": c.t02})
    ctx["level1"] = run
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
        ("model-joint", ["-s", TURTLE, "-p", NEX, "--model-joint", c.model, "-nparam", "2"], {}),
        ("link-model", ["-s", TURTLE, "-p", NEX, "-m", c.model, "--link-model", "-nparam", "2"], {}),
        ("tree-mixture", ["-s", AA, "-m", c.model + "+T", "-te", "two.trees", "-nparam", "2"], {"two.trees": two}),
        ("ghost", aa + ["-m", c.model + "+H4"], tree),
        ("site-freq", aa + ["-m", c.model, "-ft", "input.tree"], tree),
        # IQ-TREE itself refuses the flag unless GTR or GTR20 is in the model
        ("link-exchange", aa + ["-m", f"MIX{{NQC+F{NAMED},GTR20}}", "--mdef", "target.nex", "--link-exchange-rates"],
         {**named, **tree}),
        ("eigen", aa + ["-m", c.model, "--eigen"], tree),
        ("modelfinder", ["-s", AA, "-m", "MFP", "-mset", f"NQC+F{NAMED}", "--mdef", "target.nex", "-te", "input.tree"],
         {**named, **tree}),
        ("alisim", ["--alisim", "sim", "-m", c.model, "-t", "input.tree", "--seqtype", "AA"], tree),
        ("mix", aa + ["-m", f"MIX{{NQC+F{NAMED},LG}}", "--mdef", "target.nex"], {**named, **tree}),
        ("fmix", aa + ["-m", f"NQC+FMIX{{{NAMED},{NAMED}2}}", "--mdef", "target.nex"], {**fmix, **tree}),
        ("mdef-mix", aa + ["-m", "NQCMIX", "--mdef", "target.nex"], {**mixdef, **tree}),
    ]
    out = []
    for name, args, files in runs:
        run = r.run(f"refuse-{name}", args, files)
        message = next((l for l in run.stdout().splitlines() if "NQC" in l and "ERROR" in l), None)
        out.append(check(f"refused: {name}", run.code not in (0, "timeout") and message is not None,
                         exit_code=run.code, message=message))
    return out


START_LNL = re.compile(r"Init pinv, alpha: \S+ \S+ / Estimate: \S+ \S+ / LogL: (\S+)")
BEST_LNL = re.compile(r"Optimal pinv,alpha: \S+ \S+ / LogL: (\S+)")


def restore_check(label, text):
    """The +I+G path keeps the best of its starts through the model's checkpoint (decision 033): the
    model it ends with must score as the best start did, at the 3 decimals IQ-TREE prints."""
    starts, best = [float(v) for v in START_LNL.findall(text)], BEST_LNL.search(text)
    ok = bool(starts and best and abs(float(best.group(1)) - max(starts)) <= 5e-4)
    return check(f"{label}: the restored model scores as the best +I+G start", ok,
                 best_start=max(starts) if starts else None, final=best.group(1) if best else None)


def fit_checks(label, run, c, expected_parameters, level, start=None):
    if run.code != 0:
        return failed_run(run)
    text, ck = run.stdout(), nqc_checkpoint(run)
    out = []
    pi_ck = floats(ck["target"]) if "target" in ck else None
    out.append(check(f"{label}: target in the checkpoint bitwise equal to the parsed target",
                     pi_ck is not None and np.array_equal(pi_ck, c.pi)))
    res = last(RESIDUALS, text)
    worst = max(float(v) for v in res[:3]) if res else None
    out.append(check(f"{label}: residuals at most {RESIDUAL:g} (decision 026)", worst is not None and worst <= RESIDUAL,
                     stationarity_row_sums_mean_rate=res[:3] if res else None, nu=res[3] if res else None,
                     minimum_rate=res[4] if res else None))
    m = re.search(SCALARS["free_parameters"], run.report(), re.M)
    out.append(check(f"{label}: free parameters", bool(m) and int(m.group(1)) == expected_parameters,
                     reported=m.group(1) if m else None, expected=expected_parameters))
    comp = tree_comparison(c.t02.strip(), run.trees()[0]) if run.trees() else None
    out.append(check(f"{label}: root split unchanged", bool(comp and comp["root_split_kept"])))
    if level == 1:
        out.append(check(f"{label}: root-adjacent lengths unchanged", bool(comp and comp["root_lengths_kept"])))
    fits = FIT.findall(text)
    falls = [f for f in fits if float(f[3]) < float(f[1])]
    out.append(check(f"{label}: fits ending below their start (report; decision 029)", None,
                     fits=len(fits), below_start=falls))
    gaps = [abs(float(f[2]) - float(f[3])) for f in fits]
    out.append(check(f"{label}: optimizer score against the re-scored point, largest gap (report)", None,
                     largest=max(gaps) if gaps else None))
    bounds = last(BOUNDS, text)
    out.append(check(f"{label}: coordinates at the bounds (report)", None,
                     lower_upper_outside=bounds))
    if "+I+G" in label:
        out.append(restore_check(label, text))
    out.append(check(f"{label}: scaling-squaring switches (report; risk 12)", None, count=text.count(SWITCH)))
    out.append(check(f"{label}: ASSERT reached (report; risk 7)", None, found="Assertion" in text or "ASSERT" in text))
    return out


def case_fits(r, c, ctx):
    out = []
    level1 = ctx.get("level1") or r.run("fit-level1", ["-s", AA, "-m", c.model + f"+G4{{{ALPHA}}}", "-te", "input.tree",
                                                        "-blfix", "-vv"], {"input.tree": c.t02})
    ctx["level1"] = level1
    out += fit_checks("Level 1", level1, c, 360, 1)
    if level1.code == 0:
        ck = nqc_checkpoint(level1)
        z = floats(ck["coordinates"])
        ref_ok = [int(v) for v in ck["references"].split(",")] == [int(v) for v in c.ref]
        out.append(check("Level 1: references in the checkpoint equal the oracle's", ref_ok))
        oracle = -nqc.Objective(c.aa, c.tree02, c.pi, c.ref, shape=ALPHA)(z)
        accepted = float(ck["accepted_score"])
        out.append(check("Level 1: accepted score against the oracle's re-evaluation (decisions 021, 030)",
                         abs(accepted - oracle) <= TOL_LNL * abs(oracle), accepted=accepted, oracle=oracle))
        printed = re.search(SCALARS["log_likelihood"], level1.report(), re.M)
        out.append(check("Level 1: printed final score against the accepted score (report)", None,
                         printed=printed.group(1) if printed else None,
                         gap=abs(float(printed.group(1)) - accepted) if printed else None))
    for label, suffix, extra in (("Level 2 +G4", "+G4", 1), ("Level 2 +I+G4", "+I+G4", 2)):
        run = r.run("fit-" + label.lower().replace(" ", "-").replace("+", ""),
                    ["-s", AA, "-m", c.model + suffix, "-te", "input.tree", "-v"], {"input.tree": c.t02})
        out += fit_checks(label, run, c, 360 + extra + BRANCHES_AA, 2)
    run = r.run("fit-domain", ["-s", AA, "-m", c.model + f"+G4{{{ALPHA}}}", "-te", "input.tree", "-blfix",
                               "-nparam", "3", "--nqc-domain", "1e-6,100"], {"input.tree": c.t02})
    if run.code != 0:
        out += failed_run(run)
    else:
        domain = floats(nqc_checkpoint(run).get("domain", "nan,nan"))
        want = np.array([math.log(1e-6), math.log(100.0)])
        out.append(check("--nqc-domain 1e-6,100 sets the domain in the checkpoint", bool(np.array_equal(domain, want)),
                         domain=domain.tolist()))
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
    ck = nqc_checkpoint(run)
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
        ck = nqc_checkpoint(run)
        out.append(check(f"{label}: target in the checkpoint bitwise equal to the parsed target",
                         "target" in ck and np.array_equal(floats(ck["target"]), c.pi)))
        res = last(RESIDUALS, run.stdout())
        worst = max(float(v) for v in res[:3]) if res else None
        out.append(check(f"{label}: residuals at most {RESIDUAL:g} (decision 026)",
                         worst is not None and worst <= RESIDUAL, worst=worst))
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
                     "alpha": ALPHA, "p_invar": PINV, "near_zero": NEAR_ZERO},
    }


def markdown(prov, results):
    lines = ["# S2 checks of the compiled NQC model", "",
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


CASES = ("root", "seed", "gradient", "target", "refusals", "fits", "toy", "short")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--binary", required=True)
    ap.add_argument("--alisim", help="the frozen unmodified binary, which simulates the toy data")
    ap.add_argument("--calibration", help="results.json of explore/toy_calibration.py")
    ap.add_argument("--out", required=True, help="new directory outside the repository")
    ap.add_argument("--cases", help="comma-separated, from " + ",".join(CASES) + " (default: all but short)")
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
    keys = a.cases.split(",") if a.cases else [k for k in CASES if k != "short"]
    if set(keys) - set(CASES):
        sys.exit(f"unknown cases: {sorted(set(keys) - set(CASES))}")
    if {"toy", "short"} & set(keys) and not a.alisim or "toy" in keys and not a.calibration:
        sys.exit("the toy case needs --alisim and --calibration, and the short case --alisim")
    out.mkdir(parents=True)
    alisim = Path(a.alisim).expanduser().resolve() if a.alisim else None
    calibration = Path(a.calibration).expanduser().resolve() if a.calibration else None
    if calibration:
        prov["calibration"] = {"path": str(calibration), "sha256": sha256(calibration)}
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
            elif key == "toy":
                results[key] = case_toy(r, c, calibration, alisim)
            elif key == "short":
                results[key] = case_short(r, c, alisim)
        except Exception:  # a crash in one case is recorded, and the other cases still run
            results[key] = [check("case ran", False, exception=traceback.format_exc())]
        for row in results[key]:
            print(f"  {row['check']}: {({True: 'pass', False: 'FAIL', None: 'report'})[row['passed']]}", flush=True)
    prov["finished_utc"] = datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")
    write(out / "results.json", json.dumps({"provenance": prov, "runs": r.runs, "results": results}, indent=1,
                                           default=str))
    write(out / "results.md", markdown(prov, results))
    failed = [row["check"] for rows in results.values() for row in rows if row["passed"] is False]
    print(f"{len(failed)} asserted checks failed; wrote {out / 'results.md'}")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
