#!/usr/bin/env python3
"""Write the S1 unit-test fixtures from the oracle (docs/agent/CODE_PLAN.md, section 3.2).

  python test_scripts/fixedeq/make_fixtures.py [--out unittest/fixtures] [--git git.exe] [--allow-dirty]

One file per chart case: '#' header lines "key: value" (provenance, the case, kappa, and the
oracle's own error against the pass line), then one line per array, a label followed by its
values at 17 significant digits, matrices row-major. targets.txt holds the target-parser cases.

Every case is checked against decision 026's pass lines in the oracle before anything is written,
and if one fails nothing is written. The output directory must be new or empty: fixtures are
never regenerated to make a test pass (CODE_PLAN.md, section 1.1).
"""
import argparse
import datetime
import json
import math
import platform
import re
import sys
from fractions import Fraction as F
from pathlib import Path

import numpy as np
import scipy

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent / "regression"))
from oracle import builtin, cases, chart, target  # noqa: E402
from regress import REPO, capture, write  # noqa: E402

SCRIPT = Path(__file__).resolve().relative_to(REPO).as_posix()
SEED = 20261019
LO, HI = math.log(1e-5), math.log(10.0)   # the default box (decisions 007 and 027)
RESIDUAL = 1e-10                           # decision 026
MIN_FREQ = 1e-4                            # IQ-TREE's default min_state_freq (utils/tools.cpp:7263)
MANIFEST = REPO / "test_scripts/fixedeq/manifest/g0_manifest.json"
FORMAT = ("'# key: value' header lines, then one line per array: a label, then its values "
          "at 17 significant digits, matrices row-major")


def q_limit(kappa):
    """Decision 026: every entry of Q within a relative max(1e-12, 1e-14 kappa)."""
    return max(1e-12, 1e-14 * kappa)


def rel_err(a, b):
    return float(np.max(np.abs(np.asarray(a) - b) / np.abs(b)))


def g17(v):
    return f"{float(v):.17g}"


def provenance(git):
    return {
        "script": SCRIPT,
        "git_commit": capture(git, "rev-parse", "HEAD", cwd=REPO),
        "git_clean": capture(git, "status", "--porcelain", cwd=REPO) == "",
        "generated_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "seed": SEED,
    }


def case(name, kind, truth, ref_rule, pi, ref, Q, z=None, nu=None, **extra):
    """ref_rule: 'row_maxima' (ref = references(Q)), 'lg' (ref = references of LG rebuilt at pi,
    with LG's exchangeabilities in R_ref) or 'given'."""
    pi, Q = np.asarray(pi, dtype=float), np.asarray(Q, dtype=float)
    ref = np.asarray(ref, dtype=int)
    c = {"name": name, "kind": kind, "truth": truth, "ref_rule": ref_rule, "pi": pi, "ref": ref,
         "z": chart.inverse(Q, ref) if z is None else np.asarray(z, dtype=float),
         "nu": pi * -np.diag(Q) if nu is None else np.asarray(nu, dtype=float), "Q": Q}
    c.update({k: np.asarray(v, dtype=float) for k, v in extra.items()})
    return c


def exact_case(name, kind, K, nu, pi):
    """A 3-state case in exact arithmetic: nu K = nu and pi Q = 0 are checked on fractions, then
    every array is rounded once to double."""
    n = len(pi)
    assert all(sum(row) == 1 for row in K) and sum(nu) == 1 and sum(pi) == 1
    assert all(sum(nu[i] * K[i][j] for i in range(n)) == nu[j] for j in range(n))
    Q = [[nu[i] * K[i][j] / pi[i] if i != j else -nu[i] / pi[i] for j in range(n)] for i in range(n)]
    assert all(sum(pi[i] * Q[i][j] for i in range(n)) == 0 for j in range(n))
    Qf = np.array([[float(x) for x in row] for row in Q])
    ref = chart.references(Qf)
    z = [math.log(float(K[i][j] / K[i][ref[i]])) for i, free in enumerate(chart.free_destinations(ref))
         for j in free]
    return case(name, kind, "independent, exact", "row_maxima", [float(p) for p in pi], ref, Qf,
                z=z, nu=[float(v) for v in nu])


def oracle_case(name, kind, ref_rule, pi, ref, z, **extra):
    Q, nu = chart.build(z, pi, ref)
    return case(name, kind, "oracle", ref_rule, pi, ref, Q, z=z, nu=nu, **extra)


def build_cases():
    rng = np.random.default_rng(SEED)
    out = []
    for n in (3, 4, 20):
        for k in (1, 2, 3):
            pi = cases.random_target(n, rng)
            Q = cases.generator_from_flux(cases.cycle_flux(n, rng), pi)
            out.append(case(f"cycle_n{n}_{k}", "balanced flux from directed cycles", "independent",
                            "row_maxima", pi, chart.references(Q), Q))
        for k in (1, 2):
            pi = cases.random_target(n, rng)
            R = cases.random_exchangeabilities(n, rng)
            Q = chart.reversible_seed(R, pi)
            out.append(case(f"reversible_n{n}_{k}", "reversible seed from random exchangeabilities",
                            "independent", "row_maxima", pi, chart.references(Q), Q, R=R))

    third, quarter = F(1, 3), F(1, 4)
    out.append(exact_case("hand_cyclic_n3", "the 3-state cyclic example: A->B->C->A preferred, uniform nu",
                          [[0, 2 * third, third], [third, 0, 2 * third], [2 * third, third, 0]],
                          [third, third, third], [F(1, 2), quarter, quarter]))
    out.append(exact_case("hand_nonuniform_nu_n3", "3-state, nu = (4, 3, 4)/11, so a transposed solve "
                          "fails; row A is a reference tie",
                          [[0, F(1, 2), F(1, 2)], [third, 0, 2 * third], [3 * quarter, quarter, 0]],
                          [F(4, 11), F(3, 11), F(4, 11)], [F(2, 5), F(1, 5), F(2, 5)]))

    man = json.loads(MANIFEST.read_text())["manifest"]
    turtle = np.array([float(v) for v in man["target"]["values_17_digits"]])
    R_lg = builtin.exchangeabilities("LG")[0]
    lg = chart.reversible_seed(R_lg, turtle)
    ref_lg = chart.references(lg)
    assert list(ref_lg) == man["references"]["indices"], "LG references differ from the G0 manifest"
    out.append(case("lg_turtle", "LG rebuilt at the decision 023 target (the S2 default seed)",
                    "independent", "row_maxima", turtle, ref_lg, lg, R=R_lg))
    Q0 = builtin.full_matrix("NQ.PFAM")[0].copy()
    np.fill_diagonal(Q0, 0.0)
    np.fill_diagonal(Q0, -Q0.sum(axis=1))
    out.append(case("nqpfam_turtle", "NQ.pfam's jump chain transferred to the decision 023 target",
                    "oracle", "lg", turtle, ref_lg, chart.transfer_seed(Q0, turtle), Q0=Q0, R_ref=R_lg))

    def lg_ref(pi):
        return chart.references(chart.reversible_seed(R_lg, pi))

    skewed = np.full(20, 1e-4)
    skewed[0] = 1 - skewed[1:].sum()
    for k in (1, 2):
        out.append(oracle_case(f"skewed_n20_{k}", "skewed target, 19 entries at 1e-4", "lg", skewed,
                               lg_ref(skewed), rng.normal(0, 1, 360), R_ref=R_lg))

    ref4 = np.array([1, 0, 3, 2])
    out.append(oracle_case("near_reducible_n4", "near-reducible: blocks {0,1} and {2,3} joined by logits "
                           "of -11.5", "given", cases.random_target(4, rng), ref4, np.full(8, -11.5)))
    ref20 = np.array([1] + [0] * 9 + [11] + [10] * 9)
    blocks = [[0.0 if (j < 10) == (i < 10) else -11.5 for j in free]
              for i, free in enumerate(chart.free_destinations(ref20))]
    out.append(oracle_case("near_reducible_n20", "near-reducible: two blocks of 10 joined by logits of "
                           "-11.5 (explore/thresholds.py)", "given", cases.random_target(20, rng), ref20,
                           np.ravel(blocks)))

    for k in (1, 2):
        pi = cases.random_target(20, rng)
        out.append(oracle_case(f"box_spread_n20_{k}", "coordinates uniform over the default box", "lg",
                               pi, lg_ref(pi), rng.uniform(LO, HI, 360), R_ref=R_lg))
    pi = cases.random_target(20, rng)
    out.append(oracle_case("outside_box_n20", "coordinates uniform over [-16, 12], beyond the box "
                           "(decision 027: the chart is total)", "lg", pi, lg_ref(pi),
                           rng.uniform(-16.0, 12.0, 360), R_ref=R_lg))
    pi = cases.random_target(20, rng)
    z = rng.normal(0, 1, 360)
    z[0] = -690.0
    out.append(oracle_case("near_underflow_n20", "one logit at -690, near the exp underflow but every "
                           "rate a normal double", "lg", pi, lg_ref(pi), z, R_ref=R_lg))
    return out


def check(c):
    """Decision 026's lines in the oracle itself: the fixture's Q is a family member, the oracle's
    build from the fixture's z reproduces Q, and the stated references and seeds hold."""
    kappa = chart.nu_condition(chart.jump_of(c["Q"]))
    limit = q_limit(kappa)
    err = rel_err(chart.build(c["z"], c["pi"], c["ref"])[0], c["Q"])
    d = chart.diagnostics(c["Q"], c["pi"])
    nums = {"kappa": kappa, "q_limit": limit, "oracle_q_error": err,
            **{k: d[k] for k in ("stationarity", "row_sums", "mean_rate_error", "min_rate")}}
    fails = []
    if err > limit:
        fails.append(f"oracle build error {err:.3e} exceeds {limit:.3e}")
    fails += [f"{k} {d[k]:.3e} exceeds {RESIDUAL}" for k in ("stationarity", "row_sums", "mean_rate_error")
              if d[k] > RESIDUAL]
    if not d["min_rate"] > 0:
        fails.append("an off-diagonal rate is not positive")
    if c["ref_rule"] == "row_maxima" and list(chart.references(c["Q"])) != list(c["ref"]):
        fails.append("ref is not the row maxima of Q")
    if c["ref_rule"] == "lg" and list(chart.references(chart.reversible_seed(c["R_ref"], c["pi"]))) != list(c["ref"]):
        fails.append("ref is not the row maxima of LG rebuilt at pi")
    if "R" in c and rel_err(chart.reversible_seed(c["R"], c["pi"]), c["Q"]) > limit:
        fails.append("Q is not the reversible seed of R")
    if "Q0" in c and rel_err(chart.transfer_seed(c["Q0"], c["pi"]), c["Q"]) > limit:
        fails.append("Q is not the transfer of Q0")
    return nums, fails


def render(c, prov):
    nums, fails = check(c)
    head = {"note": "S1 unit-test fixture (docs/agent/CODE_PLAN.md section 3.2); never regenerated to make a test pass",
            "format": FORMAT, **prov, "name": c["name"], "kind": c["kind"], "truth": c["truth"],
            "ref_rule": c["ref_rule"], "n": len(c["pi"]), **{k: g17(v) for k, v in nums.items()}}
    lines = [f"# {k}: {v}" for k, v in head.items()]
    for label in ("pi", "ref", "z", "nu", "Q", "R", "Q0", "R_ref"):
        if label in c:
            values = np.ravel(c[label])
            lines.append(" ".join([label] + [str(int(v)) if label == "ref" else g17(v) for v in values]))
    return "\n".join(lines) + "\n", nums, fails


def target_cases(turtle_strings):
    """(name, text, min_freq). The oracle's parse_target gives each verdict."""
    turtle = [float(s) for s in turtle_strings]
    u = ["0.05"] * 20
    rest = ["0.05"] * 18
    return [
        ("turtle_17_digits_commas", ",".join(turtle_strings), MIN_FREQ),
        ("turtle_8_decimals_spaces", " ".join(f"{v:.8f}" for v in turtle), MIN_FREQ),
        ("turtle_5_decimals", ",".join(f"{v:.5f}" for v in turtle), MIN_FREQ),
        ("uniform_slashes", "/".join(u), MIN_FREQ),
        ("uniform_tabs", "\t".join(u), MIN_FREQ),
        ("uniform_doubled_and_outer_separators", ",," + ",,".join(u) + ", ", MIN_FREQ),
        ("uniform_scientific", ",".join(["5e-02"] * 20), MIN_FREQ),
        ("uniform_signs_and_bare_points", ",".join(["+0.05", ".05", "0.050", "5.e-2"] * 5), MIN_FREQ),
        ("uniform_above_a_raised_minimum", ",".join(u), 0.06),
        ("sum_plus_0.9e-6", ",".join(u[:19] + ["0.0500009"]), MIN_FREQ),
        ("sum_plus_1.1e-6", ",".join(u[:19] + ["0.0500011"]), MIN_FREQ),
        ("sum_minus_1.1e-6", ",".join(u[:19] + ["0.0499989"]), MIN_FREQ),
        ("entry_at_the_minimum", ",".join(["0.0001", "0.0999"] + rest), MIN_FREQ),
        ("entry_below_the_minimum", ",".join(["0.00005", "0.09995"] + rest), MIN_FREQ),
        ("19_entries", ",".join(rest + ["0.1"]), MIN_FREQ),
        ("21_entries", ",".join(u[:19] + ["0.025", "0.025"]), MIN_FREQ),
        ("zero_entry", ",".join(["0", "0.1"] + rest), MIN_FREQ),
        ("negative_entry", ",".join(["-0.05", "0.15"] + rest), MIN_FREQ),
        ("nan_entry", ",".join(["nan"] + u[:19]), MIN_FREQ),
        ("inf_entry", ",".join(["inf"] + u[:19]), MIN_FREQ),
        ("word_entry", ",".join(["abc"] + u[:19]), MIN_FREQ),
        ("hex_entry", ",".join(["0x1.999999999999ap-5"] + u[:19]), MIN_FREQ),
        ("overflowing_entry", ",".join(["1e400"] + u[:19]), MIN_FREQ),
        ("underflowing_entry", ",".join(["1e-400", "0.1"] + rest), MIN_FREQ),
        ("empty", "", MIN_FREQ),
    ]


def render_targets(turtle_strings, prov):
    """The parser cases with the oracle's verdicts. An accepted vector is recorded normalized as the
    compiled parser does it: summed left to right, then each entry divided by the sum."""
    lines = [f"# {k}: {v}" for k, v in {
        "note": "S1 target-parser fixture (decision 006); never regenerated to make a test pass",
        "format": "blocks of 'case NAME', 'text STRING' (the rest of the line, verbatim), 'min_freq X', "
                  "'verdict accept|reject', and when accepted 'raw_sum X' and 'values ...', ended by 'end'",
        **prov}.items()]
    fails = []
    for name, text, min_freq in target_cases(turtle_strings):
        assert "\n" not in text and "\r" not in text
        lines += [f"case {name}", f"text {text}", f"min_freq {g17(min_freq)}"]
        try:
            target.parse_target(text, min_state_freq=min_freq)
            verdict = "accept"
        except ValueError:
            verdict = "reject"
        lines.append(f"verdict {verdict}")
        if verdict == "accept":
            values = [float(t) for t in re.split(r"[,\s/]+", text.strip()) if t]   # parse_target's split
            raw = 0.0
            for v in values:
                raw += v
            if abs(raw - 1.0) > 1e-6:
                fails.append(f"{name}: the left-to-right sum {raw!r} disagrees with the oracle's verdict")
            lines += [f"raw_sum {g17(raw)}", " ".join(["values"] + [g17(v / raw) for v in values])]
        lines.append("end")
    return "\n".join(lines) + "\n", fails


def write_all(out, fixtures, turtle_strings, prov):
    """Write every fixture and targets.txt, or nothing if any check fails. Returns (summary, failures)."""
    rendered = [(c["name"], *render(c, prov)) for c in fixtures]
    targets_text, failures = render_targets(turtle_strings, prov)
    failures += [f"{name}: {f}" for name, _, _, fails in rendered for f in fails]
    if failures:
        return [], failures
    out.mkdir(parents=True, exist_ok=True)
    for name, text, _, _ in rendered:
        write(out / f"{name}.txt", text)
    write(out / "targets.txt", targets_text)
    return [(name, nums) for name, _, nums, _ in rendered], []


def read_fixture(text):
    """(header, arrays) of a fixture file; the inverse of render."""
    header, arrays = {}, {}
    for line in text.splitlines():
        if line.startswith("# "):
            key, _, value = line[2:].partition(": ")
            header[key] = value
        elif line:
            label, *values = line.split()
            arrays[label] = np.array([int(v) for v in values]) if label == "ref" else np.array([float(v) for v in values])
    return header, arrays


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", default=str(REPO / "unittest" / "fixtures"))
    ap.add_argument("--git", default="git", help="git executable (git.exe under WSL)")
    ap.add_argument("--allow-dirty", action="store_true")
    a = ap.parse_args()
    prov = provenance(a.git)   # read once, before anything is written
    if not prov["git_clean"] and not a.allow_dirty:
        sys.exit(f"the repository is not clean ({prov['git_commit']}); commit first or pass --allow-dirty")
    out = Path(a.out)
    if out.exists() and any(out.iterdir()):
        sys.exit(f"{out} is not empty: fixtures are never regenerated")
    turtle_strings = json.loads(MANIFEST.read_text())["manifest"]["target"]["values_17_digits"]
    summary, failures = write_all(out, build_cases(), turtle_strings, prov)
    if failures:
        print("nothing written; cases that fail decision 026's lines in the oracle:")
        for f in failures:
            print(f"  {f}")
        sys.exit(1)
    print(f"{len(summary)} fixtures and targets.txt written to {out} (commit {prov['git_commit']}, "
          f"clean {prov['git_clean']})")
    print(f"{'name':<24} {'kappa':>9} {'limit':>9} {'oracle':>9} {'worst residual':>14}")
    for name, s in summary:
        worst = max(s["stationarity"], s["row_sums"], s["mean_rate_error"])
        print(f"{name:<24} {s['kappa']:9.2e} {s['q_limit']:9.2e} {s['oracle_q_error']:9.2e} {worst:14.2e}")


if __name__ == "__main__":
    main()
