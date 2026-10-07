#!/usr/bin/env python3
"""The G0 run manifest (synthesis section 11, gate G0; docs/agent/CODE_PLAN.md section 4).

  python manifest.py --out test_scripts/fixedeq/manifest [--git git.exe]

Writes g0_manifest.json and g0_manifest.md. Policy entries cite the decision that sets them and do
not restate its reasoning; computed entries come from the oracle and the repository's files. Both
files are derived: regenerate them, never edit them by hand.
"""
import argparse
import datetime
import hashlib
import json
import platform
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE.parent / "probes"))
from oracle import builtin, chart, readers, target  # noqa: E402
from probes import ckp_values, read_checkpoint  # noqa: E402

REPO = HERE.parents[2]
SCRIPT = Path(__file__).resolve().relative_to(REPO).as_posix()
ALIGNMENT_CPP = REPO / "alignment" / "alignment.cpp"
TURTLE = REPO / "test_scripts" / "test_data" / "turtle_aa.fasta"
NEX = REPO / "test_scripts" / "test_data" / "turtle_aa.nex"
AA = REPO / "example" / "aa_example.phy"
BASELINE = REPO / "test_scripts" / "fixedeq" / "regression" / "baseline" / "baseline.json"
LOCK = REPO / "test_scripts" / "fixedeq" / "environment.lock.txt"
DOMAIN = (float(np.log(1e-5)), float(np.log(10.0)))   # D03, decision 007
UPSTREAM = "63c330d9 (upstream tag v3.1.4)"
FROZEN_ASAN = Path.home() / "iqtree3-baseline-asan"
# decision 025's incumbent: S0 probe (g)'s GTR20+F{pi*} fit on aa_example
INCUMBENT_CKP = Path.home() / "iqtree3-runs" / "probes" / "record-20261005T174330Z" / "g0" / "probe.ckp.gz"


def sha256_bytes(data):
    return hashlib.sha256(data).hexdigest()


def sha256(path):
    return sha256_bytes(Path(path).read_bytes())


def state_order():
    """IQ-TREE's protein state order, read from the source (alignment/alignment.cpp:34)."""
    m = re.search(r'char symbols_protein\[\] = "([A-Z]+)"', ALIGNMENT_CPP.read_text(encoding="utf-8"))
    return m.group(1).rstrip("X")


def test_target():
    """Decision 023: the turtle pooled composition under -p conventions."""
    return target.pooled_frequencies(readers.read_fasta(TURTLE), readers.read_charsets(NEX), "p")


def lg_seed(pi):
    """D05's default start: LG's exchangeabilities rebuilt at the target, unit mean rate."""
    R, _ = builtin.exchangeabilities("LG")
    return chart.reversible_seed(R, pi)


def incumbent_exchangeabilities():
    """The incumbent's R from its checkpoint; reversible rates[] is the upper triangle, row-major
    (model/modelmarkov.cpp:114-118)."""
    R = np.zeros((20, 20))
    R[np.triu_indices(20, 1)] = ckp_values(read_checkpoint(INCUMBENT_CKP)["ModelProtein!rates"])
    return R + R.T


def incumbent_placement(pi, ref):
    """Where the incumbent rebuilt at the target lies in the domain (PLAN.md risk 22)."""
    if not INCUMBENT_CKP.exists():
        return None
    z = chart.inverse(chart.reversible_seed(incumbent_exchangeabilities(), pi), ref)
    return {"checkpoint": str(INCUMBENT_CKP), "sha256": sha256(INCUMBENT_CKP),
            "min_coordinate": float(z.min()), "max_coordinate": float(z.max()),
            "outside_domain": int(np.sum((z < DOMAIN[0]) | (z > DOMAIN[1]))), "see": "PLAN.md risk 22"}


def capture(*argv, cwd=None):
    import subprocess
    try:
        return subprocess.run(argv, cwd=cwd, capture_output=True, text=True, check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError) as e:
        return f"unavailable ({e})"


def build(git):
    states = state_order()
    pi = test_target()
    Q = lg_seed(pi)
    ref = chart.references(Q)
    z = chart.inverse(Q, ref)
    seed = chart.diagnostics(Q, pi)
    base = json.loads(BASELINE.read_text(encoding="utf-8"))
    frozen = base["provenance"]["extracts"][0]["binary"]
    trees = {run: base["runs"][run]["items"]["(treefile)"]["text"] for run in ("run02", "run05", "run08", "run09")}
    san = FROZEN_ASAN / "sanitizer-build.json"
    san_rec = json.loads(san.read_text()) if san.exists() else None
    san_sum = FROZEN_ASAN / "iqtree3.sha256"
    import scipy
    return {
        "state_order": {
            "value": states,
            "source": "alignment/alignment.cpp:34 (symbols_protein, X excluded)",
            "matches_oracle": states == builtin.STATES,
        },
        "chart": {
            "value": "jump-chain log-ratio coordinates",
            "set_by": "D01",
            "coordinates": 360,
            "order": "row-major over rows i; within row i, destinations j other than i and its reference r(i), "
                     "increasing (test_scripts/fixedeq/oracle/chart.py)",
            "definition": "z_ij = log(K_ij / K_i,r(i)); K the jump matrix, nu its stationary distribution, "
                          "q_ij = nu_i K_ij / pi*_i",
            "later_slices": "none: no positive-ratio or T3 chart is compiled (decision 028)",
        },
        "references": {
            "rule": "row maxima of the off-diagonal rates of LG rebuilt at the target, ties to the lower index",
            "set_by": "D01, decision 011",
            "indices": [int(r) for r in ref],
            "destinations": "".join(states[int(r)] for r in ref),
        },
        "target": {
            "set_by": "decision 023 (provenance policy D07)",
            "values_17_digits": [f"{v:.17g}" for v in pi],
            "sum_minus_1": float(pi.sum() - 1.0),
            "minimum": float(pi.min()),
            "minimum_state": states[int(np.argmin(pi))],
            "at_least_min_state_freq_1e-4": bool(pi.min() >= 1e-4),
            "above_zero_freq_1e-10": bool(pi.min() > 1e-10),
            "estimator": "test_scripts/fixedeq/oracle/target.py, a port of IQ-TREE's estimator: per partition, "
                         "all-unknown sequences removed (at least three kept); states counted with missing taxa "
                         "padded as unknown cells (-p conventions); ambiguity codes B, Z, J and unknown "
                         "characters distributed by the 8-round fixed point of convertCountToFreq; no floor "
                         "(keep_zero_freq default); no smoothing",
            "cross_check": "reproduces IQ-TREE's printed 'Mean state frequencies' for this data under -p and -S "
                           "at the 8 decimals printed (2026-10-05, CHANGELOG.md)",
            "inputs": {str(p.relative_to(REPO)): sha256(p) for p in (TURTLE, NEX)},
        },
        "domain": {
            "value": [DOMAIN[0], DOMAIN[1]],
            "meaning": "z in [log 1e-5, log 10], the historical benchmark domain",
            "set_by": "D03, decisions 007 and 027",
            "lg_seed_min_coordinate": float(z.min()),
            "lg_seed_max_coordinate": float(z.max()),
            "lg_seed_distance_to_lower": float(z.min() - DOMAIN[0]),
            "lg_seed_distance_to_upper": float(DOMAIN[1] - z.max()),
            "lg_seed_inside": bool(DOMAIN[0] < z.min() and z.max() < DOMAIN[1]),
            "note": "no start is checked against the domain and no margin is set (decision 027)",
        },
        "lg_seed_residuals": seed,
        "derivative_policy": {
            "value": "forward difference with signed step h = eta * max(s, |x|), eta = 1e-4, s = 1, exactly "
                     "representable, backward when the forward probe would leave the domain; legacy step for "
                     "every other model",
            "set_by": "D02, decision 005",
        },
        "root_policy": {
            "set_by": "decision 016",
            "root_frequencies": "pi*",
            "input_rooting": "an unrooted input tree is rooted by IQ-TREE's conversion at the midpoint of the "
                             "longest path (per partition under -S, on the shared tree under -p and -q); a rooted "
                             "input keeps its root",
            "root_edge_search": "none: fixed topology under -te, without --root-find",
            "root_split": "fixed at Level 1 (-blfix); IQ-TREE's branch-length optimization at Level 2",
            "observed": "S0 probe (f), 2026-10-05: the root never changed branch in 13 comparisons; unrooted inputs "
                        "were rooted at the midpoint; in every Level 2 fit the root slid along its edge until one "
                        "root-adjacent length was about 2e-6",
        },
        "tree_and_rate_policy": {
            "set_by": "decisions 016 and 024",
            "level_1": {
                "one alignment": "-te TREE -blfix, rate parameters in braces",
                "-S (separate trees)": "-S PARTFILE -te TREES -blfix --model-joint, rate parameters in braces",
                "edge-proportional": "-q PARTFILE -te TREE -blfix --model-joint, with a charpartition giving each "
                                     "partition its braced rate model and its rate as {x}, taken from the preceding "
                                     "-p fit's checkpoint",
            },
            "level_2": "-p or -S with -te, as nQMaker uses them; rate parameters and branch lengths optimized",
            "fixed_rate_values_in_s0_tests": "+G4{0.5}, +I{0.2}",
            "input_trees": {f"baseline {run} tree file": sha256_bytes(text.encode("utf-8"))
                            for run, text in trees.items()},
        },
        "starts": {
            "default": "LG rebuilt at pi* (D05); its residuals are lg_seed_residuals",
            "incumbent": "GTR20+F{pi*} fitted on fixed trees and carried from its checkpoint (D05, decision 025)",
            "policy": "a start is evaluated as given and never checked against the domain; the line search "
                      "clamps each later trial point into it (decision 027)",
            "incumbent_placement": incumbent_placement(pi, ref),
            "iqtree_seed": "-seed 1 and -T 1 in every S0 run, except the thread-count runs 11 and 12",
            "oracle_seeds": "20261004 in test_scripts/fixedeq/differential/differential.py and probes/probes.py",
            "s1_fixture_seeds": "recorded in each fixture's provenance header (S1)",
        },
        "source": {
            "upstream_baseline": UPSTREAM,
            "manifest_commit": capture(git, "rev-parse", "HEAD", cwd=REPO),
            "frozen_binary": frozen,
            "sanitizer_build": (None if san_rec is None else
                                {"record": str(san), "variant": san_rec["variant"],
                                 "binary_sha256": (san_rec["binary"] or {}).get("sha256"),
                                 "built_from_commit": san_rec["git_commit"],
                                 "frozen_binary": {"path": str(FROZEN_ASAN / "iqtree3"),
                                                   "sha256": san_sum.read_text().split()[0]
                                                   if san_sum.exists() else None}}),
        },
        "thresholds": {
            "oracle": "decision 021",
            "compiled_code": "decisions 026 and 028: residuals at most 1e-10 in tests, reported and never used "
                             "to reject; Q entries and round trips within a relative max(1e-12, 1e-14 kappa), "
                             "kappa in each fixture's header; a target sum within 1e-6 of 1",
            "nesting_test": "decision 025: 1e-8 x |lnL| against the re-evaluated incumbent",
        },
        "environment": {
            "lock_file": {str(LOCK.relative_to(REPO)): sha256(LOCK)},
            "numpy": np.__version__, "scipy": scipy.__version__,
        },
    }


def markdown(m, prov):
    t = m["target"]
    d = m["domain"]
    s = m["starts"]
    place = s["incumbent_placement"]
    outside = "" if place is None else (f"; the incumbent has {place['outside_domain']} of 360 coordinates "
                                        f"outside the domain ({place['see']})")
    lines = [
        "# G0 run manifest", "",
        f"Generated by `{SCRIPT}` at commit `{prov['git_commit']}` (clean: {prov['git_clean']}) on "
        f"{prov['generated_utc']}. Do not edit by hand; `g0_manifest.json` holds every entry.", "",
        "| Item | Value | Set by |", "|---|---|---|",
        f"| State order | `{m['state_order']['value']}` | {m['state_order']['source']} |",
        f"| Chart | {m['chart']['value']}, {m['chart']['coordinates']} coordinates | {m['chart']['set_by']} |",
        f"| References | `{m['references']['destinations']}` (row maxima of LG at π*) | {m['references']['set_by']} |",
        f"| Target | turtle pooled composition, minimum {t['minimum']:.6g} ({t['minimum_state']}) | {t['set_by']} |",
        f"| Domain | [{d['value'][0]:.6f}, {d['value'][1]:.6f}]; LG seed in [{d['lg_seed_min_coordinate']:.6f}, "
        f"{d['lg_seed_max_coordinate']:.6f}] | {d['set_by']} |",
        f"| Derivative policy | {m['derivative_policy']['value']} | {m['derivative_policy']['set_by']} |",
        f"| Root policy | {m['root_policy']['root_edge_search']}; {m['root_policy']['root_split']} | "
        f"{m['root_policy']['set_by']} |",
        f"| Level 1 (edge-proportional) | {m['tree_and_rate_policy']['level_1']['edge-proportional']} | "
        f"{m['tree_and_rate_policy']['set_by']} |",
        f"| Starts | {s['default']}; {s['incumbent']}; {s['policy']}{outside} | D05, decisions 025 and 027 |",
        f"| Pass lines | oracle: {m['thresholds']['oracle']}; compiled code: {m['thresholds']['compiled_code']} | "
        f"decisions 021 and 026 |",
        f"| Upstream baseline | {m['source']['upstream_baseline']}; frozen binary `{m['source']['frozen_binary']['sha256']}` | |",
        "", "Target π* at 17 significant digits, in state order:", "",
        "| State | π* |", "|---|---|",
    ]
    lines += [f"| {s} | {v} |" for s, v in zip(m["state_order"]["value"], t["values_17_digits"])]
    lines += ["", f"LG seed residuals: {json.dumps(m['lg_seed_residuals'])}"]
    return "\n".join(lines) + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True)
    ap.add_argument("--git", default="git", help="git executable (git.exe under WSL)")
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    manifest = build(a.git)
    prov = {"script": SCRIPT, "git_commit": manifest["source"]["manifest_commit"],
            "git_clean": capture(a.git, "status", "--porcelain", cwd=REPO) == "",
            "generated_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
            "python": sys.version.split()[0], "platform": platform.platform()}
    (out / "g0_manifest.json").write_text(json.dumps({"provenance": prov, "manifest": manifest}, indent=1) + "\n",
                                          encoding="utf-8", newline="\n")
    (out / "g0_manifest.md").write_text(markdown(manifest, prov), encoding="utf-8", newline="\n")
    print(f"wrote {out / 'g0_manifest.json'} and g0_manifest.md")


if __name__ == "__main__":
    main()
