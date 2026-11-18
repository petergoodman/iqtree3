"""Exploratory: every IQ-TREE run whose ending rule an earlier driver labelled through convergence.run_state,
relabelled from IQ-TREE's printed rounds alone, against the label the earlier logic gives on the same output.

Until this script's commit, run_state appended the checkpoint's (or BEST SCORE's) final score to the printed
rounds before testing the last gain, so where IQ-TREE prints 3 decimals the last "gain" was rounding and a
fit stopped at its cap could read "round gain below the threshold". Drivers that used it: S3c's stopping study
(explore/convergence.py), S5's comparisons (reference/reference.py) and SIM-1's pilot (sim/pipeline.py, whose
-nparam cap run_state did not know; pipeline.round_ending relabels those). Each run directory given is one of
those drivers' outputs; its records stay as written.

  python test_scripts/fixedeq/explore/round_endings_recheck.py --git git.exe [--out DIR] RUN_DIR [RUN_DIR ...]
"""
import argparse
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "sim"))
import pipeline as pl  # noqa: E402  (also puts the drivers on sys.path)
import common  # noqa: E402
from convergence import CASES, ended_by, run_state  # noqa: E402
from nqcchecks import check  # noqa: E402
from regress import REPO, read  # noqa: E402

SCRIPT = Path(__file__).resolve().relative_to(REPO).as_posix()
CONVERGENCE = re.compile(r"^([ABC])-(nonrev|nqc-\S+)-(default|relaxed)$")


def settings(setting, rounds):
    """run_state's -me, -S round cap and -nparam for a setting."""
    return ((pl.ME_DEFAULT, pl.LOOP_DEFAULT, pl.NPARAM_DEFAULT) if setting == "default"
            else (pl.ME_RELAXED, rounds, rounds + 2))


def relabel(run, mode, is_nqc, setting, rounds):
    """The fixed label and the earlier one (the appended final's difference as the last gain)."""
    state = run_state(run, mode, is_nqc, setting, rounds)
    if state.get("rounds") is None:
        return None
    old = ended_by(mode, state["rounds"], state["gains"], *settings(setting, rounds))
    return {"mode": mode, "setting": setting, "rounds": state["rounds"], "recomputed": state["ended_by"],
            "earlier_logic": old, "last_gains": [round(g, 6) for g in state["gains"][-3:]],
            "final_lnl": state["final_lnl"]}


def runs_of(d):
    """(name, Run, mode, is_nqc, setting, rounds, cap) for each labelled run in a driver's output directory."""
    res = json.loads(read(d / "results.json"))
    prov = res["provenance"]
    script = prov.get("script", "")
    if script.endswith("convergence.py"):
        rounds = prov["settings"]["rounds"]
        for sub in sorted(p for p in d.iterdir() if p.is_dir() and CONVERGENCE.match(p.name)):
            key, config, setting = CONVERGENCE.match(sub.name).groups()
            yield sub.name, pl.Run(sub, sub.name, [], 0, None), CASES[key][1], config != "nonrev", setting, rounds, None
    elif script.endswith("reference.py"):
        mode = "-S" if prov["case"] == "twopart" else "single"
        for setting in ("default", "relaxed"):
            sub = d / f"nqc-{setting}"
            if sub.is_dir():
                yield sub.name, pl.Run(sub, sub.name, [], 0, None), mode, True, setting, prov["rounds"], None
    elif script.endswith("pipeline.py"):
        for name, rec in res.get("records", {}).items():
            fit = d / name / "fit"
            if rec.get("exit_code") == 0 and (fit / "stdout.txt").exists():
                yield (name, pl.Run(fit, "fit", [], 0, None), "-q", not name.startswith("variant-gtr20"),
                       "default", pl.ROUNDS, rec.get("rounds_cap") or 0)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("dirs", nargs="+")
    ap.add_argument("--out")
    ap.add_argument("--git", default="git")
    ap.add_argument("--allow-dirty", action="store_true")
    a = ap.parse_args()
    prov = common.provenance(a, SCRIPT, "round-endings-recheck", "exploratory", dirs=a.dirs)
    out = common.new_out(a.out, prov)
    rows = []
    for d in (Path(p).expanduser().resolve() for p in a.dirs):
        for name, run, mode, is_nqc, setting, rounds, cap in runs_of(d):
            r = relabel(run, mode, is_nqc, setting, rounds)
            if r is None:
                rows.append(check(f"{d.name}/{name}", None, note="no printed rounds"))
                continue
            if cap is not None:                     # the pilot: its -nparam cap, which run_state did not know
                r["cap"] = cap or pl.NPARAM_DEFAULT - 1
                r["recomputed"] = pl.round_ending(run, rounds=cap or None)
            r["changed"] = r["recomputed"] != r["earlier_logic"]
            rows.append(check(f"{d.name}/{name}", None, **r))
            print(f"{d.name}/{name}: {r['recomputed']} (earlier logic: {r['earlier_logic']}), "
                  f"{r['rounds']} rounds, last gains {r['last_gains']}", flush=True)
    changed = [r["check"] for r in rows if r.get("changed")]
    rows.append(check("labels the earlier logic got wrong", None, count=len(changed), runs=changed))
    common.finish(out, prov, "Round endings rechecked from the printed rounds", rows)


if __name__ == "__main__":
    main()
