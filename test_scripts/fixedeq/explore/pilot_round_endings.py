"""Exploratory: how each fit of a SIM-1 pilot run ended its linked rounds, relabelled by
pipeline.round_ending. The pilot at 3719111d labelled every fit "round gain below the threshold" through
convergence.run_state, which at that commit tested IQ-TREE's default cap and, on -q, a last gain that was the
checkpoint's final score minus the last printed round. The run's records stay as written.

  python test_scripts/fixedeq/explore/pilot_round_endings.py PILOT_DIR
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "sim"))
import pipeline as pl  # noqa: E402  (also puts the drivers on sys.path)

d = Path(sys.argv[1]).expanduser()
for name, rec in json.loads((d / "results.json").read_text())["records"].items():
    fit = d / name / "fit"
    if rec.get("exit_code") != 0 or not (fit / "stdout.txt").exists():
        continue
    run = pl.Run(fit, "fit", [], 0, None)
    vals, done = pl.trajectory(run, "-q")
    gains = np.diff(vals)
    cap = rec.get("rounds_cap") or pl.NPARAM_DEFAULT - 1
    print(f"{name}: {done} of {cap} rounds, ended by {pl.round_ending(run, rounds=rec.get('rounds_cap'))}; "
          f"last round gains {', '.join(f'{g:.3g}' for g in gains[-3:])}; recorded '{rec.get('ended_by')}'")
