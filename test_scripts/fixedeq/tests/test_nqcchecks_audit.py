"""The driver checks added after the 2026-10-29 release audit: refusals by their reason, the +I+G
restore line, and linked runs given --all-checkpoint that wrote no numbered checkpoint."""
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "nqc"))
import nqcchecks as nc  # noqa: E402


class FakeRun:
    """What the checks read from a Runner's run: exit code, stdout, arguments and final checkpoint."""

    def __init__(self, root, code=0, stdout="", args=(), checkpoint=None):
        self.code, self.args, self._checkpoint = code, list(args), checkpoint or {}
        self.dir = root / "run"
        self.dir.mkdir()
        (self.dir / "stdout.txt").write_text(stdout, encoding="utf-8")

    def stdout(self):
        return (self.dir / "stdout.txt").read_text(encoding="utf-8")

    def checkpoint(self):
        return self._checkpoint


def test_a_refusal_must_name_its_reason(tmp_path):
    run = FakeRun(tmp_path, code=2, stdout="ERROR: NQC is not supported with --eigen\n")
    assert nc.refused("x", run, "--eigen")["passed"]
    assert not nc.refused("x", run, "--lie-markov")["passed"]


def test_a_run_that_finishes_is_not_a_refusal(tmp_path):
    run = FakeRun(tmp_path, code=0, stdout="ERROR: NQC is not supported with --eigen\n")
    assert not nc.refused("x", run, "--eigen")["passed"]


def starts_text(starts, restored):
    lines = [f"Init pinv, alpha: 0.1, 1.0 / Estimate: 0.1, 1.0 / LogL: {s}" for s in starts]
    return "\n".join(lines + [f"Optimal pinv,alpha: 0.1, 1.0 / LogL: {restored}"]) + "\n"


def test_restore_line_is_decision_038s_capped_at_the_earlier_line():
    assert nc.restore_line(-7000.0) == 1e-8 * 7000.0
    assert nc.restore_line(-1e5) == 5e-4


def test_restore_check_holds_the_restored_score_to_the_line():
    starts = [-7000.5, -7000.0, -7000.2]
    assert nc.restore_check("x", starts_text(starts, -7000.0 + 1e-5))["passed"]
    assert not nc.restore_check("x", starts_text(starts, -7000.0 + 1e-4))["passed"]
    assert not nc.restore_check("x", starts_text(starts, "nan"))["passed"]
    assert not nc.restore_check("x", starts_text([-1e5 - 0.5, -1e5], -1e5 + 6e-4))["passed"]
    assert nc.restore_check("x", starts_text([-1e5 - 0.5, -1e5], -1e5 + 4e-4))["passed"]


def test_restore_check_reports_whether_the_best_start_was_the_last():
    assert nc.restore_check("x", starts_text([-7000.5, -7000.0, -7000.2], -7000.0))["informative"]
    assert not nc.restore_check("x", starts_text([-7000.5, -7000.0], -7000.0))["informative"]


def test_block_restore_check_judges_each_block():
    text = starts_text([-7000.5, -7000.0, -7000.2], -7000.0) + starts_text([-7000.1, -7000.3], -7000.1 + 1e-3)
    row = nc.block_restore_check("x", text)
    assert not row["passed"] and row["failing"] == [1] and row["blocks"] == 2 and row["informative"] == 2


def linked_checkpoint():
    keys = {"target": "0.5,0.5", "coordinates": "0.1,0.2", "references": "1,0", "residuals": "0,0,0,0,0,0,0"}
    return {f"PartitionModel!{name}!ModelNonrevFixedEq!{k}": v for name in ("NQC", "partA", "partB")
            for k, v in keys.items()}


def first_row(tmp_path, args):
    c = SimpleNamespace(pi=np.array([0.5, 0.5]), ref=[1, 0])
    run = FakeRun(tmp_path, args=args, checkpoint=linked_checkpoint())
    return nc.linked_invariants("x", run, c, 2)[0]


def test_a_run_without_all_checkpoint_is_judged_at_its_final_checkpoint(tmp_path):
    row = first_row(tmp_path, ["-S", "parts.nex"])
    assert row["passed"] and "the final checkpoint" in row["check"]


def test_a_run_given_all_checkpoint_must_have_written_numbered_checkpoints(tmp_path):
    row = first_row(tmp_path, ["-S", "parts.nex", "-cptime", "0", "--all-checkpoint"])
    assert not row["passed"] and row["numbered_expected"]
